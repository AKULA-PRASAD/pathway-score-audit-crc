"""Phase 8A verification: outcome-blind EDA invariants.

These tests check the things that would be silently wrong if the EDA drifted:
outcome blindness, retention of all 453 samples, absence of feature selection,
and that the frozen M2 scoring rule was applied as frozen -- including the
leakage invariant that training-fold parameters must not depend on held-out
samples.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "reports" / "phase8a"
ANALYSIS = ROOT / "data" / "processed" / "analysis"
PANELS = ROOT / "config" / "pathways"


def _load_script():
    spec = importlib.util.spec_from_file_location(
        "phase8a", ROOT / "scripts" / "32_phase8a_eda.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def mod():
    return _load_script()


@pytest.fixture(scope="module")
def summary():
    return json.loads((OUT / "phase8a_summary.json").read_text())


@pytest.fixture(scope="module")
def expr():
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[s3]
    return X.astype(np.float64), s3


# --------------------------------------------------------------- outcome blindness
def test_E01_outcome_guard_blocks_endpoint_file(mod):
    """The guard must actually raise, not merely be declared in a docstring."""
    saved = {n: getattr(pd, n) for n in ("read_csv", "read_parquet", "read_table")}
    try:
        mod._guard_outcome_files()
        with pytest.raises(PermissionError, match="OUTCOME-BLIND"):
            pd.read_csv(ANALYSIS / "endpoint.csv")
        # a non-outcome file must still read normally through the patched fn
        assert len(pd.read_csv(ANALYSIS / "cohort_membership.csv")) == 461
    finally:
        for n, f in saved.items():
            setattr(pd, n, f)


def test_E02_endpoint_csv_exists_and_was_simply_not_read():
    """Guard-by-absence would be a false pass: the file must be there."""
    assert (ANALYSIS / "endpoint.csv").exists()


def test_E03_summary_records_outcome_blindness(summary):
    assert summary["outcome_columns_opened"] is False
    assert summary["survival_modelling_performed"] is False
    assert summary["recurrence_association_tested"] is False
    assert summary["feature_selection_performed"] is False


def test_E04_no_output_file_contains_outcome_fields():
    """No Phase 8A artifact may carry a recurrence/survival column."""
    banned = {"rfs.event", "rfs.delay", "event", "time", "recurrence", "os.delay",
              "os.event", "status"}
    for f in OUT.glob("*.csv"):
        cols = {c.strip().lower() for c in
                pd.read_csv(f, nrows=0).columns}
        assert not (cols & banned), f"{f.name} exposes {cols & banned}"


# --------------------------------------------------------------- retention / selection
def test_E05_all_453_samples_retained(summary):
    assert summary["n_samples"] == 453
    assert summary["samples_excluded"] == 0
    assert len(pd.read_csv(OUT / "pca_scores_gene_level.csv")) == 453
    assert len(pd.read_csv(OUT / "sample_distribution_summary.csv")) == 453


def test_E06_no_feature_selection_all_genes_reported(summary):
    g = pd.read_csv(OUT / "gene_variability.csv")
    assert len(g) == 21815 == summary["n_genes"]
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")
    assert set(g.gene_symbol) == set(X.index), "gene table is not the full index"


def test_E07_no_normalization_or_batch_correction(summary):
    assert summary["normalization_applied"] == "none"
    assert summary["batch_correction_applied"] == "none"
    assert summary["matrix_used_as_deposited"] is True


# --------------------------------------------------------------- frozen artifacts
def test_E08_frozen_artifacts_unmodified(summary):
    assert summary["frozen_artifacts_modified"] is False
    from src.m2_amendment import resolve_hash

    for name, h in summary["frozen_artifact_hashes"].items():
        hits = list(ROOT.rglob(name))
        hits = [p for p in hits if ".venv" not in str(p)]
        assert hits, f"{name} not found"
        ok, why = resolve_hash(hits[0], h)
        assert ok, why


def test_E09_m2_panels_match_the_freeze_record():
    rec = json.loads((PANELS / "M2_FREEZE_RECORD.json").read_text())
    for role in ("primary", "sensitivity"):
        s = rec["specifications"][role]
        got = hashlib.sha256((ROOT / s["file"]).read_bytes()).hexdigest()
        assert got == s["sha256"], f"{role} panel drifted"


# --------------------------------------------------------------- frozen scoring rule
def test_E10_score_panel_refuses_implicit_training_population(mod, expr):
    """train_ids must have no default -- a full-cohort fit cannot be implicit."""
    X, _ = expr
    panel = json.loads((PANELS / "m2_panel_4set.json").read_text())
    with pytest.raises(TypeError):
        mod.score_panel(X, panel)                       # no train_ids at all
    with pytest.raises(ValueError, match="mandatory"):
        mod.score_panel(X, panel, train_ids=None)


def test_E11_training_parameters_ignore_heldout_samples(mod, expr):
    """The leakage invariant, re-derived here rather than read from the CSV."""
    X, s3 = expr
    panel = json.loads((PANELS / "m2_panel_4set.json").read_text())
    train = sorted(s3)[:300]
    held = [s for s in s3 if s not in set(train)]
    full, _ = mod.score_panel(X, panel, train_ids=train)
    trimmed, _ = mod.score_panel(X.drop(columns=held), panel, train_ids=train)
    assert np.array_equal(full.loc[train].values, trimmed.loc[train].values)


def test_E12_full_cohort_scoring_would_have_leaked(summary):
    """If the shift were 0, the fold-isolation test would prove nothing."""
    for p in summary["pathways"].values():
        assert p["fold_isolation"]["all_folds_bit_identical_when_heldout_deleted"]
        assert p["fold_isolation"]["is_phase6b_resampling_scheme"] is False
        assert p["fold_isolation"]["max_abs_score_shift_vs_full_cohort_fit"] > 0


def test_E13_full_cohort_scores_labelled_exploratory(summary):
    for p in summary["pathways"].values():
        assert p["full_cohort_scores_are_exploratory_only"] is True
        assert p["full_cohort_scores_used_as_model_input"] is False
    assert sorted(f.name for f in OUT.glob("pathway_scores_*")) == [
        "pathway_scores_4set_primary_EXPLORATORY_full_cohort.csv",
        "pathway_scores_5set_sensitivity_EXPLORATORY_full_cohort.csv"]


def test_E14_frozen_coverage_rule_met_and_unchanged(summary):
    for name, p in summary["pathways"].items():
        assert p["all_coverage_rules_met"] is True
        assert p["scores_dropped"] == []
        for gs, c in p["coverage"].items():
            assert c["frac_of_source_resolved"] >= 0.80, gs
            assert c["n_genes_used"] >= 20, gs
            assert c["n_zero_variance_dropped"] == 0, gs


def test_E15_panel_membership_is_unchanged(summary):
    rec = json.loads((PANELS / "M2_FREEZE_RECORD.json").read_text())
    assert set(summary["pathways"]["4set_primary"]["coverage"]) == set(
        rec["specifications"]["primary"]["gene_sets"])
    assert set(summary["pathways"]["5set_sensitivity"]["coverage"]) == set(
        rec["specifications"]["sensitivity"]["gene_sets"])


def test_E16_pi3k_is_sensitivity_only(summary):
    p4 = summary["pathways"]["4set_primary"]["coverage"]
    p5 = summary["pathways"]["5set_sensitivity"]["coverage"]
    assert "HALLMARK_PI3K_AKT_MTOR_SIGNALING" not in p4
    assert "HALLMARK_PI3K_AKT_MTOR_SIGNALING" in p5


def test_E17_ddof_gap_closed_by_amendment_M2_A1(summary):
    """Phase 8A found sigma's ddof unpinned (U5); amendment M2-A1 closed it.

    The constant is still verified, because the amendment's stated numerical
    consequence is exactly that factor -- pinning the convention must not
    change what the measurement says.
    """
    for p in summary["pathways"].values():
        d = p["ddof_sensitivity"]
        assert d["freeze_record_pins_ddof"] is True
        assert d["pinned_sigma_ddof"] == 1
        assert d["pinned_by_amendment"] == "M2-A1"
        assert d["read_from_specification_not_library_default"] is True
        assert d["amendment_resolves_phase8a_U5"] is True
        assert np.isclose(d["ratio_ddof0_over_ddof1_max"],
                          d["expected_constant_factor_sqrt_n_over_n_minus_1"],
                          atol=1e-9)


# --------------------------------------------------------------- PCA
def test_E18_pca_is_deterministic(mod, expr):
    """Re-running the documented method must reproduce the recorded spectrum."""
    X, _ = expr
    V = X.values
    Xc = V - V.mean(axis=1, keepdims=True)
    S = np.linalg.svd(Xc, full_matrices=False, compute_uv=False)
    ve = (S ** 2) / (S ** 2).sum()
    rec = pd.read_csv(OUT / "pca_variance_explained.csv")
    assert np.allclose(ve[:20], rec.variance_explained.values, atol=1e-10)


def test_E19_pca_rotation_not_promoted_to_model_input(summary):
    assert summary["pca"]["full_cohort_rotation_is_descriptive_only"] is True


def test_E20_gene_level_and_probe_level_pca_agree(summary):
    """A large divergence would mean median aggregation changed sample structure."""
    g = np.array(summary["pca"]["variance_explained_pc1_20"][:5])
    p = np.array(summary["pca"]["qc_phase7c_probe_level_pc1_5"])
    assert np.abs(g - p).max() < 0.02


# --------------------------------------------------------------- 2d findings
def test_E21_combat_mod_mmr_concordance_is_exact(summary):
    t = summary["technical_vs_biological"]["combat_mod_meaning"]
    assert t["s3_discordant"] == 0
    assert t["full_series"]["discordant"] == 0
    assert t["full_series"]["n_arrays_with_mmr"] == 519
    ct = pd.read_csv(OUT / "combat_mod_vs_mmr.csv", index_col=0)
    assert ct.loc["modcan1", "pMMR"] == 0 and ct.loc["modcan2", "dMMR"] == 0


def test_E22_concordance_is_not_overclaimed_as_identity(summary):
    """Perfect concordance is evidence, not proof -- the record must say so."""
    t = summary["technical_vs_biological"]["combat_mod_meaning"]
    assert t["inference_status"].startswith("HYPOTHESIS")
    assert "not_used_to_impute" in t


def test_E23_batch_pc_association_vanishes_within_stratum(summary):
    cond = pd.DataFrame(summary["technical_vs_biological"]["pc_batch_conditional"])
    pm = cond[(cond.stratum == "modcan2") & (cond.pc.isin(
        ["PC1", "PC2", "PC3", "PC4", "PC5"]))]
    assert (pm.p_value > 0.05).all(), "a batch association survived stratification"
    assert (pm.epsilon_squared < 0.01).all()


def test_E24_batch_finding_does_not_claim_absence_of_batch_effect(summary):
    interp = summary["technical_vs_biological"]["interpretation"]
    assert "does NOT prove absence" in interp
    assert "ComBat was already" in interp


def test_E25_no_exclusion_proposed_by_sex_control(summary):
    s = summary["sex_positive_control"]
    assert s["action"].startswith("REPORT ONLY")
    assert s["adjudication"].startswith("NOT ADJUDICATED")
    assert s["n_with_recorded_sex"] == 453


# --------------------------------------------------------------- deliverables
def test_E26_all_six_figures_written():
    figs = sorted(f.name for f in (OUT / "figures").glob("*.png"))
    assert figs == ["F1_expression_distributions.png", "F2_gene_variability.png",
                    "F3_sample_correlation.png", "F4_pca.png",
                    "F5_pathway_scores.png", "F6_pathway_vs_pca.png"]
    assert all((OUT / "figures" / f).stat().st_size > 10_000 for f in figs)


def test_E27_report_exists_and_separates_interpretation_tiers():
    t = (ROOT / "docs" / "PHASE8A_EXPLORATORY_ANALYSIS.md").read_text()
    for h in ("Direct observations", "Biological hypotheses",
              "Technical explanations", "Unresolved questions"):
        assert h in t, h
    assert "No biomarker is identified" in t
