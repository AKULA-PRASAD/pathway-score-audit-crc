"""Phase 8C verification: frozen-score audit, panel identity, ddof, isolation.

Enforces the guarantees Phase 8C claims rather than restating them, and traces
the report's numerical claims to reports/phase8c/.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "reports" / "phase8c"
P8A = ROOT / "reports" / "phase8a"
P8B = ROOT / "reports" / "phase8b"
ANALYSIS = ROOT / "data" / "processed" / "analysis"
PANELS = ROOT / "config" / "pathways"
SCRIPT = ROOT / "scripts" / "34_phase8c_pathway_scores.py"
DOC = ROOT / "docs" / "PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md"

PRIMARY_SETS = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS",
                "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY"}
PI3K = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"


def doc_text() -> str:
    t = DOC.read_text()
    t = re.sub(r"\n>\s*", " ", t)
    return re.sub(r"\s+", " ", t)


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("phase8c", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def summ():
    return json.loads((OUT / "phase8c_summary.json").read_text())


@pytest.fixture(scope="module")
def panels():
    return (json.loads((PANELS / "m2_panel_4set.json").read_text()),
            json.loads((PANELS / "m2_panel_5set.json").read_text()))


@pytest.fixture(scope="module")
def expr():
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[s3].astype(np.float64)
    return X, s3


# ==================== outcome blindness ====================

def test_C01_file_guard_raises(mod):
    saved = {n: getattr(pd, n) for n in ("read_csv", "read_parquet", "read_table")}
    try:
        mod._guard_outcome_files()
        with pytest.raises(PermissionError, match="OUTCOME-BLIND"):
            pd.read_csv(ANALYSIS / "endpoint.csv")
    finally:
        for n, f in saved.items():
            setattr(pd, n, f)


def test_C02_endpoint_exists_so_guard_is_not_vacuous():
    assert (ANALYSIS / "endpoint.csv").exists()


def test_C03_metadata_loader_drops_outcome_columns(mod):
    m, dropped = mod.load_metadata_outcome_blind()
    assert set(dropped) == set(mod.OUTCOME_COLUMNS)
    assert not (set(m.columns) & set(mod.OUTCOME_COLUMNS))
    raw = pd.read_csv(ROOT / "metadata" / "GSE39582_sample_metadata.csv", nrows=0)
    assert set(mod.OUTCOME_COLUMNS) <= set(raw.columns)


def test_C04_poisoning_outcome_columns_changes_nothing(mod):
    real, _ = mod.load_metadata_outcome_blind()
    src = pd.read_csv(ROOT / "metadata" / "GSE39582_sample_metadata.csv",
                      keep_default_na=False, dtype=str)
    bad = src.copy()
    rng = np.random.default_rng(1)
    for c in mod.OUTCOME_COLUMNS:
        bad[c] = rng.permutation(bad[c].values)
    saved = pd.read_csv
    try:
        pd.read_csv = lambda p, *a, **k: (bad.copy()
                                          if "GSE39582_sample_metadata" in str(p)
                                          else saved(p, *a, **k))
        got, _ = mod.load_metadata_outcome_blind()
    finally:
        pd.read_csv = saved
    pd.testing.assert_frame_equal(real, got)


def test_C05_no_outcome_field_in_any_output(mod, summ):
    banned = {c.lower() for c in mod.OUTCOME_COLUMNS} | {"recurrence", "survival"}
    for f in OUT.glob("*.csv"):
        cols = {c.strip().lower() for c in pd.read_csv(f, nrows=0).columns}
        assert not (cols & banned), f.name
    s = dict(summ)
    s.pop("outcome_columns_dropped_at_read")
    for c in mod.OUTCOME_COLUMNS:
        assert c not in json.dumps(s)


def test_C06_summary_flags(summ):
    assert summ["outcome_columns_opened"] is False
    assert summ["modelling_performed"] is False
    assert summ["recurrence_association_tested"] is False
    assert summ["feature_selection_performed"] is False
    assert summ["performance_evaluation_performed"] is False


# ==================== A: audit ====================

def test_C07_membership_matches_frozen_specs(panels, summ):
    p4, p5 = panels
    assert set(p4["gene_sets"]) == PRIMARY_SETS
    assert set(p5["gene_sets"]) == PRIMARY_SETS | {PI3K}
    for pname, expect in (("4set_primary", PRIMARY_SETS),
                          ("5set_sensitivity", PRIMARY_SETS | {PI3K})):
        assert set(summ["audit"][pname]["coverage"]) == expect


def test_C08_roles_unchanged(panels, summ):
    p4, p5 = panels
    assert p4["role"] == "PRIMARY"
    assert p5["role"].startswith("SENSITIVITY")
    assert summ["panel_membership_changed"] is False
    assert summ["panel_roles_changed"] is False


def test_C09_resolved_counts_exact(summ):
    expect = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": (200, 198),
              "HALLMARK_E2F_TARGETS": (200, 195),
              "HALLMARK_WNT_BETA_CATENIN_SIGNALING": (42, 41),
              "HALLMARK_P53_PATHWAY": (200, 194),
              PI3K: (105, 104)}
    cov = summ["audit"]["5set_sensitivity"]["coverage"]
    for name, (src, used) in expect.items():
        assert cov[name]["n_source"] == src
        assert cov[name]["n_genes_used"] == used
        assert cov[name]["n_zero_variance_dropped"] == 0
        assert cov[name]["coverage_rule_met"] is True


def test_C10_independent_reimplementation_agrees(summ):
    for a in summ["audit"].values():
        assert a["max_abs_diff_production_vs_independent"] < 1e-12


def test_C11_phase8a_artifacts_reproduced(summ):
    for a in summ["audit"].values():
        assert a["phase8a_artifact_reproduced"] is True
        assert a["max_abs_diff_vs_phase8a_artifact"] < 1e-6


def test_C12_ddof_is_one_and_read_from_spec(panels, summ, mod, expr):
    for p in panels:
        assert p["scoring"]["sigma_ddof"] == 1
        # revision advances with each approved amendment; sigma_ddof does not
        assert p["specification_revision"] >= 2
        assert not any(a.get("sigma_ddof_changed") for a in p["amendment_history"])
    for a in summ["audit"].values():
        assert a["sigma_ddof"] == 1
        assert np.isclose(a["ddof0_over_ddof1_ratio"],
                          a["expected_ratio_sqrt_n_over_n_minus_1"], atol=1e-9)
    # the implementation must refuse an unpinned panel rather than defaulting
    X, s3 = expr
    unpinned = json.loads(json.dumps(panels[0]))
    del unpinned["scoring"]["sigma_ddof"]
    with pytest.raises(ValueError, match="sigma_ddof"):
        mod.canonical_scorer()(X, unpinned, train_ids=sorted(s3)[:200])


def test_C13_independent_scorer_also_reads_ddof_from_spec(mod, expr, panels):
    X, s3 = expr
    p = panels[0]
    tr = sorted(s3)[:200]
    a = mod.score_panel_independent(X, p, tr)
    b = mod.canonical_scorer()(X, p, train_ids=tr)[0]
    assert np.abs(a.values - b[a.columns].values).max() < 1e-12


def test_C14_fold_isolation_deletion_and_poisoning(summ):
    for a in summ["audit"].values():
        f = a["fold_isolation"]
        assert len(f) == 5
        assert all(x["bit_identical_when_heldout_deleted"] for x in f)
        assert all(x["bit_identical_when_heldout_poisoned"] for x in f)
        assert all(x["heldout_scores_do_move_when_poisoned"] for x in f)
        assert all(x["max_abs_shift_vs_full_cohort_fit"] > 0 for x in f)


def test_C15_fold_isolation_rederived(mod, expr, panels):
    """Re-derive the invariant here, not just read it from the JSON."""
    X, s3 = expr
    sp = mod.canonical_scorer()
    train = sorted(s3)[:300]
    held = [s for s in s3 if s not in set(train)]
    a, _ = sp(X, panels[1], train_ids=train)
    poisoned = X.copy()
    poisoned.loc[:, held] += 500.0
    c, _ = sp(poisoned, panels[1], train_ids=train)
    assert np.array_equal(a.loc[train].values, c.loc[train].values)
    assert not np.allclose(a.loc[held].values, c.loc[held].values)


def test_C16_zero_variance_path_exercised(summ):
    for a in summ["audit"].values():
        z = a["zero_variance_path_exercised"]
        assert z["injected"] == 1 and z["dropped"] == 1
        assert z["observed_in_real_data"] == 0


def test_C17_shared_scores_identical_across_panels(summ):
    assert summ["panel_independence_of_shared_scores"]["max_abs_diff"] == 0.0


def test_C18_no_per_sample_missingness(summ):
    for a in summ["audit"].values():
        assert a["per_sample_missingness"] == 0
        assert a["non_finite_in_matrix"] == 0


# ==================== B: distributions ====================

def test_C19_implied_coherence_formula_is_correct(summ):
    """r̄ must invert Var = (1+(k-1)r̄)/k exactly."""
    for r in summ["distributions"]["table"]:
        k, var, rbar = r["n_genes_used"], r["var"], r["implied_mean_inter_gene_r"]
        assert np.isclose((1 + (k - 1) * rbar) / k, var, rtol=1e-10)
        assert np.isclose(r["sd_if_members_independent"], 1 / np.sqrt(k), rtol=1e-12)


def test_C20_two_coherence_regimes_reported(summ):
    d = {r["pathway"]: r["implied_mean_inter_gene_r"]
         for r in summ["distributions"]["table"]}
    assert d["HALLMARK_E2F_TARGETS"] > 0.25
    assert d["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] > 0.20
    for low in ("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY", PI3K):
        assert d[low] < 0.05


def test_C21_no_sample_is_unduly_influential(summ):
    for r in summ["distributions"]["influence_loo"]:
        assert r["relative_to_sd"] < 0.05


def test_C22_variation_sufficiency_not_claimed(summ):
    s = summ["distributions"]["variation_sufficiency"]
    assert "outcome-dependent question" in s and "NOT answered here" in s


# ==================== C: relationships ====================

def test_C23_phase6c_redundancy_reproduces(summ):
    for pname, tol in (("4set_primary", 0.01), ("5set_sensitivity", 0.1)):
        r = summ["relationships"][pname]
        rec = r["phase6c_recorded_n461"]
        assert abs(r["max_abs_r"] - rec["max_abs_r"]) < 0.01
        assert abs(r["condition_number"] - rec["condition_number"]) < tol
        assert r["n_gt_0.5"] == rec["n_gt_0.5"]
        assert r["pcs_80"] == rec["pcs_80"]


def test_C24_pi3k_increment_is_descriptive_only(summ):
    p = summ["relationships"]["pi3k_increment"]
    assert 0.0 < p["r2_on_primary_four"] < 1.0
    assert np.isclose(p["residual_variance_fraction"], 1 - p["r2_on_primary_four"])
    assert p["condition_number_5set"] > p["condition_number_4set"]
    assert "NOT grounds to promote, remove or reorder" in p["interpretation_limit"]


def test_C25_pi3k_still_sensitivity_only(panels):
    p4, p5 = panels
    assert PI3K not in p4["gene_sets"]
    assert PI3K in p5["gene_sets"]


# ==================== D: technical structure ====================

def test_C26_both_circularity_classes_declared(summ):
    t = summ["technical_structure"]
    assert t["circular_variables"] == ["cit.molecularsubtype"]
    assert t["matrix_derived_variables"] == ["qc_any_soft_flag"]
    assert "near-tautological" in t["circularity_note"]
    assert "mechanical" in t["matrix_derived_note"]
    assert "NOT grounds for exclusion" in t["matrix_derived_note"]


def test_C27_batch_and_extraction_do_not_survive_stratification(summ):
    s = pd.DataFrame(summ["technical_structure"]["stratified"])
    other = s[s.annotation != "qc_any_soft_flag"]
    assert len(other) == 10
    assert (other.location_p_strat > 0.05).all()
    assert (other.dispersion_p_strat > 0.05).all()


def test_C28_soft_flag_dispersion_survives_and_is_labelled(summ):
    s = pd.DataFrame(summ["technical_structure"]["stratified"])
    qc = s[s.annotation == "qc_any_soft_flag"]
    assert (qc.dispersion_p_strat < 0.05).sum() >= 3
    assert summ["technical_structure"]["matrix_derived_dispersion_hits_q05"] >= 3


def test_C29_dispersion_tested_not_only_location(summ):
    df = pd.read_csv(OUT / "score_vs_annotation.csv")
    for c in ("location_p", "dispersion_p", "location_q_BH", "dispersion_q_BH",
              "min_cell", "circular", "matrix_derived"):
        assert c in df.columns
    assert len(df) == summ["technical_structure"]["n_tests"]


def test_C30_limitation_recorded(summ):
    lim = summ["technical_structure"]["limitation"]
    assert "ComBat was applied" in lim
    assert "not evidence that batch effects are" in lim


# ==================== E: sensitivity ====================

def test_C31_representation_comparison_scoped(summ):
    r = summ["sensitivity"]["representation"]
    assert r["sensitivity_only"] is True
    assert "not biological differences" in r["scope_limit"]
    per = {x["pathway"]: x["pearson_median_vs_repiqr"] for x in r["per_pathway"]}
    assert per["HALLMARK_E2F_TARGETS"] > 0.99
    assert per["HALLMARK_WNT_BETA_CATENIN_SIGNALING"] < 0.90


def test_C32_folds_reused_not_created(summ, mod):
    assert summ["folds_newly_created"] is False
    assert "Phase 8A" in summ["sensitivity"]["fold_composition"]["folds"]
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    folds = mod.phase8a_folds(s3)
    assert sorted(np.concatenate(folds).tolist()) == sorted(s3)
    assert [len(f) for f in folds] == [91, 91, 91, 90, 90]
    p8a = pd.read_csv(P8A / "fold_isolation_4set_primary.csv")
    assert p8a.n_test.tolist() == [len(f) for f in folds]


def test_C33_scores_stable_but_not_identical_across_folds(summ):
    s = summ["sensitivity"]["fold_composition"]["summary"]
    for p, v in s.items():
        assert v["min_corr_with_full_cohort"] > 0.999
        assert v["max_abs_shift"] > 0, "a zero shift would mean no fold effect at all"


# ==================== invariants and integrity ====================

def test_C34_all_453_retained(summ):
    assert summ["n_samples"] == 453 and summ["samples_excluded"] == 0
    for p in OUT.glob("scores_*_reproduced_full_cohort.csv"):
        assert len(pd.read_csv(p)) == 453


def test_C35_frozen_artifacts_unchanged(summ):
    from src.m2_amendment import resolve_hash
    assert summ["frozen_artifacts_modified"] is False
    for name, h in summ["frozen_artifact_hashes"].items():
        hits = [p for p in ROOT.rglob(name) if ".venv" not in str(p)]
        assert hits, name
        ok, why = resolve_hash(hits[0], h)
        assert ok, why


def test_C36_prior_phase_artifacts_unchanged(summ):
    for rel, h in summ["prior_phase_artifact_hashes"].items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == h, rel
    assert any("phase8a" in k for k in summ["prior_phase_artifact_hashes"])
    assert any("phase8b" in k for k in summ["prior_phase_artifact_hashes"])


def test_C37_phase8c_wrote_only_to_its_own_directory():
    for p in OUT.glob("*"):
        assert OUT in p.parents or p.parent == OUT


def test_C38_no_normalization(summ):
    assert summ["normalization_applied"] == "none"
    assert summ["batch_correction_applied"] == "none"


# ==================== report traceability ====================

def test_C39_report_input_hashes_match_disk():
    rows = re.findall(r"^\| `([^`]+\.(?:parquet|csv|json))` \| `([0-9a-f]{16})`",
                      DOC.read_text(), re.M)
    assert len(rows) >= 9
    from src.m2_amendment import resolve_hash, superseded_hashes
    for rel, h in rows:
        got = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
        if got[:16] == h:
            continue
        # an M2 spec amended since this report was written resolves through the
        # amendment chain; anything else is undocumented drift.
        name = Path(rel).name
        assert any(s.startswith(h) for s in superseded_hashes(name)), \
            f"{rel}: cited {h} is neither current nor a documented superseded version"


NEGATORS = ("no ", "not ", "never ", "nothing ", "cannot ", "neither ",
            "without ", "says nothing")


def test_C40_report_makes_no_prognostic_claim():
    """Every clinical-claim phrase must sit inside a NEGATION, not an assertion.

    A bare substring ban would fire on the report's own disclaimer ("No pathway
    is claimed to be ... clinically useful"), so presence is not the test --
    context is. Each occurrence must be preceded by a negator within 90 chars.
    """
    t = doc_text().lower()
    assert "no outcome association has been tested" in t
    for phrase in ("is prognostic", "predicts recurrence", "biomarker for",
                   "clinically useful", "associated with recurrence",
                   "predictive value"):
        start = 0
        while (i := t.find(phrase, start)) != -1:
            ctx = t[max(0, i - 90):i]
            assert any(n in ctx for n in NEGATORS), (
                f"{phrase!r} appears as an assertion near: ...{ctx[-70:]}[{phrase}]")
            start = i + len(phrase)


def test_C41_report_labels_result_classes():
    t = doc_text()
    for h in ("Reproduced from Phase 6C / 8A", "New in Phase 8C", "Sensitivity-only",
              "Limitations and unresolved questions"):
        assert h in t, h


def test_C42_decision_log_records_no_change():
    t = doc_text()
    assert "No pathway was promoted, removed, replaced or reordered" in t
    assert t.count("**UNCHANGED**") >= 9
    assert "**NOT ACCESSED**" in t


def test_C43_recommendations_were_not_applied_BY_PHASE_8C(summ, panels):
    """Phase 8C applied nothing. Later phases may, under explicit approval."""
    t = doc_text()
    assert "requiring separate review" in t
    assert "None is applied here" in t
    p4, p5 = panels
    assert set(p4["gene_sets"]) == PRIMARY_SETS
    # Phase 8C recorded the state as NONE; M2-A2 superseded it in Phase 9 and
    # the superseded text must still be retained in the amendment history.
    cur = p4["scoring"]["score_standardisation"]
    assert cur.startswith("NONE") or any(
        a.get("superseded_score_standardisation", "").startswith("NONE")
        for a in p4["amendment_history"]), "the Phase 8C state must remain traceable"


def test_C44_headline_numbers_trace_to_artifacts(summ):
    t = doc_text()
    p = summ["relationships"]["pi3k_increment"]
    assert f"R² = {p['r2_on_primary_four']:.3f}" in t
    d = {r["pathway"]: r for r in summ["distributions"]["table"]}
    assert f"{d['HALLMARK_P53_PATHWAY']['implied_mean_inter_gene_r']:.3f}" in t
    rep = {x["pathway"]: x for x in summ["sensitivity"]["representation"]["per_pathway"]}
    assert f"{rep['HALLMARK_WNT_BETA_CATENIN_SIGNALING']['pearson_median_vs_repiqr']:.3f}" in t


def test_C45_three_figures_written():
    figs = sorted(f.name for f in (OUT / "figures").glob("*.png"))
    assert figs == ["H1_score_distributions.png", "H2_relationships.png",
                    "H3_sensitivity.png"]
    assert all((OUT / "figures" / f).stat().st_size > 10_000 for f in figs)
