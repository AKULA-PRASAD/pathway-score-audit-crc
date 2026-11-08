"""Phase 8B verification: outcome blindness, reproduction, and claim integrity.

Every numerical claim in docs/PHASE8B_PCA_BATCH_INVESTIGATION.md is checked
against reports/phase8b/, and the outcome-blindness guarantees are proven
rather than asserted.
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

OUT = ROOT / "reports" / "phase8b"
P8A = ROOT / "reports" / "phase8a"
ANALYSIS = ROOT / "data" / "processed" / "analysis"
PANELS = ROOT / "config" / "pathways"
SCRIPT = ROOT / "scripts" / "33_phase8b_pca_batch.py"
DOC = ROOT / "docs" / "PHASE8B_PCA_BATCH_INVESTIGATION.md"


def doc_text() -> str:
    """Report prose with markdown wrapping normalised away.

    A blockquote line break inserts "\n> " mid-sentence, so a raw substring
    search silently fails on text that IS present. Collapse to single spaces.
    """
    t = DOC.read_text()
    t = re.sub(r"\n>\s*", " ", t)
    return re.sub(r"\s+", " ", t)


def _load():
    spec = importlib.util.spec_from_file_location("phase8b", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def mod():
    return _load()


@pytest.fixture(scope="module")
def summ():
    return json.loads((OUT / "phase8b_summary.json").read_text())


# ==================== outcome blindness, proven ====================

def test_B01_file_guard_actually_raises(mod):
    saved = {n: getattr(pd, n) for n in ("read_csv", "read_parquet", "read_table")}
    try:
        mod._guard_outcome_files()
        with pytest.raises(PermissionError, match="OUTCOME-BLIND"):
            pd.read_csv(ANALYSIS / "endpoint.csv")
        assert len(pd.read_csv(ANALYSIS / "cohort_membership.csv")) == 461
    finally:
        for n, f in saved.items():
            setattr(pd, n, f)


def test_B02_endpoint_file_exists_so_the_guard_is_not_vacuous():
    assert (ANALYSIS / "endpoint.csv").exists()


def test_B03_metadata_loader_drops_outcome_columns(mod):
    m, dropped = mod.load_metadata_outcome_blind()
    assert set(dropped) == set(mod.OUTCOME_COLUMNS)
    assert not (set(m.columns) & set(mod.OUTCOME_COLUMNS))
    raw = pd.read_csv(ROOT / "metadata" / "GSE39582_sample_metadata.csv", nrows=0)
    assert set(mod.OUTCOME_COLUMNS) <= set(raw.columns), (
        "the source file must actually contain them, or this proves nothing")


def test_B04_corrupting_outcome_columns_changes_nothing(mod):
    """The decisive test: poison every outcome value, demand identical output."""
    real, _ = mod.load_metadata_outcome_blind()
    src = pd.read_csv(ROOT / "metadata" / "GSE39582_sample_metadata.csv",
                      keep_default_na=False, dtype=str)
    poisoned = src.copy()
    rng = np.random.default_rng(0)
    for c in mod.OUTCOME_COLUMNS:
        poisoned[c] = rng.permutation(poisoned[c].values)
        poisoned.loc[poisoned.index[:50], c] = "99999"

    saved = pd.read_csv
    try:
        pd.read_csv = lambda p, *a, **k: (
            poisoned.copy() if "GSE39582_sample_metadata" in str(p)
            else saved(p, *a, **k))
        got, _ = mod.load_metadata_outcome_blind()
    finally:
        pd.read_csv = saved
    pd.testing.assert_frame_equal(real, got)


def test_B05_script_never_names_an_outcome_column_outside_the_ban_list(mod):
    src = SCRIPT.read_text()
    decl = re.search(r"OUTCOME_COLUMNS = \([^)]*\)", src).group(0)
    rest = src.replace(decl, "")
    for c in mod.OUTCOME_COLUMNS:
        assert c not in rest, f"{c} referenced outside the ban list"


def test_B06_no_phase8b_output_carries_an_outcome_field(mod):
    banned = {c.lower() for c in mod.OUTCOME_COLUMNS} | {
        "recurrence", "survival", "event", "time_to_event"}
    for f in OUT.glob("*.csv"):
        cols = {c.strip().lower() for c in pd.read_csv(f, nrows=0).columns}
        assert not (cols & banned), f"{f.name}: {cols & banned}"
    # The summary legitimately NAMES the columns it dropped. That record must
    # be the only place an outcome column name appears.
    s = json.loads((OUT / "phase8b_summary.json").read_text())
    recorded = s.pop("outcome_columns_dropped_at_read")
    assert set(recorded) == set(mod.OUTCOME_COLUMNS)
    blob = json.dumps(s)
    for c in mod.OUTCOME_COLUMNS:
        assert c not in blob, f"{c} leaked into the summary outside the drop record"


def test_B07_summary_records_blindness(summ):
    assert summ["outcome_columns_opened"] is False
    assert summ["modelling_performed"] is False
    assert summ["recurrence_association_tested"] is False
    assert summ["feature_selection_performed"] is False
    assert set(summ["outcome_columns_dropped_at_read"]) == {
        "rfs.event", "rfs.delay", "os.event", "os.delay_(months)"}


# ==================== reproduction of Phase 8A ====================

def test_B08_pca_reproduced_by_an_independent_algorithm(summ):
    r = summ["reproduction_of_phase8a"]
    assert r["reproduced"] is True
    assert r["min_signed_correlation_scores_pc1_20"] > 1 - 1e-12
    assert r["max_abs_diff_variance_explained_pc1_20"] < 1e-7
    assert r["max_relative_diff_scores_pc1_20"] < 1e-7


def test_B09_variance_spectrum_matches_phase8a():
    p8a = pd.read_csv(P8A / "pca_variance_explained.csv")
    p8b = pd.read_csv(OUT / "pca_variance_extended.csv")
    assert np.allclose(p8a.variance_explained.values,
                       p8b.variance_explained.values[:20], atol=1e-7)


def test_B10_combat_mod_mmr_reproduced(summ):
    c = summ["combat_mod_vs_mmr"]
    assert c["s3_discordant"] == 0 and c["s3_n_with_mmr"] == 410
    assert c["sdrf_discordant"] == 0 and c["sdrf_n_with_mmr"] == 519
    assert c["reproduces_phase8a"] is True
    ct = pd.read_csv(OUT / "combat_mod_vs_mmr_S3.csv", index_col=0)
    assert ct.loc["modcan1", "pMMR"] == 0 and ct.loc["modcan2", "dMMR"] == 0


def test_B11_concordance_is_not_claimed_as_identity(summ):
    c = summ["combat_mod_vs_mmr"]
    assert "concordance is NOT identity" in c["claim_limit"]
    assert c["mmr_unrecorded_expression_probe"]["used_to_impute"] is False
    assert c["mmr_unrecorded_expression_probe"]["interpretation_status"].startswith(
        "HYPOTHESIS")


# ==================== stratified analysis and its limits ====================

def test_B12_both_strata_reported_with_cell_sizes(summ):
    s = pd.DataFrame(summ["batch_within_mmr_strata"]["per_stratum"])
    assert set(s.stratum) >= {"modcan1", "modcan2", "dMMR", "pMMR"}
    assert (s.min_cell > 0).all()
    small = s[s.stratum == "modcan1"]
    assert small.n.max() == 64 and small.min_cell.min() == 2


def test_B13_dmmr_nominal_association_is_recorded_not_suppressed(summ):
    """Phase 8A tested pMMR only; the dMMR result must not be hidden."""
    s = pd.DataFrame(summ["batch_within_mmr_strata"]["per_stratum"])
    d = s[(s.stratum == "modcan1") & (s.pc == "PC3")].iloc[0]
    assert d.p_value < 0.05 and d.epsilon_squared > 0.1
    assert "nominal PC3 × batch association" in doc_text()


def test_B14_stratified_permutation_preserves_composition(mod):
    """The null must keep every stratum x group cell count fixed."""
    strata = np.array(["a"] * 30 + ["b"] * 20)
    idx = mod.stratified_perm_indices(strata, 50, np.random.default_rng(0))
    for row in idx:
        assert sorted(row.tolist()) == list(range(50))
        assert (strata[row] == strata).all(), "permutation crossed strata"


def test_B15_permutation_test_is_null_for_pc1_to_pc5(summ):
    p = pd.DataFrame(summ["batch_within_mmr_strata"]["stratified_permutation"])
    assert (p.p_permutation > 0.05).all()
    assert (p.n_permutations == 10_000).all()


def test_B16_detection_sensitivity_is_reported_and_bounded(summ):
    b = summ["batch_within_mmr_strata"]
    pw = pd.DataFrame(b["detection_sensitivity"])
    assert pw.loc[pw.delta_in_pc3_sd == 0.0, "power_at_alpha_0.05"].iloc[0] < 0.10
    assert pw.loc[pw.delta_in_pc3_sd == 1.0, "power_at_alpha_0.05"].iloc[0] > 0.90
    assert b["min_detectable_batch05_shift_pc3_sd_at_80pct_power"] == 0.3
    assert "NOT evidence that batch effects are absent" in b["limitation"]
    assert "ComBat was applied" in b["limitation"]


def test_B17_no_claim_that_combat_removed_all_batch_effects():
    t = doc_text()
    assert "No claim is made that ComBat removed all batch effects" in t
    for bad in ("ComBat successfully removed", "batch effects are absent,",
                "no batch effects remain"):
        assert bad not in t, bad


# ==================== batch05 ====================

def test_B18_batch05_decomposition_has_uncertainty(summ):
    d = pd.DataFrame(summ["batch05_pc3"]["observation"])
    assert len(d) == 4
    assert (d.ci95_hi > d.ci95_lo).all()
    b1 = d[(d.combat_mod == "modcan1") & (d.group == "batch05")].iloc[0]
    o1 = d[(d.combat_mod == "modcan1") & (d.group == "other batches")].iloc[0]
    assert b1.ci95_lo < o1.ci95_hi, "CIs must be reported as overlapping"


def test_B19_batch05_keeps_observation_and_explanation_apart(summ):
    s = summ["batch05_pc3"]["separation_of_observation_and_explanation"]
    assert s.startswith("OBSERVATION:")
    assert "EXPLANATION (not established)" in s
    assert "perfectly confounded" in s


# ==================== sex discordance ====================

def test_B20_eight_samples_unchanged_and_unadjudicated(summ):
    s = summ["sex_discordance_review"]
    assert s["n"] == 8
    assert s["status"] == "UNRESOLVED QC OBSERVATION"
    assert s["adjudication"].startswith("NOT ADJUDICATED")
    assert len(s["batches"]) == 5
    assert s["phase8b_panel_size"] > s["phase8a_panel_size"] == 6
    ids = {r["sample_id"] for r in s["samples"]}
    p8a = set(pd.read_csv(P8A / "sex_positive_control.csv")
              .query("classification == 'discordant_both_axes'").sample_id)
    assert ids == p8a, "the eight must be the same eight Phase 8A identified"


def test_B21_soft_flag_enrichment_reported_with_its_limit(summ):
    s = summ["sex_discordance_review"]
    assert s["n_with_any_phase7c_soft_flag"] == 3
    assert s["cohort_soft_flag_rate"] == {"n_flagged": 29, "n_total": 453}
    assert s["soft_flag_enrichment_fisher"]["p_value"] < 0.05
    lim = s["evidence_limit"]
    assert "QUALITY" in lim and "IDENTITY" in lim
    assert "not evidence of a mix-up either way" in lim


def test_B22_no_swap_or_relabel_inferred():
    t = doc_text()
    assert "No independent technical evidence establishes a sample swap" in t
    assert "No sample is excluded or relabelled" in t


# ==================== new structure ====================

def test_B23_new_findings_rechecked_under_stratification(summ):
    r = pd.DataFrame(summ["new_technical_structure"]["stratified_recheck_of_new_findings"])
    assert len(r) == 4
    surv = r[r.p_stratified_permutation < 0.01].test.tolist()
    assert "AFFX_PC1 ~ batch" in surv
    assert "expression PC2 |dev| ~ batch" in surv
    gone = r[r.p_stratified_permutation > 0.05].test.tolist()
    assert all("rna.extraction.method" in g for g in gone)


def test_B24_extraction_method_confounding_is_stated(summ):
    e = summ["new_technical_structure"]["rna_extraction_method"]
    assert e["levels"] == 6
    assert e["confounding_with_batch"]["p_value"] < 1e-10
    assert sum(e["counts"].values()) == 453


def test_B25_affx_limitation_recorded(summ):
    a = summ["new_technical_structure"]["affx_control_probes"]
    assert a["n_controls"] == 62
    assert a["affx_pc1_vs_batch"]["p_value"] < 0.01
    assert "SAME deposited RMA+ComBat processing" in a["note"]


def test_B26_variance_heterogeneity_is_new_and_distinct_from_location(summ):
    v = summ["new_technical_structure"]["variance_heterogeneity"]
    res = pd.DataFrame(v["results"])
    assert res[res.pc == "PC2"].p_value.iloc[0] < 1e-4
    assert "DISPERSION effect" in v["rationale"]
    assert v["new_in_phase8b"] is True


def test_B27_cit_subtype_circularity_flagged(summ):
    w = summ["pc_associations"]["circularity_warning_cit_molecularsubtype"]
    assert "EXPRESSION-DERIVED" in w and "near-tautological" in w


def test_B28_dependent_samples_checked(summ):
    d = summ["new_technical_structure"]["dependent_samples"]
    assert d["n_in_S3"] == 0


# ==================== invariants ====================

def test_B29_all_453_retained(summ):
    assert summ["n_samples"] == 453 and summ["samples_excluded"] == 0
    assert len(pd.read_csv(OUT / "pca_variance_extended.csv")) == 50


def test_B30_m2_untouched_and_pin_in_force(summ):
    m = summ["m2_status"]
    assert m["sigma_ddof_pinned"] == 1
    assert m["read_from_specification"] is True
    assert m["primary_panel_n_sets"] == 4 and m["sensitivity_panel_n_sets"] == 5
    assert m["membership_changed_by_phase8b"] is False
    assert m["roles_changed_by_phase8b"] is False


def test_B31_phase8a_score_artifacts_untouched(summ):
    for name, h in summ["phase8a_score_artifact_hashes"].items():
        assert hashlib.sha256((P8A / name).read_bytes()).hexdigest() == h, name


def test_B32_frozen_artifacts_unchanged(summ):
    from src.m2_amendment import resolve_hash
    assert summ["frozen_artifacts_modified"] is False
    for name, h in summ["frozen_artifact_hashes"].items():
        hits = [p for p in ROOT.rglob(name) if ".venv" not in str(p)]
        assert hits, name
        ok, why = resolve_hash(hits[0], h)
        assert ok, why


def test_B33_no_normalization_or_correction(summ):
    assert summ["normalization_applied"] == "none"
    assert summ["batch_correction_applied"] == "none"
    assert summ["matrix_used_as_deposited"] is True


# ==================== report traceability ====================

def test_B34_report_input_hashes_match_disk():
    rows = re.findall(r"^\| `([^`]+\.(?:parquet|csv|json))` \| `([0-9a-f]{16})`",
                      DOC.read_text(), re.M)   # table rows: raw text is correct here
    assert len(rows) >= 7
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


def test_B35_report_separates_reproduced_from_new():
    t = doc_text()
    assert "## 9. Reproduced vs new" in t
    assert t.count("**Reproduced**") >= 5
    assert t.count("**New**") >= 8


def test_B36_decision_log_records_no_change():
    t = doc_text()
    assert "**No decision was changed by Phase 8B.**" in t
    assert t.count("**UNCHANGED**") >= 8
    assert "**NOT ACCESSED**" in t


def test_B37_headline_numbers_trace_to_artifacts(summ):
    """Spot-check report figures against the machine-readable summary."""
    t = doc_text()
    v = summ["variance_structure"]
    assert f"**{v['horn_n_pcs_above_permuted_p95']}** PCs exceed" in t
    assert f"{v['n_pcs_50pct']} / {v['n_pcs_80pct']} / {v['n_pcs_90pct']}" in t
    b = summ["batch05_pc3"]
    assert f"{b['batch05_pct_modcan1']:.1f} %".replace(" ", " ") in t or \
           f"{b['batch05_pct_modcan1']:.1f}" in t
    s = summ["sex_discordance_review"]
    assert f"OR = {s['soft_flag_enrichment_fisher']['odds_ratio']:.2f}" in t


def test_B38_three_figures_written():
    figs = sorted(f.name for f in (OUT / "figures").glob("*.png"))
    assert figs == ["G1_variance_and_power.png", "G2_batch_mmr_structure.png",
                    "G3_sex_axes.png"]
    assert all((OUT / "figures" / f).stat().st_size > 10_000 for f in figs)
