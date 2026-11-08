"""Verification that M2-A2, P6B-A1 and D8D-4-R1 are APPLIED and correct.

An amendment counts as applied only when the file changes, the behaviour and
the audit chain all confirm it. These tests check all three.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.pathway_scoring import score_panel                       # noqa: E402
from src.score_standardisation import (                            # noqa: E402
    NonEstimableReplicate, ZERO_SD_TOL, standardise_scores_training_fold)

PANELS = ROOT / "config" / "pathways"
CFG = ROOT / "config" / "analysis"
P9 = ROOT / "reports" / "phase9"
P6B = ROOT / "docs" / "PHASE6B_MODELING_PROTOCOL_DRAFT.md"
P9DOC = ROOT / "docs" / "PHASE9_MODELING_PROTOCOL_DRAFT.md"
ANALYSIS = ROOT / "data" / "processed" / "analysis"

PRIMARY_SETS = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS",
                "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY"}
PI3K = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"


def norm(p: Path) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\n>\s*", " ", p.read_text()))


@pytest.fixture(scope="module")
def panels():
    return (json.loads((PANELS / "m2_panel_4set.json").read_text()),
            json.loads((PANELS / "m2_panel_5set.json").read_text()))


@pytest.fixture(scope="module")
def record():
    return json.loads((PANELS / "M2_FREEZE_RECORD.json").read_text())


@pytest.fixture(scope="module")
def scores():
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[s3].astype(np.float64)
    p4 = json.loads((PANELS / "m2_panel_4set.json").read_text())
    sc, _ = score_panel(X, p4, train_ids=sorted(s3)[:300])
    return sc, s3, sorted(s3)[:300]


# ============ M2-A2: specification ============

def test_X01_panels_at_revision_3_with_m2a2(panels, record):
    for d in panels:
        assert d["specification_revision"] == 3
        assert [a["amendment_id"] for a in d["amendment_history"]] == ["M2-A1", "M2-A2"]
    assert record["specification_revision"] == 3
    assert [a["amendment_id"] for a in record["amendment_history"]] == ["M2-A1", "M2-A2"]


def test_X02_score_standardisation_is_training_fold_only(panels, record):
    for d in (*panels, record):
        s = d["scoring"]["score_standardisation"]
        assert s.startswith("TRAINING-FOLD ONLY")
        assert "sample SD (ddof=1)" in s
        assert "applied unchanged to held-out samples" in s
        assert "1e-12" in s and "DROPPED for that fold" in s
        assert "NON-ESTIMABLE" in s and "never silently redrawn" in s.lower() or \
               "never silently redrawn" in s
        assert "unstandardised scores are preserved" in s
    for d in panels:
        assert d["scoring"]["score_standardisation_ddof"] == 1
        assert d["scoring"]["score_standardisation_zero_sd_tol"] == 1e-12
        assert d["scoring"]["score_standardisation_implementation"].startswith(
            "src/score_standardisation.py")


def test_X03_superseded_text_retained(panels):
    for d in panels:
        a = d["amendment_history"][1]
        assert a["amendment_id"] == "M2-A2"
        assert a["superseded_score_standardisation"].startswith("NONE before the Cox")
        assert a["coefficient_interpretation"].startswith(
            "a molecular pathway coefficient is a log cause-specific hazard ratio "
            "per one training-fold standard deviation")
        assert a["outcome_data_used"] is False
        assert a["membership_changed"] is False and a["roles_changed"] is False
        assert a["sigma_ddof_changed"] is False
        assert "Clinical predictor scaling" in a["clinical_predictor_scaling_out_of_scope"]


def test_X04_membership_roles_and_gene_ddof_untouched(panels, record):
    p4, p5 = panels
    assert set(p4["gene_sets"]) == PRIMARY_SETS
    assert set(p5["gene_sets"]) == PRIMARY_SETS | {PI3K}
    assert p4["role"] == "PRIMARY" and p5["role"].startswith("SENSITIVITY")
    assert p4["scoring"]["sigma_ddof"] == p5["scoring"]["sigma_ddof"] == 1
    expect = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": (200, 198),
              "HALLMARK_E2F_TARGETS": (200, 195),
              "HALLMARK_WNT_BETA_CATENIN_SIGNALING": (42, 41),
              "HALLMARK_P53_PATHWAY": (200, 194), PI3K: (105, 104)}
    for d in panels:
        for name, gs in d["gene_sets"].items():
            n_src, n_res = expect[name]
            assert gs["n_source"] == n_src
            assert len(dict.fromkeys(gs["members_resolved"])) == n_res


def test_X05_canonical_scorer_unmodified():
    import inspect
    ref = json.loads((ROOT / "tests/fixtures/premove_scorer_reference.json").read_text())
    got = hashlib.sha256(inspect.getsource(score_panel).encode()).hexdigest()
    assert got == ref["source_sha256"]


# ============ M2-A2: behaviour ============

def test_X06_training_fold_parameters_only(scores):
    sc, s3, tr = scores
    std, d = standardise_scores_training_fold(sc, tr)
    assert d["ddof"] == 1 and d["estimated_on"] == "TRAINING FOLD SAMPLES ONLY"
    for col in std.columns:
        assert abs(std.loc[tr, col].mean()) < 1e-10
        assert abs(std.loc[tr, col].std(ddof=1) - 1.0) < 1e-10
        # held-out is NOT forced to unit variance -- that would be leakage
        held = [s for s in s3 if s not in set(tr)]
        assert abs(std.loc[held, col].std(ddof=1) - 1.0) > 1e-6


def test_X07_train_ids_mandatory(scores):
    sc, _, tr = scores
    with pytest.raises(TypeError):
        standardise_scores_training_fold(sc)
    with pytest.raises(ValueError, match="mandatory"):
        standardise_scores_training_fold(sc, None)


def test_X08_heldout_deletion_invariance(scores):
    sc, s3, tr = scores
    held = [s for s in s3 if s not in set(tr)]
    a, _ = standardise_scores_training_fold(sc, tr)
    b, _ = standardise_scores_training_fold(sc.drop(index=held), tr)
    assert np.array_equal(a.loc[tr].values, b.loc[tr].values)


def test_X09_heldout_poisoning_invariance(scores):
    sc, s3, tr = scores
    held = [s for s in s3 if s not in set(tr)]
    clean, dc = standardise_scores_training_fold(sc, tr)
    pois = sc.copy()
    pois.loc[held] += 1000.0
    dirty, dd = standardise_scores_training_fold(pois, tr)
    assert np.array_equal(clean.loc[tr].values, dirty.loc[tr].values)
    for col in clean.columns:
        assert dc["per_pathway"][col]["train_mean"] == dd["per_pathway"][col]["train_mean"]
        assert dc["per_pathway"][col]["train_sd"] == dd["per_pathway"][col]["train_sd"]
    # held-out MUST move, or the test proves nothing
    assert not np.allclose(clean.loc[held].values, dirty.loc[held].values)


def test_X10_zero_variance_score_dropped_and_recorded(scores):
    sc, _, tr = scores
    victim = sc.columns[0]
    z = sc.copy()
    z.loc[tr, victim] = 5.0
    std, d = standardise_scores_training_fold(z, tr)
    assert victim not in std.columns
    assert d["dropped_pathways"] == [victim] and d["n_dropped"] == 1
    rec = d["per_pathway"][victim]
    assert rec["dropped"] is True and rec["train_sd"] <= ZERO_SD_TOL
    assert rec["zero_sd_tol"] == 1e-12
    assert len(std.columns) == len(sc.columns) - 1


def test_X11_below_tolerance_also_dropped(scores):
    sc, _, tr = scores
    victim = sc.columns[1]
    z = sc.copy()
    z.loc[tr, victim] = 5.0
    z.loc[tr[0], victim] = 5.0 + 1e-15      # SD below ZERO_SD_TOL, not exactly 0
    std, d = standardise_scores_training_fold(z, tr)
    assert 0 < d["per_pathway"][victim]["train_sd"] <= ZERO_SD_TOL
    assert victim not in std.columns


def test_X12_non_finite_marks_replicate_non_estimable(scores):
    sc, _, tr = scores
    for bad in (np.inf, np.nan):
        z = sc.copy()
        z.loc[tr[0], sc.columns[2]] = bad
        with pytest.raises(NonEstimableReplicate) as ei:
            standardise_scores_training_fold(z, tr)
        assert "non-finite" in ei.value.reason
        assert ei.value.detail["pathway"] == sc.columns[2]


def test_X13_unstandardised_scores_preserved(scores):
    sc, _, tr = scores
    before = sc.copy(deep=True)
    standardise_scores_training_fold(sc, tr)
    pd.testing.assert_frame_equal(sc, before)


def test_X14_coefficient_interpretation_recorded(scores):
    _, d = standardise_scores_training_fold(*[scores[0], scores[2]])
    assert d["coefficient_interpretation"] == (
        "log cause-specific hazard ratio per 1 training-fold SD of the score")
    assert d["amendment"] == "M2-A2"


# ============ P6B-A1 ============

def test_X15_phase6b_amended_original_preserved():
    t = P6B.read_text()
    assert "AMENDED by P6B-A1" in t
    assert "#### 3.1.1 Amendment P6B-A1" in t
    # the superseded wording must still be findable, verbatim
    assert "ridge, λ tuned **inside** resampling" in t


def test_X16_phase6b_amendment_record():
    r = json.loads((CFG / "PHASE6B_AMENDMENT_RECORD.json").read_text())
    a = r["amendments"][0]
    assert a["amendment_id"] == "P6B-A1"
    assert a["superseded_text"].startswith("ridge, lambda tuned")
    assert a["register_A_item"] is False and a["register_A_changed"] is False
    assert a["clinical_block_identical_across_models"] is True
    assert a["clinical_block_separately_penalised_or_retuned"] is False
    assert a["m1_nested_in_m2_when_molecular_coefficients_zero"] is True
    assert "inner cross-validation inside each training fold" in a["molecular_tuning_location"]
    assert a["outcome_data_used"] is False
    assert hashlib.sha256(P6B.read_bytes()).hexdigest() == a["document_sha256_after"]


def test_X17_clinical_block_unpenalised_and_identical():
    t = norm(P6B)
    assert "force-included and **UNPENALISED** in M1, M2 and M3" in t or \
           "force-included and unpenalised** in M1, M2 and M3" in t
    assert "the clinical block is identical in all three" in t.lower()
    assert "never separately penalised or retuned" in t


def test_X18_molecular_tuning_nested():
    t = norm(P6B)
    assert "tuned by inner cross-validation inside each training fold" in t
    v = json.loads((P9 / "validation_specification.json").read_text())
    assert v["tuning_nesting"]["inner_folds"] == 5
    assert "ever selected using" in v["tuning_nesting"]["what_is_never_tuned_outside"]


def test_X19_delta_c_interpretation_binding():
    r = json.loads((CFG / "PHASE6B_AMENDMENT_RECORD.json").read_text())
    s = r["amendments"][0]["required_delta_c_interpretation"]
    assert "CONSERVATIVE" in s
    assert "MUST NOT be described as an unbiased estimate" in s
    t = norm(P6B)
    assert "must not be described as an unbiased estimate" in t.lower()


def test_X20_m1_nested_in_m2():
    t = norm(P6B)
    assert "M1 is exactly M2 with all molecular coefficients set to zero" in t
    assert "same data and the same risk sets" in t


# ============ D8D-4-R1 ============

@pytest.fixture(scope="module")
def d4():
    reg = json.loads((P9 / "decision_register.json").read_text())
    d = next(x for x in reg["decisions"] if x["id"] == "D8D-4")
    return d["sparse_level_handling"]["prespecified_failure_response"]


def test_X21_estimability_criteria_only(d4):
    crit = " ".join(d4["declared_non_estimable_if_and_only_if"])
    assert "converge" in crit
    assert "singular or not positive definite" in crit
    assert "separation sentinel" in crit
    assert "50%" not in crit and "%" not in crit


def test_X22_no_event_cutoff_or_pooling_rule_anywhere():
    reg = json.loads((P9 / "decision_register.json").read_text())
    d4 = next(x for x in reg["decisions"] if x["id"] == "D8D-4")
    f = d4["sparse_level_handling"]["prespecified_failure_response"]
    # the withdrawn rule may be QUOTED in the supersedes record -- that is the
    # audit trail. It must appear nowhere else, and never as an active rule.
    quoting = [k for k, v in f.items()
               if isinstance(v, str) and "fewer than 50%" in v]
    assert quoting == ["supersedes"], quoting
    rest = json.dumps({k: v for k, v in f.items() if k != "supersedes"})
    for banned in ("fewer than 50%", "events per stratum threshold",
                   "pool strata", "pooling threshold", "% of events"):
        assert banned not in rest, banned
    s = next(x for x in reg["decisions"] if x["id"] == "D8D-4")["sparse_level_handling"]
    assert s["RECOMMENDED_RULE"].startswith("All nine")
    assert sum(v for k, v in s["stratum_sizes_complete_case"].items()
               if k != "total") == 452
    assert sum(v for k, v in s["stratum_sizes_s3_cohort_for_contrast"].items()
               if k != "total") == 453


def test_X23_squared_se_ratio_rejected(d4):
    r = d4["REJECTED_the_squared_SE_ratio"]
    assert r["verdict"].startswith("REJECTED")
    assert r["not_forced_into_the_protocol"] is True
    joined = " ".join(r["reasons"])
    assert "ESTIMAND MISMATCH" in joined
    assert "PENALISED SEs ARE NOT NOMINAL" in joined
    assert "NO GUARANTEED DIRECTION" in joined


def test_X24_replacement_diagnostic_well_defined(d4):
    rep = d4["REPLACEMENT_information_diagnostic"]
    assert rep["name"] == "risk-set comparison retention"
    assert "N_strat / N_unstrat" in rep["formula"]
    assert "|R(t_i)| - 1" in rep["formula"]
    assert "MODEL-FREE" in " ".join(rep["why_it_is_defensible"])
    assert "BOUNDED BY CONSTRUCTION" in " ".join(rep["why_it_is_defensible"])
    lim = " ".join(rep["limitations_stated_not_hidden"])
    assert "COMPARISONS, not Fisher information" in lim
    assert rep["no_threshold"].startswith("DESCRIPTIVE ONLY")


def test_X25_diagnostics_necessary_not_sufficient(d4):
    n = d4["diagnostics_are_necessary_not_sufficient"]
    assert "NOT sufficient" in n["statement"]
    joined = " ".join(n["explicitly"])
    assert "positive-definite information matrix does NOT prove" in joined
    assert "numerical failure is NOT the same thing as a wide confidence interval" in joined
    assert "zero-event stratum" in joined and "MUST NOT be treated as a numerical failure" in joined
    assert "one-event stratum is not automatically informative" in joined


def test_X26_risk_set_retention_is_bounded_by_construction():
    """The claim 'bounded by [0,1]' must be true, not just asserted."""
    rng = np.random.default_rng(0)
    for _ in range(50):
        n = int(rng.integers(20, 120))
        strata = rng.integers(0, 5, size=n)
        order = rng.permutation(n)
        n_un = n_st = 0
        for pos, i in enumerate(order[: n // 3]):          # arbitrary "events"
            at_risk = order[pos:]
            n_un += len(at_risk) - 1
            n_st += int((strata[at_risk] == strata[i]).sum()) - 1
        assert 0 <= n_st <= n_un


# ============ chain + frozen state ============

def test_X27_amendment_chain_resolves_all_revisions():
    from src.m2_amendment import resolve_hash, superseded_hashes
    for f, rev1 in (("m2_panel_4set.json", "bab99da9f139db34"),
                    ("m2_panel_5set.json", "6058f733b2acda1c")):
        sup = superseded_hashes(f)
        assert len(sup) == 2, f
        assert any(h.startswith(rev1) for h in sup)
        p = PANELS / f
        for h in sup:
            ok, why = resolve_hash(p, h)
            assert ok and "superseded by" in why, (f, h)
        ok, why = resolve_hash(p, hashlib.sha256(p.read_bytes()).hexdigest())
        assert ok and why == "current"
        assert not resolve_hash(p, "f" * 64)[0]


def test_X28_freeze_record_own_chain_resolves():
    from src.m2_amendment import resolve_hash, superseded_hashes
    sup = superseded_hashes("M2_FREEZE_RECORD.json")
    assert len(sup) >= 2
    p = PANELS / "M2_FREEZE_RECORD.json"
    for h in sup:
        assert resolve_hash(p, h)[0], h


def test_X29_untouched_records_byte_identical():
    expect = {
        "config/analysis/A19_FREEZE_RECORD.json": "a77de7b37fe09eeb",
        "config/analysis/QC_DECISIONS_RECORD.json": "4e026fc547b801ba",
        "config/analysis/QC_SPECIFICATION.json": "d5d6f5b08ac6ab12",
        "src/pathway_scoring.py": "38cac035c0b2d450",
        "data/processed/frozen/frozen_gene_to_probeset_median.csv": "08f8715fb7eca5a0",
        "data/processed/frozen/frozen_gene_to_probe_discovery443.csv": "c386d5f11cc511dc",
        "data/processed/analysis/expression_median_S0.parquet": "a1ca41807c534584",
        "data/processed/analysis/cohort_membership.csv": "bd3671f8df64c594",
        "data/processed/analysis/endpoint.csv": "802a68575a34a3b6",
    }
    for rel, pre in expect.items():
        got = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
        assert got.startswith(pre), rel


def test_X30_phase8_artifacts_unchanged():
    s = json.loads((ROOT / "reports/phase8d/phase8d_summary.json").read_text())
    bad = [k for k, h in s["read_only_artifact_hashes"].items()
           if "phase8" in k and
           hashlib.sha256((ROOT / k).read_bytes()).hexdigest() != h]
    assert not bad, bad


def test_X31_protocol_v3_and_history_preserved():
    s = json.loads((P9 / "phase9_summary.json").read_text())
    assert s["amendments_applied"] == ["M2-A2", "P6B-A1", "D8D-4-R1"]
    log = json.loads((P9 / "protocol_change_log.json").read_text())
    vers = {r["version"]: r["status"] for r in log["revisions"]}
    # v3.0 applied the amendments; later revisions may supersede it, but the
    # three amendments and the full history must survive unchanged.
    assert {"1.0", "2.0", "2.1", "3.0"} <= set(vers)
    assert vers["3.0"] in ("current", "superseded")
    assert [v for v, st in vers.items() if st == "current"] == [s["version"]]


def test_X32_still_zero_models_and_no_outcome():
    s = json.loads((P9 / "phase9_summary.json").read_text())
    assert s["models_fitted"] == 0
    assert s["outcome_columns_opened"] is False
    assert s["event_counts_recomputed"] is False


def test_X33_protocol_states_the_frozen_invariants():
    t = norm(P9DOC)
    for claim in ("4-set PRIMARY, 5-set SENSITIVITY-only",
                  "median aggregation primary; representative-IQR sensitivity",
                  "S3, 453 samples, 0 exclusions"):
        assert claim in t, claim
    assert "Explicitly prohibited while U1 is unresolved" in t
    assert "no suitable cohort exists" in t.lower()
    assert "no analysis here is external validation" in t.lower()


def test_X34_no_stale_pending_claim_in_the_protocol():
    """v3.0 must not still assert that nothing was applied."""
    t = norm(P9DOC)
    assert "All four recommendations are `PENDING_USER_APPROVAL`" not in t
    assert "APPROVED" in t and "APPLIED" in t
    assert "M2-A2" in t and "P6B-A1" in t and "D8D-4-R1" in t
    # historical notes must be scoped to their revision, not read as current
    assert "As of v2.1, no frozen record had been modified" in t
    assert "That changed in v3.0" in t
