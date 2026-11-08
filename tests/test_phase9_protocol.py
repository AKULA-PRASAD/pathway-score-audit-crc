"""Phase 9 planning verification: consistency with every frozen decision.

The protocol is a PLAN. These tests check it does not contradict a frozen
record, does not apply anything, and does not smuggle in an outcome-dependent
or U1-blocked quantity.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "reports" / "phase9"
PANELS = ROOT / "config" / "pathways"
CFG = ROOT / "config" / "analysis"
DOC = ROOT / "docs" / "PHASE9_MODELING_PROTOCOL_DRAFT.md"
SCRIPT = ROOT / "scripts" / "36_phase9_protocol_plan.py"

PRIMARY_SETS = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS",
                "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY"}
PI3K = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
PENDING = "PENDING_USER_APPROVAL"


def doc_text() -> str:
    t = DOC.read_text()
    t = re.sub(r"\n>\s*", " ", t)
    return re.sub(r"\s+", " ", t)


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("phase9", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def summ():
    return json.loads((OUT / "phase9_summary.json").read_text())


@pytest.fixture(scope="module")
def reg():
    return json.loads((OUT / "decision_register.json").read_text())


@pytest.fixture(scope="module")
def panels():
    return (json.loads((PANELS / "m2_panel_4set.json").read_text()),
            json.loads((PANELS / "m2_panel_5set.json").read_text()))


# ==================== nothing executed, nothing applied ====================

def test_P01_no_model_fitted_no_outcome_read(summ):
    assert summ["models_fitted"] == 0
    assert summ["outcome_columns_opened"] is False
    assert summ["event_counts_recomputed"] is False
    assert summ["stage"] == "PLANNING"


def test_P02_outcome_file_guard_raises(mod):
    saved = {n: getattr(pd, n) for n in ("read_csv", "read_parquet", "read_table")}
    try:
        mod._guard_outcome_files()
        with pytest.raises(PermissionError, match="OUTCOME-BLIND"):
            pd.read_csv(ROOT / "data/processed/analysis/endpoint.csv")
    finally:
        for n, f in saved.items():
            setattr(pd, n, f)


def test_P03_only_authorised_amendments_applied(summ, reg):
    assert summ["models_fitted"] == 0
    assert summ["amendments_applied"] == ["M2-A2", "P6B-A1", "D8D-4-R1"]
    # D8D-2 / D8D-3 were approved in v3.1 as REPORTING protocol, no amendment
    assert summ["decisions_still_pending"] == []
    assert sorted(summ["decisions_approved"]) == [
        "D8D-1", "D8D-2", "D8D-3", "D8D-4", "D9-1"]
    assert all(d["status"].startswith("APPROVED") for d in reg["decisions"])


def test_P04_score_standardisation_now_training_fold_only(summ, panels):
    """M2-A2 is applied; the superseded NONE text must remain traceable."""
    for d in panels:
        assert d["scoring"]["score_standardisation"].startswith("TRAINING-FOLD ONLY")
        assert any(a.get("superseded_score_standardisation", "").startswith("NONE")
                   for a in d["amendment_history"])


def test_P05_m2a2_applied_with_full_provenance(reg):
    d1 = next(d for d in reg["decisions"] if d["id"] == "D8D-1")
    amd = d1["applied_amendment"]
    assert d1["proposed_amendment_id"] == "M2-A2"
    assert d1["amendment_applied"] is True
    assert amd["from"].startswith("NONE")
    assert "TRAINING-FOLD ONLY" in amd["to"]
    assert amd["membership_changed"] is False and amd["roles_changed"] is False
    assert amd["outcome_data_used"] is False
    rec = json.loads((PANELS / "M2_FREEZE_RECORD.json").read_text())
    assert [a["amendment_id"] for a in rec["amendment_history"]] == ["M2-A1", "M2-A2"]


# ==================== consistency with frozen decisions ====================

def test_P06_cohort_unchanged(summ):
    assert summ["cohort"]["n"] == 453
    assert summ["cohort"]["exclusions"] == 0
    assert summ["cohort"]["complete_case_n_M1"] == 452
    assert summ["cohort"]["complete_case_n_M2"] == 452
    assert summ["cohort"]["a19_dropped"] == ["GSM972262"]


def test_P07_complete_case_event_count_not_determined(summ):
    cc = summ["cohort"]["complete_case_events"]
    assert "not determined here" in cc
    assert "[135, 136]" in cc
    lo, hi = summ["epv_M2_primary_range"]
    assert 10 < lo <= hi


def test_P08_panel_membership_and_roles(summ, panels):
    p4, p5 = panels
    assert set(p4["gene_sets"]) == PRIMARY_SETS
    assert set(p5["gene_sets"]) == PRIMARY_SETS | {PI3K}
    assert summ["frozen_state"]["panel_primary_role"] == "PRIMARY"
    assert summ["frozen_state"]["panel_sensitivity_role"].startswith("SENSITIVITY")
    assert summ["frozen_state"]["sigma_ddof"] == 1


def test_P09_canonical_scorer(summ):
    from src.pathway_scoring import score_panel
    assert summ["frozen_state"]["canonical_scorer_module"] == "src.pathway_scoring"
    assert score_panel.__module__ == "src.pathway_scoring"


def test_P10_representation_unchanged():
    fr = json.loads((ROOT / "data/processed/frozen/frozen_representation.json").read_text())
    assert fr["approved"]["primary"]["strategy"] == "median"
    assert fr["approved"]["sensitivity"]["strategy"] == "representative_iqr"


def test_P11_missingness_policy_is_complete_case(summ):
    a19 = json.loads((CFG / "A19_FREEZE_RECORD.json").read_text())
    assert "COMPLETE" in json.dumps(a19["decision"]).upper()
    leak = pd.read_csv(OUT / "leakage_control_matrix.csv")
    imp = leak[leak.operation == "imputation"]
    assert len(imp) == 1
    # "N/A" would be read back as NaN; the token must survive the round-trip
    assert imp.permitted_location.iloc[0] == "NOT_APPLICABLE"
    assert leak.permitted_location.notna().all(), "a permitted_location was lost to NaN"


def test_P12_no_real_conflicts(summ):
    assert summ["real_conflicts"] == 0
    assert not any(c["type"] == "REAL_CONFLICT" for c in summ["conflicts"])


# ==================== documented figures are cited, not recomputed ====================

def test_P13_cited_figures_resolve_in_their_named_sources(summ):
    for c in summ["documented_figures_cited"]:
        txt = (ROOT / c["source"]).read_text()
        assert c["matched_text"][:60] in txt, c["label"]
        assert c["value"] in c["matched_text"]


def test_P14_event_counts_match_phase6a_register(summ):
    cited = {c["label"]: c["value"] for c in summ["documented_figures_cited"]}
    assert cited["s3_n"] == "453"
    assert cited["s3_events"] == "136"
    assert cited["competing_events"] == "46"
    assert cited["composite_dfs_events"] == "182"


# ==================== leakage matrix ====================

def test_P15_leakage_matrix_places_estimation_inside():
    leak = pd.read_csv(OUT / "leakage_control_matrix.csv")
    inside = set(leak[leak.permitted_location.str.startswith("INSIDE")].operation)
    for op in ("gene-level mu/sigma for pathway scoring",
               "zero-variance gene exclusion",
               "expression/variance filter cut-point (M3)",
               "univariate screening (M3)",
               "feature selection (M3)",
               "penalty tuning (lambda, alpha)"):
        assert op in inside, op


def test_P16_no_full_cohort_screening_permitted():
    leak = pd.read_csv(OUT / "leakage_control_matrix.csv")
    row = leak[leak.operation.str.contains("univariate screening")].iloc[0]
    assert "NO full-cohort screening" in row.constraint
    t = doc_text()
    assert "No full-cohort univariate screening occurs at any point" in t


def test_P17_score_standardisation_row_is_conditional():
    leak = pd.read_csv(OUT / "leakage_control_matrix.csv")
    row = leak[leak.operation.str.contains("score standardisation")].iloc[0]
    assert "IF D8D-1 APPROVED" in row.permitted_location
    assert PENDING in row.authority


# ==================== validation specification ====================

def test_P18_validation_quantities_are_u1_compatible():
    v = json.loads((OUT / "validation_specification.json").read_text())
    est = " ".join(v["quantities_estimated"]).lower()
    for blocked in ("t-year", "time-dependent auc", "brier", "cif curve",
                    "absolute risk"):
        assert blocked not in est, blocked
    notest = " ".join(v["quantities_NOT_estimated"]).lower()
    for b in ("t-year", "time-dependent auc", "brier", "cloglog"):
        assert b in notest, b


def test_P19_monte_carlo_and_sampling_kept_separate():
    v = json.loads((OUT / "validation_specification.json").read_text())
    mc = v["monte_carlo_vs_sampling_uncertainty"]
    assert "does NOT shrink" in mc["sampling_SE"]
    assert "DOES shrink" in mc["monte_carlo_SE"]
    assert "never summed" in mc["reporting_rule"]
    assert "10%" in mc["criterion"]


def test_P20_five_seeds_declared():
    v = json.loads((OUT / "validation_specification.json").read_text())
    s = v["seed_policy"]
    assert s["n_seeds"] >= 5 and len(set(s["declared_seeds"])) == s["n_seeds"]
    assert 20261102 in s["declared_seeds"]


def test_P21_B_is_an_output_not_an_input():
    v = json.loads((OUT / "validation_specification.json").read_text())
    m = v["monte_carlo_stability"]
    assert m["B_is_an_output"] is True
    assert "NEVER lowered" in m["B_rule"]


def test_P22_event_stratification_flagged_as_outcome_dependent():
    v = json.loads((OUT / "validation_specification.json").read_text())
    note = v["stratification"]["outcome_blindness_note"]
    assert "REQUIRES the outcome" in note
    assert "MODEL-FITTING-TIME" in note


def test_P23_no_external_validation_claim():
    v = json.loads((OUT / "validation_specification.json").read_text())
    assert "NONE" in v["internal_validation_design"]["external_validation"]
    assert "APPROXIMATION" in v["phase5_status"]["statement"]
    assert "NOT a faithful implementation" in v["phase5_status"]["statement"]


# ==================== decision recommendations ====================

def test_P24_five_decisions_each_with_rationale(reg):
    ids = [d["id"] for d in reg["decisions"]]
    assert ids == ["D8D-1", "D8D-2", "D8D-3", "D8D-4", "D9-1"]
    for d in reg["decisions"]:
        assert d["status"] == PENDING or d["status"].startswith("APPROVED")
        assert {"current_frozen_state", "current_protocol_text"} & set(d), d["id"]
        assert "recommendation" in d


def test_P25_d1_states_its_own_limitations(reg):
    d = next(x for x in reg["decisions"] if x["id"] == "D8D-1")
    lim = " ".join(d["limitations_of_the_recommendation"])
    assert "Low dispersion is NOT low biological relevance" in lim
    assert "not like-for-like" in lim
    assert any("training" in c.lower() for c in d["leakage_constraints_if_approved"])


def test_P26_d2_preserves_roles_and_forbids_switching(reg):
    d = next(x for x in reg["decisions"] if x["id"] == "D8D-2")
    rules = " ".join(d["interpretation_rules_fixed_in_advance"])
    assert "reported as primary ALWAYS" in rules
    assert "NOT grounds to switch panels" in rules
    assert "0.631" in rules or "r = 0.631" in rules
    assert "NO role change" in d["recommendation"]


def test_P27_d3_uses_existing_estimand_s7(reg):
    d = next(x for x in reg["decisions"] if x["id"] == "D8D-3")
    assert "S-7" in d["recommendation"]
    assert "Do NOT create a new policy" in d["recommendation"]
    guard = d["comparison_specification"]["selection_guard"]
    assert "REGARDLESS" in guard
    assert "NOT independent confirmation" in \
        " ".join(d["comparison_specification"]["interpretation"])


def test_P28_d4_excludes_confounded_batch_variables(reg):
    d = next(x for x in reg["decisions"] if x["id"] == "D8D-4")
    e = d["eligible_batch_variables"]
    assert e["combat_batch"].startswith("ELIGIBLE")
    assert e["combat_mod"].startswith("NOT ELIGIBLE")
    assert e["rna.extraction.method"].startswith("NOT ELIGIBLE")
    assert "NO batch covariate is added" in d["recommendation"]
    assert len(d["why_never_primary"]) >= 3
    s = d["sparse_level_handling"]
    assert s["RECOMMENDED_RULE"].startswith("All nine")
    assert "WITHDRAWN" in s["v1_proposal_REVERSED"]


# ==================== estimand separation ====================

def test_P29_estimands_kept_distinct():
    e = pd.read_csv(OUT / "estimand_separation.csv")
    assert set(e.id) >= {"PRIMARY", "S-1", "S-2", "S-3", "S-7"}
    s1 = e[e.id == "S-1"].iloc[0]
    assert "NOT comparable" in s1.distinctness
    s2 = e[e.id == "S-2"].iloc[0]
    assert "not more power" in s2.distinctness
    s3 = e[e.id == "S-3"].iloc[0]
    assert "CONDITIONAL on U1" in s3.distinctness
    t = doc_text()
    assert "never pooled, averaged or substituted" in t


# ==================== report discipline ====================

def test_P30_no_prognostic_or_external_claim():
    t = doc_text().lower()
    for phrase in ("externally validated", "is prognostic", "biomarker",
                   "predictive performance", "clinically useful"):
        start = 0
        while (i := t.find(phrase, start)) != -1:
            # a negation may PRECEDE ("not clinically useful") or FOLLOW
            # ("predictive performance ... none has been computed"), so the
            # window spans both sides of the phrase.
            ctx = t[max(0, i - 95):i + len(phrase) + 95]
            assert any(n in ctx for n in ("no ", "not ", "never ", "nothing ",
                                          "none ", "cannot ", "without ")), \
                f"{phrase!r} asserted near: ...{ctx}"
            start = i + len(phrase)


def test_P31_u1_blocked_quantities_absent_from_estimated_set():
    t = doc_text()
    assert "all blocked by U1" in t
    assert "The Cox partial likelihood uses only the **ordering** of times" in t \
        or "uses only the ordering" in t


def test_P32_c_index_limitation_stated():
    t = doc_text()
    assert "censoring-distribution dependent" in t
    assert "not as a transportable performance figure" in t


def test_P33_filename_correction_recorded():
    t = doc_text()
    assert "PHASE6C_HALLMARK_PANEL_REVIEW.md" in t
    assert "does not exist" in t.lower()
    assert "PHASE6C_PATHWAY_PANEL_REVIEW.md" in t
    assert not (ROOT / "docs/PHASE6C_HALLMARK_PANEL_REVIEW.md").exists()


def test_P34_chemotherapy_not_causal():
    t = doc_text()
    assert "never interpreted as a treatment effect" in t


def test_P35_m3_is_exploratory_only():
    cmp = pd.read_csv(OUT / "model_comparison_plan.csv")
    m3 = cmp[cmp.model == "M3"].iloc[0]
    assert m3.status == "EXPLORATORY only"
    assert "no confirmatory claim" in m3.role_constraint
    t = doc_text()
    assert "Exploratory status is not a disclaimer" in t


def test_P36_readiness_checklist_marks_five_decisions():
    c = pd.read_csv(OUT / "phase9_readiness_checklist.csv")
    assert (c.status == "DECISION REQUIRED").sum() == 0
    assert (c.status == "APPROVED - APPLIED").sum() == 3
    assert (c.status == "APPROVED - FORMALISED").sum() == 2
    assert (c.status == "READY").sum() >= 12
    assert "AWAITING USER" in set(c.status)


# ==================== integrity ====================

def test_P37_read_only_artifacts_unchanged(summ):
    bad = [k for k, h in summ["read_only_artifact_hashes"].items()
           if hashlib.sha256((ROOT / k).read_bytes()).hexdigest() != h]
    assert not bad, bad
    assert summ["read_only_artifacts_modified"] is False
    keys = summ["read_only_artifact_hashes"]
    for p in ("phase8a", "phase8b", "phase8c", "phase8d"):
        assert any(p in k for k in keys), p


def test_P38_amendment_chain_resolves():
    from src.m2_amendment import resolve_hash
    for n in ("m2_panel_4set.json", "m2_panel_5set.json"):
        p = PANELS / n
        ok, why = resolve_hash(p, hashlib.sha256(p.read_bytes()).hexdigest())
        assert ok and why == "current", n


def test_P39_phase9_wrote_only_to_its_own_directory():
    for p in OUT.glob("*"):
        assert p.parent == OUT


def test_P40_report_cites_machine_readable_artifacts():
    t = DOC.read_text()
    for p in ("reports/phase9/decision_register.json",
              "reports/phase9/leakage_control_matrix.csv",
              "reports/phase9/validation_specification.json",
              "reports/phase9/estimand_separation.csv",
              "reports/phase6d/definitive_audit.json",
              "reports/phase8c/representation_sensitivity.csv",
              "src/pathway_scoring.py"):
        assert p in t, p


# ==================== v2.0 finalization ====================

def test_P41_current_version_with_full_history_preserved(summ):
    """v2.1 is current; every earlier revision survives as superseded."""
    log = json.loads((OUT / "protocol_change_log.json").read_text())
    vers = {r["version"]: r for r in log["revisions"]}
    assert summ["version"] == max(vers)
    current = [v for v, r in vers.items() if r["status"] == "current"]
    assert current == [summ["version"]], current
    assert all(r["status"] == "superseded" for v, r in vers.items()
               if v != summ["version"])
    # revisions that changed nothing say NONE; a revision that applied an
    # approved amendment must say exactly what moved.
    # exactly one revision applied a specification amendment and must name it;
    # every other revision -- including reporting-only v3.1 -- records NONE.
    naming = {v for v, r in vers.items()
              if "M2-A2" in r["frozen_records_modified"]}
    assert naming == {"3.0"}, naming
    for v, r in vers.items():
        if v not in naming:
            assert r["frozen_records_modified"] == "NONE", (v, r["frozen_records_modified"])
    assert set(vers) >= {"1.0", "2.0", "2.1", "3.0", "3.1"}
    t = doc_text()
    assert "Revision history" in t and "**1.0**" in t and "superseded" in t


def test_P42_d8d4_pooling_rule_reversed(reg):
    d = next(x for x in reg["decisions"] if x["id"] == "D8D-4")
    s = d["sparse_level_handling"]
    assert "WITHDRAWN" in s["v1_proposal_REVERSED"]
    assert "estimates NO parameters for its strata" in s["why_the_threshold_was_wrong"]
    assert s["RECOMMENDED_RULE"].startswith("All nine")
    assert sum(v for k, v in s["stratum_sizes_complete_case"].items()
               if k != "total") == s["stratum_sizes_complete_case"]["total"] == 452
    f = s["prespecified_failure_response"]
    # the v2.0 "50% of events" threshold was WITHDRAWN by the Phase 9 final
    # audit; estimability is now mathematical/numerical, not a percentage.
    assert f["amendment_id"] == "D8D-4-R1"
    assert "WITHDRAWN" in f["supersedes"]
    assert "threshold_fixed_now" not in f
    assert "50%" not in " ".join(f["declared_non_estimable_if_and_only_if"])
    assert f["REJECTED_the_squared_SE_ratio"]["verdict"].startswith("REJECTED")
    assert f["REPLACEMENT_information_diagnostic"]["name"] == \
        "risk-set comparison retention"
    assert "We do NOT fall back to pooling" in f["response_if_declared_not_estimable"]


def test_P43_estimability_depends_on_outcome_is_stated(reg):
    d = next(x for x in reg["decisions"] if x["id"] == "D8D-4")
    s = d["sparse_level_handling"]["estimability_depends_on_outcome_information"]
    assert "outcome data and was NOT read" in s


def test_P44_d9_1_penalisation_scope(reg):
    d = next(x for x in reg["decisions"] if x["id"] == "D9-1")
    assert d["recommendation"] == "B"
    assert d["is_this_a_frozen_item"] is False
    assert "UNPENALISED" in d["recommendation_text"]
    assert "scale-equivariant" in " ".join(d["rationale"])
    assert "change from the 6B draft" in " ".join(d["limitations_of_the_recommendation"])
    assert d["amendment_applied"] is True
    assert d["applied_amendment"]["amendment_id"] == "P6B-A1"


def test_P45_d8d1_scaling_guards_and_clinical_answer(reg):
    d = next(x for x in reg["decisions"] if x["id"] == "D8D-1")
    g = d["scaling_guard_rules"]
    assert "DROPPED for that fold" in g["zero_variance_score_in_a_training_fold"]
    assert "non-estimable replicate" in g["non_finite_score"]
    c = d["clinical_predictor_scaling"]
    assert c["answer"].startswith("NO")
    assert "D9-1" in c["resolution"]
    i = d["coefficient_interpretation_after_standardisation"]
    assert "per 1 training-fold SD" in i["scale"]
    assert "does NOT mean" in i["what_it_does_NOT_mean"]
    assert "NOT a statement about biological relevance" in i["dispersion_vs_relevance"]


def test_P46_age_spline_specification():
    s = json.loads((OUT / "age_spline_specification.json").read_text())
    assert s["form"].startswith("restricted cubic spline, 3 knots")
    k = s["knot_rule"]
    assert "TRAINING SAMPLES OF THAT FOLD" in k["placement"]
    assert k["recomputed_per_fold"] is True
    r = s["reference_age_for_reporting"]
    assert r["value"] == 69 and r["fixed_now"] is True
    assert "non-comparable across folds" in r["why_fixed_rather_than_per_fold"]
    d = s["degenerate_cases"]
    assert "LINEAR AGE" in d["duplicate_knots"]
    assert "linear age" in d["rank_deficient_basis"]
    assert d["fallback_is_recorded_not_silent"] is True
    assert s["all_transformations_inside_the_training_fold"] is True
    assert len(s["covariate_set_unchanged"]) == 5


def test_P47_optimism_and_nonestimable_policy():
    v = json.loads((OUT / "validation_specification.json").read_text())
    o = v["optimism_correction"]
    assert len(o["procedure"]) == 4
    assert "never reported alone" in o["reported_together"]
    n = v["convergence_and_non_estimable_replicates"]
    assert n["threshold_fixed_now"] is True
    assert "5%" in n["pre_registered_threshold"]
    assert "never silently re-drawn" in n["pre_registered_threshold"]
    t = v["tuning_nesting"]
    assert t["inner_folds"] == 5
    assert "ever selected using" in t["what_is_never_tuned_outside"]


def test_P48_explicit_u1_prohibition_block():
    t = doc_text()
    assert "Explicitly prohibited while U1 is unresolved" in t
    for item in ("five-year", "named-horizon calibration", "time-dependent AUC",
                 "cumulative-incidence"):
        assert item in t, item
    assert "prohibited**, not merely omitted" in t


def test_P49_only_m2_specs_changed_everything_else_frozen(summ):
    """M2-A2 touched the M2 specs only; no other frozen record moved."""
    import hashlib
    for rel, pre in (("config/analysis/A19_FREEZE_RECORD.json", "a77de7b37fe09eeb"),
                     ("config/analysis/QC_DECISIONS_RECORD.json", "4e026fc547b801ba"),
                     ("config/analysis/QC_SPECIFICATION.json", "d5d6f5b08ac6ab12"),
                     ("src/pathway_scoring.py", "38cac035c0b2d450")):
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest().startswith(pre), rel
    p4 = json.loads((PANELS / "m2_panel_4set.json").read_text())
    assert set(p4["gene_sets"]) == PRIMARY_SETS
    assert p4["role"] == "PRIMARY" and p4["scoring"]["sigma_ddof"] == 1


def test_P50_change_log_records_the_reversal():
    log = json.loads((OUT / "protocol_change_log.json").read_text())
    v2 = next(r for r in log["revisions"] if r["version"] == "2.0")
    joined = " ".join(v2["changes"])
    assert "REVERSED" in joined
    assert "D9-1 ADDED" in joined
    assert "AGE SPLINE FINALISED" in joined
