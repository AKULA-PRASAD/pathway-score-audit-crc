"""Phase 9 Final Protocol Audit verification (v2.0 -> v2.1)."""

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

P9 = ROOT / "reports" / "phase9"
PANELS = ROOT / "config" / "pathways"
AUDIT_DOC = ROOT / "docs" / "PHASE9_FINAL_AUDIT.md"
PROTO_DOC = ROOT / "docs" / "PHASE9_MODELING_PROTOCOL_DRAFT.md"
SCRIPT = ROOT / "scripts" / "37_phase9_final_audit.py"
PENDING = "PENDING_USER_APPROVAL"


def norm(p: Path) -> str:
    t = re.sub(r"\n>\s*", " ", p.read_text())
    return re.sub(r"\s+", " ", t)


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("audit", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def aud():
    return json.loads((P9 / "final_audit_summary.json").read_text())


@pytest.fixture(scope="module")
def summ():
    return json.loads((P9 / "phase9_summary.json").read_text())


@pytest.fixture(scope="module")
def reg():
    return json.loads((P9 / "decision_register.json").read_text())


# ==================== boundary ====================

def test_A01_outcome_guard_raises(mod):
    saved = {n: getattr(pd, n) for n in ("read_csv", "read_parquet", "read_table")}
    try:
        mod._guard_outcome_files()
        with pytest.raises(PermissionError, match="OUTCOME-BLIND"):
            pd.read_csv(ROOT / "data/processed/analysis/endpoint.csv")
    finally:
        for n, f in saved.items():
            setattr(pd, n, f)


def test_A02_nothing_executed(aud, summ):
    """The audit itself fits nothing. Amendments were applied later, by approval."""
    assert aud["models_fitted"] == 0
    assert aud["outcome_columns_opened"] is False
    assert aud["event_counts_read_or_inferred"] is False
    assert summ["models_fitted"] == 0
    fr = aud["frozen_records_modified"]
    assert "M2-A2" in fr and "NONE" in fr, fr


def test_A03_453_retained(aud):
    assert aud["samples_retained"] == 453 and aud["exclusions"] == 0


# ==================== finding 1 ====================

def test_A04_batch_reconciliation_independently_rederived():
    """Re-derive the complete-case set and the counts here, not from the JSON."""
    mem = pd.read_csv(ROOT / "data/processed/analysis/cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    cov = pd.read_csv(ROOT / "data/processed/analysis/covariates.csv",
                      keep_default_na=False, dtype=str).set_index("sample_id")
    pb = pd.read_csv(ROOT / "data/processed/analysis/partition_and_batch.csv"
                     ).set_index("sample_id")
    a19 = json.loads((ROOT / "config/analysis/A19_FREEZE_RECORD.json").read_text())
    dropped = [e["sample_id"] for e in a19["excluded_observations"]]
    M1 = ["age.at.diagnosis_(year)", "sex", "tnm.stage", "tumor.location",
          "chemotherapy.adjuvant"]
    cc = [s for s in s3 if all(cov.loc[s, v] != "" for v in M1)]
    assert len(s3) == 453 and len(cc) == 452
    assert sorted(set(s3) - set(cc)) == sorted(dropped) == ["GSM972262"]
    b453 = pb.loc[s3, "combat_batch"].value_counts().sort_index()
    b452 = pb.loc[cc, "combat_batch"].value_counts().sort_index()
    assert b453.sum() == 453 and b452.sum() == 452
    assert [k for k in b453.index if b453[k] != b452.get(k, 0)] == ["batch05"]
    assert pb.loc["GSM972262", "combat_batch"] == "batch05"


def test_A05_batch_join_has_no_gaps(aud):
    j = aud["finding_1_batch_reconciliation"]["join_integrity"]
    assert j["s3_ids_absent_from_batch_metadata"] == []
    assert j["duplicate_sample_ids_in_batch_metadata"] == 0
    assert j["s3_ids_with_blank_or_NA_batch"] == []
    assert j["s3_ids_with_multiple_batch_rows"] == []
    assert j["distinct_batch_labels_over_s3"] == 9
    assert aud["finding_1_batch_reconciliation"][
        "bundle_batch_equals_sdrf_batch_for_all_453"] is True


def test_A06_premise_recorded_as_incorrect_not_silently_fixed(aud):
    f = aud["finding_1_batch_reconciliation"]
    assert f["premise_status"].startswith("INCORRECT")
    assert "NOT A DATA ERROR" in f["verdict"]
    assert f["correction_issued"].startswith("CLARITY, not fact")
    assert f["event_status_of_dropped_sample"] == "NOT INSPECTED -- outcome data"


def test_A07_both_count_tables_published(summ):
    b = summ["batch_counts_both_sets"]
    assert sum(b["s3_cohort_453"].values()) == 453
    assert sum(b["complete_case_452"].values()) == 452
    assert b["s3_cohort_453"]["batch05"] - b["complete_case_452"]["batch05"] == 1
    t = norm(PROTO_DOC)
    assert "S3 cohort (453)" in t and "complete-case (452)" in t


# ==================== finding 2 ====================

def test_A08_fifty_percent_rule_withdrawn(aud, reg):
    f = aud["finding_2_fifty_percent_rule"]
    assert f["verdict"].startswith("WITHDRAW")
    assert "VACUOUS" in f["defect_1_logical"]
    assert "no statistical basis" in f["defect_2_statistical"]
    blob = json.dumps(reg)
    assert "fewer than 50% of events" not in blob or "WITHDRAWN" in blob
    d4 = next(d for d in reg["decisions"] if d["id"] == "D8D-4")
    r = d4["sparse_level_handling"]["prespecified_failure_response"]
    assert r["amendment_id"] == "D8D-4-R1"
    assert "WITHDRAWN" in r["supersedes"]
    assert r["status"].startswith("APPROVED")


def test_A09_replacement_is_estimability_not_precision(reg):
    d4 = next(d for d in reg["decisions"] if d["id"] == "D8D-4")
    r = d4["sparse_level_handling"]["prespecified_failure_response"]
    crit = " ".join(r["declared_non_estimable_if_and_only_if"])
    assert "converge" in crit and "singular or not positive definite" in crit
    assert "separation sentinel" in crit
    assert "50%" not in crit
    assert len(r["reported_always_never_gated_on"]) >= 6
    rej = r["REJECTED_the_squared_SE_ratio"]
    assert rej["verdict"].startswith("REJECTED")
    assert rej["not_forced_into_the_protocol"] is True
    rep = r["REPLACEMENT_information_diagnostic"]
    assert rep["name"] == "risk-set comparison retention"
    assert rep["no_threshold"].startswith("DESCRIPTIVE ONLY")
    nec = r["diagnostics_are_necessary_not_sufficient"]
    assert "NOT sufficient" in nec["statement"]


def test_A10_three_concepts_distinguished(aud):
    c = aud["finding_2_fifty_percent_rule"]["three_concepts_conflated"]
    assert set(c) == {"mathematical_estimability", "numerical_convergence",
                      "useful_information"}
    assert "identified iff" in c["mathematical_estimability"]
    assert "CONTINUOUS" in c["useful_information"]


def test_A11_partial_likelihood_contributions_explained(aud):
    h = aud["finding_2_fifty_percent_rule"]["how_strata_contribute"]
    assert "factor of 1" in h["zero_event_stratum"]
    assert "log(1) = 0" in h["event_with_a_degenerate_risk_set"]
    assert "non-degenerate" in h["consequence"]


def test_A12_no_pooling_and_no_post_outcome_redefinition(reg, aud):
    d4 = next(d for d in reg["decisions"] if d["id"] == "D8D-4")
    assert d4["sparse_level_handling"]["RECOMMENDED_RULE"].startswith("All nine")
    unchanged = " ".join(aud["finding_2_fifty_percent_rule"]["proposed_replacement"]
                         ["unchanged_from_v2_0"])
    assert "NO pooling" in unchanged
    assert "no post-outcome pooling" in unchanged


# ==================== finding 3 ====================

def test_A13_d9_1_audited_sound_with_amendment(aud, reg):
    f = aud["finding_3_penalisation_scope"]
    assert f["verdict"].startswith("SOUND")
    assert f["phase6b_permission"]["does_register_A_cover_it"] is False
    assert f["phase6b_permission"]["is_an_amendment_required"] is True
    assert f["phase6b_permission"]["amendment_id"] == "P6B-A1"
    d9 = next(d for d in reg["decisions"] if d["id"] == "D9-1")
    assert d9["applied_amendment"]["amendment_id"] == "P6B-A1"
    assert d9["applied_amendment"]["applied"] is True
    assert d9["applied_amendment"]["original_wording_preserved"] is True
    assert d9["status"].startswith("APPROVED")


def test_A14_nesting_and_comparability_findings(aud):
    f = aud["finding_3_penalisation_scope"]
    assert "properly nested" in f["m1_m2_interpretability"]
    assert "recovers M1 exactly" in f["nesting_coherence"]
    assert f["clinical_coefficient_comparability"].startswith("IMPROVED")


def test_A15_conservative_estimand_statement_is_mandatory(aud):
    e = aud["finding_3_penalisation_scope"]["effect_on_the_molecular_estimand"]
    assert e["changed"] is True
    s = e["statement_required_whenever_deltaC_is_reported"]
    assert "CONSERVATIVE" in s and "biased against" in s
    assert "conservative" in norm(PROTO_DOC).lower()


def test_A16_phase6b_original_wording_preserved_after_amendment():
    """P6B-A1 is applied, but the superseded wording must survive verbatim."""
    t = (ROOT / "docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md").read_text()
    assert "ridge, λ tuned **inside** resampling" in t, "original wording lost"
    assert "AMENDED by P6B-A1" in t
    assert "3.1.1 Amendment P6B-A1" in t


# ==================== finding 4 ====================

def test_A17_all_consistency_checks_pass(aud):
    f = aud["finding_4_consistency"]
    assert f["all_consistent"] is True
    assert f["n_inconsistent"] == 0
    assert f["n_checks"] >= 30
    areas = {c["area"] for c in f["checks"]}
    for a in ("Register A", "Model roles", "Leakage", "Tuning", "M2-A1 vs M2-A2",
              "Age spline", "Batch", "Bootstrap/MC", "Replicates", "U1/U2",
              "Estimands", "External validation", "Counts"):
        assert a in areas, a


# ==================== version / change log ====================

def test_A18_audit_issued_2_1_and_history_is_complete(aud, summ):
    """The audit issued v2.1; v3.0 applied the amendments. History is complete."""
    assert aud["audited_version"] == "2.0" and aud["issued_version"] == "2.1"
    log = json.loads((P9 / "protocol_change_log.json").read_text())
    vers = {r["version"]: r for r in log["revisions"]}
    assert {"1.0", "2.0", "2.1", "3.0"} <= set(vers)
    assert [v for v, r in vers.items() if r["status"] == "current"] == [summ["version"]]
    for v in ("1.0", "2.0", "2.1"):
        assert vers[v]["status"] == "superseded", v
        assert vers[v]["frozen_records_modified"] == "NONE", v
    # the revision that APPLIED M2-A2 must name it; later reporting-only
    # revisions correctly record NONE.
    assert "M2-A2" in vers["3.0"]["frozen_records_modified"]
    t = norm(PROTO_DOC)
    assert "Revision history" in t and "**1.0**" in t and "superseded" in t


def test_A19_three_amendments_now_applied(aud):
    ids = {a["id"]: a for a in aud["new_amendments_proposed"]}
    assert set(ids) == {"M2-A2", "P6B-A1", "D8D-4-R1"}
    assert all(a["applied"] is True and a["status"].startswith("APPROVED")
               for a in ids.values())
    ps = aud["post_amendment_state"]
    assert ps["amendment_chain"] == ["M2-A1", "M2-A2"]
    assert ps["m2_specification_revision"] == 3
    # D8D-2 / D8D-3 were approved in v3.1 as reporting protocol
    assert ps["decisions_still_pending"] == []
    assert sorted(ps["decisions_approved"]) == [
        "D8D-1", "D8D-2", "D8D-3", "D8D-4", "D9-1"]
    assert ps["reporting_protocols_formalised"] == ["D8D-2", "D8D-3"]


# ==================== frozen state ====================

def test_A20_untouched_records_unchanged_amended_ones_traceable(aud):
    bad = [k for k, h in aud["read_only_artifact_hashes"].items()
           if hashlib.sha256((ROOT / k).read_bytes()).hexdigest() != h]
    assert not bad, bad
    p4 = json.loads((PANELS / "m2_panel_4set.json").read_text())
    rec = json.loads((PANELS / "M2_FREEZE_RECORD.json").read_text())
    assert p4["scoring"]["score_standardisation"].startswith("TRAINING-FOLD ONLY")
    assert any(a.get("superseded_score_standardisation", "").startswith("NONE")
               for a in p4["amendment_history"])
    assert [a["amendment_id"] for a in rec["amendment_history"]] == ["M2-A1", "M2-A2"]
    assert aud["canonical_scorer_module"] == "src.pathway_scoring"


def test_A21_amendment_chain_intact():
    from src.m2_amendment import resolve_hash
    for n in ("m2_panel_4set.json", "m2_panel_5set.json"):
        p = PANELS / n
        ok, why = resolve_hash(p, hashlib.sha256(p.read_bytes()).hexdigest())
        assert ok and why == "current", n


def test_A22_audit_doc_states_the_four_findings():
    t = norm(AUDIT_DOC)
    assert "Finding 1" in t and "Finding 2" in t and "Finding 3" in t and "Finding 4" in t
    assert "no error of fact" in t.lower()
    assert "WITHDRAWN" in t
    assert "33 / 33 checks pass" in t or "33/33" in t
    assert "Frozen records modified: NONE" in t
