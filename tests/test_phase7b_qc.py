"""Phase 7B: QC specification tests.

Validate the QC specification's structure and its outcome-blind restrictions.
These do NOT execute QC and do NOT validate QC results -- no QC has been run.
"""

from __future__ import annotations

import hashlib
import json
import re

import pytest

from src import paths

SPEC = paths.ROOT / "config" / "analysis" / "QC_SPECIFICATION.json"
DOC = paths.ROOT / "docs" / "PHASE7B_QC_PROTOCOL.md"
TIERS = {"HARD_EXCLUSION", "SOFT_INVESTIGATION", "DESCRIPTIVE"}


@pytest.fixture(scope="module")
def spec():
    return json.loads(SPEC.read_text())


@pytest.fixture(scope="module")
def crit(spec):
    return {c["id"]: c for c in spec["criteria"]}


# ---------- structure ----------

def test_Q01_spec_approved_qc_not_executed(spec):
    assert spec["status"].startswith("APPROVED")
    assert "thresholds immutable" in spec["status"]
    assert spec["qc_executed"] is False
    assert spec["approval"]["criteria_count_preserved"] == 20
    assert spec["approval"]["thresholds_preserved"] is True


def test_Q02_every_criterion_fully_specified(spec):
    required = {"id", "domain", "tier", "metric", "rule", "justification",
                "action", "level"}
    for c in spec["criteria"]:
        missing = required - set(c)
        assert not missing, f"{c.get('id')}: missing {missing}"
        assert c["tier"] in TIERS, c["id"]
        assert c["level"] in {"sample", "gene", "cohort"}, c["id"]
        assert len(c["justification"]) > 40, c["id"]


def test_Q03_ids_unique_and_tier_prefix_consistent(crit):
    assert len(crit) == 20
    for cid, c in crit.items():
        pre = cid.split("-")[0]
        assert pre in {"H", "S", "D"}, cid
        expect = {"H": "HARD_EXCLUSION", "S": "SOFT_INVESTIGATION",
                  "D": "DESCRIPTIVE"}[pre]
        assert c["tier"] == expect, f"{cid} prefix/tier mismatch"


def test_Q04_tier_counts(spec):
    import collections
    c = collections.Counter(x["tier"] for x in spec["criteria"])
    assert c["HARD_EXCLUSION"] == 9
    assert c["SOFT_INVESTIGATION"] == 3
    assert c["DESCRIPTIVE"] == 8


# ---------- the core safety property ----------

def test_Q05_only_deterministic_criteria_exclude_automatically(crit):
    """No judgement-based criterion may exclude without human approval."""
    auto = [c for c in crit.values()
            if c["tier"] == "HARD_EXCLUSION" and "exclude" in c["action"].lower()
            and "approval" not in c["action"].lower()]
    for c in auto:
        j = c["justification"].lower()
        assert any(k in j for k in
                   ["not a measurement", "deterministic", "frozen", "no threshold"]), c["id"]


def test_Q06_soft_flags_require_written_approval(crit):
    for c in crit.values():
        if c["tier"] == "SOFT_INVESTIGATION":
            a = c["action"].lower()
            assert "investigate" in a, c["id"]
            assert "approval" in a, c["id"]


def test_Q07_descriptive_never_triggers_exclusion(crit):
    for c in crit.values():
        if c["tier"] == "DESCRIPTIVE":
            assert "report only" in c["action"].lower(), c["id"]
            assert "exclude" not in c["action"].lower().replace("never", ""), c["id"]


# ---------- threshold discipline ----------

def test_Q08_distributions_not_inspected(spec):
    mp = spec["methodological_posture"]
    assert mp["distributions_not_inspected"] is True
    assert "circular" in mp["rationale"].lower()


def test_Q09_tukey_rule_cited_not_invented(crit):
    for cid in ("S-S1", "S-S2", "S-S3"):
        r = crit[cid]["rule"]
        assert "Tukey" in r
        assert "arrayQualityMetrics" in r
        assert "Kauffmann" in r
        assert "DATA-ADAPTIVE" in r


def test_Q10_no_invented_pca_threshold(crit):
    for cid in ("D-S1", "D-S2"):
        assert "NO numeric outlier rule" in crit[cid]["rule"] or \
               "no cut-height" in crit[cid]["rule"], cid


def test_Q11_affymetrix_ratio_threshold_not_applied(crit):
    c = crit["D-S4"]
    assert "NO threshold applied" in c["rule"]
    assert "does NOT transfer" in c["justification"]
    assert "MAS5" in c["justification"]
    assert "ComBat" in c["justification"]
    assert "NOT applied" in c["action"]


def test_Q12_gnuse_excluded_with_reason(spec):
    inap = spec["inapplicable_metrics"]
    assert "GNUSE" in inap
    assert "fRMA" in inap["GNUSE"]
    for m in ("NUSE / RLE", "percent present / MAS5 calls"):
        assert m in inap


# ---------- frozen-artifact protection ----------

def test_Q13_no_gene_criterion_removes_from_frozen(crit):
    for c in crit.values():
        if c["level"] == "gene":
            eff = c.get("effect_on_frozen", "")
            act = c["action"].lower()
            assert ("none" in eff.lower() or "escalate" in act
                    or "approval" in eff.lower()), c["id"]
            # catch an action that silently removes; "do NOT silently drop" is
            # a prohibition, not a violation, so match the verb phrase itself
            assert not re.search(r"\bsilently (drop|remove|delete)", act) or \
                   re.search(r"(not|never|no)\s+silently", act), c["id"]


def test_Q14_zero_variance_defers_to_frozen_m2_rule(crit):
    j = crit["D-G1"]["justification"]
    assert "FROZEN" in j and "PER FOLD" in j
    panel = json.loads((paths.ROOT / "config/pathways/m2_panel_4set.json").read_text())
    assert "EXCLUDED" in panel["scoring"]["zero_variance_genes"]


def test_Q15_affx_not_removed(crit, spec):
    c = crit["D-G3"]
    assert "no removal" in c["rule"].lower()
    assert "NONE" in c["effect_on_frozen"]
    s = spec["affx_sensitivity_analysis"]
    assert s["frozen_specs_unchanged"] is True
    assert s["affected_in_M2_primary_panel"] is False
    assert set(s["affected_genes"]) == {"ACTB", "GAPDH", "STAT1"}
    assert s["pre_specified_now"] is True


# ---------- cohort integrity ----------

def test_Q16_cohort_checks_halt_not_reconcile(crit):
    for cid in ("H-C1", "H-C2", "H-C3"):
        assert "HALT" in crit[cid]["action"], cid


def test_Q17_gsm972262_exclusion_verified(crit):
    c = crit["H-C4"]
    assert "GSM972262" in c["metric"]
    assert c["expected_triggers"] == 1
    assert "chemotherapy.adjuvant" in c["justification"]
    rec = json.loads((paths.ROOT / "config/analysis/A19_FREEZE_RECORD.json").read_text())
    assert rec["excluded_observations"][0]["sample_id"] == "GSM972262"


def test_Q18_normals_excluded_from_tumour_cohort(crit):
    c = crit["H-S4"]
    assert "non-tumoural" in c["metric"]
    assert "19" in c["rule"]


# ---------- prohibitions ----------

def test_Q19_prohibited_actions_complete(spec):
    p = " ".join(spec["prohibited_actions"]).lower()
    for item in ["outcome-guided sample exclusion", "outcome-guided threshold",
                 "frozen phase 2", "frozen m2", "batch correction",
                 "weak outcome association", "expression filter"]:
        assert item in p, item


def test_Q20_outcome_blind(spec):
    assert spec["outcome_data_used"] is False
    blob = json.dumps(spec).lower()
    for pat in ["rfs_event", "rfs_time", "os_event", "survival association",
                "recurrence association"]:
        assert pat not in blob, pat


def test_Q21_doc_and_spec_agree(spec):
    doc = DOC.read_text()
    assert "APPROVED WITH CLARIFICATIONS" in doc
    assert "QC NOT EXECUTED" in doc
    assert "deliberately" in doc.lower()
    for cid in ("H-S1", "S-S1", "D-S4", "H-C4", "D-G3"):
        assert cid in doc, cid


def test_Q22_frozen_artifacts_unchanged():
    """Phase 2/3 artifacts are byte-frozen; M2 specs follow the amendment chain."""
    byte_frozen = {
        "data/processed/frozen/frozen_gene_to_probeset_median.csv": "08f8715fb7eca5a0",
        "data/processed/frozen/frozen_gene_to_probe_discovery443.csv": "c386d5f11cc511dc",
        "data/processed/analysis/expression_median_S0.parquet": "a1ca41807c534584",
    }
    for rel, pre in byte_frozen.items():
        h = hashlib.sha256((paths.ROOT / rel).read_bytes()).hexdigest()
        assert h.startswith(pre), rel

    # The hashes Phase 7B recorded for the M2 specs are a historical
    # attestation; amendment M2-A1 supersedes them with an approved record.
    from src.m2_amendment import resolve_hash

    for rel, pre_amendment in (
            ("config/pathways/m2_panel_4set.json",
             "bab99da9f139db34ccc2d5f5f954cda9e92fc6cf33f41a076024e5146b065dfe"),
            ("config/pathways/m2_panel_5set.json",
             "6058f733b2acda1c7a77746a83bb5d72dd8a9ae64522c81958299af9da79cc27")):
        ok, why = resolve_hash(paths.ROOT / rel, pre_amendment)
        assert ok, why


# ===================== review clarifications (C1-C6) =====================

def test_Q23_thresholds_immutable_after_inspection(spec):
    """C1: thresholds may not be changed once distributions are seen."""
    ti = spec["threshold_immutability"]
    assert ti["applies_to_all_tiers"] is True
    assert "must NOT be changed after observing QC distributions" in ti["rule"]
    for req in ("separate documented decision", "explicit approval"):
        assert any(req in x for x in ti["change_requires"]), req
    assert "no longer pre-specified" in ti["prohibited"]
    assert len(spec["criteria"]) == 20


def test_Q24_hard_exclusions_not_broadened(spec):
    """C2: the hard list is fixed and both conditions are required."""
    h = spec["hard_exclusion_scope"]
    assert h["may_not_be_broadened"] is True
    assert len(h["conditions_both_required"]) == 2
    assert any("deterministic" in c for c in h["conditions_both_required"])
    assert any("demonstrably invalid" in c for c in h["conditions_both_required"])
    fixed = set(h["fixed_criteria"])
    actual = {c["id"] for c in spec["criteria"] if c["tier"] == "HARD_EXCLUSION"}
    assert fixed == actual, fixed ^ actual
    assert len(fixed) == 9
    assert "demonstrably" in spec["tiers"]["HARD_EXCLUSION"].lower()


def test_Q25_no_automatic_soft_exclusion(spec, crit):
    """C3: THE core rule -- a fence crossing can never exclude on its own."""
    sf = spec["soft_flag_policy"]
    assert sf["investigation_only"] is True
    assert sf["automatic_exclusion"] == "PROHIBITED"
    p = sf["explicit_prohibition"]
    for term in ("Tukey", "PCA", "clustering", "median", "IQR",
                 "between-array distance"):
        assert term in p, term
    assert "never a reason to exclude" in p
    # every non-hard criterion must carry the prohibition flag
    for c in spec["criteria"]:
        if c["tier"] in ("SOFT_INVESTIGATION", "DESCRIPTIVE"):
            assert c.get("automatic_exclusion") == "PROHIBITED", c["id"]


def test_Q26_soft_actions_state_prohibition_and_approval(crit):
    for cid in ("S-S1", "S-S2", "S-S3"):
        a = crit[cid]["action"]
        assert "INVESTIGATE ONLY" in a, cid
        assert "PROHIBITED" in a, cid
        assert "approval" in a.lower(), cid


def test_Q27_exclusion_pathway_terminates_in_approval(spec):
    """C4: the only route from flag to exclusion passes through approval."""
    path = spec["soft_flag_policy"]["exclusion_pathway"]
    for step in ("investigation record", "request explicit approval",
                 "only if approved"):
        assert step in path, step


def test_Q28_investigation_record_has_seven_fields(spec):
    tpl = spec["investigation_record_template"]
    fields = [f["field"] for f in tpl["required_fields"]]
    assert len(fields) == 7, fields
    for f in ("sample_id", "flagged_metric_and_value", "metadata_and_batch",
              "independent_technical_evidence", "biological_plausibility",
              "effect_on_cohort", "proposed_action_and_justification"):
        assert f in fields, f
    assert tpl["one_record_per_sample"] is True
    assert "STOP" in tpl["terminal_step"]


def test_Q29_cohort_effect_never_guides_exclusion(spec):
    """C4 field 6: reported for transparency, never used to decide."""
    tpl = spec["investigation_record_template"]
    assert "PROHIBITED" in tpl["outcome_guided_selection"]
    assert "never inform the decision" in tpl["outcome_guided_selection"]
    eff = next(f for f in tpl["required_fields"] if f["field"] == "effect_on_cohort")
    assert "NEVER used to choose" in eff["description"]


def test_Q30_independent_evidence_required_or_stated_absent(spec):
    tpl = spec["investigation_record_template"]
    f = next(x for x in tpl["required_fields"]
             if x["field"] == "independent_technical_evidence")
    assert "independent of the flag itself" in f["description"]
    assert "state explicitly if NONE" in f["description"]


def test_Q31_affx_approved_conditional_and_non_rescuing(spec):
    """C5: approved, conditional, and cannot rescue an unfavourable result."""
    s = spec["affx_sensitivity_analysis"]
    assert s["status"] == "APPROVED - CONDITIONAL"
    assert "ONLY if M3 materially implicates" in s["trigger"]
    assert s["primary_representation"].startswith("PRESERVED")
    assert "BOTH analyses reported" in s["reporting_requirement"]
    a = s["anti_rescue_clause"]
    assert "must NOT be used to rescue" in a
    assert "primary remains primary" in a
    assert "prohibited" in a.lower()


def test_Q32_preprocessing_decision_preserved(spec):
    """C6: deposited matrix as-is; no normalization/ComBat/fRMA/CEL."""
    pp = spec["preprocessing_decision"]
    assert pp["preserved"] is True
    assert "AS-IS" in pp["decision"]
    joined = " ".join(pp["prohibited"]).lower()
    for item in ("normalization", "combat", "frma", "cel"):
        assert item in joined, item


def test_Q33_doc_records_all_clarifications():
    doc = DOC.read_text()
    for tag in ("C1", "C2", "C3", "C4", "C5", "C6"):
        assert f"### {tag} —" in doc, tag
    assert "never a reason to exclude" in doc
    assert "Anti-rescue clause" in doc
    assert "QC NOT EXECUTED" in doc


def test_Q34_spec_remains_the_preregistration(spec):
    """The specification is a PRE-REGISTRATION artifact.

    It is not mutated by execution: `qc_executed` stays False in the spec, and
    the fact of execution is recorded in the Phase 7C outputs instead. This
    keeps the approved, threshold-immutable document immutable.
    """
    assert spec["qc_executed"] is False
    assert "qc_execution" not in spec, (
        "execution facts must live in the decisions record and reports/, "
        "not in the pre-registration document")
    assert spec["preregistration_invariant"].startswith("This document is the pre-registration")
    assert spec["status"].startswith("APPROVED")


def test_Q35_no_sample_was_excluded_by_qc():
    """THE durable invariant: QC execution may never exclude a sample."""
    import json as _j
    summ = paths.ROOT / "reports" / "phase7c" / "qc_summary.json"
    if not summ.exists():
        pytest.skip("Phase 7C not executed")
    d = _j.loads(summ.read_text())
    assert d["samples_excluded"] == 0
    assert d["exclusions_performed"] is False
    assert d["outcome_columns_opened"] is False
    assert d["event_counts_computed"] is False
    assert d["normalization_applied"] == "none"
    assert d["batch_correction_applied"] == "none"
    assert d["frozen_artifacts_modified"] is False


def test_Q36_every_soft_flagged_sample_has_a_full_record():
    """Each flagged sample must carry all seven required fields."""
    import json as _j
    recs_p = paths.ROOT / "reports" / "phase7c" / "investigation_records.json"
    summ_p = paths.ROOT / "reports" / "phase7c" / "qc_summary.json"
    if not recs_p.exists():
        pytest.skip("Phase 7C not executed")
    recs = _j.loads(recs_p.read_text())
    summ = _j.loads(summ_p.read_text())
    assert len(recs) == summ["samples_soft_flagged"]
    assert {r["sample_id"] for r in recs} == set(summ["sample_ids_requiring_investigation"])
    required = {"sample_id", "flagged_metric_and_value", "metadata_and_batch",
                "independent_technical_evidence", "biological_plausibility",
                "effect_on_cohort", "proposed_action_and_justification"}
    for r in recs:
        assert required <= set(r), r["sample_id"]
        assert r["proposed_action_and_justification"]["proposed_action"] == \
            "NO EXCLUSION - investigation only"
        assert "NOT COMPUTED" in r["effect_on_cohort"]["event_counts"]
        assert "HYPOTHESIS, not evidence" in r["biological_plausibility"]["note"]
