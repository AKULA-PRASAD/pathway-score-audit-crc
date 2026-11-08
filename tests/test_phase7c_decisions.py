"""Phase 7C: approved-decision tests.

Verify the three approved decisions are recorded correctly and that the
no-exclusion / 453-retention invariants hold. These do not re-run QC.
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd
import pytest

from src import paths

SPEC = paths.ROOT / "config" / "analysis" / "QC_SPECIFICATION.json"
DEC = paths.ROOT / "config" / "analysis" / "QC_DECISIONS_RECORD.json"
OUT = paths.ROOT / "reports" / "phase7c"
DOC = paths.ROOT / "docs" / "PHASE7C_QC_REPORT.md"


@pytest.fixture(scope="module")
def spec():
    return json.loads(SPEC.read_text())


@pytest.fixture(scope="module")
def dec():
    return json.loads(DEC.read_text())


@pytest.fixture(scope="module")
def summ():
    return json.loads((OUT / "qc_summary.json").read_text())


# ---------- D1: the amended H-C1 rule ----------

def test_R01_spec_at_revision_3_with_history(spec):
    assert spec["specification_revision"] == 3
    hist = spec["revision_history"]
    assert [h["revision"] for h in hist] == [1, 2, 3]
    assert "H-C1 AMENDED" in hist[2]["status"]


def test_R02_hc1_amended_to_subset_rule(spec):
    c = next(x for x in spec["criteria"] if x["id"] == "H-C1")
    assert c["amended_in_revision"] == 3
    assert "SUBSET" in c["metric"]
    assert "zero bundle IDs missing" in c["metric"]
    assert "HALTS QC" in c["rule"]
    assert "are NOT a failure" in c["rule"]
    # still a hard criterion, still one of nine
    assert c["tier"] == "HARD_EXCLUSION"
    assert len([x for x in spec["criteria"]
                if x["tier"] == "HARD_EXCLUSION"]) == 9


def test_R03_original_hc1_preserved_verbatim(spec):
    c = next(x for x in spec["criteria"] if x["id"] == "H-C1")
    a = spec["revision_history"][2]["amendment"]
    assert "equality" in c["superseded_metric"].lower()
    assert c["superseded_metric"] == a["original_metric"]
    assert c["superseded_rule"] == a["original_rule"]


def test_R04_amendment_records_defect_rationale_approval(spec):
    a = spec["revision_history"][2]["amendment"]
    for k in ("defect", "rationale", "approved_by", "approved_utc",
              "corrected_metric", "corrected_rule"):
        assert a.get(k), k
    assert "never achievable" in a["defect"]
    assert "469 - 461 = 8" in a["rationale"]


def test_R05_original_execution_result_preserved(spec):
    r = spec["revision_history"][2]["amendment"]["original_execution_result"]
    assert r["as_written_set_equality"] is False
    assert r["matrix_samples"] == 585
    assert r["frozen_bundle_samples"] == 461
    assert "FAIL" in r["status_as_written"]
    assert r["preserved_for_auditability"] is True


def test_R06_amendment_neither_broadens_nor_rethresholds(spec, dec):
    a = spec["revision_history"][2]["amendment"]
    assert a["broadens_hard_list"] is False
    assert a["changes_any_threshold"] is False
    d1 = dec["decisions"]["D1_hc1_amendment"]
    assert d1["broadens_hard_criteria_list"] is False
    assert d1["changes_any_threshold"] is False
    assert len(spec["criteria"]) == 20


def test_R07_amended_rule_holds_on_the_real_artifacts():
    """Independently re-evaluate the corrected H-C1 rule."""
    import pyarrow.parquet as pq
    cols = set(pq.ParquetFile(paths.PROBE_MATRIX).schema.names) - {"probe_id", "__index_level_0__"}
    bundle = set(pd.read_csv(OUT.parent.parent / "data/processed/analysis/cohort_membership.csv").sample_id)
    assert bundle <= cols
    assert not (bundle - cols), "a frozen-bundle sample is missing from the matrix"


# ---------- D2: retention ----------

def test_R08_all_453_retained(dec, summ):
    d2 = dec["decisions"]["D2_retain_all_453"]
    assert d2["samples_retained"] == 453
    assert d2["samples_excluded"] == 0
    assert d2["no_soft_flagged_sample_authorized_for_exclusion"] is True
    assert summ["samples_excluded"] == 0
    assert summ["exclusions_performed"] is False
    mem = pd.read_csv(paths.ROOT / "data/processed/analysis/cohort_membership.csv")
    assert int(mem.in_S3.astype(bool).sum()) == 453


def test_R09_all_29_records_preserved(dec, summ):
    recs = json.loads((OUT / "investigation_records.json").read_text())
    assert len(recs) == 29 == dec["decisions"]["D2_retain_all_453"]["investigation_records_preserved"]
    assert {r["sample_id"] for r in recs} == set(dec["decisions"]["D2_retain_all_453"]["flagged_sample_ids"])
    assert {r["sample_id"] for r in recs} == set(summ["sample_ids_requiring_investigation"])
    for r in recs:
        assert r["proposed_action_and_justification"]["proposed_action"] == \
            "NO EXCLUSION - investigation only"


def test_R10_explicit_no_evidence_statement(dec):
    s = dec["decisions"]["D2_retain_all_453"]["explicit_statement"]
    assert "No independent evidence demonstrated technical invalidity" in s
    assert "did NOT discriminate" in s
    inv = pd.read_csv(OUT / "investigation_summary.csv")
    assert int((inv.independent_technical_evidence == "NONE FOUND").sum()) == 26
    assert len(inv) == 29


# ---------- D3: observations, not justifications ----------

def test_R11_observations_are_not_exclusion_justifications(dec):
    obs = dec["decisions"]["D3_observations_not_justifications"]["observations"]
    assert len(obs) == 2
    ids = {o["id"] for o in obs}
    assert ids == {"OBS-1", "OBS-2"}
    for o in obs:
        assert o["is_exclusion_justification"] is False
        assert o["status"] == "OBSERVATION ONLY"
        assert "none" in o["action_taken"].lower()
    b = next(o for o in obs if o["id"] == "OBS-1")
    assert "batch05" in b["finding"] and "0.226" in b["finding"]
    g = next(o for o in obs if o["id"] == "OBS-2")
    assert "GSM972311" in g["finding"]
    assert "inside their fences" in g["interpretation"] or \
           "WITHOUT being" in g["interpretation"]


# ---------- preservation ----------

def test_R12_thresholds_unchanged_by_the_amendment(spec):
    crit = {c["id"]: c for c in spec["criteria"]}
    assert "Tukey" in crit["S-S1"]["rule"]
    assert "Q3 + 1.5*IQR" in crit["S-S1"]["rule"]
    for cid in ("S-S1", "S-S2", "S-S3"):
        assert crit[cid].get("amended_in_revision") is None, cid
    assert spec["threshold_immutability"]["applies_to_all_tiers"] is True
    assert dec_unchanged_hashes()


def dec_unchanged_hashes():
    """The Phase 7C record is a HISTORICAL attestation and is never rewritten.

    A recorded hash must therefore resolve either to the current file or to a
    version an approved amendment documents as superseded (see
    src/m2_amendment). Undocumented drift still fails.
    """
    from src.m2_amendment import resolve_hash

    d = json.loads(DEC.read_text())["frozen_artifacts_unmodified"]
    for rel, h in d.items():
        cands = [paths.ROOT / "data/processed/frozen" / rel,
                 paths.ROOT / "config/pathways" / rel,
                 paths.ROOT / "config/analysis" / rel,
                 paths.ROOT / "data/processed/analysis" / rel]
        p = next(c for c in cands if c.exists())
        ok, why = resolve_hash(p, h)
        assert ok, why
    return True


def test_R13_no_prohibited_action_performed(dec, summ):
    na = [x.lower() for x in dec["actions_not_performed"]]
    for item in ("normalization", "batch correction", "modeling",
                 "feature selection", "sample exclusion", "outcome analysis"):
        assert item in na, item
    assert summ["normalization_applied"] == "none"
    assert summ["batch_correction_applied"] == "none"
    assert summ["outcome_columns_opened"] is False


def test_R14_report_records_all_three_decisions():
    doc = DOC.read_text()
    for tag in ("### D1 —", "### D2 —", "### D3 —"):
        assert tag in doc, tag
    assert "OBSERVATION ONLY" in doc
    assert "no independent evidence demonstrated technical invalidity" in doc.lower()
    assert "All 453 S3 samples retained" in doc
