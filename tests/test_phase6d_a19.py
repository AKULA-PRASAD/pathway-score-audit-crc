"""Phase 6D: A-19 missing-data strategy — policy tests.

These are POLICY assertions on the protocol text and the preliminary audit
artifacts. They establish that A-19 is recorded as PROPOSED (not frozen), that
the never-impute list is complete, and that the frozen M2 scoring rule is
unaltered by the missing-data proposal.

They do NOT establish that imputation is implemented correctly — no imputation
code exists. Behavioural fold-isolation and outcome-leakage tests are
specified in the protocol (§6.4.10) and will be added with the implementation.
"""

from __future__ import annotations

import json
import re

import pandas as pd
import pytest

from src import paths

PROTOCOL = paths.ROOT / "docs" / "PHASE6B_MODELING_PROTOCOL_DRAFT.md"
AUDIT_CSV = paths.ROOT / "reports" / "phase6d" / "preliminary_missingness_S3.csv"
AUDIT_JSON = paths.ROOT / "reports" / "phase6d" / "preliminary_missingness_summary.json"
PANEL = paths.ROOT / "config" / "pathways" / "m2_panel_4set.json"

NEVER_IMPUTE = ["rfs_event_source", "rfs_time_source", "os_",
                "Sample identifiers", "Cohort membership", "Provenance fields"]


@pytest.fixture(scope="module")
def protocol():
    return PROTOCOL.read_text()


@pytest.fixture(scope="module")
def audit():
    return json.loads(AUDIT_JSON.read_text())


# ---------- status ----------

RECORD = paths.ROOT / "config" / "analysis" / "A19_FREEZE_RECORD.json"
EXCL = paths.ROOT / "reports" / "phase6d" / "complete_case_exclusions.json"
PLANNED = ["M1 primary (clinical only)", "M2 primary (frozen 4-set pathway panel)",
           "M2 sensitivity (frozen 5-set pathway panel)",
           "M3 exploratory (transcriptome-wide)"]


@pytest.fixture(scope="module")
def record():
    return json.loads(RECORD.read_text())


def test_A01_a19_is_frozen(protocol, record):
    assert "## 6.4 A-19 — Missing-data strategy ✅ **FROZEN**" in protocol
    assert "✅ **YES — FROZEN" in protocol
    assert record["status"] == "FROZEN"
    assert record["record_type"] == "A19_MISSING_DATA_STRATEGY_FREEZE"


def test_A01b_strategy_is_complete_case(record):
    assert record["decision"]["strategy"] == "COMPLETE-CASE ANALYSIS"
    assert record["decision"]["multiple_imputation_required"] is False


def test_A01c_all_four_models_covered(record):
    assert set(record["decision"]["applies_to"]) == set(PLANNED)
    per = {m["model"]: m for m in record["evidence"]["per_model"]}
    assert len(per) == 4
    for m in per.values():
        assert m["complete_cases"] == 452
        assert m["lost"] == 1


def test_A01d_excluded_observation_identified(record):
    ex = record["excluded_observations"]
    assert len(ex) == 1
    e = ex[0]
    assert e["sample_id"] == "GSM972262"
    assert e["title"] == "CIT306"
    assert e["missing_variables"] == ["chemotherapy.adjuvant"]
    assert e["raw_values"]["chemotherapy.adjuvant"] == "N/A"
    assert e["dataset"] == "discovery"
    assert e["combat_batch"] == "batch05"


def test_A01e_exclusion_also_in_protocol_and_artifact(protocol):
    assert "GSM972262" in protocol and "CIT306" in protocol
    d = json.loads(EXCL.read_text())
    assert d["n_excluded"] == 1 and d["n_retained"] == 452
    assert d["excluded"][0]["sample_id"] == "GSM972262"


def test_A01f_prohibited_imputation_recorded(record):
    pr = [x.lower() for x in record["prohibited_imputation"]]
    for item in ["expression values", "pathway members", "endpoint fields",
                 "sample identifiers", "cohort membership", "provenance fields"]:
        assert item in pr, item


def test_A01g_m2_coverage_rules_unchanged(record):
    cv = record["m2_coverage_rules"]
    assert cv["unchanged"] is True
    assert ">=80%" in cv["rule"] and ">=20 genes" in cv["rule"]
    assert "never substituted" in cv["missing_members"]
    # The A-19 record attests to the bytes Phase 6D read. It is never
    # rewritten by a later amendment, so the recorded hash must resolve to the
    # current file OR to an approved superseded version; undocumented drift
    # still fails. Membership and the coverage rules are asserted separately
    # above and are what A-19 actually depends on.
    from src.m2_amendment import resolve_hash

    for f in ("m2_panel_4set.json", "m2_panel_5set.json"):
        path = paths.ROOT / "config" / "pathways" / f
        ok, why = resolve_hash(path, record["upstream_frozen_hashes"][f])
        assert ok, why


def test_A01h_contingency_requires_approval(record):
    c = record["contingency"]
    assert c["mi_auto_activation"] is False
    assert "tp53.mutation" in c["optional_covariates_excluded_from_M1"]
    joined = " ".join(c["if_proposed_in_future_require"]).lower()
    for req in ["rationale", "complete-case retention", "multiple imputation",
                "batch", "explicit approval"]:
        assert req in joined, req


def test_A01i_pooling_deferred_not_invented(record):
    d = record["deferred"]["pooling_rubin_vs_averaged_predictions"]
    assert d["status"] == "DEFERRED — NOT APPLICABLE"
    assert "none is invented" in d["reason"]
    assert "approved" in d["reopens_if"]


def test_A01j_phase6d_distinctions_preserved(record):
    f = record["preserved_phase6d_findings"]
    assert "ND is used 0 times" in f["nd_na_distinction"]
    assert "must NOT be treated as ordinary missing" in f["structural_na"]
    assert "batch" in f["batch_association"]
    assert "not proven" in f["mechanism_qualified"]
    assert "GSM972262" in f["single_missing_chemotherapy"]


def test_A02_register_a_marks_a19_frozen(protocol):
    m = re.search(r"\|\s*A-19\s*\|[^\r\n]*", protocol)
    assert m, "Register A must contain an A-19 row"
    row = m.group(0)
    assert "FROZEN" in row, row
    assert "complete-case" in row.lower(), row
    assert "NOT FROZEN" not in row, row


def test_A03_all_register_a_items_settled(protocol):
    assert "Every Register A item is now settled" in protocol


# ---------- never-impute list ----------

def test_A04_never_impute_list_complete(protocol):
    s = protocol[protocol.index("### 6.4.3"):protocol.index("### 6.4.4")]
    for item in NEVER_IMPUTE:
        assert item in s, f"never-impute list missing: {item}"


def test_A05_endpoint_ambiguities_not_imputable(protocol):
    s = protocol[protocol.index("### 6.4.3"):protocol.index("### 6.4.4")]
    assert "must NOT be used to resolve the documented endpoint ambiguities" in s
    assert "U1" in s and "U2" in s


def test_A06_expression_imputation_not_authorised(protocol):
    s = protocol[protocol.index("### 6.4.3"):protocol.index("### 6.4.4")]
    assert "separate scientific justification and approval" in s
    assert "not authorised by this proposal" in s


# ---------- M2 scoring rule preserved ----------

def test_A07_m2_scoring_unchanged_by_a19():
    """The frozen M2 scoring spec must be untouched by the missing-data proposal."""
    d = json.loads(PANEL.read_text())
    sc = d["scoring"]
    assert d["status"] == "FROZEN"
    assert "EXCLUDED" in sc["zero_variance_genes"]
    assert "never imputed" in sc["missing_gene_members"]
    assert ">=80%" in sc["minimum_set_coverage"]
    assert ">=20 genes" in sc["minimum_set_coverage"]


def test_A08_protocol_restates_frozen_coverage_rule(protocol):
    s = protocol[protocol.index("### 6.4.3"):protocol.index("### 6.4.4")]
    assert "80%" in s and "20 genes" in s
    assert "no silent substitution" in s.lower()


# ---------- fold isolation policy ----------

def test_A09_fold_isolation_rules_present(protocol):
    s = protocol[protocol.index("### 6.4.4"):protocol.index("### 6.4.5")]
    assert "fitted on training data only" in s
    assert "Never" in s and "complete dataset before cross-validation" in s
    for leak in ["imputation parameters", "feature filtering", "standardisation",
                 "pathway-scoring parameters", "feature selection", "penalty tuning"]:
        assert leak in s, leak


def test_A10_heldout_not_conditioned_on_outcome(protocol):
    s = protocol[protocol.index("### 6.4.6"):protocol.index("### 6.4.7")]
    assert "must NOT condition on held-out event status or survival time" in s
    assert "unavailable at prediction time" in s


def test_A11_pooling_not_assumed_interchangeable(protocol):
    s = protocol[protocol.index("### 6.4.7"):protocol.index("### 6.4.8")]
    assert "NOT interchangeable" in s
    assert "non-linear" in s
    assert "Not decided here" in s


def test_A12_sensitivity_not_proof(protocol):
    s = protocol[protocol.index("### 6.4.8"):protocol.index("### 6.4.9")]
    assert "never as proof that one method is correct" in s


def test_A13_methods_deliberately_unspecified(protocol):
    s = protocol[protocol.index("### 6.4.5"):protocol.index("### 6.4.6")]
    assert "DELIBERATELY UNSPECIFIED" in protocol or "left **unspecified**" in s
    assert "NOT the primary strategy" in s


# ---------- preliminary audit artifacts ----------

def test_A14_preliminary_audit_preserved(audit):
    """The preliminary audit is retained, still labelled preliminary."""
    assert "PRELIMINARY" in audit["status"]
    assert "pre-QC" in audit["status"]


def test_A14b_definitive_audit_present():
    """Phase 6D definitive audit artifacts exist and are self-consistent."""
    d = json.loads((paths.ROOT / "reports" / "phase6d"
                    / "definitive_audit.json").read_text())
    assert d["outcome_columns_opened"] is False
    assert d["cohort"]["n"] == 453
    assert d["encoding"]["ND_occurrences_total"] == 0
    assert d["encoding"]["ND_documented_but_unused"] is True
    assert d["expression"]["nan_values"] == 0
    assert d["expression"]["zero_variance_genes"] == 0


def test_A14c_planned_models_complete_case():
    """All four planned models lose exactly one patient."""
    df = pd.read_csv(paths.ROOT / "reports" / "phase6d" / "model_complete_case.csv")
    planned = df[df.model.str.contains("M1 primary|M2 primary|M2 sensitivity|M3 explor")]
    assert len(planned) == 4
    assert (planned.complete_cases == 452).all()
    assert (planned.lost == 1).all()


def test_A14d_molecular_covariates_not_in_M1():
    """High-missingness molecular covariates must not enter M1."""
    v = pd.read_csv(paths.ROOT / "reports" / "phase6d" / "variable_missingness.csv")
    m1 = set(v.loc[v.group == "M1_primary", "variable"])
    mol = set(v.loc[v.group == "molecular_optional", "variable"])
    assert not (m1 & mol)
    assert "tp53.mutation" in mol
    assert v.loc[v.variable == "tp53.mutation", "pct"].iloc[0] > 30


def test_A15_audit_numbers_consistent(audit):
    df = pd.read_csv(AUDIT_CSV)
    assert audit["n"] == 453
    assert len(df) == 15
    assert set(df.loc[df.complete, "variable"]) == set(audit["complete_variables"])
    assert audit["endpoint_missing_in_S3"] == 0


def test_A16_m1_nearly_complete(audit):
    """The finding that drives the proposal: M1 barely needs imputation."""
    assert audit["M1_complete_case_retained"] == 452
    assert audit["M1_complete_case_pct"] >= 99.0
    assert audit["all_variable_complete_case_pct"] < 50.0


def test_A17_definitive_audit_recorded(protocol):
    s = protocol[protocol.index("### 6.4.2"):protocol.index("### 6.4.3")]
    assert "COMPLETED (Phase 6D)" in s
    assert "by model design" in s


def test_A18_planned_tests_documented(protocol):
    s = protocol[protocol.index("### 6.4.10"):]
    for tid in ("T-A19-1", "T-A19-2", "T-A19-3", "T-A19-4"):
        assert tid in s, tid
    assert "Not implementable yet" in s
