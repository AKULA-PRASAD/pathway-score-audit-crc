"""D8D-2 and D8D-3: approved reporting protocols, and what they must NOT change.

These decisions formalise REPORTING. They carry no specification amendment, so
the tests check both that the protocol now states them and that nothing frozen
moved as a side effect.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

P9 = ROOT / "reports" / "phase9"
P8C = ROOT / "reports" / "phase8c"
PANELS = ROOT / "config" / "pathways"
DOC = ROOT / "docs" / "PHASE9_MODELING_PROTOCOL_DRAFT.md"

PRIMARY_SETS = ["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS",
                "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY"]
PI3K = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"


def norm() -> str:
    return re.sub(r"\s+", " ", re.sub(r"\n>\s*", " ", DOC.read_text()))


@pytest.fixture(scope="module")
def reg():
    return json.loads((P9 / "decision_register.json").read_text())


@pytest.fixture(scope="module")
def d2(reg):
    return next(d for d in reg["decisions"] if d["id"] == "D8D-2")


@pytest.fixture(scope="module")
def d3(reg):
    return next(d for d in reg["decisions"] if d["id"] == "D8D-3")


@pytest.fixture(scope="module")
def panels():
    return (json.loads((PANELS / "m2_panel_4set.json").read_text()),
            json.loads((PANELS / "m2_panel_5set.json").read_text()))


# ==================== D8D-2 ====================

def test_Y01_d8d2_approved_with_provenance(d2):
    a = d2["approval"]
    assert a["approved"] is True and a["date"] == "2026-11-03"
    assert a["approved_by"] == "project owner"
    assert a["original_proposal_preserved_above"] is True
    assert d2["status"].startswith("APPROVED")
    # the original proposal must still be there, unrewritten
    assert "interpretation_rules_fixed_in_advance" in d2
    assert "prespecified_divergence_metric" in d2


def test_Y02_four_set_is_primary(d2, panels):
    f = d2["formalised_decision"]["primary_panel"]
    assert f["panel"] == "4-set" and f["role"].startswith("PRIMARY")
    assert sorted(f["members"]) == sorted(PRIMARY_SETS)
    p4, _ = panels
    assert p4["role"] == "PRIMARY"
    assert sorted(p4["gene_sets"]) == sorted(PRIMARY_SETS)


def test_Y03_five_set_is_sensitivity_only(d2, panels):
    s = d2["formalised_decision"]["sensitivity_panel"]
    assert s["panel"] == "5-set" and s["role"] == "SENSITIVITY-ONLY"
    assert s["adds"] == PI3K
    _, p5 = panels
    assert p5["role"].startswith("SENSITIVITY")
    assert sorted(p5["gene_sets"]) == sorted(PRIMARY_SETS + [PI3K])


def test_Y04_pi3k_not_promoted(d2, panels):
    p4, p5 = panels
    assert PI3K not in p4["gene_sets"], "PI3K silently promoted to primary"
    assert PI3K in p5["gene_sets"]
    proh = " ".join(d2["formalised_decision"]["prohibitions"])
    assert "MUST NOT be promoted to primary" in proh
    assert "model performance or statistical significance" in proh
    assert "REGARDLESS of which panel performs better" in proh
    rec = json.loads((PANELS / "M2_FREEZE_RECORD.json").read_text())
    clause = rec["approval"]["immutability_clause"]
    assert "must NOT be promoted to the primary panel" in clause


def test_Y05_collinearity_reporting_mandatory_and_correct(d2):
    m = d2["formalised_decision"]["mandatory_reporting_with_the_5_set_analysis"]
    c = pd.read_csv(P8C / "score_correlation_pearson_5set_sensitivity.csv", index_col=0)
    actual = float(c.loc["HALLMARK_P53_PATHWAY", PI3K])
    assert abs(m["p53_pi3k_pearson_r"] - actual) < 5e-5
    s = json.loads((P8C / "phase8c_summary.json").read_text())["relationships"]
    d = m["collinearity_diagnostics"]
    assert abs(d["condition_number_4set"] - s["4set_primary"]["condition_number"]) < 1e-3
    assert abs(d["condition_number_5set"] - s["5set_sensitivity"]["condition_number"]) < 1e-3
    assert abs(d["min_eigenvalue_4set"] - s["4set_primary"]["min_eigenvalue"]) < 1e-4
    assert abs(d["min_eigenvalue_5set"] - s["5set_sensitivity"]["min_eigenvalue"]) < 1e-4
    assert d["n_pairs_abs_r_gt_0.5_4set"] == s["4set_primary"]["n_gt_0.5"] == 0
    assert d["n_pairs_abs_r_gt_0.5_5set"] == s["5set_sensitivity"]["n_gt_0.5"] == 1
    assert abs(d["pi3k_r2_on_the_primary_four"]
               - s["pi3k_increment"]["r2_on_primary_four"]) < 1e-4
    assert "P53 and PI3K TOGETHER" in d["also_report"]


def test_Y06_pi3k_correlation_is_not_biological_dismissal(d2):
    w = d2["formalised_decision"]["why_pi3k_is_sensitivity_only"]
    assert "CONDITIONING" in w["explicitly_NOT_the_reason"]
    assert "NOT evidence that PI3K" in w["explicitly_NOT_the_reason"]
    assert "biologically irrelevant" in w["explicitly_NOT_the_reason"]
    assert "No outcome association has been tested" in w["explicitly_NOT_the_reason"]


def test_Y07_no_biological_claim_from_dispersion_or_correlation(d2):
    proh = " ".join(d2["formalised_decision"]["prohibitions"])
    assert "dispersion or correlation MUST NOT be interpreted as proof of " \
           "biological activity or clinical relevance" in proh


def test_Y08_hierarchy_distinguished_from_performance(d2):
    h = d2["formalised_decision"]["analysis_hierarchy_vs_observed_performance"]
    assert "PRE-SPECIFIED" in h["rule"] and "OBSERVED MODEL PERFORMANCE" in h["rule"]
    ops = " ".join(h["operationally"])
    assert "4-set, always" in ops
    assert "no evidence of a difference" in ops


# ==================== D8D-3 ====================

def test_Y09_d8d3_approved_with_provenance(d3):
    a = d3["approval"]
    assert a["approved"] is True and a["date"] == "2026-11-03"
    assert a["original_proposal_preserved_above"] is True
    assert "comparison_specification" in d3        # original proposal intact
    assert "observed_evidence" in d3
    assert d3["status"].startswith("APPROVED")


def test_Y10_median_primary_repiqr_sensitivity(d3):
    f = d3["formalised_decision"]
    assert f["primary_representation"]["role"] == "PRIMARY"
    assert "median" in f["primary_representation"]["name"]
    assert f["sensitivity_representation"]["role"] == "SENSITIVITY-ONLY"
    assert f["sensitivity_representation"]["registered_estimand"] == "S-7"
    fr = json.loads((ROOT / "data/processed/frozen/frozen_representation.json").read_text())
    assert fr["approved"]["primary"]["strategy"] == "median"
    assert fr["approved"]["sensitivity"]["strategy"] == "representative_iqr"
    for key in ("primary_representation", "sensitivity_representation"):
        assert (ROOT / f[key]["artifact"]).exists()


def test_Y11_no_replacement_on_performance(d3):
    proh = " ".join(d3["formalised_decision"]["prohibitions"])
    assert "MUST NOT be replaced on the basis of observed model performance" in proh
    assert "universally superior" in proh
    assert "chosen after seeing outcomes" in proh


def test_Y12_wnt_weakest_agreement_documented_and_accurate(d3):
    r = d3["formalised_decision"]["representation_agreement_phase8c"]
    src = pd.read_csv(P8C / "representation_sensitivity.csv").set_index("pathway")
    for path in PRIMARY_SETS + [PI3K]:
        assert abs(r[path] - src.loc[path, "pearson_median_vs_repiqr"]) < 5e-5, path
    assert r["weakest_agreement"] == "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
    assert abs(r["weakest_value"] - 0.8577) < 5e-5
    # WNT must genuinely be the minimum, and the only one below 0.9
    vals = {k: r[k] for k in PRIMARY_SETS + [PI3K]}
    assert min(vals, key=vals.get) == "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
    assert [k for k, v in vals.items() if v < 0.9] == \
        ["HALLMARK_WNT_BETA_CATENIN_SIGNALING"]
    assert "WEAKEST representation agreement" in r["acknowledgement"]
    assert "not_evidence_of_superiority" in r


def test_Y13_robustness_is_not_performance(d3):
    rb = d3["formalised_decision"]["robustness_is_not_performance"]
    assert "never substituted for one another" in rb["statement"]
    assert "NONE has been computed" in rb["performance"]
    assert "NOT independent replication" in rb["and_neither_is_replication"]


def test_Y14_reporting_requirements_listed(d3):
    rep = d3["formalised_decision"]["mandatory_reporting_when_the_analyses_are_authorised"]
    joined = " ".join(rep)
    for k in ("DIRECTIONALLY", "QUANTITATIVELY", "pathway scores", "coefficients",
              "uncertainty", "discrimination"):
        assert k in joined, k


# ==================== nothing frozen moved ====================

def test_Y15_prior_amendments_intact():
    rec = json.loads((PANELS / "M2_FREEZE_RECORD.json").read_text())
    assert [a["amendment_id"] for a in rec["amendment_history"]] == ["M2-A1", "M2-A2"]
    assert rec["specification_revision"] == 3
    p4 = json.loads((PANELS / "m2_panel_4set.json").read_text())
    assert p4["scoring"]["sigma_ddof"] == 1                       # M2-A1
    assert p4["scoring"]["score_standardisation"].startswith("TRAINING-FOLD ONLY")  # M2-A2
    p6b = json.loads((ROOT / "config/analysis/PHASE6B_AMENDMENT_RECORD.json").read_text())
    assert p6b["amendments"][0]["amendment_id"] == "P6B-A1"       # P6B-A1
    d4 = next(d for d in json.loads((P9 / "decision_register.json").read_text())
              ["decisions"] if d["id"] == "D8D-4")
    f = d4["sparse_level_handling"]["prespecified_failure_response"]
    assert f["amendment_id"] == "D8D-4-R1"                        # D8D-4-R1
    assert f["REJECTED_the_squared_SE_ratio"]["verdict"].startswith("REJECTED")


def test_Y16_m2_panel_hashes_unchanged_by_this_approval():
    """D8D-2/D8D-3 are reporting protocol: the specs must not have moved."""
    expect = {"m2_panel_4set.json": "43ce15b2f0e06f75",
              "m2_panel_5set.json": "d792368a85b1f3c0",
              "M2_FREEZE_RECORD.json": "64a86fa8ed0d6a8a"}
    for f, pre in expect.items():
        got = hashlib.sha256((PANELS / f).read_bytes()).hexdigest()
        assert got.startswith(pre), f


def test_Y17_protected_files_unchanged():
    expect = {
        "config/analysis/A19_FREEZE_RECORD.json": "a77de7b37fe09eeb",
        "config/analysis/QC_DECISIONS_RECORD.json": "4e026fc547b801ba",
        "config/analysis/QC_SPECIFICATION.json": "d5d6f5b08ac6ab12",
        "src/pathway_scoring.py": "38cac035c0b2d450",
        "src/score_standardisation.py": "41b5c091125e1b11",
        "data/processed/frozen/frozen_representation.json": "d5ec63736d874258",
        "data/processed/frozen/frozen_gene_to_probeset_median.csv": "08f8715fb7eca5a0",
        "data/processed/frozen/frozen_gene_to_probe_discovery443.csv": "c386d5f11cc511dc",
        "data/processed/analysis/expression_median_S0.parquet": "a1ca41807c534584",
        "data/processed/analysis/expression_repiqr_disc443_S0.parquet": "e636e259b4feabfa",
        "data/processed/analysis/cohort_membership.csv": "bd3671f8df64c594",
        "data/processed/analysis/endpoint.csv": "802a68575a34a3b6",
    }
    for rel, pre in expect.items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest().startswith(pre), rel


def test_Y18_amendment_chain_still_resolves():
    from src.m2_amendment import resolve_hash, superseded_hashes
    for f in ("m2_panel_4set.json", "m2_panel_5set.json", "M2_FREEZE_RECORD.json"):
        q = PANELS / f
        for h in superseded_hashes(f) + [hashlib.sha256(q.read_bytes()).hexdigest()]:
            ok, _ = resolve_hash(q, h)
            assert ok, (f, h)
        assert not resolve_hash(q, "f" * 64)[0]


def test_Y19_historical_attestations_still_resolve():
    from src.m2_amendment import resolve_hash
    bad = []
    for s in ("reports/phase8d/phase8d_summary.json",
              "reports/phase9/phase9_summary.json"):
        for k, h in json.loads((ROOT / s).read_text()).get(
                "read_only_artifact_hashes", {}).items():
            ok, why = resolve_hash(ROOT / k, h)
            if not ok:
                bad.append((s, k, why))
    assert not bad, bad


# ==================== protocol state ====================

def test_Y20_protocol_v3_1_with_full_history():
    s = json.loads((P9 / "phase9_summary.json").read_text())
    assert s["version"] == "3.1"
    assert s["decisions_still_pending"] == []
    assert sorted(s["decisions_approved"]) == ["D8D-1", "D8D-2", "D8D-3", "D8D-4", "D9-1"]
    assert s["reporting_protocols_formalised"] == ["D8D-2", "D8D-3"]
    log = json.loads((P9 / "protocol_change_log.json").read_text())
    vers = {r["version"]: r["status"] for r in log["revisions"]}
    assert vers == {"1.0": "superseded", "2.0": "superseded", "2.1": "superseded",
                    "3.0": "superseded", "3.1": "current"}
    v31 = next(r for r in log["revisions"] if r["version"] == "3.1")
    assert v31["frozen_records_modified"] == "NONE"


def test_Y21_document_states_both_decisions():
    t = norm()
    assert "**Version:** 3.1" in t
    assert "D8D-2 and D8D-3 APPROVED" in t
    assert "four-set panel = **pre-specified primary**" in t.replace("Four", "four")
    assert "0.8577" in t and "0.9975" in t and "0.6314" in t
    low = t.lower()
    assert "not evidence that pi3k" in low        # needle must also be lowercased
    assert "biologically irrelevant" in low
    assert "not independent replication" in low
    assert "No Phase 9 decision remains pending" in t


def test_Y22_still_zero_models_and_no_outcome():
    s = json.loads((P9 / "phase9_summary.json").read_text())
    assert s["models_fitted"] == 0
    assert s["outcome_columns_opened"] is False
    assert s["event_counts_recomputed"] is False
    assert s["cohort"]["n"] == 453 and s["cohort"]["exclusions"] == 0
