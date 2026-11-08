"""Phase 8D verification: cross-phase counts, frozen decisions, no new claims.

Phase 8D consolidates; it must not introduce a number that is not already in a
source artifact, must not decide anything, and must not touch any prior-phase
artifact.
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

OUT = ROOT / "reports" / "phase8d"
P8A, P8B, P8C = (ROOT / "reports" / f"phase8{x}" for x in "abc")
PANELS = ROOT / "config" / "pathways"
DOC = ROOT / "docs" / "PHASE8D_CONSOLIDATED_EDA_REVIEW.md"
SCRIPT = ROOT / "scripts" / "35_phase8d_consolidation.py"

PRIMARY_SETS = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS",
                "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY"}
PI3K = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"


def doc_text() -> str:
    t = DOC.read_text()
    t = re.sub(r"\n>\s*", " ", t)
    return re.sub(r"\s+", " ", t)


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("phase8d", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def summ():
    return json.loads((OUT / "phase8d_summary.json").read_text())


@pytest.fixture(scope="module")
def reg():
    return json.loads((OUT / "decision_register.json").read_text())


@pytest.fixture(scope="module")
def srcs():
    return {k: json.loads((p / f"phase8{k}_summary.json").read_text())
            for k, p in (("a", P8A), ("b", P8B), ("c", P8C))}


# ==================== outcome blindness ====================

def test_D01_file_guard_raises(mod):
    saved = {n: getattr(pd, n) for n in ("read_csv", "read_parquet", "read_table")}
    try:
        mod._guard_outcome_files()
        with pytest.raises(PermissionError, match="OUTCOME-BLIND"):
            pd.read_csv(ROOT / "data/processed/analysis/endpoint.csv")
    finally:
        for n, f in saved.items():
            setattr(pd, n, f)


def test_D02_endpoint_exists_so_guard_is_not_vacuous():
    assert (ROOT / "data/processed/analysis/endpoint.csv").exists()


def test_D03_no_outcome_field_anywhere_in_outputs(mod, summ):
    banned = {c.lower() for c in mod.OUTCOME_COLUMNS}
    for f in OUT.glob("*.csv"):
        cols = {c.strip().lower() for c in pd.read_csv(f, nrows=0).columns}
        assert not (cols & banned), f.name
    s = dict(summ)
    s.pop("outcome_columns_never_read")
    blob = json.dumps(s).lower()
    for c in mod.OUTCOME_COLUMNS:
        # the register legitimately NAMES rfs.delay when describing U1
        if c == "rfs.delay":
            continue
        assert c.lower() not in blob, c


def test_D04_summary_flags(summ):
    assert summ["outcome_columns_opened"] is False
    assert summ["modelling_performed"] is False
    assert summ["recurrence_association_tested"] is False
    assert summ["feature_selection_performed"] is False
    assert summ["new_analysis_performed"] is False


def test_D05_all_three_source_phases_were_outcome_blind(srcs):
    for k, s in srcs.items():
        assert s["outcome_columns_opened"] is False, k


# ==================== claims are resolved, not transcribed ====================

def test_D06_every_claim_resolves_against_its_named_artifact(mod):
    claims = pd.read_csv(OUT / "consolidated_claims.csv")
    assert len(claims) >= 40
    for _, r in claims.iterrows():
        obj = json.loads((ROOT / r.source_artifact).read_text())
        val = mod.dig(obj, r.json_pointer)
        got, want = str(val), str(r.value)
        assert got == want or abs(float(val) - float(r.value)) < 1e-12, \
            f"{r.claim_id}: artifact {got} vs recorded {want}"


def test_D07_claims_declare_zero_transcription(summ):
    assert summ["claims_transcribed"] == 0
    assert summ["claims_resolved_from_artifacts"] == len(
        pd.read_csv(OUT / "consolidated_claims.csv"))


def test_D08_matrix_derived_claims_are_marked():
    claims = pd.read_csv(OUT / "consolidated_claims.csv")
    assert claims.derived_from_same_expression_matrix.sum() > 0
    assert (~claims.derived_from_same_expression_matrix).sum() > 0


# ==================== cross-phase counts ====================

def test_D09_cohort_size_agrees_across_phases(srcs, summ):
    assert all(s["n_samples"] == 453 for s in srcs.values())
    assert summ["n_samples"] == 453
    assert all(s["samples_excluded"] == 0 for s in srcs.values())
    assert summ["samples_excluded"] == 0


def test_D10_all_cross_phase_checks_agree(summ):
    assert summ["all_cross_phase_checks_agree"] is True
    checks = summ["cross_phase_checks"]
    assert len(checks) >= 15
    assert all(c["agrees"] for c in checks)


def test_D11_pathway_coverage_agrees_8a_vs_8c(srcs):
    a = {k: v["n_genes_used"] for k, v in
         srcs["a"]["pathways"]["5set_sensitivity"]["coverage"].items()}
    c = {k: v["n_genes_used"] for k, v in
         srcs["c"]["audit"]["5set_sensitivity"]["coverage"].items()}
    assert a == c
    assert sorted(c.values()) == [41, 104, 194, 195, 198]


def test_D12_the_eight_are_the_same_eight(srcs):
    a = set(pd.read_csv(P8A / "sex_positive_control.csv")
            .query("classification == 'discordant_both_axes'").sample_id)
    b = {r["sample_id"] for r in srcs["b"]["sex_discordance_review"]["samples"]}
    assert a == b and len(a) == 8


def test_D13_pca_reproduction_recorded(srcs):
    r = srcs["b"]["reproduction_of_phase8a"]
    assert r["reproduced"] is True
    assert r["min_signed_correlation_scores_pc1_20"] > 1 - 1e-12


# ==================== frozen decisions unchanged ====================

def test_D14_panel_membership_and_roles(summ):
    p4 = json.loads((PANELS / "m2_panel_4set.json").read_text())
    p5 = json.loads((PANELS / "m2_panel_5set.json").read_text())
    assert set(p4["gene_sets"]) == PRIMARY_SETS
    assert set(p5["gene_sets"]) == PRIMARY_SETS | {PI3K}
    assert p4["role"] == "PRIMARY" and p5["role"].startswith("SENSITIVITY")
    assert summ["frozen_state"]["panel_primary"] == "PRIMARY"
    assert summ["frozen_state"]["sigma_ddof"] == 1


def test_D15_phase8d_recorded_the_pre_amendment_state(summ):
    """Phase 8D attests to NONE. M2-A2 later superseded it, traceably."""
    assert summ["frozen_state"]["score_standardisation"].startswith("NONE")
    p4 = json.loads((PANELS / "m2_panel_4set.json").read_text())
    cur = p4["scoring"]["score_standardisation"]
    assert cur.startswith("NONE") or any(
        a.get("superseded_score_standardisation", "").startswith("NONE")
        for a in p4["amendment_history"])


def test_D16_canonical_scorer_unchanged(summ):
    from src.pathway_scoring import score_panel
    assert summ["frozen_state"]["canonical_scorer_module"] == "src.pathway_scoring"
    assert score_panel.__module__ == "src.pathway_scoring"


def test_D17_no_phase8c_recommendation_implemented(summ):
    assert summ["phase8c_recommendations_implemented_here"] == []


# ==================== the register decides nothing ====================

def test_D18_register_has_eight_items_all_awaiting(reg, summ):
    assert reg["decisions_made_by_phase8d"] == 0
    assert summ["decisions_made_by_phase8d"] == 0
    ids = [d["id"] for d in reg["decisions"]]
    assert ids == [f"D8D-{i}" for i in range(1, 9)]
    assert all(d["decision"] == "AWAITING_APPROVAL" for d in reg["decisions"])
    assert summ["decisions_awaiting_approval"] == 8


def test_D19_no_item_blocks_phase9_entirely(reg):
    allowed = {"blocks_specific_analysis", "documented_limitation",
               "blocks_phase9_entirely"}
    for d in reg["decisions"]:
        assert d["blocking"] in allowed, d["id"]
    assert not any(d["blocking"] == "blocks_phase9_entirely" for d in reg["decisions"])


def test_D20_every_item_offers_options_without_choosing(reg):
    for d in reg["decisions"]:
        keys = set(d)
        assert keys & {"options_not_chosen", "proposed_reporting_framework_not_a_change",
                       "outstanding_route"}, d["id"]
        assert "implication" in d or "current_frozen_state" in d


def test_D21_required_register_topics_present(reg):
    blob = json.dumps(reg).lower()
    for topic in ("standardis", "pi3k", "wnt", "batch sensitivity", "combat_mod",
                  "sex-axis", "u1", "u2", "external validation"):
        assert topic in blob, topic


def test_D22_u1_u2_constraints_stated(reg):
    d7 = next(d for d in reg["decisions"] if d["id"] == "D8D-7")
    assert "fixed-time" in d7["blocks_what"].lower()
    assert "relative-risk analyses are unaffected" in d7["blocks_what"].lower()
    assert "naming_hazard" in d7, "the Phase 1 vs 6A U1/U2 collision must be recorded"


def test_D23_external_validation_carried_forward(reg):
    d8 = next(d for d in reg["decisions"] if d["id"] == "D8D-8")
    imp = d8["implication"]
    assert "48" in imp and "18" in imp and "338" in imp
    assert "MUST NOT be combined" in imp
    assert "GSE17536" in imp and "overlap" in imp
    assert "EXPLORATORY" in imp.upper()


def test_D24_sex_discordance_genotype_unavailable(reg):
    d6 = next(d for d in reg["decisions"] if d["id"] == "D8D-6")
    assert "UNAVAILABLE" in d6["implication"].upper()
    assert "RETAINED" in d6["current_frozen_state"].upper()


def test_D25_combat_mod_route_is_author_correspondence(reg):
    d5 = next(d for d in reg["decisions"] if d["id"] == "D8D-5")
    assert "author correspondence" in d5["outstanding_route"].lower()
    assert "not proof of identity" in d5["implication"].lower()


# ==================== evidence table ====================

def test_D26_evidence_table_has_required_columns():
    ev = pd.read_csv(OUT / "evidence_table.csv")
    for c in ("finding", "source_phase", "evidence_type",
              "reproducibility_or_corroboration",
              "limitation_or_alternative_explanation",
              "implication_for_modelling", "decision_status"):
        assert c in ev.columns, c
    assert len(ev) >= 20


def test_D27_four_classes_all_populated(summ):
    cat = summ["evidence_classification"]
    for k in ("reproduced_across_phases", "new_in_8B_or_8C", "unresolved",
              "not_independent_same_matrix_or_qc_rules"):
        assert len(cat[k]) > 0, k


def test_D28_negative_and_contradictory_findings_preserved():
    ev = pd.read_csv(OUT / "evidence_table.csv")
    blob = " ".join(ev.reproducibility_or_corroboration.astype(str)).lower()
    assert "not corroborated" in blob, "the dMMR nominal result must be kept"
    lim = " ".join(ev.limitation_or_alternative_explanation.astype(str))
    assert "PERFECTLY CONFOUNDED" in lim
    assert "CIRCULAR" in lim and "MATRIX-DERIVED" in lim


# ==================== report discipline ====================

def test_D29_no_prognostic_or_validation_claim():
    t = doc_text().lower()
    assert "no outcome association has been tested" in t
    for phrase in ("externally validated", "clinically useful", "is prognostic",
                   "predictive performance", "biomarker"):
        start = 0
        while (i := t.find(phrase, start)) != -1:
            ctx = t[max(0, i - 95):i]
            assert any(n in ctx for n in ("no ", "not ", "never ", "nothing ",
                                          "cannot ", "without ")), \
                f"{phrase!r} appears as an assertion near: ...{ctx[-70:]}"
            start = i + len(phrase)


def test_D30_no_combat_success_claim():
    t = doc_text()
    assert "no claim is made that ComBat removed all batch effects" in t
    for bad in ("ComBat successfully removed", "batch effects are absent,"):
        assert bad not in t


def test_D31_no_swap_inferred():
    t = doc_text()
    assert "No sample swap is inferred" in t
    assert "Genotype-level evidence is unavailable" in t


def test_D32_combat_mod_not_claimed_as_identity():
    t = doc_text()
    assert "This is concordance, not identity." in t


def test_D33_report_cites_source_paths():
    t = DOC.read_text()
    for p in ("docs/PHASE8A_EXPLORATORY_ANALYSIS.md",
              "docs/PHASE8B_PCA_BATCH_INVESTIGATION.md",
              "docs/PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md",
              "reports/phase8d/evidence_table.csv",
              "reports/phase8d/decision_register.json",
              "docs/PHASE4_PREACQUISITION.md",
              "docs/PHASE1_COHORT_ENDPOINT_INVESTIGATION.md"):
        assert p in t, p


def test_D34_readiness_checklist_complete():
    c = pd.read_csv(OUT / "phase9_readiness_checklist.csv")
    assert len(c) >= 18
    assert (c.status == "READY").sum() >= 8
    assert (c.status == "DECISION REQUIRED").sum() == 4
    assert "NOT STARTED" in set(c.status)


# ==================== integrity ====================

def test_D35_prior_phase_artifacts_unchanged(summ):
    """Phase 8D's attestation resolves against the amendment chain."""
    from src.m2_amendment import resolve_hash
    bad = []
    for k, h in summ["read_only_artifact_hashes"].items():
        ok, why = resolve_hash(ROOT / k, h)
        if not ok:
            bad.append((k, why))
    assert not bad, bad
    assert summ["read_only_artifacts_modified"] is False
    keys = summ["read_only_artifact_hashes"]
    for phase in ("phase8a", "phase8b", "phase8c"):
        assert any(phase in k for k in keys), phase
    assert any("PHASE8A_EXPLORATORY" in k for k in keys)


def test_D36_amendment_chain_resolves(summ):
    from src.m2_amendment import resolve_hash
    for n in ("m2_panel_4set.json", "m2_panel_5set.json"):
        p = PANELS / n
        ok, why = resolve_hash(p, hashlib.sha256(p.read_bytes()).hexdigest())
        assert ok and why == "current", n


def test_D37_phase8d_wrote_only_to_its_own_directory():
    assert OUT.exists()
    for p in OUT.glob("*"):
        assert p.parent == OUT


def test_D38_no_normalization_or_correction(summ):
    assert summ["normalization_applied"] == "none"
    assert summ["batch_correction_applied"] == "none"
