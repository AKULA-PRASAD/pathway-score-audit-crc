"""Phase 6C: M2 pathway specification tests — 22 tests in five groups.

A: freeze integrity   B: specification validity   C: membership reproducibility
D: outcome-blindness  E: immutability policy

Establishes that both specifications are valid, reproducible from the
authoritative MSigDB source, and unchanged since the freeze. Does not
establish biological validity of the panel choice.
"""

from __future__ import annotations

import hashlib
import json
import re

import pandas as pd
import pytest

from src import paths

PANEL_DIR = paths.ROOT / "config" / "pathways"
PRIMARY = PANEL_DIR / "m2_panel_4set.json"
SENSITIVITY = PANEL_DIR / "m2_panel_5set.json"
RECORD = PANEL_DIR / "M2_FREEZE_RECORD.json"
GMT_SYM = paths.ROOT / "reports" / "phase6c" / "h.all.v2026.2.Hs.symbols.gmt"
GMT_ENT = paths.ROOT / "reports" / "phase6c" / "h.all.v2026.2.Hs.entrez.gmt"

PRIMARY_SETS = {
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY",
}
PI3K = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
EXPECTED_N = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": 198,
              "HALLMARK_E2F_TARGETS": 195,
              "HALLMARK_WNT_BETA_CATENIN_SIGNALING": 41,
              "HALLMARK_P53_PATHWAY": 194, PI3K: 104}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


@pytest.fixture(scope="module")
def primary():
    return json.loads(PRIMARY.read_text())


@pytest.fixture(scope="module")
def sensitivity():
    return json.loads(SENSITIVITY.read_text())


@pytest.fixture(scope="module")
def record():
    return json.loads(RECORD.read_text())


@pytest.fixture(scope="module")
def gmt_symbols():
    d = {}
    for line in GMT_SYM.read_text().splitlines():
        p = line.split("\t")
        d[p[0]] = [x for x in p[2:] if x]
    return d


# ===================== A. freeze integrity =====================

def test_M01_freeze_record_exists(record):
    assert record["record_type"] == "M2_PATHWAY_SPECIFICATION_FREEZE"
    assert record["frozen_utc"]


def test_M02_specs_match_recorded_hashes(record):
    for role, path in (("primary", PRIMARY), ("sensitivity", SENSITIVITY)):
        assert sha256(path) == record["specifications"][role]["sha256"], role


def test_M03_both_specs_marked_frozen(primary, sensitivity):
    for d in (primary, sensitivity):
        assert d["status"] == "FROZEN"
        assert d["immutable"] is True
        assert d["frozen_utc"]


def test_M04_roles_assigned(primary, sensitivity):
    assert primary["role"] == "PRIMARY"
    assert sensitivity["role"].startswith("SENSITIVITY")


def test_M05_gmt_checksums_preserved(record, primary):
    db = record["source_database"]
    assert db["version"] == "v2026.2.Hs"
    assert db["release_date"] == "2026-10-15"
    assert db["sha256_symbols"] == sha256(GMT_SYM)
    assert db["sha256_entrez"] == sha256(GMT_ENT)
    assert primary["database"]["sha256_symbols"] == db["sha256_symbols"]


# ===================== B. specification validity =====================

def test_M06_primary_is_four_sets(primary):
    assert primary["n_sets"] == 4
    assert set(primary["gene_sets"]) == PRIMARY_SETS


def test_M07_sensitivity_is_five_sets(sensitivity):
    assert sensitivity["n_sets"] == 5
    assert set(sensitivity["gene_sets"]) == PRIMARY_SETS | {PI3K}


def test_M08_pi3k_is_sensitivity_only(primary, sensitivity, record):
    assert PI3K not in primary["gene_sets"], "PI3K must NOT be in the primary panel"
    assert PI3K in sensitivity["gene_sets"]
    assert record["approval"]["sensitivity_only_addition"] == PI3K
    clause = record["approval"]["immutability_clause"]
    assert "must NOT be promoted" in clause
    assert "must NOT be removed" in clause
    assert "observed model performance" in clause


def test_M09_sensitivity_is_primary_plus_pi3k(primary, sensitivity):
    """The sensitivity panel must be a strict superset differing only by PI3K."""
    assert set(sensitivity["gene_sets"]) - set(primary["gene_sets"]) == {PI3K}
    for s in primary["gene_sets"]:
        assert (primary["gene_sets"][s]["members_resolved"]
                == sensitivity["gene_sets"][s]["members_resolved"]), s


def test_M10_member_counts_as_specified(primary, sensitivity):
    for d in (primary, sensitivity):
        for s, g in d["gene_sets"].items():
            assert g["n_resolved_unique"] == EXPECTED_N[s], s
            assert len(g["members_resolved"]) == EXPECTED_N[s], s


def test_M11_no_duplicate_members(primary, sensitivity):
    for d in (primary, sensitivity):
        for s, g in d["gene_sets"].items():
            m = g["members_resolved"]
            assert len(m) == len(set(m)), s
            assert g["duplicates_within_set"] == 0, s


def test_M12_unmatched_reconciles(primary, sensitivity):
    """source - resolved must equal the true-unmatched count exactly."""
    for d in (primary, sensitivity):
        for s, g in d["gene_sets"].items():
            assert g["n_source"] - g["n_resolved_unique"] == g["true_unmatched_n"], s
            assert len(g["true_unmatched"]) == g["true_unmatched_n"], s


def test_M13_union_counts_correct(primary, sensitivity):
    for d in (primary, sensitivity):
        u = set()
        for g in d["gene_sets"].values():
            u |= set(g["members_resolved"])
        assert d["union_unique_genes"] == len(u)
    assert primary["union_unique_genes"] == 607
    assert sensitivity["union_unique_genes"] == 702


# ===================== C. membership reproducibility =====================

def test_M14_sets_exist_in_source_gmt(primary, sensitivity, gmt_symbols):
    for d in (primary, sensitivity):
        for s, g in d["gene_sets"].items():
            assert s in gmt_symbols, s
            assert g["n_source"] == len(gmt_symbols[s]), s


def test_M15_members_are_subset_of_frozen_index(primary, sensitivity):
    med = pd.read_csv(paths.PROCESSED_DIR / "frozen"
                      / "frozen_gene_to_probeset_median.csv")
    idx = set(med.gene_symbol.unique())
    assert len(idx) == 21815
    for d in (primary, sensitivity):
        for s, g in d["gene_sets"].items():
            missing = set(g["members_resolved"]) - idx
            assert not missing, f"{s}: {sorted(missing)[:5]}"


def test_M16_each_gene_exactly_one_matrix_row(primary):
    """Scores use the gene-level matrix: one row per gene, no duplication."""
    X = pd.read_parquet(paths.PROCESSED_DIR / "analysis"
                        / "expression_median_S0.parquet")
    assert X.index.is_unique
    for s, g in primary["gene_sets"].items():
        present = [x for x in g["members_resolved"] if x in X.index]
        assert len(present) == len(set(present)), s
        assert len(present) == g["n_resolved_unique"], s


def test_M17_true_unmatched_genuinely_absent(primary, gmt_symbols):
    """Every true-unmatched entry names a gene absent from our index."""
    med = pd.read_csv(paths.PROCESSED_DIR / "frozen"
                      / "frozen_gene_to_probeset_median.csv")
    idx = set(med.gene_symbol.unique())
    for s, g in primary["gene_sets"].items():
        for entry in g["true_unmatched"]:
            sym = entry.split(" (")[0]
            assert sym not in idx, f"{s}: {sym} IS in the index"


# ===================== D. outcome-blindness =====================

def test_M18_declared_outcome_blind(primary, sensitivity, record):
    for d in (primary, sensitivity):
        assert d["outcome_data_used"] is False
        assert "No recurrence, survival or outcome" in d["outcome_blind_statement"]
    a = record["approval"]
    assert a["outcome_blind"] is True
    assert a["frozen_before_outcome_analysis"] is True
    assert any("outcome" in c for c in a["selection_criteria_NOT_used"])


def test_M19_no_outcome_field_in_specs(primary, sensitivity):
    pat = re.compile(r"rfs[._]|os[._]event|os[._]delay|recurrence_event|survival_time",
                     re.I)
    for d in (primary, sensitivity):
        blob = json.dumps(d)
        hits = [m.group(0) for m in pat.finditer(blob)]
        assert not hits, hits


def test_M20_scoring_is_fold_safe(primary):
    sc = primary["scoring"]
    assert "TRAINING FOLD SAMPLES ONLY" in sc["mu_sigma_estimated_on"]
    assert "training-derived" in sc["applied_to_heldout"]
    assert "EXCLUDED" in sc["zero_variance_genes"]
    assert "never imputed" in sc["missing_gene_members"]
    assert sc["no_new_packages_required"] is True


def test_M21_minimum_coverage_satisfied(primary, sensitivity):
    """All sets meet the >=80% resolved and >=20 gene coverage rule."""
    for d in (primary, sensitivity):
        assert ">=80%" in d["scoring"]["minimum_set_coverage"]
        for s, g in d["gene_sets"].items():
            assert g["n_resolved_unique"] / g["n_source"] >= 0.80, s
            assert g["n_resolved_unique"] >= 20, s


# ===================== E. immutability =====================

def test_M22_modification_policy_recorded(primary, sensitivity):
    for d in (primary, sensitivity):
        pol = d["modification_policy"]
        assert "Immutable" in pol
        assert "M3-exploratory" in pol
        assert "non-outcome" in pol
