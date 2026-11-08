"""Phase 7A: preprocessing compatibility tests.

Verify the identifier compatibility facts underpinning the fRMA hard gate and
the deposited-matrix assessment. These are compatibility checks, NOT a claim
that fRMA has been validated -- the loaded-object check remains outstanding.
"""

from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
import pytest

from src import paths

REPORT = paths.ROOT / "reports" / "phase7" / "preprocessing_compatibility.json"
DOC = paths.ROOT / "docs" / "PHASE7A_PREPROCESSING_FEASIBILITY.md"
AFFY = re.compile(r"^(AFFX-)?[A-Za-z0-9._\-]+_(at|s_at|x_at|a_at|i_at|g_at|r_at|b_at|f_at)$")


@pytest.fixture(scope="module")
def rep():
    return json.loads(REPORT.read_text())


@pytest.fixture(scope="module")
def expr():
    return pd.read_parquet(paths.PROBE_MATRIX, engine="pyarrow")


@pytest.fixture(scope="module")
def frozen():
    med = pd.read_csv(paths.PROCESSED_DIR / "frozen/frozen_gene_to_probeset_median.csv")
    rip = pd.read_csv(paths.PROCESSED_DIR / "frozen/frozen_gene_to_probe_discovery443.csv")
    return set(med.probe_id), set(rip.probe_id)


# ---------- deposited matrix ----------

def test_P01_dimensions_and_scale(rep, expr):
    d = rep["deposited_matrix"]
    assert d["dimensions"] == {"probes": 54675, "samples": 585}
    assert expr.shape == (54675, 585)
    assert d["scale"] == "log2"
    assert d["non_finite"] == 0


def test_P02_fitting_population_undocumented(rep):
    assert rep["deposited_matrix"]["fitting_population_documented"] is False


def test_P03_quantile_normalisation_evidence(rep):
    f = rep["deposited_matrix"]["forensic_evidence"]
    assert f["per_array_median_sd"] < 0.1, "per-array medians should be near-constant"
    assert f["sorted_vectors_exactly_identical"] is False
    assert "quantile" in f["interpretation"]


def test_P04_nontumour_run_inference_is_qualified(rep):
    f = rep["deposited_matrix"]["forensic_evidence"]
    assert f["cross_vs_within_ratio"] > 2.0
    doc = DOC.read_text()
    assert "Inference (not proof)" in doc
    assert "probably not" in doc


# ---------- fRMA identifier gate ----------

def test_P05_probe_sets_identical_to_platform(rep, expr):
    g = rep["frma_gate"]
    assert g["identical_sets"] is True
    assert g["expression_probe_ids"] == g["annotation_probe_ids"] == 54675


def test_P06_identifiers_are_standard_affy(rep, expr):
    nonmatching = [p for p in expr.index if not AFFY.match(p)]
    assert len(nonmatching) == 13
    assert all(p.startswith("AFFX-") for p in nonmatching), nonmatching
    assert rep["frma_gate"]["pct_standard_grammar"] > 99.9


def test_P07_identifier_set_is_standard_cdf_not_entrezg(rep, expr):
    """Discriminate standard CDF from a Brainarray ENTREZG CDF.

    Counting '<digits>_at' is NOT a discriminator -- standard Affymetrix IDs
    have that form too (202763_at). The valid discriminators are the probe-set
    COUNT, the presence of non-plain suffixes, and AFFX controls: an ENTREZG
    CDF has ~20k probe sets, all plain '_at', and no AFFX controls.
    """
    g = rep["frma_gate"]
    assert g["expression_probe_ids"] == 54675
    assert g["non_plain_suffix_probesets"] == 17261
    assert g["suffix_distribution"]["s_at"] == 11295
    assert g["suffix_distribution"]["x_at"] == 4207
    assert g["affx_control_probesets"] == 62
    # independently recomputed from the matrix
    non_plain = sum(any(p.endswith("_" + s) for s in
                        ("s_at", "x_at", "a_at", "g_at", "r_at", "i_at", "f_at"))
                    for p in expr.index)
    assert non_plain == 17261
    assert sum(p.startswith("AFFX-") for p in expr.index) == 62


def test_P08_frozen_sets_are_platform_subsets(rep, frozen, expr):
    med, rip = frozen
    assert med <= set(expr.index)
    assert rip <= set(expr.index)
    assert len(med) == 43163 and len(rip) == 21815
    assert rep["frma_gate"]["frozen_median_subset_of_platform"] is True
    assert rep["frma_gate"]["frozen_repiqr_subset_of_platform"] is True


def test_P09_cdf_ambiguity_documented_resolution(rep):
    d = rep["frma_gate"]["documented_vector_objects"]
    assert d["standard_cdf_object"] == "hgu133plus2frmavecs"
    assert d["entrez_altcdf_object"] == "hgu133plus2hsentrezgfrmavecs"
    assert "if applicable" in " ".join(d["format"]).lower()
    assert "DOCUMENTATION" in d["verified_level"]


def test_P10_gate_is_conditional_not_claimed_complete(rep):
    v = rep["frma_gate"]["verdict"]
    assert v.startswith("CONDITIONAL GO")
    assert "remains outstanding" in v
    assert "requires installation" in v


# ---------- AFFX finding ----------

def test_P11_affx_controls_in_frozen_representations(frozen):
    med, rip = frozen
    affx_med = {p for p in med if p.startswith("AFFX-")}
    affx_rip = {p for p in rip if p.startswith("AFFX-")}
    assert len(affx_med) == 10
    assert len(affx_rip) == 2
    assert "AFFX-HSAC07/X00351_5_at" in affx_rip
    assert "AFFX-HUMGAPDH/M33197_5_at" in affx_rip


def test_P12_affx_genes_absent_from_frozen_m2_panel():
    """The frozen M2 primary panel must contain no AFFX-derived gene."""
    med = pd.read_csv(paths.PROCESSED_DIR / "frozen/frozen_gene_to_probeset_median.csv")
    affx_genes = set(med.loc[med.probe_id.str.startswith("AFFX-"), "gene_symbol"])
    assert affx_genes == {"ACTB", "GAPDH", "STAT1"}
    panel = json.loads((paths.ROOT / "config/pathways/m2_panel_4set.json").read_text())
    members = set()
    for g in panel["gene_sets"].values():
        members |= set(g["members_resolved"])
    assert not (affx_genes & members), affx_genes & members


# ---------- restrictions ----------

def test_P13_no_prohibited_action_taken(rep):
    assert rep["normalization_performed"] == "none"
    assert rep["combat_performed"] == "none"
    assert rep["frma_executed"] is False
    assert rep["packages_installed"] == "none"
    assert rep["cel_files_downloaded"] == 0
    assert rep["outcome_data_used"] is False


def test_P14_frozen_artifacts_untouched():
    import hashlib
    expected = {
        "frozen_gene_to_probeset_median.csv": "08f8715fb7eca5a0",
        "frozen_gene_to_probe_discovery443.csv": "c386d5f11cc511dc",
    }
    for name, pre in expected.items():
        p = paths.PROCESSED_DIR / "frozen" / name
        assert hashlib.sha256(p.read_bytes()).hexdigest().startswith(pre), name
