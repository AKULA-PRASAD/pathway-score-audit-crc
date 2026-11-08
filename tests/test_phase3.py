"""Phase 3 verification suite -- 36 tests in seven groups.

A: cohort membership    B: alignment          C: endpoint
D: leakage prevention   E: integrity          F: reproducibility
G: additional integrity checks (Correction 4)

Establishes computational correctness and leakage prevention. It does NOT
establish biological validity: the deposited values' processing provenance is
incompletely documented (see src/provenance.py) and no normalization or batch
correction has been applied.
"""

from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
import pytest

from src import paths
from src.aggregation import aggregate
from src.bundle import ANALYSIS_DIR, align_by_id, sha256_file
from src.cache import load_probe_matrix
from src.cohort import COHORTS, exclusion_table, load_manifest, materialise
from src.endpoint import HORIZON_UNITS, build_endpoint
from src.gene_mapping import build_mapping
from src.provenance import FORBIDDEN_PROVENANCE_PHRASES

FROZEN = paths.PROCESSED_DIR / "frozen"
CAND = paths.PROCESSED_DIR / "candidates"

EXPECTED = {"S3": (453, 136, 130), "S0": (461, 139, 133),
            "S1": (458, 136, 130), "S2": (456, 139, 133)}
EXPECTED_PARTITION = {
    "S3": {"discovery": (354, 107), "validation": (99, 29)},
    "S0": {"discovery": (359, 109), "validation": (102, 30)},
    "S1": {"discovery": (357, 107), "validation": (101, 29)},
    "S2": {"discovery": (356, 109), "validation": (100, 30)},
}
C1_SAMPLES = {"GSM972062", "GSM972206", "GSM972402"}
C2_SAMPLES = {"GSM972229", "GSM972246", "GSM972258", "GSM972401", "GSM972456"}


@pytest.fixture(scope="module")
def manifest():
    return json.loads((ANALYSIS_DIR / "analysis_bundle_manifest.json").read_text())


@pytest.fixture(scope="module")
def membership():
    return pd.read_csv(ANALYSIS_DIR / "cohort_membership.csv")


@pytest.fixture(scope="module")
def endpoint():
    return pd.read_csv(ANALYSIS_DIR / "endpoint.csv")


@pytest.fixture(scope="module")
def partition():
    return pd.read_csv(ANALYSIS_DIR / "partition_and_batch.csv")


@pytest.fixture(scope="module")
def med():
    return pd.read_parquet(ANALYSIS_DIR / "expression_median_S0.parquet")


@pytest.fixture(scope="module")
def rep():
    return pd.read_parquet(ANALYSIS_DIR / "expression_repiqr_disc443_S0.parquet")


@pytest.fixture(scope="module")
def disc443():
    return (FROZEN / "selection_population_discovery443.txt").read_text().split()


def ids_of(membership, cohort):
    return membership.loc[membership[f"in_{cohort}"].astype(bool), "sample_id"].tolist()


# ===================== A. cohort membership =====================

def test_T01_s3_exact_counts(membership, endpoint):
    ids = set(ids_of(membership, "S3"))
    e = endpoint[endpoint.sample_id.isin(ids)]
    assert len(e) == 453
    assert int((e.rfs_event_full == 1).sum()) == 136
    assert int((e.rfs_event_60units == 1).sum()) == 130


def test_T02_all_cohort_counts(membership, endpoint):
    for c, (n, ev, ev60) in EXPECTED.items():
        ids = set(ids_of(membership, c))
        e = endpoint[endpoint.sample_id.isin(ids)]
        assert (len(e), int((e.rfs_event_full == 1).sum()),
                int((e.rfs_event_60units == 1).sum())) == (n, ev, ev60), c


def test_T03_partition_counts(membership, endpoint, partition):
    pm = partition.set_index("sample_id")["dataset"]
    for c, exp in EXPECTED_PARTITION.items():
        ids = ids_of(membership, c)
        e = endpoint[endpoint.sample_id.isin(set(ids))].copy()
        e["dataset"] = e.sample_id.map(pm)
        for ds, (n, ev) in exp.items():
            g = e[e.dataset == ds]
            assert (len(g), int((g.rfs_event_full == 1).sum())) == (n, ev), (c, ds)


def test_T04_subset_relations(membership):
    s = {c: set(ids_of(membership, c)) for c in COHORTS}
    assert s["S3"] < s["S1"] < s["S0"]
    assert s["S3"] < s["S2"] < s["S0"]


def test_T05_exclusions_named_with_criterion(membership):
    m = membership.set_index("sample_id")["exclusion_criterion"].fillna("")
    c1 = set(m[m == "fails_criterion_1"].index)
    c2 = set(m[m == "fails_criterion_2"].index)
    assert c1 == C1_SAMPLES
    assert c2 == C2_SAMPLES


def test_T06_membership_matches_phase1_manifest(membership):
    man = load_manifest()
    ref = load_probe_matrix().columns.tolist()
    for c in COHORTS:
        assert set(ids_of(membership, c)) == set(
            materialise(c, man, ref).sample_ids), c


# ===================== B. alignment =====================

def test_T07_alignment_is_label_based_not_positional(med, membership):
    ids = ids_of(membership, "S3")
    straight = align_by_id(med, ids)
    shuffled = med.sample(axis=1, frac=1.0, random_state=7)
    assert list(shuffled.columns) != list(med.columns)
    from_shuffled = align_by_id(shuffled, ids)
    assert list(from_shuffled.columns) == list(straight.columns)
    assert np.array_equal(from_shuffled.to_numpy(), straight.to_numpy())


def test_T08_columns_equal_manifest_in_order(med, rep, membership):
    for c in COHORTS:
        ids = ids_of(membership, c)
        assert list(align_by_id(med, ids).columns) == ids, c
        assert list(align_by_id(rep, ids).columns) == ids, c


def test_T09_reindex_introduces_no_nan(med, rep, membership):
    for c in COHORTS:
        ids = ids_of(membership, c)
        assert not align_by_id(med, ids).isna().to_numpy().any()
        assert not align_by_id(rep, ids).isna().to_numpy().any()


def test_T10_representations_share_axes(med, rep):
    assert list(med.index) == list(rep.index)
    assert list(med.columns) == list(rep.columns)


def test_T11_missing_sample_raises(med):
    with pytest.raises(KeyError):
        align_by_id(med, [med.columns[0], "GSM_DOES_NOT_EXIST"])


# ===================== C. endpoint =====================

def test_T12_source_columns_byte_identical(endpoint):
    meta = pd.read_csv(paths.SAMPLE_METADATA).set_index("sample_id")
    e = endpoint.set_index("sample_id")
    for src, col in [("rfs.event", "rfs_event_source"),
                     ("rfs.delay", "rfs_time_source"),
                     ("os.event", "os_event_source"),
                     ("os.delay_(months)", "os_time_source")]:
        a = meta[src].reindex(e.index).to_numpy()
        b = e[col].to_numpy()
        assert np.allclose(a, b, equal_nan=True), col


def test_T13_full_equals_source(endpoint):
    assert np.allclose(endpoint.rfs_event_full, endpoint.rfs_event_source,
                       equal_nan=True)
    assert np.allclose(endpoint.rfs_time_full, endpoint.rfs_time_source,
                       equal_nan=True)


def test_T14_horizon_rule_rowwise(endpoint, partition):
    over = endpoint.rfs_time_source > HORIZON_UNITS
    assert (endpoint.loc[over, "rfs_time_60units"] == HORIZON_UNITS).all()
    assert (endpoint.loc[over, "rfs_event_60units"] == 0).all()
    assert np.allclose(endpoint.loc[~over, "rfs_time_60units"],
                       endpoint.loc[~over, "rfs_time_source"])
    pm = partition.set_index("sample_id")["dataset"]
    recl = endpoint[endpoint.reclassified_by_horizon]
    assert len(recl) == 6
    assert set(recl.sample_id.map(pm)) == {"discovery"}


def test_T15_validation_reclassifies_zero(endpoint, partition, membership):
    pm = partition.set_index("sample_id")["dataset"]
    e = endpoint.copy()
    e["dataset"] = e.sample_id.map(pm)
    val = e[e.dataset == "validation"]
    assert int(val.reclassified_by_horizon.sum()) == 0
    assert val.loc[val.rfs_event_full == 1, "rfs_time_full"].max() <= HORIZON_UNITS


def test_T16_assumption_fields_present(manifest):
    ep = manifest["endpoint"]
    assert ep["time_unit_assumption"].strip()
    assert "PROVISIONAL" in ep["time_unit_assumption"]
    assert ep["censoring_interpretation"].strip()
    assert "PROVISIONAL" in ep["censoring_interpretation"]


def test_T17_no_year_named_columns(endpoint, med, rep):
    bad = re.compile(r"(_5y|_5yr|_years?$|five_?year)", re.I)
    for df in (endpoint,):
        assert not [c for c in df.columns if bad.search(str(c))]
    assert not [c for c in med.columns if bad.search(str(c))]


# ===================== D. leakage prevention =====================

def test_T18_no_validation_sample_in_selection_population(disc443):
    meta = pd.read_csv(paths.SAMPLE_METADATA)
    val = set(meta.loc[meta.dataset == "validation", "sample_id"])
    assert len(disc443) == 443
    assert not (set(disc443) & val)


def test_T19_frozen_map_reproduces_from_discovery(disc443):
    expr = load_probe_matrix()
    um = build_mapping().unique_map
    got = aggregate(expr, um, "representative_iqr",
                    selection_columns=disc443).representative_probe
    frozen = pd.read_csv(FROZEN / "frozen_gene_to_probe_discovery443.csv",
                         index_col=0)["probe_id"]
    assert got.reindex(frozen.index).equals(frozen)


def test_T20_map_unchanged_when_validation_deleted(disc443):
    """Physically delete every non-discovery column, then re-select."""
    expr = load_probe_matrix()
    um = build_mapping().unique_map
    stripped = expr[disc443].copy()
    assert stripped.shape[1] == 443
    got = aggregate(stripped, um, "representative_iqr",
                    selection_columns=disc443).representative_probe
    frozen = pd.read_csv(FROZEN / "frozen_gene_to_probe_discovery443.csv",
                         index_col=0)["probe_id"]
    assert got.reindex(frozen.index).equals(frozen)


def test_T21_no_transformation_fitted(manifest):
    assert manifest["transformation_fitted_in_phase3"] == "none"
    assert manifest["normalization_applied"] == "none"
    assert manifest["batch_correction_applied"] == "none"


def test_T22_no_endpoint_column_in_expression_modules():
    pat = re.compile(r"rfs\.|os\.event|os\.delay|rfs_event|rfs_time")
    for name in ("gene_mapping.py", "aggregation.py", "external_annotation.py"):
        assert not pat.search((paths.ROOT / "src" / name).read_text()), name


def test_T23_gene_set_identical_across_cohorts(med, rep, membership):
    ref = list(med.index)
    for c in COHORTS:
        ids = ids_of(membership, c)
        assert list(align_by_id(med, ids).index) == ref, c
        assert list(align_by_id(rep, ids).index) == ref, c


# ===================== E. integrity =====================

def test_T24_input_hashes_match(manifest):
    paths_map = {
        "probe_matrix": paths.PROBE_MATRIX,
        "frozen_median_spec": FROZEN / "frozen_gene_to_probeset_median.csv",
        "frozen_repiqr_spec": FROZEN / "frozen_gene_to_probe_discovery443.csv",
        "selection_population": FROZEN / "selection_population_discovery443.txt",
        "cohort_manifest": paths.ROOT / "reports/phase1/cohort_manifest.csv",
        "sample_metadata": paths.SAMPLE_METADATA,
        "sdrf": paths.ROOT / "data/external/E-GEOD-39582.sdrf.txt",
    }
    for k, p in paths_map.items():
        assert manifest["input_hashes"][k] == sha256_file(p), k


def test_T25_matrices_finite_float64(med, rep):
    for m in (med, rep):
        assert str(m.dtypes.iloc[0]) == "float64"
        assert np.isfinite(m.to_numpy()).all()
        assert m.shape == (21815, 461)


def test_T26_upstream_artifacts_unchanged():
    from src.cache import cache_is_current
    ok, _ = cache_is_current()
    assert ok
    assert len(list(CAND.glob("*.parquet"))) == 6


# ===================== F. reproducibility =====================

def test_T27_outputs_hash_as_recorded(manifest):
    for name, rec in manifest["outputs"].items():
        p = ANALYSIS_DIR / name
        assert p.exists(), name
        assert sha256_file(p) == rec["sha256"], name


def test_T28_endpoint_rebuild_deterministic(membership):
    meta = pd.read_csv(paths.SAMPLE_METADATA)
    ids = ids_of(membership, "S0")
    a = build_endpoint(meta, ids)
    b = build_endpoint(meta, ids)
    assert a.equals(b)
    stored = pd.read_csv(ANALYSIS_DIR / "endpoint.csv")
    assert list(a.columns) == list(stored.columns)
    assert np.allclose(a.drop(columns=["sample_id", "reclassified_by_horizon"]),
                       stored.drop(columns=["sample_id",
                                            "reclassified_by_horizon"]),
                       equal_nan=True)


# ===================== G. additional integrity (Correction 4) =====

def test_T29_selection_cannot_access_validation_columns(disc443):
    """The selection procedure runs on a matrix that HAS no validation data."""
    expr = load_probe_matrix()
    um = build_mapping().unique_map
    discovery_only = expr[disc443].copy()
    meta = pd.read_csv(paths.SAMPLE_METADATA)
    val = set(meta.loc[meta.dataset == "validation", "sample_id"])
    assert not (set(discovery_only.columns) & val)
    res = aggregate(discovery_only, um, "representative_iqr",
                    selection_columns=disc443)
    frozen = pd.read_csv(FROZEN / "frozen_gene_to_probe_discovery443.csv",
                         index_col=0)["probe_id"]
    got = res.representative_probe.reindex(frozen.index)
    assert got.equals(frozen)
    assert res.diagnostics["selection_samples"] == 443


def test_T30_provenance_guard_on_selection_population(disc443):
    """A selection labelled discovery_443 must contain no validation sample."""
    meta = pd.read_csv(paths.SAMPLE_METADATA)
    val = sorted(meta.loc[meta.dataset == "validation", "sample_id"])
    contaminated = disc443 + [val[0]]
    assert len(contaminated) == 444
    assert set(contaminated) & set(val)
    # The frozen artifact must never have been built from such a list.
    stored = (FROZEN / "selection_population_discovery443.txt").read_text().split()
    assert len(stored) == 443
    assert not (set(stored) & set(val))


def test_T31_batch_preserved_not_injected(partition, med):
    assert partition["combat_batch"].notna().all()
    assert partition["combat_batch"].nunique() == 9
    cov = pd.read_csv(ANALYSIS_DIR / "covariates.csv")
    assert not [c for c in cov.columns if "batch" in str(c).lower()]
    assert not [g for g in med.index if "batch" in str(g).lower()]


def test_T32_no_batch_correction_artifact(manifest):
    assert manifest["batch_correction_applied"] == "none"
    assert manifest["batch_covariate_auto_included"] is False
    names = " ".join(manifest["outputs"].keys()).lower()
    for bad in ("combat", "corrected", "adjusted", "residual", "batchfit"):
        assert bad not in names, bad


def test_T33_os_fields_preserved_separately(endpoint):
    meta = pd.read_csv(paths.SAMPLE_METADATA).set_index("sample_id")
    e = endpoint.set_index("sample_id")
    assert np.allclose(meta["os.event"].reindex(e.index),
                       e["os_event_source"], equal_nan=True)
    assert np.allclose(meta["os.delay_(months)"].reindex(e.index),
                       e["os_time_source"], equal_nan=True)
    # OS must not be folded into any RFS column.
    assert not np.allclose(e["rfs_event_full"], e["os_event_source"],
                           equal_nan=True)


def test_T34_endpoint_hierarchy_recorded(manifest):
    h = manifest["endpoint"]["hierarchy"]
    assert h["rfs_full"] == "primary"
    assert h["rfs_60units"].startswith("conditional")
    assert h["os"].startswith("separate_preserved")


def test_T35_processing_provenance_wording(manifest):
    prov = manifest["processing_provenance"].lower()
    assert "does not document the fitting population" in prov
    for bad in FORBIDDEN_PROVENANCE_PHRASES:
        assert bad not in prov, bad
    blob = json.dumps(manifest).lower()
    assert "jointly rma" not in blob
    assert "across all 585" not in blob


def test_T36_batch_imbalance_recorded(manifest):
    comp = manifest["cohorts"]["S3"]["batch_composition_by_partition"]
    assert comp["batch01"] == {"discovery": 33, "validation": 0}
    assert comp["batch02"] == {"discovery": 14, "validation": 0}
    assert comp["batch09"] == {"discovery": 1, "validation": 14}
    assert "NOT imply independence" in manifest["partition_independence"]
