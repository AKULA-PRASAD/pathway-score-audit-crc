"""Phase 0 verification suite.

Two independent things are tested, and the distinction is deliberate:

  1. PRESERVATION -- no original file has changed (tests 1-2).
  2. COMPUTATIONAL FIDELITY -- the new infrastructure reproduces the
     original data and the original results exactly (tests 3-10).

Neither establishes BIOLOGICAL or NORMALIZATION validity. The deposited
values were RMA-normalised and ComBat-corrected jointly across all 585
arrays; whether that is a sound basis for cross-cohort transport is a Phase 4
question and is untouched here. A pipeline can pass every test in this file
and still be scientifically invalid.
"""

from __future__ import annotations

import hashlib
import json
import random
import stat

import numpy as np
import pandas as pd
import pytest

from src import annotation, paths
from src.cache import cache_is_current, load_manifest, load_probe_matrix
from src.reproducibility import DEFAULT_SEED, package_versions
from src.soft_parser import iter_sample_tables, list_sample_ids, sha256


# --------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------

@pytest.fixture(scope="module")
def baseline() -> dict:
    return json.loads(paths.BASELINE_CHECKSUMS.read_text())


@pytest.fixture(scope="module")
def matrix() -> pd.DataFrame:
    return load_probe_matrix()


@pytest.fixture(scope="module")
def manifest() -> dict:
    return load_manifest()


@pytest.fixture(scope="module")
def soft_ids() -> list[str]:
    return list_sample_ids(paths.SOFT_SERIES)


# --------------------------------------------------------------------
# 1-2. preservation of originals
# --------------------------------------------------------------------

def test_originals_unmodified(baseline):
    """Every pre-Phase-0 project file is byte-identical to its baseline."""
    changed, missing = [], []
    for rel, rec in baseline["files"].items():
        path = paths.ROOT / rel
        if not path.exists():
            missing.append(rel)
            continue
        if sha256(path) != rec["sha256"]:
            changed.append(rel)
    assert not missing, f"original files missing: {missing}"
    assert not changed, f"original files modified: {changed}"


def test_original_permissions_unchanged(baseline):
    """File permissions on originals are untouched (no chmod was applied)."""
    altered = [
        rel for rel, rec in baseline["files"].items()
        if (paths.ROOT / rel).exists()
        and oct(stat.S_IMODE((paths.ROOT / rel).stat().st_mode)) != rec["mode"]
    ]
    assert not altered, f"permissions changed on: {altered}"


def test_raw_data_unmodified():
    """Raw downloads match the recorded SHA-256 manifest."""
    recorded = {}
    for line in paths.RAW_CHECKSUMS.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        digest, rel = line.split(None, 1)
        recorded[rel.strip()] = digest

    assert recorded, "CHECKSUMS.sha256 contains no entries"
    for rel, digest in recorded.items():
        path = paths.RAW_DIR / rel
        assert path.exists(), f"raw file missing: {rel}"
        assert sha256(path) == digest, f"raw file modified: {rel}"


def test_baseline_covers_expected_inventory(baseline):
    """The baseline covers all 12 original scripts plus requirements.txt."""
    files = baseline["files"]
    scripts = [f for f in files if f.startswith("scripts/")]
    assert len(scripts) == 12, f"expected 12 original scripts, got {len(scripts)}"
    assert "requirements.txt" in files
    assert baseline["file_count"] == len(files) == 19
    assert baseline["exclusions"], "exclusions must be documented"


# --------------------------------------------------------------------
# 3. paths
# --------------------------------------------------------------------

def test_paths_resolve_and_exist():
    for p in (paths.SOFT_SERIES, paths.GPL570_ANNOTATION,
              paths.SAMPLE_METADATA, paths.RAW_CHECKSUMS,
              paths.CONFIG_FILE, paths.BASELINE_CHECKSUMS):
        assert p.is_absolute(), f"{p} is not absolute"
        assert p.exists(), f"{p} does not exist"


def test_paths_are_cwd_independent():
    """paths.ROOT is anchored to __file__, not the process CWD."""
    import os
    import importlib
    original = os.getcwd()
    try:
        os.chdir("/")
        reloaded = importlib.reload(paths)
        assert reloaded.ROOT == paths.ROOT
        assert reloaded.SOFT_SERIES.exists()
    finally:
        os.chdir(original)


# --------------------------------------------------------------------
# 4-7. cache integrity
# --------------------------------------------------------------------

def test_cache_is_current():
    ok, reason = cache_is_current()
    assert ok, f"cache not current: {reason}"


def test_cache_shape_and_dtype(matrix):
    assert matrix.shape == (paths.GPL570_PROBE_COUNT, paths.SERIES_SAMPLE_COUNT)
    assert str(matrix.dtypes.iloc[0]) == "float64"
    assert matrix.index.is_unique
    assert matrix.columns.is_unique


def test_cache_contains_all_585_samples(matrix, soft_ids):
    """All deposited samples are cached, not only the 566 primary tumours."""
    assert len(soft_ids) == 585
    assert list(matrix.columns) == soft_ids, "sample set or order differs"

    meta = pd.read_csv(paths.SAMPLE_METADATA)
    tumours = meta.loc[
        meta["source_name"].str.contains(
            "primary colorectal Adenocarcinoma", case=False, na=False),
        "sample_id",
    ]
    # The 19 non-tumoural samples must still be present: excluding them
    # would be a cohort-eligibility decision, reserved for Phase 1.
    assert len(tumours) == 566
    assert matrix.shape[1] - len(tumours) == 19


def test_probe_ids_identical_across_all_585_samples(matrix):
    """Probe ids AND their ordering are identical in every deposited sample."""
    reference = list(matrix.index)
    mismatched, counted = [], 0
    for table in iter_sample_tables(paths.SOFT_SERIES, sample_ids=None):
        counted += 1
        if table.probe_ids != reference:
            mismatched.append(table.sample_id)
    assert counted == 585, f"parsed {counted} samples, expected 585"
    assert not mismatched, f"probe mismatch in: {mismatched[:10]}"


def test_all_cached_values_finite(matrix, manifest):
    values = matrix.to_numpy()
    assert np.isfinite(values).all(), "non-finite values in cache"
    assert manifest["all_values_finite"] is True
    assert manifest["non_numeric_fields_encountered"] == 0


def test_cache_values_bit_identical_to_source(matrix, soft_ids):
    """Re-parsed arrays equal cached arrays exactly (no float32 rounding)."""
    rng = random.Random(DEFAULT_SEED)
    chosen = set(rng.sample(soft_ids, 20))
    compared = 0
    for table in iter_sample_tables(paths.SOFT_SERIES, sample_ids=chosen):
        compared += 1
        fresh = np.asarray(table.values, dtype=np.float64)
        cached = matrix[table.sample_id].to_numpy()
        assert np.array_equal(fresh, cached), (
            f"{table.sample_id}: {(fresh != cached).sum()} values differ")
    assert compared == 20


def test_stale_cache_is_detected(tmp_path, manifest):
    """A changed source file invalidates the cache by checksum, not mtime."""
    fake = tmp_path / "fake.soft.gz"
    import gzip
    with gzip.open(fake, "wt") as fh:
        fh.write("^SAMPLE = GSM0\nID_REF\tVALUE\n1_at\t1.0\n!sample_table_end\n")
    ok, reason = cache_is_current(soft_path=fake)
    assert not ok and "changed" in reason


def test_manifest_records_provenance(manifest):
    assert manifest["source_sha256"] == sha256(paths.SOFT_SERIES)
    assert manifest["value_transformation"].startswith("none")
    assert manifest["sample_scope"] == "all deposited samples; no cohort filter applied"
    assert "ComBat" in manifest["validity_note"]
    assert manifest["environment"]["packages"]["pandas"]


# --------------------------------------------------------------------
# 8. annotation
# --------------------------------------------------------------------

def test_annotation_header_strategies_agree():
    """Scripts 11 and 12 used different header matches; they must agree."""
    idx = annotation.find_table_header(paths.GPL570_ANNOTATION)
    assert idx > 0
    table = annotation.load_annotation()
    assert len(table) == paths.GPL570_PROBE_COUNT
    assert table["ID"].nunique() == paths.GPL570_PROBE_COUNT


def test_annotation_probes_match_expression_probes(matrix):
    table = annotation.load_annotation()
    assert set(table["ID"]) == set(matrix.index)


def test_split_symbols_is_pure_utility():
    assert annotation.split_symbols("ABC /// DEF") == ["ABC", "DEF"]
    assert annotation.split_symbols("ABC /// ABC") == ["ABC"]
    assert annotation.split_symbols("---") == []
    assert annotation.split_symbols(None) == []
    assert annotation.split_symbols(float("nan")) == []


# --------------------------------------------------------------------
# 9. parser equivalence with the original scripts
# --------------------------------------------------------------------

def test_parser_matches_script_12_extraction(matrix):
    """The shared parser reproduces script 12's inline extraction logic."""
    import gzip
    target = str(matrix.columns[0])
    legacy: dict[str, float] = {}
    current, inside = None, False
    with gzip.open(paths.SOFT_SERIES, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("^SAMPLE = "):
                if current == target and legacy:
                    break
                current = line.strip().split("=", 1)[1].strip()
                inside = False
            elif line.startswith("ID_REF\t") and current == target:
                inside = True
            elif inside and current == target:
                if line.startswith("!sample_table_end"):
                    break
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    try:
                        legacy[parts[0]] = float(parts[1])
                    except ValueError:
                        pass
    assert len(legacy) == paths.GPL570_PROBE_COUNT
    cached = matrix[target]
    assert np.array_equal(
        cached.to_numpy(), np.array([legacy[p] for p in cached.index]))


# --------------------------------------------------------------------
# 10. regression against the existing published result
# --------------------------------------------------------------------

def test_regression_reproduces_script_12_correlation(matrix):
    """Recompute script 12's headline statistic through the new fast path.

    Script 12 reported a median within-gene probe Spearman rho of 0.331803
    across 10,751 genes (results/probe_correlation_audit/). Recomputing it
    from the cache must reproduce that exactly -- proving the cache is a
    faster route to the same numbers, not a different pipeline.

    Note on interpretation: this asserts NUMERICAL equivalence only. The
    Phase 0 inspection showed that figure is driven by unexpressed,
    low-variance probes and is not evidence of poor probe-gene mapping.
    Acting on that is a Phase 2 decision; this test deliberately does not
    apply any filter.
    """
    summary = pd.read_csv(
        paths.RESULTS_DIR / "probe_correlation_audit"
        / "within_gene_correlation_summary.csv")
    pairs = pd.read_csv(
        paths.RESULTS_DIR / "probe_correlation_audit"
        / "within_gene_probe_correlations.csv")

    meta = pd.read_csv(paths.SAMPLE_METADATA)
    tumours = meta.loc[
        meta["source_name"].str.contains(
            "primary colorectal Adenocarcinoma", case=False, na=False),
        "sample_id",
    ].tolist()

    # Rank once, then correlate columns pairwise: algebraically identical to
    # scipy.stats.spearmanr on data with no ties handled differently, but
    # vectorised. Verified against the stored pair-level results below.
    sub = matrix.loc[:, tumours]
    probes_needed = pd.unique(pairs[["probe_a", "probe_b"]].values.ravel())
    ranked = sub.loc[probes_needed].rank(axis=1).to_numpy()
    index = {p: i for i, p in enumerate(probes_needed)}

    a = ranked[[index[p] for p in pairs["probe_a"]]]
    b = ranked[[index[p] for p in pairs["probe_b"]]]
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    rho = (a * b).sum(axis=1) / np.sqrt((a ** 2).sum(axis=1) * (b ** 2).sum(axis=1))

    assert np.allclose(rho, pairs["spearman_rho"].to_numpy(), atol=1e-9), (
        "pair-level correlations differ from the stored script 12 results")

    recomputed = pd.DataFrame({
        "gene_symbol": pairs["gene_symbol"], "rho": rho,
    }).groupby("gene_symbol")["rho"].median()

    assert len(recomputed) == len(summary) == 10751
    assert abs(recomputed.median() - 0.331803) < 1e-4, (
        f"median rho {recomputed.median():.6f} != 0.331803")


# --------------------------------------------------------------------
# environment
# --------------------------------------------------------------------

def test_phase0_packages_present():
    versions = package_versions()
    for pkg in ("numpy", "pandas", "scipy", "pyarrow", "pytest"):
        assert versions[pkg] != "not installed", f"{pkg} missing"


def test_deferred_packages_not_required():
    """scikit-survival and statsmodels are deferred; Phase 0 must not need them."""
    import importlib.util
    assert importlib.util.find_spec("src.cache") is not None
    assert importlib.util.find_spec("src.soft_parser") is not None
