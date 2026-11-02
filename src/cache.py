"""Probe-level expression matrix cache.

Turns a ~60 second SOFT re-parse into a ~2 second Parquet load.

Design decisions and their reasons
----------------------------------
*All 585 deposited samples are cached, not only the 566 primary tumours.*
Restricting the cache would bake a cohort-eligibility decision into
infrastructure. Eligibility is a Phase 1 decision made in the open, so the
cache stays complete and filtering happens downstream.

*float64, not float32.* SOFT values carry 9-10 significant digits
(e.g. "9.561278312"); float32 holds about 7. float32 would halve the cache to
128 MB, but would make bit-exact round-trip verification impossible. The
exactness guarantee is worth 128 MB of disk.

*No value transformation.* Values are stored exactly as deposited. The cache
is a faster container for the same bytes, nothing more. It carries no claim
about whether those values are appropriate for any analysis.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from src import paths
from src.reproducibility import environment_report
from src.soft_parser import iter_sample_tables, sha256

CACHE_FORMAT_VERSION = 1


def build_cache(
    soft_path: Path | None = None,
    out_path: Path | None = None,
    manifest_path: Path | None = None,
    verbose: bool = True,
) -> dict:
    """Parse the SOFT series once and write the matrix plus its manifest.

    Returns the manifest dict. Raises on any probe-order or probe-count
    inconsistency rather than silently reconciling it.
    """
    soft_path = Path(soft_path or paths.SOFT_SERIES)
    out_path = Path(out_path or paths.PROBE_MATRIX)
    manifest_path = Path(manifest_path or paths.PROBE_MATRIX_MANIFEST)
    paths.ensure_processed_dir()

    reference_probes: list[str] | None = None
    sample_ids: list[str] = []
    columns: list[np.ndarray] = []
    non_numeric_total = 0

    for table in iter_sample_tables(soft_path, sample_ids=None):
        if reference_probes is None:
            reference_probes = list(table.probe_ids)
            if len(set(reference_probes)) != len(reference_probes):
                raise ValueError(
                    f"Duplicate probe ids within sample {table.sample_id}."
                )
        elif table.probe_ids != reference_probes:
            raise ValueError(
                f"Probe ids or ordering differ in sample {table.sample_id}; "
                "the series is not a rectangular matrix."
            )

        sample_ids.append(table.sample_id)
        columns.append(np.asarray(table.values, dtype=np.float64))
        non_numeric_total += table.non_numeric

        if verbose and len(sample_ids) % 100 == 0:
            print(f"  parsed {len(sample_ids)} samples...", flush=True)

    if reference_probes is None:
        raise ValueError(f"No sample tables found in {soft_path}.")

    matrix = pd.DataFrame(
        np.column_stack(columns),
        index=pd.Index(reference_probes, name="probe_id"),
        columns=sample_ids,
    )

    matrix.to_parquet(out_path, engine="pyarrow", compression="snappy")

    finite = bool(np.isfinite(matrix.to_numpy()).all())
    manifest = {
        "cache_format_version": CACHE_FORMAT_VERSION,
        "source_file": str(soft_path.relative_to(paths.ROOT)),
        "source_sha256": sha256(soft_path),
        "n_probes": int(matrix.shape[0]),
        "n_samples": int(matrix.shape[1]),
        "dtype": str(matrix.dtypes.iloc[0]),
        "all_values_finite": finite,
        "non_numeric_fields_encountered": non_numeric_total,
        "value_transformation": "none (values stored exactly as deposited)",
        "sample_scope": "all deposited samples; no cohort filter applied",
        "validity_note": (
            "Computational fidelity only. The deposited values were produced "
            "by RMA + ComBat applied jointly across the whole series; this "
            "cache makes no claim about their suitability for cross-cohort "
            "analysis."
        ),
        "environment": environment_report(),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    if verbose:
        print(f"  wrote {out_path.name}: {matrix.shape[0]} x {matrix.shape[1]}")
    return manifest


def load_manifest(manifest_path: Path | None = None) -> dict:
    """Read the cache manifest."""
    manifest_path = Path(manifest_path or paths.PROBE_MATRIX_MANIFEST)
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"No cache manifest at {manifest_path}. "
            "Run scripts/13_build_probe_cache.py first."
        )
    return json.loads(manifest_path.read_text())


def cache_is_current(soft_path: Path | None = None) -> tuple[bool, str]:
    """Check the cache against the live raw file.

    Returns ``(ok, reason)``. A stale cache is detected by SHA-256, not by
    timestamp, so touching the raw file does not produce a false alarm and
    editing it cannot go unnoticed.
    """
    soft_path = Path(soft_path or paths.SOFT_SERIES)
    if not paths.PROBE_MATRIX.exists():
        return False, "cache file missing"
    try:
        manifest = load_manifest()
    except FileNotFoundError as exc:
        return False, str(exc)

    actual = sha256(soft_path)
    if manifest.get("source_sha256") != actual:
        return False, (
            "source file changed since the cache was built "
            f"(manifest {manifest.get('source_sha256', '?')[:12]}..., "
            f"actual {actual[:12]}...)"
        )
    if manifest.get("cache_format_version") != CACHE_FORMAT_VERSION:
        return False, "cache format version mismatch"
    return True, "current"


def load_probe_matrix(validate: bool = True) -> pd.DataFrame:
    """Load the cached probe x sample matrix.

    With ``validate=True`` (default) a stale or missing cache raises rather
    than returning data that no longer matches the raw file.
    """
    if validate:
        ok, reason = cache_is_current()
        if not ok:
            raise RuntimeError(
                f"Probe matrix cache is not usable: {reason}. "
                "Rebuild with scripts/13_build_probe_cache.py."
            )
    return pd.read_parquet(paths.PROBE_MATRIX, engine="pyarrow")
