"""Phase 0: independently verify the probe matrix cache.

Re-parses randomly chosen arrays directly from the SOFT file and compares
them against the cache bit-for-bit, then checks structural invariants across
ALL 585 deposited samples.

This verifies COMPUTATIONAL fidelity -- that the cache is a faithful
container for the deposited values. It does not and cannot verify that those
values are biologically or normalizationally appropriate for the research
question; that is addressed in Phases 1 and 4.

Usage:
    .venv/bin/python scripts/14_verify_cache.py [--n 20]
"""

import argparse
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                       # noqa: E402
from src.cache import cache_is_current, load_manifest, load_probe_matrix  # noqa: E402
from src.reproducibility import DEFAULT_SEED                # noqa: E402
from src.soft_parser import iter_sample_tables, list_sample_ids, sha256  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=20,
                        help="Arrays to re-parse and compare exactly.")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()

    print("=" * 62)
    print("PHASE 0 -- CACHE VERIFICATION")
    print("=" * 62)

    failures: list[str] = []

    def check(label: str, passed: bool, detail: str = "") -> None:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {label}" + (f" -- {detail}" if detail else ""))
        if not passed:
            failures.append(label)

    ok, reason = cache_is_current()
    check("cache is current vs raw SHA-256", ok, reason)
    if not ok:
        print("\nCannot verify further. Rebuild the cache first.")
        return 1

    manifest = load_manifest()
    matrix = load_probe_matrix()

    print("\n-- structural checks (all deposited samples) --")
    soft_ids = list_sample_ids(paths.SOFT_SERIES)
    check("sample count == 585",
          len(soft_ids) == paths.SERIES_SAMPLE_COUNT, f"{len(soft_ids)}")
    check("cache has every deposited sample",
          list(matrix.columns) == soft_ids,
          f"{matrix.shape[1]} columns, order preserved")
    check("probe count == 54675",
          matrix.shape[0] == paths.GPL570_PROBE_COUNT, f"{matrix.shape[0]}")
    check("probe ids unique", matrix.index.is_unique)
    check("dtype is float64", str(matrix.dtypes.iloc[0]) == "float64")
    check("all values finite", bool(np.isfinite(matrix.to_numpy()).all()))
    check("no non-numeric fields encountered",
          manifest["non_numeric_fields_encountered"] == 0)

    print("\n-- probe identity and ordering across ALL 585 samples --")
    reference = list(matrix.index)
    mismatched, counted = [], 0
    for table in iter_sample_tables(paths.SOFT_SERIES, sample_ids=None):
        counted += 1
        if table.probe_ids != reference:
            mismatched.append(table.sample_id)
    check("every sample parsed", counted == paths.SERIES_SAMPLE_COUNT,
          f"{counted} samples")
    check("identical probe ids and ordering in all samples",
          not mismatched,
          "all 585 match" if not mismatched else f"{len(mismatched)} differ")

    print(f"\n-- exact value comparison ({args.n} random arrays) --")
    rng = random.Random(args.seed)
    chosen = set(rng.sample(soft_ids, min(args.n, len(soft_ids))))
    compared, exact = 0, 0
    for table in iter_sample_tables(paths.SOFT_SERIES, sample_ids=chosen):
        compared += 1
        fresh = np.asarray(table.values, dtype=np.float64)
        cached = matrix[table.sample_id].to_numpy()
        if np.array_equal(fresh, cached):
            exact += 1
        else:
            n_diff = int((fresh != cached).sum())
            print(f"    {table.sample_id}: {n_diff} differing values")
    check(f"all {compared} arrays bit-identical to source",
          compared == len(chosen) and exact == compared,
          f"{exact}/{compared} exact")

    print("\n-- raw data integrity --")
    actual = sha256(paths.SOFT_SERIES)
    check("series SOFT SHA-256 unchanged",
          actual == manifest["source_sha256"], f"{actual[:16]}...")

    print("\n" + "=" * 62)
    if failures:
        print(f"VERIFICATION FAILED -- {len(failures)} check(s):")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("VERIFICATION PASSED -- cache is computationally faithful.")
    print("Scope: computational fidelity only; see Phases 1 and 4 for")
    print("normalization and endpoint validity.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
