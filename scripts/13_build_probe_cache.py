"""Phase 0: build the probe-level expression matrix cache.

Parses data/raw/gse39582/GSE39582_family.soft.gz once and writes a Parquet
matrix plus a provenance manifest to data/processed/.

Caches ALL deposited samples. No cohort filter, no value transformation, no
probe filtering, no aggregation. Original files are opened read-only; the
only write target is data/processed/.

Usage:
    .venv/bin/python scripts/13_build_probe_cache.py [--force]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                    # noqa: E402
from src.cache import build_cache, cache_is_current      # noqa: E402
from src.soft_parser import series_field                 # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force", action="store_true",
        help="Rebuild even if the existing cache is current.",
    )
    args = parser.parse_args()

    print("=" * 62)
    print("PHASE 0 -- PROBE MATRIX CACHE BUILD")
    print("=" * 62)

    if not paths.SOFT_SERIES.exists():
        print(f"ERROR: source not found: {paths.SOFT_SERIES}")
        return 1

    ok, reason = cache_is_current()
    if ok and not args.force:
        print(f"Cache is already current ({reason}). Use --force to rebuild.")
        return 0
    if not ok:
        print(f"Rebuilding: {reason}")

    # Surface the deposited preprocessing provenance at build time so that
    # downstream users see it rather than discovering it later.
    processing = series_field(paths.SOFT_SERIES, "Sample_data_processing")
    distinct = sorted(set(processing))
    print("\nDeposited preprocessing (from the SOFT record):")
    for entry in distinct:
        print(f"  {entry}")
    print(
        "\nNOTE: this cache guarantees computational fidelity to the\n"
        "      deposited values only. It makes no claim that those values\n"
        "      are suitable for cross-cohort analysis -- see Phase 4."
    )

    print(f"\nParsing {paths.SOFT_SERIES.name} ...")
    manifest = build_cache(verbose=True)

    print("\n" + "=" * 62)
    print("BUILD SUMMARY")
    print("=" * 62)
    print(f"  probes           : {manifest['n_probes']}")
    print(f"  samples          : {manifest['n_samples']}")
    print(f"  dtype            : {manifest['dtype']}")
    print(f"  all finite       : {manifest['all_values_finite']}")
    print(f"  non-numeric seen : {manifest['non_numeric_fields_encountered']}")
    print(f"  source sha256    : {manifest['source_sha256'][:32]}...")
    print(f"\n  matrix   -> {paths.PROBE_MATRIX}")
    print(f"  manifest -> {paths.PROBE_MATRIX_MANIFEST}")

    expected_probes = paths.GPL570_PROBE_COUNT
    expected_samples = paths.SERIES_SAMPLE_COUNT
    problems = []
    if manifest["n_probes"] != expected_probes:
        problems.append(f"probe count {manifest['n_probes']} != {expected_probes}")
    if manifest["n_samples"] != expected_samples:
        problems.append(f"sample count {manifest['n_samples']} != {expected_samples}")
    if not manifest["all_values_finite"]:
        problems.append("non-finite values present")

    if problems:
        print("\nFAILED:")
        for p in problems:
            print(f"  - {p}")
        return 1

    print("\nOK. Original data were not modified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
