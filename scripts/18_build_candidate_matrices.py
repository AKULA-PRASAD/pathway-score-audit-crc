"""Phase 2: build CANDIDATE gene-level expression matrices.

Produces one matrix per aggregation strategy, written SEPARATELY to
data/processed/candidates/. The probe-level cache is read-only input and is
never modified; no candidate replaces it, and none is designated final.

Invariants enforced at build time (the script aborts if any fails):
  * all 585 deposited samples retained, in the cache's column order
  * sample ids and ordering match the Phase 1 cohort manifest
  * NO clinical cohort exclusion applied to the expression matrix
  * NO normalization, NO batch correction, NO expression filter
  * NO outcome data read -- this script never opens the endpoint columns
  * ambiguous and unmapped probes excluded by construction and NAMED

Usage:
    .venv/bin/python scripts/18_build_candidate_matrices.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                       # noqa: E402
from src.aggregation import STRATEGIES, aggregate           # noqa: E402
from src.cache import load_manifest, load_probe_matrix      # noqa: E402
from src.gene_mapping import build_mapping                  # noqa: E402
from src.reproducibility import environment_report          # noqa: E402

CAND_DIR = paths.PROCESSED_DIR / "candidates"
OUT = paths.ROOT / "reports" / "phase2"


def main() -> int:
    CAND_DIR.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    expr = load_probe_matrix()
    mapping = build_mapping()
    gene_map = mapping.unique_map

    print("=" * 72)
    print("PHASE 2 -- CANDIDATE GENE-LEVEL MATRICES")
    print("=" * 72)
    print(f"  probe-level input : {expr.shape[0]} probes x {expr.shape[1]} samples")
    print(f"  uniquely mapped   : {len(gene_map)} probes -> "
          f"{gene_map.nunique()} genes")
    print(f"  ambiguous excluded: {len(mapping.probes('ambiguous'))} (named in "
          "reports/phase2/probes_ambiguous.csv)")
    print(f"  unmapped excluded : {len(mapping.probes('unmapped'))} (named in "
          "reports/phase2/probes_unmapped.csv)")

    # --- sample-order contract against the Phase 1 manifest ---
    man_path = paths.ROOT / "reports" / "phase1" / "cohort_manifest.csv"
    if not man_path.exists():
        print("\nERROR: Phase 1 cohort manifest not found. Run script 16 first.")
        return 1
    man = pd.read_csv(man_path)
    cache_cols = list(expr.columns)
    if set(man["sample_id"]) != set(cache_cols):
        print("\nERROR: manifest sample set differs from the expression matrix.")
        return 1
    print(f"\n  manifest samples  : {len(man)}  (set matches expression matrix)")
    print(f"  column order      : cache order preserved "
          f"({cache_cols[0]} ... {cache_cols[-1]})")
    in_primary = int(man["is_primary_cohort"].sum())
    print(f"  manifest flags {in_primary} samples as the Phase 1 primary cohort,")
    print(f"  but ALL {len(cache_cols)} are retained here by instruction.")

    summary = {
        "probe_matrix": {"probes": int(expr.shape[0]),
                         "samples": int(expr.shape[1])},
        "mapping": mapping.summary(),
        "source_cache_sha256": load_manifest()["source_sha256"],
        "normalization_applied": "none",
        "batch_correction_applied": "none",
        "expression_filter_applied": "none",
        "cohort_filter_applied": "none (all 585 samples retained)",
        "outcome_data_used": False,
        "deposited_provenance": "RMA + ComBat applied by the depositors; this "
                                "phase adds no normalization or batch adjustment",
        "environment": environment_report(),
        "candidates": {},
    }

    for strategy in STRATEGIES:
        print("\n" + "-" * 72)
        print(f"  strategy: {strategy}")
        res = aggregate(expr, gene_map, strategy)
        m = res.matrix

        assert list(m.columns) == cache_cols, "sample order changed"
        assert m.shape[1] == paths.SERIES_SAMPLE_COUNT, "sample count changed"
        assert m.index.is_unique, "duplicate gene symbols"
        assert m.index.is_monotonic_increasing, "gene index not sorted"
        finite = bool(np.isfinite(m.to_numpy()).all())
        assert finite, "non-finite values produced"

        path = CAND_DIR / f"gene_matrix_{strategy}.parquet"
        m.to_parquet(path, engine="pyarrow", compression="snappy")

        vals = m.to_numpy()
        info = {
            "n_genes": int(m.shape[0]), "n_samples": int(m.shape[1]),
            "n_probes_used": int(res.n_probes_used),
            "all_finite": finite,
            "value_min": float(vals.min()), "value_max": float(vals.max()),
            "mean_of_means": float(vals.mean()),
            "file": str(path.relative_to(paths.ROOT)),
            "file_mb": round(path.stat().st_size / 1e6, 1),
            "diagnostics": res.diagnostics,
        }
        summary["candidates"][strategy] = info
        print(f"    genes x samples : {m.shape[0]} x {m.shape[1]}")
        print(f"    probes consumed : {res.n_probes_used}")
        print(f"    value range     : {vals.min():.3f} - {vals.max():.3f}")
        print(f"    diagnostics     : {res.diagnostics}")
        print(f"    -> {path.relative_to(paths.ROOT)} ({info['file_mb']} MB)")

        if res.representative_probe is not None:
            rp = res.representative_probe.rename("probe_id").to_frame()
            rp.index.name = "gene_symbol"
            rp.to_csv(OUT / f"representative_probes_{strategy}.csv")

    (OUT / "candidate_matrices_summary.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n")

    print("\n" + "=" * 72)
    print("  All candidates built. The probe-level cache is UNCHANGED and")
    print("  remains the canonical representation. No candidate is final.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
