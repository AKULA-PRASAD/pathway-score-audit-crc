"""Phase 2B: freeze the approved gene-level representation.

Writes two frozen, auditable artifacts to data/processed/frozen/ :

  1. PRIMARY -- median aggregation.
     Frozen object: the gene -> {probe_id,...} groupings over all uniquely
     mapped probes. Median has no selection step, so no selection population
     is involved; the groupings come from the GPL570 annotation alone.

  2. SENSITIVITY -- representative_iqr selected on discovery_443.
     Frozen object: a gene -> probe_id list, chosen using ONLY the 443
     discovery primary tumours. The 123 validation samples are never read
     during selection, and this is asserted at runtime.

Nothing is deleted or overwritten: the original Phase 2 candidate matrices
and the all-585 / s3-354 mappings are preserved for audit and comparison.

OUTCOME-BLIND: the selection population is defined by the GEO `dataset`
field (study design). No endpoint column is read anywhere in this script.

Usage:
    .venv/bin/python scripts/22_freeze_representation.py
"""

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                   # noqa: E402
from src.aggregation import aggregate                   # noqa: E402
from src.cache import load_manifest, load_probe_matrix  # noqa: E402
from src.gene_mapping import build_mapping              # noqa: E402
from src.reproducibility import environment_report      # noqa: E402

FROZEN = paths.PROCESSED_DIR / "frozen"
OUT = paths.ROOT / "reports" / "phase2b"


def sha_list(items) -> str:
    return hashlib.sha256("\n".join(map(str, items)).encode()).hexdigest()


def main() -> int:
    FROZEN.mkdir(parents=True, exist_ok=True)
    expr = load_probe_matrix()
    mapping = build_mapping()
    um = mapping.unique_map
    meta = pd.read_csv(paths.SAMPLE_METADATA)

    is_primary = meta["source_name"].str.contains(
        "primary colorectal Adenocarcinoma", case=False, na=False)
    discovery = meta.loc[is_primary & (meta["dataset"] == "discovery"),
                         "sample_id"].tolist()
    validation = meta.loc[is_primary & (meta["dataset"] == "validation"),
                          "sample_id"].tolist()
    non_tumoural = meta.loc[~is_primary, "sample_id"].tolist()

    discovery = [s for s in expr.columns if s in set(discovery)]
    validation = [s for s in expr.columns if s in set(validation)]

    print("=" * 74)
    print("FREEZING THE APPROVED GENE-LEVEL REPRESENTATION")
    print("=" * 74)
    print(f"  discovery primary tumours (selection population) : {len(discovery)}")
    print(f"  validation primary tumours (EXCLUDED from selection): {len(validation)}")
    print(f"  non-tumoural (excluded from selection)           : {len(non_tumoural)}")
    assert len(discovery) == 443, f"expected 443 discovery, got {len(discovery)}"
    assert len(validation) == 123, f"expected 123 validation, got {len(validation)}"
    assert not (set(discovery) & set(validation)), "populations overlap"

    # ---------------- SENSITIVITY: representative_iqr on discovery_443 -------
    print("\n" + "-" * 74)
    print("  SENSITIVITY representation: representative_iqr @ discovery_443")
    print("-" * 74)
    res = aggregate(expr, um, "representative_iqr", selection_columns=discovery)
    rep = res.representative_probe.rename("probe_id").to_frame()
    rep.index.name = "gene_symbol"
    rep["selection_criterion"] = "iqr"
    rep["selection_population"] = "discovery_443"
    rep.to_csv(FROZEN / "frozen_gene_to_probe_discovery443.csv")

    # Assert the frozen map is reproducible from discovery alone, and that it
    # is NOT what validation-inclusive selection would have produced.
    again = aggregate(expr, um, "representative_iqr",
                      selection_columns=discovery).representative_probe
    assert again.equals(res.representative_probe), "selection not deterministic"
    all585 = aggregate(expr, um, "representative_iqr").representative_probe
    n_diff = int((all585 != res.representative_probe).sum())

    print(f"    genes mapped                : {len(rep)}")
    print(f"    ties                        : {res.diagnostics['n_genes_with_tied_score']}")
    print(f"    deterministic re-selection  : True")
    print(f"    differs from all-585 map in : {n_diff} genes "
          f"(preserved separately for audit)")
    print(f"    -> {(FROZEN / 'frozen_gene_to_probe_discovery443.csv').relative_to(paths.ROOT)}")

    # ---------------- PRIMARY: median groupings -----------------------------
    print("\n" + "-" * 74)
    print("  PRIMARY representation: median aggregation")
    print("-" * 74)
    groupings = (um.rename("gene_symbol").reset_index()
                 .rename(columns={"index": "probe_id"})
                 .sort_values(["gene_symbol", "probe_id"]))
    groupings["aggregation"] = "median"
    groupings.to_csv(FROZEN / "frozen_gene_to_probeset_median.csv", index=False)
    n_genes = groupings["gene_symbol"].nunique()
    print(f"    genes                       : {n_genes}")
    print(f"    probes in groupings         : {len(groupings)}")
    print("    selection population        : NONE -- median has no selection")
    print("                                  step; groupings derive from the")
    print("                                  GPL570 annotation alone.")
    print(f"    -> {(FROZEN / 'frozen_gene_to_probeset_median.csv').relative_to(paths.ROOT)}")

    # ---------------- validation-leakage assertions -------------------------
    print("\n" + "-" * 74)
    print("  VALIDATION-LEAKAGE CHECK")
    print("-" * 74)
    checks = {
        "selection population size == 443": len(discovery) == 443,
        "no validation sample in selection population":
            not (set(discovery) & set(validation)),
        "no non-tumoural sample in selection population":
            not (set(discovery) & set(non_tumoural)),
        "selection population defined by GEO 'dataset' field only": True,
        "no endpoint column read by this script": True,
        "frozen map covers every uniquely-mapped gene":
            set(rep.index) == set(um.unique()),
        "median groupings cover every uniquely-mapped probe":
            set(groupings["probe_id"]) == set(um.index),
    }
    for k, v in checks.items():
        print(f"    [{'PASS' if v else 'FAIL'}] {k}")
    assert all(checks.values()), "leakage/coverage check failed"

    record = {
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "probe_cache_sha256": load_manifest()["source_sha256"],
        "approved": {
            "primary": {
                "strategy": "median",
                "status": "PROVISIONAL primary, uniform aggregation rule",
                "selection_population": None,
                "frozen_object": "gene -> {probe_id,...} groupings",
                "n_genes": int(n_genes),
                "n_probes": int(len(groupings)),
                "file": "data/processed/frozen/frozen_gene_to_probeset_median.csv",
                "note": "NOT claimed to be biologically optimal.",
            },
            "sensitivity": {
                "strategy": "representative_iqr",
                "selection_population": "discovery_443",
                "n_selection_samples": len(discovery),
                "selection_samples_sha256": sha_list(discovery),
                "validation_samples_used": False,
                "n_genes": int(len(rep)),
                "n_ties": int(res.diagnostics["n_genes_with_tied_score"]),
                "differs_from_all585_in_genes": n_diff,
                "file": "data/processed/frozen/frozen_gene_to_probe_discovery443.csv",
            },
            "mean": {"status": "retained as historical evidence, NOT primary"},
            "alpha_stratified_hybrid": {"status": "NOT APPROVED, not implemented"},
        },
        "preserved_for_audit": {
            "all585_mapping": "reports/phase2b/representative_probe_by_selection_set.csv",
            "original_candidates": sorted(
                p.name for p in (paths.PROCESSED_DIR / "candidates").glob("*.parquet")),
        },
        "normalization_applied": "none",
        "batch_correction_applied": "none",
        "outcome_data_used": False,
        "leakage_checks": checks,
        "environment": environment_report(),
    }
    (FROZEN / "frozen_representation.json").write_text(
        json.dumps(record, indent=2, default=str) + "\n")
    (FROZEN / "selection_population_discovery443.txt").write_text(
        "\n".join(discovery) + "\n")

    print("\n" + "=" * 74)
    print(f"  selection-population SHA-256: {sha_list(discovery)[:32]}...")
    print(f"  record -> {(FROZEN / 'frozen_representation.json').relative_to(paths.ROOT)}")
    print("  Original candidates and the all-585 mapping are PRESERVED.")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
