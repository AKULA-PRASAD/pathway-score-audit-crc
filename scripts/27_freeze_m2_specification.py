"""Phase 6C: freeze the M2 pathway specification.

Marks both panel specifications FROZEN, records the approval decision, and
writes an immutable freeze record with SHA-256 of each specification.

OUTCOME-BLIND: this script reads no recurrence, survival or outcome column.
The freeze occurs BEFORE any outcome-association analysis has been performed.

Does not modify Phase 2 frozen representations or Phase 3 analysis data.

Usage:
    .venv/bin/python scripts/27_freeze_m2_specification.py [--verify-only]
"""

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                  # noqa: E402
from src.reproducibility import environment_report     # noqa: E402

PANEL_DIR = paths.ROOT / "config" / "pathways"
PRIMARY = PANEL_DIR / "m2_panel_4set.json"
SENSITIVITY = PANEL_DIR / "m2_panel_5set.json"
RECORD = PANEL_DIR / "M2_FREEZE_RECORD.json"

APPROVAL = {
    "approved_by": "project owner",
    "approval_basis": "Phase 6C rev. 2 specification audit",
    "primary_panel": "4-set",
    "sensitivity_panel": "5-set",
    "sensitivity_only_addition": "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "immutability_clause": (
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING is a SENSITIVITY-ONLY addition. It "
        "must NOT be promoted to the primary panel, and must NOT be removed "
        "from the sensitivity panel, on the basis of observed model "
        "performance. Any such change would convert a pre-specified analysis "
        "into a post-hoc one and reclassifies M2 as M3-exploratory."
    ),
    "outcome_blind": True,
    "outcome_blind_statement": (
        "Panel membership and the primary/sensitivity assignment were decided "
        "using gene-set metadata, published colorectal cancer biology, and the "
        "outcome-blind expression correlation structure. No recurrence, "
        "survival or outcome variable was accessed at any point."
    ),
    "frozen_before_outcome_analysis": True,
    "selection_criteria_used": [
        "effective dimensionality (identical: 3 PCs for 80%, Kaiser 2)",
        "score-matrix conditioning (2.95 vs 6.95)",
        "smallest eigenvalue (0.480 vs 0.280)",
        "count of |r| > 0.5 pairs (0 vs 1)",
        "parameter budget against 136 events",
        "published recurrence-prediction evidence strength",
    ],
    "selection_criteria_NOT_used": [
        "any recurrence or survival association",
        "any model fit or performance metric",
        "any outcome column whatsoever",
    ],
}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def mark_frozen(path: Path, role: str, stamp: str) -> None:
    d = json.loads(path.read_text())
    d["status"] = "FROZEN"
    d["role"] = role
    d["frozen_utc"] = stamp
    d["approval"] = APPROVAL
    d["immutable"] = True
    d["modification_policy"] = (
        "Immutable. Permitted revisions are limited to non-outcome corrections "
        "(e.g. a database version correction), which must be recorded with "
        "reason and timestamp, retain the prior version, and re-freeze. Any "
        "revision informed by outcome data reclassifies M2 as M3-exploratory, "
        "permanently and irreversibly."
    )
    path.write_text(json.dumps(d, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify-only", action="store_true")
    args = ap.parse_args()

    print("=" * 76)
    print("M2 PATHWAY SPECIFICATION FREEZE")
    print("=" * 76)

    if not args.verify_only:
        stamp = datetime.now(timezone.utc).isoformat()
        mark_frozen(PRIMARY, "PRIMARY", stamp)
        mark_frozen(SENSITIVITY, "SENSITIVITY (pre-specified)", stamp)

        prim = json.loads(PRIMARY.read_text())
        sens = json.loads(SENSITIVITY.read_text())
        record = {
            "record_type": "M2_PATHWAY_SPECIFICATION_FREEZE",
            "frozen_utc": stamp,
            "approval": APPROVAL,
            "specifications": {
                "primary": {
                    "file": str(PRIMARY.relative_to(paths.ROOT)),
                    "panel": prim["panel"], "n_sets": prim["n_sets"],
                    "gene_sets": sorted(prim["gene_sets"]),
                    "union_unique_genes": prim["union_unique_genes"],
                    "sha256": sha256(PRIMARY),
                },
                "sensitivity": {
                    "file": str(SENSITIVITY.relative_to(paths.ROOT)),
                    "panel": sens["panel"], "n_sets": sens["n_sets"],
                    "gene_sets": sorted(sens["gene_sets"]),
                    "union_unique_genes": sens["union_unique_genes"],
                    "sha256": sha256(SENSITIVITY),
                },
            },
            "source_database": prim["database"],
            "mapping_provenance": prim["mapping"],
            "scoring": prim["scoring"],
            "target_index": prim["target_index"],
            "upstream_hashes": {
                "frozen_gene_to_probeset_median.csv": sha256(
                    paths.PROCESSED_DIR / "frozen/frozen_gene_to_probeset_median.csv"),
            },
            "environment": environment_report(),
        }
        RECORD.write_text(json.dumps(record, indent=2, default=str) + "\n")
        print(f"\n  freeze record -> {RECORD.relative_to(paths.ROOT)}")

    rec = json.loads(RECORD.read_text())
    print(f"\n  frozen at : {rec['frozen_utc']}")
    print(f"  outcome-blind : {rec['approval']['outcome_blind']}")
    print(f"  frozen before outcome analysis : "
          f"{rec['approval']['frozen_before_outcome_analysis']}")
    for role in ("primary", "sensitivity"):
        s = rec["specifications"][role]
        print(f"\n  {role.upper()}: {s['panel']}  ({s['n_sets']} sets, "
              f"{s['union_unique_genes']} union genes)")
        print(f"    sha256 {s['sha256'][:32]}...")
        for g in s["gene_sets"]:
            print(f"      - {g}")
    print(f"\n  GMT checksums preserved:")
    print(f"    symbols {rec['source_database']['sha256_symbols'][:32]}...")
    print(f"    entrez  {rec['source_database']['sha256_entrez'][:32]}...")

    ok = all(sha256(paths.ROOT / rec["specifications"][r]["file"])
             == rec["specifications"][r]["sha256"] for r in ("primary", "sensitivity"))
    print(f"\n  specification files match the record: {ok}")
    print("\n" + "=" * 76)
    print("  M2 IS FROZEN. PI3K/AKT/mTOR is sensitivity-only and must not be")
    print("  promoted or removed on the basis of observed model performance.")
    print("=" * 76)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
