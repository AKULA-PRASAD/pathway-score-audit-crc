"""Phase 6D: freeze A-19 — missing-data strategy.

Records the approved decision: COMPLETE-CASE analysis for M1, M2 primary
(4-set), M2 sensitivity (5-set) and M3. Multiple imputation is NOT required
for the current model specifications.

No imputation, no modelling, no preprocessing. Phase 2 / Phase 3 / frozen M2
artifacts are read-only and unmodified.

Usage:
    .venv/bin/python scripts/29_freeze_a19_missing_data.py [--verify-only]
"""

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                  # noqa: E402
from src.reproducibility import environment_report     # noqa: E402

RECORD = paths.ROOT / "config" / "analysis" / "A19_FREEZE_RECORD.json"
EXCL = paths.ROOT / "reports" / "phase6d" / "complete_case_exclusions.json"
AUDIT = paths.ROOT / "reports" / "phase6d" / "definitive_audit.json"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def build() -> dict:
    excl = json.loads(EXCL.read_text())
    cc = pd.read_csv(paths.ROOT / "reports" / "phase6d" / "model_complete_case.csv")
    planned = cc[cc.model.str.contains("M1 primary|M2 primary|M2 sensitivity|M3 explor")]
    return {
        "record_type": "A19_MISSING_DATA_STRATEGY_FREEZE",
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "status": "FROZEN",
        "decision": {
            "strategy": "COMPLETE-CASE ANALYSIS",
            "applies_to": ["M1 primary (clinical only)",
                           "M2 primary (frozen 4-set pathway panel)",
                           "M2 sensitivity (frozen 5-set pathway panel)",
                           "M3 exploratory (transcriptome-wide)"],
            "multiple_imputation_required": False,
            "basis": "Phase 6D definitive missingness audit",
            "approved_by": "project owner",
        },
        "evidence": {
            "cohort": "S3", "n": 453,
            "complete_cases_per_planned_model": 452,
            "pct": 99.8, "excluded": 1,
            "per_model": planned[["model", "complete_cases", "pct", "lost"]]
                          .to_dict(orient="records"),
            "expression_nan_values": 0,
            "expression_zero_variance_genes": 0,
        },
        "excluded_observations": excl["excluded"],
        "prohibited_imputation": [
            "expression values", "pathway members", "endpoint fields",
            "sample identifiers", "cohort membership", "provenance fields",
        ],
        "m2_coverage_rules": {
            "unchanged": True,
            "rule": ">=80% of resolved members AND >=20 genes surviving "
                    "per-fold zero-variance exclusion",
            "missing_members": "EXCLUDED from the score and enumerated; never substituted",
        },
        "contingency": {
            "optional_covariates_excluded_from_M1": [
                "tnm.t", "tnm.n", "tnm.m", "mmr.status",
                "cimp.status", "cin.status", "kras.mutation",
                "braf.mutation", "tp53.mutation",
            ],
            "if_proposed_in_future_require": [
                "reassess the scientific rationale and missingness",
                "recalculate complete-case retention",
                "reconsider whether multiple imputation is justified",
                "account for assay/batch structure",
                "obtain explicit approval before changing the model specification",
            ],
            "mi_auto_activation": False,
            "note": "Multiple imputation must NOT be activated automatically.",
        },
        "deferred": {
            "pooling_rubin_vs_averaged_predictions": {
                "status": "DEFERRED — NOT APPLICABLE",
                "reason": "Pooling across imputations has no role under a "
                          "complete-case strategy. No pooling strategy is "
                          "specified, and none is invented.",
                "reopens_if": "multiple imputation is proposed AND approved "
                              "in a future analysis",
            }
        },
        "preserved_phase6d_findings": {
            "nd_na_distinction": "The deposit documents 'ND: not determined' "
                                 "and 'NA: Not applicable', but ND is used 0 "
                                 "times. Both meanings collapse to 'N/A'.",
            "structural_na": "Structural N/A values (e.g. chemotherapy type "
                             "when no chemotherapy was given; mutation detail "
                             "columns when WT) must NOT be treated as ordinary "
                             "missing observations.",
            "batch_association": "Molecular assay availability is strongly "
                                 "associated with processing batch: "
                                 "tp53.mutation missingness is 0% in four "
                                 "batches and 100% in three.",
            "mechanism_qualified": "MCAR is inconsistent with the evidence. "
                                   "MAR (conditional on batch) and MNAR/"
                                   "structural remain QUALIFIED ASSESSMENTS, "
                                   "not proven classifications.",
            "single_missing_chemotherapy": "One observation (GSM972262 / "
                                           "CIT306) has chemotherapy.adjuvant "
                                           "= 'N/A' and is the sole complete-"
                                           "case exclusion. It is documented, "
                                           "never silently discarded.",
        },
        "source_artifacts": {
            "definitive_audit": {"file": "reports/phase6d/definitive_audit.json",
                                 "sha256": sha256(AUDIT)},
            "exclusions": {"file": "reports/phase6d/complete_case_exclusions.json",
                           "sha256": sha256(EXCL)},
        },
        "upstream_frozen_hashes": {
            "m2_panel_4set.json": sha256(paths.ROOT / "config/pathways/m2_panel_4set.json"),
            "m2_panel_5set.json": sha256(paths.ROOT / "config/pathways/m2_panel_5set.json"),
        },
        "modification_policy": (
            "Immutable. A change of strategy requires a new approval, a "
            "recomputed complete-case assessment, and a re-freeze. Multiple "
            "imputation must never be activated automatically."
        ),
        "environment": environment_report(),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify-only", action="store_true")
    args = ap.parse_args()
    RECORD.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 76)
    print("A-19 MISSING-DATA STRATEGY FREEZE")
    print("=" * 76)
    if not args.verify_only:
        RECORD.write_text(json.dumps(build(), indent=2, default=str) + "\n")
        print(f"\n  record -> {RECORD.relative_to(paths.ROOT)}")

    r = json.loads(RECORD.read_text())
    print(f"\n  status   : {r['status']}")
    print(f"  frozen   : {r['frozen_utc']}")
    print(f"  strategy : {r['decision']['strategy']}")
    print(f"  MI required: {r['decision']['multiple_imputation_required']}")
    print("\n  applies to:")
    for m in r["decision"]["applies_to"]:
        print(f"    - {m}")
    print("\n  per-model complete cases:")
    for m in r["evidence"]["per_model"]:
        print(f"    {m['model']:40s} {m['complete_cases']}/453 "
              f"({m['pct']}%)  lost {m['lost']}")
    print("\n  EXCLUDED OBSERVATION(S) — documented, not silently discarded:")
    for e in r["excluded_observations"]:
        print(f"    {e['sample_id']} ({e['title']}) dataset={e['dataset']} "
              f"batch={e['combat_batch']}")
        print(f"      missing: {e['missing_variables']} = "
              f"{list(e['raw_values'].values())}")
    print(f"\n  pooling question: {r['deferred']['pooling_rubin_vs_averaged_predictions']['status']}")
    print(f"  MI auto-activation: {r['contingency']['mi_auto_activation']}")

    # The A-19 record attests to the M2 bytes Phase 6D read. It is a historical
    # attestation and is never rewritten, so a later APPROVED amendment makes a
    # bare byte comparison fail even though nothing is wrong. Resolve through
    # the amendment chain instead; undocumented drift still fails.
    from src.m2_amendment import resolve_hash

    ok = True
    print()
    for f in ("m2_panel_4set.json", "m2_panel_5set.json"):
        good, why = resolve_hash(paths.ROOT / "config" / "pathways" / f,
                                 r["upstream_frozen_hashes"][f])
        ok &= good
        print(f"  M2 spec {f:22s} {'OK  ' if good else 'FAIL'} {why}")
    print(f"  frozen M2 specs verified against the amendment chain: {ok}")
    print("\n" + "=" * 76)
    print("  A-19 IS FROZEN: COMPLETE-CASE. MI is NOT auto-activated.")
    print("=" * 76)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
