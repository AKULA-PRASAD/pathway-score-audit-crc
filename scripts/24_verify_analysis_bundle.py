"""Phase 3: standalone analysis-bundle verification (no pytest required).

Reports the bundle's state and checks every recorded invariant. Read-only.
Exit 0 if all checks pass, 1 otherwise.

Usage:
    .venv/bin/python scripts/24_verify_analysis_bundle.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                   # noqa: E402
from src.bundle import ANALYSIS_DIR, sha256_file        # noqa: E402
from src.cohort import COHORTS                          # noqa: E402
from src.provenance import FORBIDDEN_PROVENANCE_PHRASES  # noqa: E402

EXPECTED = {"S3": (453, 136, 130), "S0": (461, 139, 133),
            "S1": (458, 136, 130), "S2": (456, 139, 133)}


def main() -> int:
    mf_path = ANALYSIS_DIR / "analysis_bundle_manifest.json"
    if not mf_path.exists():
        print("ERROR: no bundle. Run scripts/23_build_analysis_bundle.py.")
        return 1
    mf = json.loads(mf_path.read_text())
    fails: list[str] = []

    def check(label: str, ok: bool, detail: str = "") -> None:
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}"
              + (f" -- {detail}" if detail else ""))
        if not ok:
            fails.append(label)

    print("=" * 74); print("PHASE 3 BUNDLE VERIFICATION"); print("=" * 74)

    print("\n-- input hashes --")
    for k, h in mf["input_hashes"].items():
        p = {"probe_matrix": paths.PROBE_MATRIX,
             "frozen_median_spec": paths.PROCESSED_DIR / "frozen/frozen_gene_to_probeset_median.csv",
             "frozen_repiqr_spec": paths.PROCESSED_DIR / "frozen/frozen_gene_to_probe_discovery443.csv",
             "selection_population": paths.PROCESSED_DIR / "frozen/selection_population_discovery443.txt",
             "cohort_manifest": paths.ROOT / "reports/phase1/cohort_manifest.csv",
             "sample_metadata": paths.SAMPLE_METADATA,
             "sdrf": paths.ROOT / "data/external/E-GEOD-39582.sdrf.txt"}[k]
        check(f"input unchanged: {k}", sha256_file(p) == h)

    print("\n-- output hashes --")
    for name, rec in mf["outputs"].items():
        p = ANALYSIS_DIR / name
        check(f"output intact: {name}",
              p.exists() and sha256_file(p) == rec["sha256"])

    print("\n-- cohort counts --")
    memb = pd.read_csv(ANALYSIS_DIR / "cohort_membership.csv")
    ep = pd.read_csv(ANALYSIS_DIR / "endpoint.csv")
    for c in COHORTS:
        ids = set(memb.loc[memb[f"in_{c}"].astype(bool), "sample_id"])
        e = ep[ep.sample_id.isin(ids)]
        got = (len(e), int((e.rfs_event_full == 1).sum()),
               int((e.rfs_event_60units == 1).sum()))
        check(f"{c} counts", got == EXPECTED[c],
              f"n/ev/ev60 = {got} vs expected {EXPECTED[c]}")

    print("\n-- matrices --")
    for f in ("expression_median_S0.parquet",
              "expression_repiqr_disc443_S0.parquet"):
        m = pd.read_parquet(ANALYSIS_DIR / f)
        check(f"{f} shape/dtype/finite",
              m.shape == (21815, 461) and str(m.dtypes.iloc[0]) == "float64"
              and bool(np.isfinite(m.to_numpy()).all()),
              f"{m.shape}")

    print("\n-- no transformation fitted, no correction applied --")
    check("normalization_applied == none",
          mf["normalization_applied"] == "none")
    check("batch_correction_applied == none",
          mf["batch_correction_applied"] == "none")
    check("batch covariate not auto-included",
          mf["batch_covariate_auto_included"] is False)
    check("no transformation fitted",
          mf["transformation_fitted_in_phase3"] == "none")

    print("\n-- provenance wording --")
    prov = mf["processing_provenance"].lower()
    check("states the fitting population is undocumented",
          "does not document the fitting population" in prov)
    for bad in FORBIDDEN_PROVENANCE_PHRASES:
        check(f"no forbidden phrase: {bad!r}", bad not in prov)

    print("\n-- leakage --")
    check("validation samples not used in probe selection",
          mf["representations"]["sensitivity"]
          ["validation_samples_used_in_selection"] is False)

    print("\n" + "=" * 74)
    if fails:
        print(f"VERIFICATION FAILED -- {len(fails)} check(s)")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("VERIFICATION PASSED")
    print("Scope: computational correctness and leakage prevention only.")
    print("Biological validity is NOT established; processing provenance is")
    print("incompletely documented and no normalization has been applied.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
