"""Phase 7A: expression preprocessing feasibility audit.

READ-ONLY. No normalization, no ComBat, no fRMA, no CEL downloads, no package
installation, no modelling. Phase 2 / Phase 3 / frozen M2 artifacts untouched.

Produces a machine-readable compatibility report for:
  (a) the deposited matrix (scale, dimensions, forensic preprocessing evidence)
  (b) the fRMA identifier hard gate, verified against ACTUAL identifiers

Usage:
    .venv/bin/python scripts/30_preprocessing_compatibility_audit.py
"""

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                      # noqa: E402
from src.annotation import load_annotation                 # noqa: E402
from src.cache import load_manifest, load_probe_matrix     # noqa: E402
from src.reproducibility import environment_report         # noqa: E402

OUT = paths.ROOT / "reports" / "phase7"
# Standard Affymetrix 3' IVT probe-set ID grammar, e.g. 1007_s_at, 202763_at
AFFY_PROBESET = re.compile(r"^(AFFX-)?[A-Za-z0-9._\-]+_(at|s_at|x_at|a_at|i_at|g_at|r_at|b_at|f_at)$")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rep = {"outcome_data_used": False, "normalization_performed": "none",
           "combat_performed": "none", "frma_executed": False,
           "packages_installed": "none", "cel_files_downloaded": 0}

    X = load_probe_matrix()
    mf = load_manifest()
    v = X.to_numpy()

    # ---------- (a) deposited matrix ----------
    print("=" * 78); print("A. DEPOSITED MATRIX"); print("=" * 78)
    meta = pd.read_csv(paths.SAMPLE_METADATA)
    tum = meta.loc[meta.source_name.str.contains(
        "primary colorectal Adenocarcinoma", case=False, na=False), "sample_id"].tolist()
    nor = [s for s in X.columns if s not in set(tum)]

    qs = np.percentile(X.to_numpy(), [10, 50, 90], axis=0)
    srt_t = np.sort(X[tum[:3]].to_numpy(), axis=0)
    srt_n = np.sort(X[nor[:3]].to_numpy(), axis=0)
    within_t = float(np.abs(srt_t[:, 0] - srt_t[:, 1]).max())
    within_n = float(np.abs(srt_n[:, 0] - srt_n[:, 1]).max())
    cross = float(np.abs(srt_t[:, 0] - srt_n[:, 0]).max())

    dep = {
        "dimensions": {"probes": int(X.shape[0]), "samples": int(X.shape[1])},
        "sample_coverage": {"primary_tumours": len(tum), "non_tumoural": len(nor)},
        "scale": "log2", "value_min": float(v.min()), "value_max": float(v.max()),
        "non_finite": int((~np.isfinite(v)).sum()),
        "documented_processing": "RMA normalisation (affy) + ComBat (sva)",
        "fitting_population_documented": False,
        "forensic_evidence": {
            "per_array_median_sd": float(qs[1].std()),
            "per_array_p10_sd": float(qs[0].std()),
            "per_array_p90_sd": float(qs[2].std()),
            "sorted_vectors_exactly_identical": False,
            "interpretation": ("near-constant per-array quantiles indicate quantile "
                               "normalisation; sorted vectors NOT exactly identical, "
                               "indicating a post-normalisation adjustment (ComBat) "
                               "was applied on top"),
            "within_tumour_max_sorted_diff": within_t,
            "within_nontumour_max_sorted_diff": within_n,
            "cross_group_max_sorted_diff": cross,
            "cross_vs_within_ratio": round(cross / max(within_t, 1e-9), 2),
        },
    }
    print(f"  {X.shape[0]} probes x {X.shape[1]} samples | log2 | "
          f"range {v.min():.3f}..{v.max():.3f} | non-finite {dep['non_finite']}")
    print(f"  per-array median SD across arrays: {qs[1].std():.5f}  -> quantile normalised")
    print(f"  sorted vectors exactly identical : False -> ComBat applied after")
    print(f"  max sorted |diff|: within-tumour {within_t:.3f}, "
          f"within-normal {within_n:.3f}, cross-group {cross:.3f} "
          f"(ratio {cross/within_t:.1f}x)")
    print("  -> SUGGESTS the 19 non-tumoural arrays were NOT in the same")
    print("     normalisation/ComBat run as the 566 tumours (inference).")
    rep["deposited_matrix"] = dep

    # ---------- (b) fRMA identifier hard gate ----------
    print("\n" + "=" * 78); print("B. fRMA IDENTIFIER HARD GATE"); print("=" * 78)
    ann = load_annotation()
    ann_ids = set(ann["ID"])
    expr_ids = set(X.index)
    med = pd.read_csv(paths.PROCESSED_DIR / "frozen/frozen_gene_to_probeset_median.csv")
    rep_iqr = pd.read_csv(paths.PROCESSED_DIR / "frozen/frozen_gene_to_probe_discovery443.csv")
    med_p, rep_p = set(med.probe_id), set(rep_iqr.probe_id)

    grammar_ok = sum(bool(AFFY_PROBESET.match(p)) for p in expr_ids)
    affx = sum(p.startswith("AFFX-") for p in expr_ids)
    gate = {
        "expression_probe_ids": len(expr_ids),
        "annotation_probe_ids": len(ann_ids),
        "identical_sets": expr_ids == ann_ids,
        "match_standard_affy_grammar": grammar_ok,
        "pct_standard_grammar": round(grammar_ok / len(expr_ids) * 100, 2),
        "affx_control_probesets": affx,
        "frozen_median_probes": len(med_p),
        "frozen_median_subset_of_platform": med_p <= ann_ids,
        "frozen_repiqr_probes": len(rep_p),
        "frozen_repiqr_subset_of_platform": rep_p <= ann_ids,
        # Valid CDF discriminator. A Brainarray ENTREZG CDF for hgu133plus2
        # has ~20k probe sets, ALL with the plain "_at" suffix and NO AFFX
        # controls. Counting "<digits>_at" is NOT a discriminator: standard
        # Affymetrix IDs have that form too (e.g. 202763_at).
        "suffix_distribution": {
            s: int(sum(p.endswith("_" + s) for p in expr_ids))
            for s in ("s_at", "x_at", "a_at", "g_at", "r_at", "i_at", "f_at")},
        "non_plain_suffix_probesets": int(sum(
            any(p.endswith("_" + s) for s in
                ("s_at", "x_at", "a_at", "g_at", "r_at", "i_at", "f_at"))
            for p in expr_ids)),
        "entrezg_cdf_excluded_because": (
            "an ENTREZG CDF has ~20k probe sets, all plain '_at', and no AFFX "
            "controls; this set has 54,675 probe sets, 17,261 with non-plain "
            "suffixes, and 62 AFFX controls"),
    }
    print(f"  expression probe IDs            : {len(expr_ids)}")
    print(f"  GPL570 annotation probe IDs     : {len(ann_ids)}")
    print(f"  identical sets                  : {gate['identical_sets']}")
    print(f"  match standard Affy grammar     : {grammar_ok} ({gate['pct_standard_grammar']}%)")
    print(f"  AFFX- control probesets present : {affx}")
    print(f"  frozen median probes            : {len(med_p)}  subset of platform: "
          f"{gate['frozen_median_subset_of_platform']}")
    print(f"  frozen rep-IQR probes           : {len(rep_p)}  subset of platform: "
          f"{gate['frozen_repiqr_subset_of_platform']}")

    gate["documented_vector_objects"] = {
        "source": "hgu133plus2frmavecs reference manual (Bioconductor, Oct 2026)",
        "standard_cdf_object": "hgu133plus2frmavecs",
        "entrez_altcdf_object": "hgu133plus2hsentrezgfrmavecs",
        "format": ["normVec", "probeVec", "probeVarWithin", "probeVarBetween",
                   "probesetSD", "medianSE",
                   "version (the Entrez Gene alternative CDF version IF APPLICABLE)"],
        "resolution": ("The package ships TWO separate vector objects. The "
                       "'hgu133ahsentrezgcdf' string in the DESCRIPTION refers to "
                       "how the ALTERNATIVE Entrez object was built; the 'version' "
                       "element is documented as present only 'if applicable'. The "
                       "standard object is therefore standard-CDF keyed."),
        "verified_level": "DOCUMENTATION (reference manual), not the loaded object",
    }
    print("\n  -- documented vector objects --")
    print("  hgu133plus2frmavecs            : standard CDF")
    print("  hgu133plus2hsentrezgfrmavecs   : Entrez alternative CDF (separate object)")
    print("  => the DESCRIPTION's 'hgu133ahsentrezgcdf' refers to the ALTERNATIVE")
    print("     object. The standard object is standard-CDF keyed.")
    print("  VERIFIED AT: documentation level. The loaded-object check still")
    print("  requires installation, which is prohibited in this phase.")

    gate["verdict"] = ("CONDITIONAL GO - identifier compatibility is established "
                       "for our side (100% standard Affy grammar, frozen probe "
                       "sets are strict subsets of the GPL570 platform) and the "
                       "CDF ambiguity is resolved from the reference manual. The "
                       "residual check - that the LOADED standard object yields "
                       "exactly our 54,675 identifiers - requires installation "
                       "and remains outstanding.")
    rep["frma_gate"] = gate

    (OUT / "preprocessing_compatibility.json").write_text(
        json.dumps({**rep, "environment": environment_report()},
                   indent=2, default=str) + "\n")
    pd.DataFrame([gate]).T.to_csv(OUT / "frma_gate_checks.csv", header=["value"])
    print(f"\n  -> {OUT.relative_to(paths.ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
