"""Phase 3: build the analysis bundle.

Assembles the approved expression representations, cohort membership,
endpoint, covariates, partition and batch labels into a verified,
hash-recorded bundle at data/processed/analysis/.

Performs NO survival modelling, NO feature selection, NO normalization and
NO batch correction. Fits nothing.

Usage:
    .venv/bin/python scripts/23_build_analysis_bundle.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                       # noqa: E402
from src.aggregation import aggregate                       # noqa: E402
from src.bundle import (ANALYSIS_DIR, EXCLUDED_COVARIATES, align_by_id,
                        base_manifest, batch_composition,
                        build_covariates, build_partition_and_batch,
                        write_manifest)                      # noqa: E402
from src.cache import load_probe_matrix                      # noqa: E402
from src.cohort import (COHORTS, PRIMARY_COHORT, exclusion_table,
                        load_manifest, materialise)          # noqa: E402
from src.endpoint import build_endpoint, endpoint_metadata   # noqa: E402
from src.external_annotation import load_batch_assignments   # noqa: E402
from src.gene_mapping import build_mapping                   # noqa: E402

FROZEN = paths.PROCESSED_DIR / "frozen"
CAND = paths.PROCESSED_DIR / "candidates"

EXPECTED = {"S3": (453, 136, 130), "S0": (461, 139, 133),
            "S1": (458, 136, 130), "S2": (456, 139, 133)}


def section(t: str) -> None:
    print("\n" + "=" * 74); print(t); print("=" * 74)


def main() -> int:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    # ---------------- 3.1 verify frozen representations ----------------
    section("3.1  VERIFY FROZEN REPRESENTATIONS")
    expr = load_probe_matrix()
    um = build_mapping().unique_map
    disc443 = (FROZEN / "selection_population_discovery443.txt"
               ).read_text().split()
    print(f"  selection population: {len(disc443)} discovery samples")

    med_existing = pd.read_parquet(CAND / "gene_matrix_median.parquet")
    med_rebuilt = aggregate(expr, um, "median").matrix
    ok_med = np.array_equal(med_existing.to_numpy(), med_rebuilt.to_numpy()) \
        and list(med_existing.index) == list(med_rebuilt.index) \
        and list(med_existing.columns) == list(med_rebuilt.columns)
    print(f"  median        rebuild bit-identical: {ok_med}")

    rep_existing = pd.read_parquet(
        CAND / "gene_matrix_representative_iqr_discovery_443.parquet")
    rep_rebuilt = aggregate(expr, um, "representative_iqr",
                            selection_columns=disc443).matrix
    ok_rep = np.array_equal(rep_existing.to_numpy(), rep_rebuilt.to_numpy()) \
        and list(rep_existing.index) == list(rep_rebuilt.index) \
        and list(rep_existing.columns) == list(rep_rebuilt.columns)
    print(f"  rep-IQR @443  rebuild bit-identical: {ok_rep}")
    assert ok_med and ok_rep, "frozen representation did not reproduce"

    # ---------------- 3.2 cohort membership ----------------
    section("3.2  COHORT MEMBERSHIP (read from the Phase 1 manifest)")
    man = load_manifest()
    ref_cols = list(expr.columns)
    cohorts = {c: materialise(c, man, ref_cols) for c in COHORTS}
    s0 = cohorts["S0"]
    membership = pd.DataFrame({"sample_id": s0.sample_ids})
    for c in COHORTS:
        s = set(cohorts[c].sample_ids)
        membership[f"in_{c}"] = membership["sample_id"].isin(s)
    excl = exclusion_table(man)
    crit = excl.set_index("sample_id")["criterion"]
    membership["exclusion_criterion"] = membership["sample_id"].map(
        crit).fillna("")
    print(f"  S0 (superset stored): {len(s0)}")
    for c in COHORTS:
        print(f"    {c}: {len(cohorts[c])}")

    # ---------------- 3.3-3.5 tables ----------------
    section("3.3-3.5  ENDPOINT, COVARIATES, PARTITION AND BATCH")
    meta = pd.read_csv(paths.SAMPLE_METADATA)
    endpoint = build_endpoint(meta, s0.sample_ids)
    covars = build_covariates(meta, s0.sample_ids)
    batch = load_batch_assignments(allow_fetch=False)
    part = build_partition_and_batch(meta, batch, s0.sample_ids)
    print(f"  endpoint rows   : {len(endpoint)}")
    print(f"  covariate rows  : {len(covars)}  ({len(covars.columns)-1} columns)")
    print(f"  batch labels    : {int(part['combat_batch'].notna().sum())}/{len(part)}")
    print(f"  flagged stage/M : {int(covars['flag_stage_m_conflict'].sum())}"
          f" ({covars.loc[covars.flag_stage_m_conflict,'sample_id'].tolist()})")

    # ---------------- 3.6 write bundle ----------------
    section("3.6  WRITE BUNDLE (S0-scoped matrices + tables)")
    med_s0 = align_by_id(med_existing, s0.sample_ids)
    rep_s0 = align_by_id(rep_existing, s0.sample_ids)
    for name, m in [("expression_median_S0", med_s0),
                    ("expression_repiqr_disc443_S0", rep_s0)]:
        p = ANALYSIS_DIR / f"{name}.parquet"
        m.to_parquet(p, engine="pyarrow", compression="snappy")
        print(f"  {name}: {m.shape[0]} x {m.shape[1]}  "
              f"-> {p.stat().st_size/1e6:.1f} MB")

    membership.to_csv(ANALYSIS_DIR / "cohort_membership.csv", index=False)
    endpoint.to_csv(ANALYSIS_DIR / "endpoint.csv", index=False)
    covars.to_csv(ANALYSIS_DIR / "covariates.csv", index=False)
    pd.DataFrame([{"variable": k, "reason": v}
                  for k, v in EXCLUDED_COVARIATES.items()]).to_csv(
        ANALYSIS_DIR / "covariates_excluded.csv", index=False)
    part.to_csv(ANALYSIS_DIR / "partition_and_batch.csv", index=False)
    excl.to_csv(ANALYSIS_DIR / "exclusions.csv", index=False)

    # ---------------- 3.7 manifest ----------------
    section("3.7  MANIFEST AND COUNT VERIFICATION")
    inputs = {
        "probe_matrix": paths.PROBE_MATRIX,
        "frozen_median_spec": FROZEN / "frozen_gene_to_probeset_median.csv",
        "frozen_repiqr_spec": FROZEN / "frozen_gene_to_probe_discovery443.csv",
        "selection_population": FROZEN / "selection_population_discovery443.txt",
        "cohort_manifest": paths.ROOT / "reports/phase1/cohort_manifest.csv",
        "sample_metadata": paths.SAMPLE_METADATA,
        "sdrf": paths.ROOT / "data/external/E-GEOD-39582.sdrf.txt",
    }
    mf = base_manifest(inputs)
    mf["representations"] = {
        "primary": {"name": "median", "file": "expression_median_S0.parquet",
                    "genes": int(med_s0.shape[0]),
                    "samples": int(med_s0.shape[1]),
                    "status": "PROVISIONAL primary, uniform aggregation rule; "
                              "NOT claimed biologically optimal"},
        "sensitivity": {"name": "representative_iqr",
                        "selection_population": "discovery_443",
                        "validation_samples_used_in_selection": False,
                        "file": "expression_repiqr_disc443_S0.parquet",
                        "genes": int(rep_s0.shape[0]),
                        "samples": int(rep_s0.shape[1])},
        "historical_only": ["mean"],
        "not_approved": ["alpha_stratified_hybrid"],
    }
    mf["endpoint"] = endpoint_metadata()
    mf["cohorts"] = {}
    print("  cohort counts (actual vs expected):")
    all_ok = True
    for c in COHORTS:
        ids = set(cohorts[c].sample_ids)
        e = endpoint[endpoint.sample_id.isin(ids)]
        p = part[part.sample_id.isin(ids)]
        n, ev, ev60 = len(e), int((e.rfs_event_full == 1).sum()), \
            int((e.rfs_event_60units == 1).sum())
        exp = EXPECTED[c]
        ok = (n, ev, ev60) == exp
        all_ok &= ok
        print(f"    {c}: n={n:3d}/{exp[0]:3d}  ev={ev:3d}/{exp[1]:3d}  "
              f"ev60={ev60:3d}/{exp[2]:3d}  [{'OK' if ok else 'MISMATCH'}]")
        mf["cohorts"][c] = {
            "n": n, "events_full": ev, "events_60units": ev60,
            "expected": {"n": exp[0], "events_full": exp[1],
                         "events_60units": exp[2]},
            "matches_expected": ok,
            "role": "PRIMARY" if c == PRIMARY_COHORT else "sensitivity",
            "by_partition": {
                ds: {"n": int((p.dataset == ds).sum()),
                     "events_full": int(
                         (e[e.sample_id.isin(
                             p.loc[p.dataset == ds, "sample_id"])]
                          .rfs_event_full == 1).sum())}
                for ds in ("discovery", "validation")},
            "batch_composition_by_partition": batch_composition(
                part, membership, c),
        }
    assert all_ok, "cohort counts do not match expected values"

    mf["exclusions"] = {
        "total": int(len(excl)),
        "by_criterion": excl.groupby("criterion").size().to_dict(),
        "samples": excl[["sample_id", "criterion"]].to_dict(orient="records"),
    }
    mf["outputs"] = {}
    for p_ in sorted(ANALYSIS_DIR.glob("*")):
        if p_.name == "analysis_bundle_manifest.json":
            continue
        from src.bundle import sha256_file
        mf["outputs"][p_.name] = {"sha256": sha256_file(p_),
                                  "bytes": p_.stat().st_size}
    write_manifest(ANALYSIS_DIR / "analysis_bundle_manifest.json", mf)

    print(f"\n  bundle -> {ANALYSIS_DIR.relative_to(paths.ROOT)}/")
    print("  Original and frozen inputs unchanged. Nothing fitted.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
