"""Phase 2B: selection-set sensitivity, batch association and probe quality.

Read-only over the probe cache and the four existing candidate matrices,
which are preserved untouched. New candidates are written separately.

OUTCOME-BLIND: no endpoint column is read anywhere. Sample subsets are
defined by study design (discovery/validation) and by the Phase 1 eligibility
manifest, never by recurrence status.

NO normalization and NO batch correction is performed. ComBat batch labels
are used DESCRIPTIVELY, to ask whether a selection criterion tracks batch.

Usage:
    .venv/bin/python scripts/20_phase2b_selection_sensitivity.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                      # noqa: E402
from src.aggregation import aggregate                      # noqa: E402
from src.cache import load_probe_matrix                    # noqa: E402
from src.external_annotation import load_batch_assignments  # noqa: E402
from src.gene_mapping import build_mapping                 # noqa: E402

CAND = paths.PROCESSED_DIR / "candidates"
OUT = paths.ROOT / "reports" / "phase2b"


def section(t: str) -> None:
    print("\n" + "=" * 74); print(t); print("=" * 74)


def eta_squared_by_group(values: np.ndarray, codes: np.ndarray,
                         n_groups: int) -> np.ndarray:
    """Fraction of each probe's variance explained by a grouping factor.

    One-way ANOVA eta^2 per probe, vectorised. Purely descriptive: nothing is
    adjusted, corrected or removed.
    """
    n = values.shape[1]
    grand = values.mean(axis=1, keepdims=True)
    sst = ((values - grand) ** 2).sum(axis=1)
    ssb = np.zeros(values.shape[0])
    for g in range(n_groups):
        m = codes == g
        k = int(m.sum())
        if k == 0:
            continue
        gm = values[:, m].mean(axis=1)
        ssb += k * (gm - grand[:, 0]) ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(sst > 0, ssb / sst, np.nan)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    findings: dict = {}

    expr = load_probe_matrix()
    mapping = build_mapping()
    um = mapping.unique_map
    per_gene = um.groupby(um.values).size()
    multi_genes = sorted(per_gene[per_gene > 1].index)
    probes_multi = um[um.isin(multi_genes)].index

    meta = pd.read_csv(paths.SAMPLE_METADATA)
    man = pd.read_csv(paths.ROOT / "reports" / "phase1" / "cohort_manifest.csv")

    all_585 = list(expr.columns)
    disc_443 = meta.loc[(meta["dataset"] == "discovery")
                        & meta["source_name"].str.contains(
                            "primary colorectal Adenocarcinoma",
                            case=False, na=False), "sample_id"].tolist()
    s3_354 = man.loc[man["in_S3"] & (man["dataset"] == "discovery"),
                     "sample_id"].tolist()
    disc_443 = [s for s in all_585 if s in set(disc_443)]
    s3_354 = [s for s in all_585 if s in set(s3_354)]

    SETS = {"all_585": all_585, "discovery_443": disc_443,
            "s3_discovery_354": s3_354}

    section("0. SELECTION SETS (defined by design and eligibility, not outcome)")
    for k, v in SETS.items():
        print(f"  {k:20s} n = {len(v)}")
    assert len(disc_443) == 443 and len(s3_354) == 354, "selection set size"
    findings["selection_sets"] = {k: len(v) for k, v in SETS.items()}

    # ---------------------------------------------------------------
    section("1. REPRESENTATIVE PROBE CHOICE vs SELECTION SET (IQR)")
    # ---------------------------------------------------------------
    reps = {}
    for name, cols in SETS.items():
        res = aggregate(expr, um, "representative_iqr", selection_columns=cols)
        reps[name] = res.representative_probe
        if name != "all_585":
            path = CAND / f"gene_matrix_representative_iqr_{name}.parquet"
            res.matrix.to_parquet(path, engine="pyarrow", compression="snappy")
            print(f"  new candidate -> {path.relative_to(paths.ROOT)}")
        print(f"  {name:20s} ties={res.diagnostics['n_genes_with_tied_score']}")

    rep_df = pd.DataFrame(reps)
    rep_df.index.name = "gene_symbol"
    rep_df.to_csv(OUT / "representative_probe_by_selection_set.csv")

    print("\n  pairwise differences in SELECTED PROBE:")
    rows = []
    keys = list(SETS)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = keys[i], keys[j]
            diff = rep_df[a] != rep_df[b]
            diff_multi = diff.reindex(multi_genes)
            rows.append({
                "comparison": f"{a} vs {b}",
                "genes_differing": int(diff.sum()),
                "pct_all_genes": round(diff.mean() * 100, 2),
                "pct_multi_probe": round(diff_multi.mean() * 100, 2),
            })
    cmp_df = pd.DataFrame(rows)
    print(cmp_df.to_string(index=False))
    cmp_df.to_csv(OUT / "selection_set_differences.csv", index=False)
    findings["selection_differences"] = cmp_df.to_dict(orient="records")
    print("\n  (single-probe genes can never differ; pct_multi_probe is the")
    print("   meaningful figure.)")

    # ---- which genes are most affected ----
    vals = expr.to_numpy()
    q1, q3 = np.percentile(vals, [25, 75], axis=1)
    pstat = pd.DataFrame({"mean": vals.mean(axis=1),
                          "sd": vals.std(axis=1, ddof=0),
                          "iqr": q3 - q1}, index=expr.index)

    wg = pd.read_csv(paths.ROOT / "reports" / "phase2"
                     / "within_gene_correlation_summary.csv",
                     index_col=0)["median_spearman_rho"]

    changed = rep_df.index[(rep_df["all_585"] != rep_df["s3_discovery_354"])]
    changed = [g for g in changed if g in set(multi_genes)]
    aff = pd.DataFrame({
        "gene_symbol": changed,
        "n_probes": per_gene.reindex(changed).values,
        "within_gene_rho": wg.reindex(changed).values,
        "probe_all585": rep_df.loc[changed, "all_585"].values,
        "probe_s3_354": rep_df.loc[changed, "s3_discovery_354"].values,
    })
    aff["mean_all585"] = pstat["mean"].reindex(aff.probe_all585).values
    aff["mean_s3"] = pstat["mean"].reindex(aff.probe_s3_354).values
    aff["iqr_all585"] = pstat["iqr"].reindex(aff.probe_all585).values
    aff["iqr_s3"] = pstat["iqr"].reindex(aff.probe_s3_354).values
    aff["abs_iqr_gap"] = (aff.iqr_all585 - aff.iqr_s3).abs()
    aff = aff.sort_values("abs_iqr_gap", ascending=False)
    aff.to_csv(OUT / "genes_affected_by_selection_set.csv", index=False)

    print(f"\n  genes changing probe (all_585 -> s3_discovery_354): {len(aff)}")
    print(f"    median within-gene rho of AFFECTED genes   : "
          f"{aff['within_gene_rho'].median():.3f}")
    print(f"    median within-gene rho of UNAFFECTED genes : "
          f"{wg.reindex([g for g in multi_genes if g not in set(changed)]).median():.3f}")
    print("\n  most affected genes (largest IQR gap between the two choices):")
    print(aff.head(12)[["gene_symbol", "n_probes", "within_gene_rho",
                        "probe_all585", "probe_s3_354", "iqr_all585",
                        "iqr_s3"]].round(3).to_string(index=False))
    findings["affected_genes"] = {
        "n_changed_585_vs_354": int(len(aff)),
        "median_rho_affected": float(aff["within_gene_rho"].median()),
        "median_rho_unaffected": float(
            wg.reindex([g for g in multi_genes if g not in set(changed)]).median()),
        "top_affected": aff.head(20)["gene_symbol"].tolist(),
    }

    # ---------------------------------------------------------------
    section("2. BATCH ASSOCIATION OF SELECTED PROBES (descriptive only)")
    # ---------------------------------------------------------------
    batch = load_batch_assignments()
    if batch is None:
        print("  batch assignments unavailable; section skipped.")
        findings["batch"] = {"status": "unavailable"}
    else:
        bmap = batch.set_index("sample_id")["combat_batch"]
        cols = [s for s in all_585 if s in bmap.index]
        print(f"  samples with a recovered batch label: {len(cols)} / 585")
        print(f"  samples WITHOUT one (non-tumoural)  : {585 - len(cols)}")
        codes_s = pd.Categorical(bmap.reindex(cols))
        codes = codes_s.codes.astype(int)
        ng = len(codes_s.categories)
        print(f"  batches: {ng} -> {list(codes_s.categories)}")

        sub = expr.loc[probes_multi, cols]
        eta = eta_squared_by_group(sub.to_numpy(), codes, ng)
        eta_s = pd.Series(eta, index=sub.index, name="batch_eta2")
        eta_s.to_frame().join(pstat).to_csv(OUT / "probe_batch_eta2.csv")

        sel585 = set(rep_df.loc[multi_genes, "all_585"])
        is_sel = pd.Series(sub.index.isin(sel585), index=sub.index)
        print("\n  variance fraction explained by ComBat batch (eta^2):")
        print(f"    SELECTED probes      n={int(is_sel.sum()):5d}  "
              f"median={eta_s[is_sel].median():.4f}  "
              f"mean={eta_s[is_sel].mean():.4f}  "
              f"p90={eta_s[is_sel].quantile(.9):.4f}")
        print(f"    NOT-selected probes  n={int((~is_sel).sum()):5d}  "
              f"median={eta_s[~is_sel].median():.4f}  "
              f"mean={eta_s[~is_sel].mean():.4f}  "
              f"p90={eta_s[~is_sel].quantile(.9):.4f}")

        print("\n  eta^2 by IQR decile (is high IQR batch-driven?):")
        dec = pd.qcut(pstat["iqr"].reindex(sub.index), 10, labels=False)
        tbl = pd.DataFrame({"eta2": eta_s, "dec": dec}).groupby("dec").agg(
            n=("eta2", "size"), median_eta2=("eta2", "median"),
            mean_eta2=("eta2", "mean")).round(4)
        print(tbl.to_string())
        r_iqr_eta = float(pd.Series(eta_s).corr(
            pstat["iqr"].reindex(sub.index), method="spearman"))
        r_mean_eta = float(pd.Series(eta_s).corr(
            pstat["mean"].reindex(sub.index), method="spearman"))
        print(f"\n    Spearman(IQR,  batch eta^2) = {r_iqr_eta:+.4f}")
        print(f"    Spearman(mean, batch eta^2) = {r_mean_eta:+.4f}")
        hi = eta_s[is_sel].gt(0.25).mean()
        hi2 = eta_s[~is_sel].gt(0.25).mean()
        print(f"\n    selected probes with eta^2 > 0.25    : {hi*100:.1f}%")
        print(f"    not-selected probes with eta^2 > 0.25: {hi2*100:.1f}%")
        findings["batch"] = {
            "status": "ok", "n_samples_with_batch": len(cols),
            "n_batches": int(ng),
            "median_eta2_selected": float(eta_s[is_sel].median()),
            "median_eta2_not_selected": float(eta_s[~is_sel].median()),
            "spearman_iqr_eta2": r_iqr_eta,
            "spearman_mean_eta2": r_mean_eta,
            "frac_selected_eta2_gt_025": float(hi),
            "frac_not_selected_eta2_gt_025": float(hi2),
        }

    (OUT / "phase2b_part1.json").write_text(
        json.dumps(findings, indent=2, default=str) + "\n")
    print("\n  part 1 written.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
