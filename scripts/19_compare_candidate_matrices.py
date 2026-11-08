"""Phase 2: compare candidate gene-level matrices and run QC.

Read-only over the candidates built by script 18. Uses no outcome data.

Usage:
    .venv/bin/python scripts/19_compare_candidate_matrices.py
"""

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                   # noqa: E402
from src.aggregation import STRATEGIES, aggregate       # noqa: E402
from src.cache import load_probe_matrix                 # noqa: E402
from src.gene_mapping import build_mapping              # noqa: E402

CAND = paths.PROCESSED_DIR / "candidates"
OUT = paths.ROOT / "reports" / "phase2"


def load(strategy: str) -> pd.DataFrame:
    return pd.read_parquet(CAND / f"gene_matrix_{strategy}.parquet",
                           engine="pyarrow")


def section(t: str) -> None:
    print("\n" + "=" * 72); print(t); print("=" * 72)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    findings: dict = {}
    mapping = build_mapping()
    um = mapping.unique_map
    per_gene = um.groupby(um.values).size()
    single = set(per_gene[per_gene == 1].index)
    multi = sorted(per_gene[per_gene > 1].index)

    section("1. DIMENSIONS AND QC")
    rows = []
    mats = {}
    for s in STRATEGIES:
        m = load(s); mats[s] = m
        v = m.to_numpy()
        rows.append({
            "strategy": s, "genes": m.shape[0], "samples": m.shape[1],
            "non_finite": int((~np.isfinite(v)).sum()),
            "n_missing": int(np.isnan(v).sum()),
            "min": round(float(v.min()), 3), "max": round(float(v.max()), 3),
            "mean": round(float(v.mean()), 4),
            "median": round(float(np.median(v)), 4),
            "sd": round(float(v.std()), 4),
        })
    qc = pd.DataFrame(rows)
    print(qc.to_string(index=False))
    qc.to_csv(OUT / "candidate_qc.csv", index=False)
    findings["qc"] = qc.to_dict(orient="records")

    ref_cols = list(mats[STRATEGIES[0]].columns)
    ref_idx = list(mats[STRATEGIES[0]].index)
    print(f"\n  identical sample ids and order across all candidates: "
          f"{all(list(mats[s].columns) == ref_cols for s in STRATEGIES)}")
    print(f"  identical gene index across all candidates: "
          f"{all(list(mats[s].index) == ref_idx for s in STRATEGIES)}")
    print(f"  all 585 samples retained: {len(ref_cols) == 585}")

    section("2. WHERE STRATEGIES CAN DIFFER AT ALL")
    print(f"  genes total                  : {len(ref_idx)}")
    print(f"  single-probe genes           : {len(single)}  "
          f"({len(single)/len(ref_idx)*100:.1f}%)")
    print(f"  multi-probe genes            : {len(multi)}  "
          f"({len(multi)/len(ref_idx)*100:.1f}%)")
    a = mats["representative_iqr"].loc[sorted(single)].to_numpy()
    b = mats["mean"].loc[sorted(single)].to_numpy()
    print(f"\n  single-probe genes identical across strategies: "
          f"{bool(np.allclose(a, b, atol=0, rtol=0))}")
    print("  => all disagreement is confined to the multi-probe genes.")
    findings["genes"] = {"total": len(ref_idx), "single_probe": len(single),
                         "multi_probe": len(multi)}

    section("3. PAIRWISE AGREEMENT ON MULTI-PROBE GENES")
    sub = {s: mats[s].loc[multi].to_numpy() for s in STRATEGIES}

    def per_gene_pearson(x, y):
        xc = x - x.mean(axis=1, keepdims=True)
        yc = y - y.mean(axis=1, keepdims=True)
        den = np.sqrt((xc ** 2).sum(1) * (yc ** 2).sum(1))
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den > 0, (xc * yc).sum(1) / den, np.nan)

    pair_rows = []
    for s1, s2 in itertools.combinations(STRATEGIES, 2):
        r = per_gene_pearson(sub[s1], sub[s2])
        d = np.abs(sub[s1] - sub[s2])
        pair_rows.append({
            "pair": f"{s1} vs {s2}",
            "median_r": round(float(np.nanmedian(r)), 4),
            "frac_r<0.9": round(float(np.nanmean(r < 0.9)), 4),
            "frac_r<0.7": round(float(np.nanmean(r < 0.7)), 4),
            "median_abs_diff": round(float(np.median(d)), 4),
            "p95_abs_diff": round(float(np.percentile(d, 95)), 4),
            "identical_genes": int((d.max(axis=1) == 0).sum()),
        })
    pairs = pd.DataFrame(pair_rows)
    print(pairs.to_string(index=False))
    pairs.to_csv(OUT / "candidate_pairwise_agreement.csv", index=False)
    findings["pairwise"] = pairs.to_dict(orient="records")

    section("4. DOES DISAGREEMENT TRACK WITHIN-GENE PROBE AGREEMENT?")
    wg = OUT / "within_gene_correlation_summary.csv"
    if wg.exists():
        gm = pd.read_csv(wg, index_col=0)["median_spearman_rho"]
        gm = gm.reindex(multi)
        r_ir = per_gene_pearson(sub["representative_iqr"], sub["mean"])
        comp = pd.DataFrame({"within_gene_rho": gm.values,
                             "rep_iqr_vs_mean_r": r_ir}, index=multi)
        comp["rho_decile"] = pd.qcut(comp["within_gene_rho"], 10, labels=False)
        agg = comp.groupby("rho_decile").agg(
            n=("rep_iqr_vs_mean_r", "size"),
            median_within_gene_rho=("within_gene_rho", "median"),
            median_strategy_agreement=("rep_iqr_vs_mean_r", "median")).round(3)
        print(agg.to_string())
        print("\n  => where a gene's probes agree, the strategy choice is almost")
        print("     irrelevant. Where they disagree, the choice decides what the")
        print("     gene-level value means.")
        comp.to_csv(OUT / "strategy_disagreement_vs_probe_agreement.csv")
        lo = comp[comp["within_gene_rho"] < 0.3]
        hi = comp[comp["within_gene_rho"] >= 0.7]
        print(f"\n  genes with probe rho < 0.3 : {len(lo):5d}, "
              f"median strategy agreement {lo['rep_iqr_vs_mean_r'].median():.3f}")
        print(f"  genes with probe rho >= 0.7: {len(hi):5d}, "
              f"median strategy agreement {hi['rep_iqr_vs_mean_r'].median():.3f}")
        findings["disagreement_vs_agreement"] = {
            "low_rho_genes": int(len(lo)),
            "low_rho_median_strategy_agreement": float(lo["rep_iqr_vs_mean_r"].median()),
            "high_rho_genes": int(len(hi)),
            "high_rho_median_strategy_agreement": float(hi["rep_iqr_vs_mean_r"].median()),
        }

    section("5. REPRESENTATIVE-PROBE CHOICE: IQR vs MEAN CRITERION")
    r1 = pd.read_csv(OUT / "representative_probes_representative_iqr.csv",
                     index_col=0)["probe_id"]
    r2 = pd.read_csv(OUT / "representative_probes_representative_mean.csv",
                     index_col=0)["probe_id"]
    same = (r1 == r2.reindex(r1.index))
    print(f"  genes where both criteria pick the SAME probe : {int(same.sum())} "
          f"({same.mean()*100:.1f}%)")
    print(f"  genes where they differ                       : {int((~same).sum())}")
    mp = [g for g in multi]
    same_m = same.reindex(mp)
    print(f"  among multi-probe genes only                  : "
          f"{int(same_m.sum())}/{len(mp)} ({same_m.mean()*100:.1f}%) agree")
    findings["representative_agreement"] = {
        "same_overall": int(same.sum()), "differ": int((~same).sum()),
        "same_multi_probe": int(same_m.sum()), "multi_probe": len(mp)}

    section("6. EXPRESSION DISTRIBUTION BY STRATEGY")
    dist = pd.DataFrame({s: pd.Series(mats[s].to_numpy().ravel()).describe(
        percentiles=[.05, .25, .5, .75, .95]) for s in STRATEGIES}).round(3)
    print(dist.to_string())
    dist.to_csv(OUT / "candidate_distributions.csv")
    print("\n  gene-level mean expression, multi-probe genes only:")
    gm_means = pd.DataFrame(
        {s: mats[s].loc[multi].to_numpy().mean(axis=1) for s in STRATEGIES},
        index=multi)
    print(gm_means.describe(percentiles=[.25, .5, .75]).round(3).to_string())
    print("\n  representative_* sit higher than median/mean because selecting")
    print("  one probe avoids averaging a bright probe with dim siblings.")

    section("7. REPRODUCIBILITY -- deterministic rebuild")
    expr = load_probe_matrix()
    ok = {}
    for s in STRATEGIES:
        again = aggregate(expr, um, s).matrix
        ident = bool(np.array_equal(again.to_numpy(), mats[s].to_numpy()))
        same_ax = (list(again.index) == list(mats[s].index)
                   and list(again.columns) == list(mats[s].columns))
        ok[s] = ident and same_ax
        print(f"  {s:22s} rebuild bit-identical: {ok[s]}")
    findings["reproducible"] = ok
    print(f"\n  all strategies deterministic: {all(ok.values())}")

    (OUT / "candidate_comparison.json").write_text(
        json.dumps(findings, indent=2, default=str) + "\n")
    section("DONE")
    print(f"  evidence -> {OUT.relative_to(paths.ROOT)}/")
    print("  No outcome data used. Probe cache unchanged. No strategy selected.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
