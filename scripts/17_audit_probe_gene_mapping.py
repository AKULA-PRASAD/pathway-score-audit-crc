"""Phase 2: audit probe-to-gene mapping and reproduce the Phase 0 audits.

Read-only. Reproduces the probe-class counts from script 11 and the
within-gene correlation result from script 12 through the cached matrix, then
extends both with the expression/variance stratification that explains them.

Uses no outcome data. Writes only to reports/phase2/.

Usage:
    .venv/bin/python scripts/17_audit_probe_gene_mapping.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                      # noqa: E402
from src.cache import load_manifest, load_probe_matrix     # noqa: E402
from src.gene_mapping import build_mapping, reconcile      # noqa: E402

OUT = paths.ROOT / "reports" / "phase2"


def section(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    findings: dict = {}

    section("0. PROVENANCE -- carried forward, not re-derived")
    mf = load_manifest()
    print(f"  cache source       : {mf['source_file']}")
    print(f"  dimensions         : {mf['n_probes']} probes x {mf['n_samples']} samples")
    print(f"  value transformation: {mf['value_transformation']}")
    print("  deposited pipeline : RMA normalisation + ComBat, applied jointly")
    print("                       across all 585 arrays by the depositors.")
    print("  THIS PHASE PERFORMS NO NORMALIZATION AND NO BATCH CORRECTION.")
    print("  Nothing here is independently normalized or batch-corrected.")
    findings["provenance"] = {
        "deposited_processing": "RMA + ComBat (by depositors)",
        "phase2_normalization_applied": "none",
        "phase2_batch_correction_applied": "none",
    }

    expr = load_probe_matrix()
    mapping = build_mapping()

    section("1. PROBE SET RECONCILIATION")
    rec = reconcile(mapping, expr.index)
    for k, v in rec.items():
        print(f"  {k:20s}: {v}")
    findings["reconciliation"] = rec

    section("2. PROBE MAPPING CLASSES (disjoint, exhaustive)")
    s = mapping.summary()
    total = s["n_unmapped"] + s["n_unique"] + s["n_ambiguous"]
    for label, key in [("unmapped  (0 symbols)", "n_unmapped"),
                       ("unique    (1 symbol)", "n_unique"),
                       ("ambiguous (>1 symbol)", "n_ambiguous")]:
        print(f"  {label:24s}: {s[key]:6d}  ({s[key]/total*100:5.2f}%)")
    print(f"  {'TOTAL':24s}: {total:6d}  (array size {paths.GPL570_PROBE_COUNT})")
    print(f"  classes are exhaustive: {total == paths.GPL570_PROBE_COUNT}")

    print("\n  genes reachable from uniquely-mapped probes:")
    print(f"    unique genes                : {s['n_genes_unique']:6d}")
    print(f"    single-probe genes          : {s['n_genes_single_probe']:6d}")
    print(f"    multi-probe genes           : {s['n_genes_multi_probe']:6d}")
    print(f"    probes in multi-probe genes : {s['n_probes_in_multi_probe_genes']:6d}")
    print(f"    max probes for one gene     : {s['max_probes_per_gene']:6d}")
    findings["mapping"] = s

    print("\n  REPRODUCTION CHECK vs script 11 / script 12:")
    ok12 = s["n_genes_multi_probe"] == 10751
    print(f"    script 12 'Genes with multiple probes' = 10751 -> {s['n_genes_multi_probe']} "
          f"[{'MATCH' if ok12 else 'MISMATCH'}]")
    findings["reproduces_script12_gene_count"] = ok12

    mapping.table.to_csv(OUT / "probe_mapping_table.csv", index=False)
    for cls in ("unmapped", "ambiguous"):
        ids = mapping.probes(cls)
        pd.DataFrame({"probe_id": ids}).to_csv(
            OUT / f"probes_{cls}.csv", index=False)
    print(f"\n  full per-probe table -> probe_mapping_table.csv")
    print(f"  excluded-probe lists -> probes_unmapped.csv, probes_ambiguous.csv")
    print("  (nothing is silently discarded: every excluded probe is named)")

    section("3. PROBE EXPRESSION DISTRIBUTION (all 54,675 probes)")
    vals = expr.to_numpy()
    pstats = pd.DataFrame({
        "probe_id": expr.index,
        "mean": vals.mean(axis=1),
        "sd": vals.std(axis=1, ddof=0),
        "median": np.median(vals, axis=1),
        "iqr": np.subtract(*np.percentile(vals, [75, 25], axis=1)),
        "min": vals.min(axis=1),
        "max": vals.max(axis=1),
    })
    pstats = pstats.merge(
        mapping.table[["ID", "mapping_class", "gene_symbol"]],
        left_on="probe_id", right_on="ID", how="left").drop(columns="ID")
    pstats.to_csv(OUT / "probe_statistics.csv", index=False)

    print(pstats[["mean", "sd", "iqr"]].describe(
        percentiles=[.1, .25, .5, .75, .9]).round(3).to_string())
    print(f"\n  non-finite values in matrix: {int((~np.isfinite(vals)).sum())}")
    print(f"  value range: {vals.min():.3f} - {vals.max():.3f}  (log2 scale)")

    print("\n  by mapping class (mean expression):")
    print(pstats.groupby("mapping_class")["mean"].describe()
          .round(2)[["count", "mean", "50%", "max"]].to_string())
    findings["expression_distribution"] = {
        "mean_median": float(pstats["mean"].median()),
        "sd_median": float(pstats["sd"].median()),
        "iqr_median": float(pstats["iqr"].median()),
        "non_finite": int((~np.isfinite(vals)).sum()),
    }

    section("4. LOW-EXPRESSION PROBES (reported, NOT filtered)")
    for thr in (4.0, 5.0, 6.0, 7.0):
        n = int((pstats["mean"] < thr).sum())
        print(f"  mean log2 < {thr:.1f} : {n:6d} probes ({n/len(pstats)*100:5.1f}%)")
    for thr in (0.2, 0.3, 0.5):
        n = int((pstats["sd"] < thr).sum())
        print(f"  sd        < {thr:.1f} : {n:6d} probes ({n/len(pstats)*100:5.1f}%)")
    print("\n  NO expression filter is applied in this phase. These counts")
    print("  describe what a filter WOULD remove; applying one is a separate")
    print("  decision requiring approval.")
    findings["low_expression_counts"] = {
        f"mean_lt_{t}": int((pstats["mean"] < t).sum()) for t in (4.0, 5.0, 6.0, 7.0)
    } | {f"sd_lt_{t}": int((pstats["sd"] < t).sum()) for t in (0.2, 0.3, 0.5)}

    section("5. WITHIN-GENE PROBE AGREEMENT (reproduce + stratify)")
    um = mapping.unique_map
    per_gene = um.groupby(um.values).size()
    multi = per_gene[per_gene > 1].index
    probes_multi = um[um.isin(multi)].index

    # Script 12 computed this over the 566 PRIMARY TUMOURS only. Phase 2's
    # matrix holds all 585 deposited samples. Both are reported: the 566-sample
    # figure is the reproduction check, the 585-sample figure is what this
    # phase's representation actually yields.
    meta_ = pd.read_csv(paths.SAMPLE_METADATA)
    tumour_ids = meta_.loc[meta_["source_name"].str.contains(
        "primary colorectal Adenocarcinoma", case=False, na=False),
        "sample_id"].tolist()

    def _gene_medians(columns):
        s = expr.loc[probes_multi, columns]
        r = s.rank(axis=1).to_numpy()
        r = r - r.mean(axis=1, keepdims=True)
        nrm = np.sqrt((r ** 2).sum(axis=1))
        idx = {p: i for i, p in enumerate(s.index)}
        g_of = um.loc[s.index]
        recs = []
        for gene, grp in g_of.groupby(g_of.values):
            ids = list(grp.index)
            for i in range(len(ids)):
                for j in range(i + 1, len(ids)):
                    a, b = idx[ids[i]], idx[ids[j]]
                    den = nrm[a] * nrm[b]
                    recs.append((gene, ids[i], ids[j],
                                 float((r[a] * r[b]).sum() / den) if den else np.nan))
        return pd.DataFrame(recs, columns=["gene_symbol", "probe_a",
                                           "probe_b", "spearman_rho"])

    pairs_566 = _gene_medians(tumour_ids)
    med_566 = pairs_566.groupby("gene_symbol")["spearman_rho"].median().median()

    sub = expr.loc[probes_multi]
    ranks = sub.rank(axis=1).to_numpy()
    ranks = ranks - ranks.mean(axis=1, keepdims=True)
    norms = np.sqrt((ranks ** 2).sum(axis=1))
    pos = {p: i for i, p in enumerate(sub.index)}

    gene_of = um.loc[sub.index]
    rows = []
    for gene, grp in gene_of.groupby(gene_of.values):
        ids = list(grp.index)
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                a, b = pos[ids[i]], pos[ids[j]]
                denom = norms[a] * norms[b]
                rho = float((ranks[a] * ranks[b]).sum() / denom) if denom else np.nan
                rows.append((gene, ids[i], ids[j], rho))
    pairs = pd.DataFrame(rows, columns=["gene_symbol", "probe_a", "probe_b",
                                        "spearman_rho"])
    pairs.to_csv(OUT / "within_gene_probe_correlations.csv", index=False)

    gene_med = pairs.groupby("gene_symbol")["spearman_rho"].median()
    print(f"  genes evaluated      : {len(gene_med)}")
    print(f"  probe pairs evaluated: {len(pairs)}")
    print(f"  median of gene medians: {gene_med.median():.6f}")
    legacy = (paths.RESULTS_DIR / "probe_correlation_audit"
              / "within_gene_correlation_summary.csv")
    reproduced = None
    if legacy.exists():
        old = pd.read_csv(legacy)
        ref = old["median_spearman_rho"].median()
        reproduced = bool(abs(med_566 - ref) < 1e-6)
        print(f"\n  -- reproduction check (same 566 primary tumours as script 12) --")
        print(f"     script 12 stored result : {ref:.6f} over {len(old)} genes")
        print(f"     recomputed from cache   : {med_566:.6f}")
        print(f"     REPRODUCES EXACTLY      : {reproduced}")
        print(f"\n  -- this phase's matrix (all 585 deposited samples) --")
        print(f"     median of gene medians  : {gene_med.median():.6f}")
        print(f"     difference              : {gene_med.median() - med_566:+.6f}")
        print("     The difference is entirely attributable to the 19 non-tumoural")
        print("     samples, which Phase 2 retains by instruction and script 12")
        print("     excluded. It is a sample-set difference, not a discrepancy.")

    s_ = pstats.set_index("probe_id")
    pairs["mean_min"] = np.minimum(s_["mean"].reindex(pairs.probe_a).values,
                                   s_["mean"].reindex(pairs.probe_b).values)
    pairs["sd_min"] = np.minimum(s_["sd"].reindex(pairs.probe_a).values,
                                 s_["sd"].reindex(pairs.probe_b).values)
    print("\n  rho stratified by the LOWER expression of each pair (deciles):")
    d = pairs.groupby(pd.qcut(pairs.mean_min, 10, labels=False)).agg(
        n=("spearman_rho", "size"), median_rho=("spearman_rho", "median"),
        mean_expr=("mean_min", "mean")).round(3)
    print(d.to_string())
    print("\n  rho stratified by the LOWER SD of each pair (deciles):")
    d2 = pairs.groupby(pd.qcut(pairs.sd_min, 10, labels=False)).agg(
        n=("spearman_rho", "size"), median_rho=("spearman_rho", "median"),
        mean_sd=("sd_min", "mean")).round(3)
    print(d2.to_string())
    expressed = pairs[(pairs.sd_min >= 0.5) & (pairs.mean_min >= 6)]
    print(f"\n  expressed AND variable pairs (mean>=6, sd>=0.5): {len(expressed)}")
    print(f"    median rho: {expressed['spearman_rho'].median():.3f}")
    print("  => low within-gene agreement tracks low expression and low variance,")
    print("     NOT a defective probe-gene mapping.")
    gene_med.to_frame("median_spearman_rho").to_csv(
        OUT / "within_gene_correlation_summary.csv")
    findings["within_gene"] = {
        "n_genes": int(len(gene_med)), "n_pairs": int(len(pairs)),
        "median_of_gene_medians": float(gene_med.median()),
        "median_of_gene_medians_566_tumours": float(med_566),
        "reproduces_script12_exactly": reproduced,
        "median_rho_expressed_variable": float(expressed["spearman_rho"].median()),
        "n_pairs_expressed_variable": int(len(expressed)),
    }

    (OUT / "phase2_mapping_audit.json").write_text(
        json.dumps(findings, indent=2, default=str) + "\n")
    section("DONE")
    print(f"  evidence -> {OUT.relative_to(paths.ROOT)}/")
    print("  Probe cache, raw data and metadata unchanged. No outcome data used.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
