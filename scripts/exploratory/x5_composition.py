"""X5. Computational composition estimates: ESTIMATE (re-implemented, validated) and MCP-counter. Registered plan §X5.

These are enrichment / abundance SCORES estimated from bulk expression. They are not measured cell fractions.
"""
import numpy as np
import pandas as pd

import estimate_py as E
import xcommon as X

MCP_POPS = ["T cells", "CD8 T cells", "Cytotoxic lymphocytes", "B lineage", "NK cells", "Monocytic lineage",
            "Myeloid dendritic cells", "Neutrophils", "Endothelial cells", "Fibroblasts"]


def validate_estimate():
    inp = pd.read_csv(E.PKG / "sample_input.gct", sep="\t", skiprows=2, index_col=0).drop(columns="Description")
    ref = pd.read_csv(E.PKG / "sample_estimate.gct", sep="\t", skiprows=2, index_col=0).drop(columns="Description").T
    pur = pd.read_csv(E.PKG / "sample_purity.gct", sep="\t", skiprows=2, index_col=0).drop(columns="Description").T
    r = E.estimate_score(E.filter_common(inp))
    rel = np.abs((r[["StromalScore", "ImmuneScore", "ESTIMATEScore"]].values - ref.values) / ref.values).max()
    rel_p = np.abs((r["TumorPurity"].values - pur["ESTIMATE-predicted tumor purity"].values) / pur["ESTIMATE-predicted tumor purity"].values).max()
    ok = rel < 1e-6   # registered gate: sample_input.gct -> sample_estimate.gct (Stromal, Immune, ESTIMATE scores)
    return {"registered_gate": "sample_estimate.gct scores within 1e-6 relative error", "max_rel_err_scores": rel, "passed": bool(ok),
            "observation_purity_file": {"max_rel_diff_vs_bundled_sample_purity_gct": rel_p,
                                        "note": "The bundled sample_purity.gct disagrees with the package's own purity formula cos(0.6049872018 + 0.0001467884*ESTIMATEScore), "
                                                "which is identical in estimateScore.R and plotPurity.R of v1.0.13; the file appears to predate that formula. "
                                                "Purity here follows the package code. Purity is a monotone transform of ESTIMATEScore, so rank-based results are unaffected."}}


def mcp(expr, markers):
    res, cov = {}, {}
    for pop in MCP_POPS:
        feats = [f for f in markers.get(pop, []) if f in expr.index]
        cov[pop] = {"n_markers": len(markers.get(pop, [])), "n_present": len(feats)}
        if feats:
            res[pop] = expr.loc[feats].mean(axis=0)
    return pd.DataFrame(res), cov


def main():
    out = X.outdir("x5_composition")
    gate = validate_estimate()
    if not gate["passed"]:
        raise SystemExit(f"ESTIMATE validation gate failed: {gate}")
    ctx = X.context()
    est = E.estimate_score(E.filter_common(ctx.X))
    est_cov = {"genes_after_common_filter": int(E.filter_common(ctx.X).shape[0]), "overlap": est.attrs["overlap"]}
    # MCP-counter, native probe-set mode
    pm = pd.read_parquet(X.ROOT / "data/processed/probe_matrix.parquet")[ctx.ids]
    ps = pd.read_csv(X.EXT / "MCPcounter_probesets.txt", sep="\t", header=None, names=["probe", "pop"])
    mk = ps.groupby("pop").probe.apply(list).to_dict()
    mcp_p, cov_p = mcp(pm, mk)
    gn = pd.read_csv(X.EXT / "MCPcounter_genes.txt", sep="\t")
    mkg = gn.groupby("Cell population")["HUGO symbols"].apply(list).to_dict()
    mcp_g, cov_g = mcp(ctx.X, mkg)
    agree = {pop: float(np.corrcoef(mcp_p[pop], mcp_g[pop])[0, 1]) for pop in mcp_p.columns if pop in mcp_g.columns}
    comp = est.join(mcp_p.add_prefix("MCP_"))
    comp.index.name = "sample_id"
    comp.to_csv(out / "composition_estimates.csv")
    mcp_g.add_prefix("MCPgene_").to_csv(out / "mcp_gene_mode_check.csv", index_label="sample_id")
    desc = comp.describe().T[["mean", "std", "min", "50%", "max"]]
    X.write_json(out / "x5_summary.json", {"plan_section": "X5", "label": "computational estimates, not measurements",
                                          "estimate_validation_gate": gate, "estimate_coverage": est_cov,
                                          "mcp_probe_mode_coverage": cov_p, "mcp_gene_mode_coverage": cov_g,
                                          "mcp_probe_vs_gene_mode_pearson": agree, "n": len(ctx.ids),
                                          "purity_out_of_bounds_n": int(est.TumorPurity.isna().sum()),
                                          "descriptives": desc.round(4).to_dict("index")})
    print(gate, est_cov)
    print(desc.round(3).to_string())
    print({k: round(v, 3) for k, v in agree.items()})
    print(cov_p)


if __name__ == "__main__":
    main()
