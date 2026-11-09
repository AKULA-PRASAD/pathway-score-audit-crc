"""Independent-challenge tests (post-review, declared deviations; outcome-free).

C1  partial-correlation null: partial rho(EMT, random gene set | stromal, immune) for 1,000 random 194-gene sets
C2  R6 gene-disjoint: EMT minus MCP fibroblast markers vs fibroblast markers, per cohort (GSE39582, GSE41258, GSE143985)
C3  GSE41258 de-duplicated by patient id: R1, R5, R6 re-estimated (one primary tumour per patient)
C4  Wnt module BRAF effect adjusted for dMMR (and stromal, immune), bootstrap CI; BRAF effect within pMMR
"""
import numpy as np
import pandas as pd
from scipy import stats

import rcommon as C
from rcommon import X
from e2_conditioned_structure import partial_rho
from e3_replication import COHORTS, gene_matrix, labels, platform_map, read_matrix


def main():
    rng = np.random.default_rng(C.SEED)
    ctx = X.context()
    Z = C.zrows(ctx.X)
    z, _ = X.cohort_scores(ctx)
    comp = pd.read_csv(X.OUT / "x5_composition" / "composition_estimates.csv", index_col=0).loc[ctx.ids]
    W = comp[["StromalScore", "ImmuneScore"]].to_numpy()
    out = {}
    # C1
    expressed = Z.index[ctx.X.loc[Z.index].mean(axis=1) > ctx.X.mean(axis=1).median()]
    null = [partial_rho(z["EMT"].to_numpy(), C.score(Z, list(rng.choice(expressed, 194, replace=False))).to_numpy(), W) for _ in range(1000)]
    obs = pd.read_csv(C.RES / "e2_partial_correlations.csv").set_index("pair").partial_rho.to_dict()
    out["C1_partial_null"] = {"null_median": float(np.median(null)), "null_95": [float(np.percentile(null, 2.5)), float(np.percentile(null, 97.5))],
                              "observed": obs, "exceeds_null_97.5": {k: bool(v > np.percentile(null, 97.5)) for k, v in obs.items()},
                              "null_genes": "random 194-gene sets from above-median-expressed genes; 1,000 draws"}
    # C2
    mcp = pd.read_csv(X.EXT / "MCPcounter_genes.txt", sep="\t")
    fib = mcp.loc[mcp["Cell population"] == "Fibroblasts", "HUGO symbols"].tolist()
    emt = ctx.panel["gene_sets"]["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]["members_resolved"]
    emt_dis = [g for g in emt if g not in fib]
    c2 = {"n_fibroblast_markers": len(fib), "n_shared_with_EMT": len(set(fib) & set(emt)), "shared": sorted(set(fib) & set(emt))}

    def disjoint_rho(G):
        Zg = C.zrows(G)
        e = C.score(Zg, emt_dis)
        f = G.loc[[g for g in fib if g in G.index]].mean(axis=0)
        return C.boot_ci(lambda a, b: stats.spearmanr(a, b)[0], [e.to_numpy(), f.to_numpy()], rng)
    c2["GSE39582"] = disjoint_rho(ctx.X)
    gm = {}
    for acc in ("GSE41258", "GSE143985"):
        ann, expr = read_matrix(acc)
        G = gene_matrix(expr, platform_map(COHORTS[acc]["gpl"]), COHORTS[acc]["linear"])
        ids, L = labels(acc, ann)
        gm[acc] = (ann, G[[i for i in ids if i in G.columns]], L)
        c2[acc] = disjoint_rho(gm[acc][1])
    out["C2_R6_gene_disjoint"] = c2
    # C3
    ann, G, L = gm["GSE41258"]
    pid = ann.loc[G.columns, "patient id"]
    keep = ~pid.duplicated(keep="first")
    Gd = G.loc[:, keep.values]
    Zd = C.zrows(Gd)
    sc = {lab: C.score(Zd, [g for g in ctx.panel["gene_sets"][p]["members_resolved"] if g in Zd.index]) for p, lab in C.SETS}
    c3 = {"n_samples": int(G.shape[1]), "n_patients": int(keep.sum())}
    for hid, ann_name, s in (("R1", "TP53", "p53"), ("R5", "MSI", "p53")):
        lab = L[ann_name].reindex(Gd.columns)
        pos, neg = sc[s][(lab == True).values].to_numpy(), sc[s][(lab == False).values].to_numpy()  # noqa: E712
        v, g = np.r_[pos, neg], np.r_[np.ones(len(pos)), np.zeros(len(neg))]
        c3[hid] = (len(pos), len(neg)) + C.boot_ci(lambda vv, gg: C.rank_biserial(vv[gg == 1], vv[gg == 0])[0], [v, g], rng, strat=g)
    f = Gd.loc[[x for x in fib if x in Gd.index]].mean(axis=0)
    c3["R6"] = C.boot_ci(lambda a, b: stats.spearmanr(a, b)[0], [sc["EMT"].to_numpy(), f.to_numpy()], rng)
    out["C3_GSE41258_dedup"] = c3
    # C4
    genes = pd.read_csv(C.RES / "e1_gene_contributions.csv")
    wm = C.score(Z, list(genes.gene[(genes.pathway == "Wnt") & (genes.module == 1)])).to_numpy()
    cov = ctx.cov.loc[ctx.ids]
    braf, mmr = cov["braf.mutation"].astype(str).to_numpy(), cov["mmr.status"].astype(str).to_numpy()
    sel = np.isin(braf, ["M", "WT"]) & np.isin(mmr, ["dMMR", "pMMR"])
    x = (braf[sel] == "M").astype(float)
    d = (mmr[sel] == "dMMR").astype(float)

    def coef(y, xx, *adj):
        A = np.column_stack([np.ones(len(y)), xx, *adj])
        return np.linalg.lstsq(A, y, rcond=None)[0][1]
    c4 = {"n": int(sel.sum()), "braf_mut": int(x.sum()), "braf_mut_pMMR": int(((x == 1) & (d == 0)).sum()),
          "braf_only": C.boot_ci(lambda y, a: coef(y, a), [wm[sel], x], rng),
          "braf_plus_dMMR": C.boot_ci(lambda y, a, b: coef(y, a, b), [wm[sel], x, d], rng),
          "braf_plus_dMMR_composition": C.boot_ci(lambda y, a, b, s1, s2: coef(y, a, b, s1, s2), [wm[sel], x, d, W[sel, 0], W[sel, 1]], rng)}
    p = (d == 0)
    c4["within_pMMR_rank_biserial"] = C.rank_biserial(wm[sel][p & (x == 1)], wm[sel][p & (x == 0)])
    out["C4_wnt_module_braf_vs_dMMR"] = c4
    X.write_json(C.REP / "challenge_tests.json", {"label": "post-review challenge tests (declared deviations); outcome-free", **out})
    import json
    print(json.dumps(out, indent=1, default=float)[:4000])


if __name__ == "__main__":
    main()
