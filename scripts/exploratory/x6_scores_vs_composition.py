"""X6. Pathway scores versus estimated composition; marker overlap; EMT-ex sensitivity score. Registered plan §X6.

Descriptive only (Spearman rho with bootstrap 95 % intervals). The original EMT score is unchanged; EMT-ex is a separate
exploratory sensitivity score.
"""
import numpy as np
import pandas as pd
from scipy import stats

import estimate_py as E
import xcommon as X

COMP = ["StromalScore", "ImmuneScore", "TumorPurity", "MCP_Fibroblasts", "MCP_Endothelial cells",
        "MCP_Cytotoxic lymphocytes", "MCP_T cells", "MCP_Monocytic lineage"]


def boot_rho(a, b, rng, B=2000):
    n = len(a)
    r0 = stats.spearmanr(a, b)[0]
    rs = np.empty(B)
    for i in range(B):
        ix = rng.integers(0, n, n)
        rs[i] = stats.spearmanr(a[ix], b[ix])[0]
    return r0, np.percentile(rs, 2.5), np.percentile(rs, 97.5)


def main():
    ctx = X.context()
    out = X.outdir("x6_scores_vs_composition")
    z, _ = X.cohort_scores(ctx)
    comp = pd.read_csv(X.OUT / "x5_composition" / "composition_estimates.csv", index_col=0).loc[ctx.ids]
    gs = E.gene_sets()
    gn = pd.read_csv(X.EXT / "MCPcounter_genes.txt", sep="\t")
    mcpg = gn.groupby("Cell population")["HUGO symbols"].apply(set).to_dict()
    markers = {"ESTIMATE stromal": set(gs["StromalSignature"]), "ESTIMATE immune": set(gs["ImmuneSignature"]),
               "MCP fibroblasts": mcpg["Fibroblasts"], "MCP endothelial": mcpg["Endothelial cells"]}
    overlap = {}
    for pid, lab in X.PATHWAYS:
        mem = set(ctx.panel["gene_sets"][pid]["members_resolved"])
        overlap[lab] = {"n_members": len(mem)}
        for k, s in markers.items():
            ov = sorted(mem & s)
            overlap[lab][k] = {"n": len(ov), "pct_of_pathway": round(100 * len(ov) / len(mem), 1), "genes": ov}
        union = mem & (markers["ESTIMATE stromal"] | markers["MCP fibroblasts"] | markers["MCP endothelial"])
        overlap[lab]["stromal_union"] = {"n": len(union), "pct_of_pathway": round(100 * len(union) / len(mem), 1)}
    # ---- EMT-ex: frozen scoring formula, stromal-marker genes removed (all 452 as estimation population)
    emt_id = X.PATHWAYS[0][0]
    excl = markers["ESTIMATE stromal"] | markers["MCP fibroblasts"] | markers["MCP endothelial"]
    keep = [g for g in ctx.panel["gene_sets"][emt_id]["members_resolved"] if g not in excl and g in ctx.X.index]
    Xk = ctx.X.loc[keep]
    sd = Xk.std(axis=1, ddof=1)
    Xk = Xk.loc[sd > 0]
    emt_ex = Xk.sub(Xk.mean(axis=1), axis=0).div(Xk.std(axis=1, ddof=1), axis=0).mean(axis=0)
    emt_ex = (emt_ex - emt_ex.mean()) / emt_ex.std(ddof=1)
    emt_ex.rename("EMT_ex").to_csv(out / "emt_ex_score.csv", index_label="sample_id")
    # ---- correlations
    rng = np.random.default_rng(X.SEED)
    rows = []
    scores = z.join(emt_ex.rename("EMT_ex"))
    for lab in list(z.columns) + ["EMT_ex"]:
        for c in COMP:
            r, lo, hi = boot_rho(scores[lab].to_numpy(), comp[c].to_numpy(), rng)
            rows.append(dict(score=lab, composition=c, rho=r, lo=lo, hi=hi, n=len(scores)))
    cor = pd.DataFrame(rows)
    cor.to_csv(out / "score_composition_spearman.csv", index=False)
    r_emt = boot_rho(scores["EMT"].to_numpy(), scores["EMT_ex"].to_numpy(), rng)
    # MCP fibroblast vs ESTIMATE stromal (agreement of two independent estimators)
    r_est = boot_rho(comp["StromalScore"].to_numpy(), comp["MCP_Fibroblasts"].to_numpy(), rng)
    stromal_emtex = cor[(cor.score == "EMT_ex") & cor.composition.isin(["StromalScore", "MCP_Fibroblasts"])]
    rule = bool((stromal_emtex.rho.abs() >= 0.7).any())
    X.write_json(out / "x6_summary.json", {"plan_section": "X6", "label": "descriptive; estimated composition, not measurement",
                                          "overlap": overlap, "emt_ex": {"n_genes_original": 198, "n_genes_removed": 198 - len(keep),
                                                                         "n_genes_used": int(Xk.shape[0]),
                                                                         "rho_with_original_EMT": r_emt},
                                          "estimators_agreement_StromalScore_vs_MCP_fibroblasts": r_est,
                                          "registered_rule_EMT_ex_abs_rho_ge_0.7_with_stromal": rule,
                                          "bootstrap": {"B": 2000, "seed": X.SEED}})
    piv = cor.pivot(index="score", columns="composition", values="rho").round(3)
    print(piv.to_string())
    print({k: {kk: (vv["n"], vv["pct_of_pathway"]) for kk, vv in v.items() if isinstance(vv, dict)} for k, v in overlap.items()})
    print("EMT_ex genes", Xk.shape[0], "rho(EMT,EMT_ex)", r_emt, "stromal vs MCP fib", r_est, "rule", rule)
    print(stromal_emtex.round(3).to_string())


if __name__ == "__main__":
    main()
