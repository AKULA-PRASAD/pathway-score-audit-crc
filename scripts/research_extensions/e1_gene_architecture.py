"""E1. Within-pathway gene architecture (registered extension E1). Outcome-free.

Per-gene item–rest correlations, leave-one-gene-out stability, fixed-rule co-expression modules (average linkage, distance
1 − r, cut 0.7, ≥10 genes; Wnt ≥5), bootstrap stability (ARI), module–annotation and module–composition associations,
and a representative-probe replication. Tests H-E1.
"""
import json

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

import rcommon as C
from rcommon import X

CUT = 0.7
ANN = {"tp53.mutation": ("M", "WT"), "kras.mutation": ("M", "WT"), "braf.mutation": ("M", "WT"),
       "cimp.status": ("+", "-"), "cin.status": ("+", "-"), "mmr.status": ("dMMR", "pMMR")}
COMP = ["StromalScore", "ImmuneScore", "TumorPurity", "MCP_Fibroblasts", "MCP_Cytotoxic lymphocytes"]


def modules(Zs, minsize):
    R = np.corrcoef(Zs.to_numpy())
    D = np.clip(1 - R, 0, 2)
    np.fill_diagonal(D, 0)
    lab = fcluster(linkage(squareform(D, checks=False), "average"), t=CUT, criterion="distance")
    s = pd.Series(lab, index=Zs.index)
    sizes = s.value_counts()
    keep = sizes[sizes >= minsize].index
    s = s.where(s.isin(keep), 0)
    # relabel modules by decreasing size: 1, 2, ...
    order = {old: i + 1 for i, old in enumerate(s[s > 0].value_counts().index)}
    return s.map(lambda v: order.get(v, 0))


def analyse(Z, genes_by_set, ctx, comp, rng, boot=True):
    out = {}
    for pid, lab in C.SETS:
        g = [x for x in genes_by_set[pid] if x in Z.index]
        Zs = Z.loc[g]
        full = C.score(Z, g)
        item = {}
        logo = {}
        for gene in g:
            rest = C.score(Z, [x for x in g if x != gene])
            item[gene] = float(np.corrcoef(Zs.loc[gene], rest)[0, 1])
            logo[gene] = float(np.corrcoef(rest, full)[0, 1])
        minsize = 5 if lab == "Wnt" else 10
        mod = modules(Zs, minsize)
        ari = []
        if boot:
            for _ in range(200):
                ix = rng.integers(0, Zs.shape[1], Zs.shape[1])
                Zb = Zs.iloc[:, ix]
                Zb = Zb.sub(Zb.mean(axis=1), axis=0).div(Zb.std(axis=1, ddof=1).replace(0, np.nan), axis=0).fillna(0)
                mb = modules(Zb, minsize)
                inmod = mod > 0
                ari.append(adjusted_rand_score(mod[inmod], mb[inmod]) if inmod.sum() > 1 else np.nan)
        mscores = {}
        for m in sorted(mod[mod > 0].unique()):
            mscores[f"{lab}_M{m}"] = C.score(Z, list(mod[mod == m].index))
        out[lab] = {"genes": g, "item_rest": item, "logo_r": logo, "modules": mod, "ari": ari, "module_scores": mscores,
                    "full": full}
    return out


def main():
    C.dirs()
    ctx = X.context()
    gsets = {pid: ctx.panel["gene_sets"][pid]["members_resolved"] for pid, _ in C.SETS}
    rng = np.random.default_rng(C.SEED)
    Z = C.zrows(ctx.X)
    comp = pd.read_csv(X.OUT / "x5_composition" / "composition_estimates.csv", index_col=0).loc[ctx.ids]
    cov = ctx.cov.loc[ctx.ids]
    A = analyse(Z, gsets, ctx, comp, rng)
    # representative-probe replication
    rep = pd.read_parquet(X.ROOT / "data/processed/candidates/gene_matrix_representative_iqr.parquet")[ctx.ids].astype(float)
    Zr = C.zrows(rep)
    Ar = analyse(Zr, gsets, ctx, comp, np.random.default_rng(C.SEED), boot=False)
    # ---- tables
    gene_rows, mod_rows, assoc_rows = [], [], []
    for lab, a in A.items():
        for gene in a["genes"]:
            gene_rows.append(dict(pathway=lab, gene=gene, item_rest_r=a["item_rest"][gene], logo_r_with_full=a["logo_r"][gene],
                                  module=int(a["modules"].get(gene, 0)), module_repIQR=int(Ar[lab]["modules"].get(gene, 0))))
        for name, s in a["module_scores"].items():
            m = int(name.split("_M")[1])
            members = list(a["modules"][a["modules"] == m].index)
            # representative-probe score of the SAME gene list
            srep = C.score(Zr, members)
            row = dict(pathway=lab, module=name, n_genes=len(members), r_with_pathway_score=float(np.corrcoef(s, a["full"])[0, 1]),
                       r_same_genes_repIQR=float(np.corrcoef(s, srep)[0, 1]),
                       canonical_p53_targets=", ".join(sorted(set(members) & set(C.CANONICAL_P53))),
                       top_genes=", ".join(sorted(members, key=lambda x: -a["item_rest"][x])[:8]))
            mod_rows.append(row)
            for ann, (g1, g2) in ANN.items():
                v = cov[ann].astype(str)
                for repname, sc in (("median", s), ("repIQR", srep)):
                    a1, a2 = sc[(v == g1).values], sc[(v == g2).values]
                    rb, p = C.rank_biserial(a1, a2)
                    assoc_rows.append(dict(module=name, annotation=ann, representation=repname, n1=len(a1), n2=len(a2), effect=rb, p=p, kind="rank-biserial"))
            for c in COMP:
                for repname, sc in (("median", s), ("repIQR", srep)):
                    r, p = stats.spearmanr(sc, comp[c])
                    assoc_rows.append(dict(module=name, annotation=c, representation=repname, n1=len(sc), n2=np.nan, effect=r, p=p, kind="spearman"))
    genes = pd.DataFrame(gene_rows)
    mods = pd.DataFrame(mod_rows)
    assoc = pd.DataFrame(assoc_rows)
    m = (assoc.representation == "median") & (assoc.kind == "rank-biserial")
    assoc["q_BH"] = np.nan
    assoc.loc[m, "q_BH"] = X.bh(assoc.loc[m, "p"].to_numpy())
    genes.to_csv(C.RES / "e1_gene_contributions.csv", index=False)
    mods.to_csv(C.RES / "e1_modules.csv", index=False)
    assoc.to_csv(C.RES / "e1_module_associations.csv", index=False)
    stab = {lab: {"n_modules": int(len(a["module_scores"])), "n_unassigned": int((a["modules"] == 0).sum()), "n_genes": len(a["genes"]),
                  "median_ARI": float(np.nanmedian(a["ari"])) if a["ari"] else None,
                  "ARI_IQR": [float(np.nanpercentile(a["ari"], 25)), float(np.nanpercentile(a["ari"], 75))] if a["ari"] else None,
                  "stable": bool(a["ari"] and np.nanmedian(a["ari"]) >= 0.7),
                  "min_logo_r": float(min(a["logo_r"].values())),
                  "partition_agreement_repIQR_ARI": float(adjusted_rand_score(a["modules"], Ar[lab]["modules"].reindex(a["modules"].index).fillna(0)))}
            for lab, a in A.items()}
    # ---- H-E1
    p53 = mods[mods.pathway == "p53"]
    he1 = {"evaluable": False}
    if len(p53) >= 2 and stab["p53"]["stable"]:
        canon = p53[p53.canonical_p53_targets != ""].sort_values("n_genes", ascending=False)
        res = {}
        for repname in ("median", "repIQR"):
            t = assoc[(assoc.annotation == "tp53.mutation") & (assoc.representation == repname) & assoc.module.str.startswith("p53_")].set_index("module").effect
            im = assoc[(assoc.annotation == "ImmuneScore") & (assoc.representation == repname) & assoc.module.str.startswith("p53_")].set_index("module").effect
            canon_mod = canon.iloc[0].module if len(canon) else None
            res[repname] = {"canonical_module": canon_mod, "tp53_effect_by_module": t.round(3).to_dict(), "immune_rho_by_module": im.round(3).to_dict(),
                            "tp53_strongest_in_canonical": bool(canon_mod and t.abs().idxmax() == canon_mod),
                            "immune_strongest_elsewhere": bool(canon_mod and im.abs().idxmax() != canon_mod)}
        he1 = {"evaluable": True, **res, "supported": all(res[k]["tp53_strongest_in_canonical"] and res[k]["immune_strongest_elsewhere"] for k in res)}
    X.write_json(C.REP / "e1_summary.json", {"extension": "E1", "label": "post-hoc extension, outcome-free", "n": len(ctx.ids),
                                             "rule": {"linkage": "average", "distance": "1 - Pearson r", "cut": CUT, "min_module": "10 (Wnt 5)", "bootstrap": 200},
                                             "stability": stab, "H_E1": he1})
    print(json.dumps(stab, indent=1))
    print(mods.to_string())
    print(json.dumps(he1, indent=1))
    sig = assoc[(assoc.representation == "median") & (assoc.kind == "rank-biserial") & (assoc.q_BH < 0.05)]
    print(sig[["module", "annotation", "n1", "n2", "effect", "q_BH"]].round(4).to_string())
    print(assoc[(assoc.representation == "median") & (assoc.kind == "spearman")].pivot(index="module", columns="annotation", values="effect").round(2).to_string())


if __name__ == "__main__":
    main()
