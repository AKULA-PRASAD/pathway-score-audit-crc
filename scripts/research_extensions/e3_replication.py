"""E3. Independent outcome-free biological replication (registered extension E3) in GSE41258, GSE143985, GSE103479.

Frozen gene sets and coverage rule; within-cohort scoring; MCP-counter gene-symbol mode; pre-registered hypotheses R1–R6
with a fixed decision rule. Cohorts are never pooled.
"""
import gzip
import io
import json
import re

import numpy as np
import pandas as pd
from scipy import stats

import rcommon as C
from rcommon import X

COHORTS = {
    "GSE41258": {"gpl": "GPL96_full.txt", "linear": True},
    "GSE143985": {"gpl": None, "linear": True},          # GPL570: frozen project platform file
    "GSE103479": {"gpl": "GPL23985_full.txt", "linear": False},
}


def read_matrix(acc):
    with gzip.open(C.EXT / f"{acc}_series_matrix.txt.gz", "rt", errors="replace") as f:
        lines = f.read().splitlines()
    meta = {}
    chars = []
    start = end = None
    for i, l in enumerate(lines):
        if l.startswith("!Sample_characteristics_ch1"):
            chars.append([v.strip('"') for v in l.split("\t")[1:]])
        elif l.startswith("!Sample_geo_accession"):
            meta["gsm"] = [v.strip('"') for v in l.split("\t")[1:]]
        elif l.startswith("!Sample_source_name_ch1"):
            meta["source"] = [v.strip('"') for v in l.split("\t")[1:]]
        elif l.startswith("!series_matrix_table_begin"):
            start = i + 1
        elif l.startswith("!series_matrix_table_end"):
            end = i
    ann = pd.DataFrame(index=meta["gsm"])
    for row in chars:
        for gsm, v in zip(meta["gsm"], row):
            if ":" in v:
                k, val = v.split(":", 1)
                ann.loc[gsm, k.strip().lower()] = val.strip()
    ann["source"] = meta.get("source")
    expr = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0, na_values=["null", "NA", ""])
    expr.columns = [c.strip('"') for c in expr.columns]
    expr.index = [str(i).strip('"') for i in expr.index]
    return ann, expr


def platform_map(fname):
    if fname is None:
        txt = gzip.open(X.ROOT / "data/raw/gpl570/GPL570_platform.soft", "rt") if False else open(X.ROOT / "data/raw/gpl570/GPL570_platform.soft", errors="replace")
        lines = txt.read().splitlines()
        s = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")) + 1
        e = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
        tab = pd.read_csv(io.StringIO("\n".join(lines[s:e])), sep="\t", dtype=str, low_memory=False)
        sym = "Gene Symbol"
    else:
        lines = open(C.EXT / fname, errors="replace").read().splitlines()
        s = next(i for i, l in enumerate(lines) if l.startswith("ID\t"))
        e = next((i for i, l in enumerate(lines) if l.startswith("!platform_table_end")), len(lines))
        tab = pd.read_csv(io.StringIO("\n".join(lines[s:e])), sep="\t", dtype=str, low_memory=False)
        sym = "Gene Symbol" if "Gene Symbol" in tab.columns else "Gene.Symbol"
    tab = tab[["ID", sym]].dropna()
    tab = tab[~tab[sym].str.contains("///")]
    return dict(zip(tab["ID"], tab[sym].str.strip()))


def gene_matrix(expr, pmap, linear):
    E = expr.copy()
    if linear:
        E = np.log2(E.clip(lower=1))
    E = E.loc[[i for i in E.index if i in pmap]]
    E["_g"] = [pmap[i] for i in E.index]
    return E.groupby("_g").median()


def labels(acc, ann):
    """Primary tumours and binary annotations: dict name -> Series(True=positive, False=negative), NaN otherwise."""
    a = ann.copy()
    L = {}
    if acc == "GSE41258":
        a = a[a["tissue"] == "Primary Tumor"]
        L["TP53"] = a["p53 mutation status"].map({"Mutant": True, "Wild-type": False})
        L["MSI"] = a["microsattelite instability"].map({"MSI-high": True, "MSS": False})
    elif acc == "GSE143985":
        L["TP53"] = a["tp53_mutation"].map({"Y": True, "N": False})
        L["BRAF"] = a["braf_mutation"].map({"Y": True, "N": False})
        L["MSI"] = a["msi_status"].map({"MSI_H": True, "MSS": False})
    elif acc == "GSE103479":
        L["TP53"] = a["p53"].map({"MT": True, "WT": False})
        L["BRAF"] = a["braf (v600e)"].map({"MT": True, "WT": False})
    return a.index.tolist(), {k: v.reindex(a.index) for k, v in L.items()}


HYP = [("R1", "TP53", "p53", -1), ("R2", "TP53", "E2F", +1), ("R3", "BRAF", "Wnt", -1), ("R4", "BRAF", "p53", +1), ("R5", "MSI", "p53", +1)]


def main():
    C.dirs()
    ctx = X.context()
    panel = ctx.panel
    mcp = pd.read_csv(X.EXT / "MCPcounter_genes.txt", sep="\t")
    fib = mcp.loc[mcp["Cell population"] == "Fibroblasts", "HUGO symbols"].tolist()
    rng = np.random.default_rng(C.SEED)
    rows, covg, cohort_info = [], [], {}
    for acc, cfg in COHORTS.items():
        ann, expr = read_matrix(acc)
        G = gene_matrix(expr, platform_map(cfg["gpl"]), cfg["linear"])
        ids, L = labels(acc, ann)
        G = G[[i for i in ids if i in G.columns]]
        Z = C.zrows(G)
        scores = {}
        for pid, lab in C.SETS:
            gs = panel["gene_sets"][pid]
            present = [g for g in gs["members_resolved"] if g in Z.index]
            ok = len(present) / int(gs["n_source"]) >= 0.8 and len(present) >= 20
            covg.append(dict(cohort=acc, pathway=lab, n_source=int(gs["n_source"]), n_present=len(present), coverage_rule_met=ok))
            if ok:
                scores[lab] = C.score(Z, present)
        fpres = [g for g in fib if g in Z.index]
        mfib = G.loc[fpres].mean(axis=0)
        cohort_info[acc] = {"n_tumours": int(G.shape[1]), "n_genes": int(G.shape[0]), "mcp_fibroblast_genes_present": f"{len(fpres)}/{len(fib)}",
                            "groups": {k: {"pos": int((v == True).sum()), "neg": int((v == False).sum())} for k, v in L.items()}}  # noqa: E712
        for hid, ann_name, sc, direction in HYP:
            if ann_name not in L or sc not in scores:
                rows.append(dict(hypothesis=hid, cohort=acc, eligible=False, reason="annotation or score unavailable"))
                continue
            lab = L[ann_name].reindex(G.columns)
            pos = scores[sc][(lab == True).values].to_numpy()  # noqa: E712
            neg = scores[sc][(lab == False).values].to_numpy()  # noqa: E712
            if len(pos) < 10 or len(neg) < 10:
                rows.append(dict(hypothesis=hid, cohort=acc, eligible=False, reason=f"groups {len(pos)}/{len(neg)} < 10", n_pos=len(pos), n_neg=len(neg)))
                continue
            vals = np.r_[pos, neg]
            grp = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
            est, lo, hi = C.boot_ci(lambda v, g: C.rank_biserial(v[g == 1], v[g == 0])[0], [vals, grp], rng, strat=grp)
            consistent = (lo > 0 if direction > 0 else hi < 0)
            opposite = (hi < 0 if direction > 0 else lo > 0)
            rows.append(dict(hypothesis=hid, cohort=acc, eligible=True, n_pos=len(pos), n_neg=len(neg), effect=est, lo=lo, hi=hi,
                             expected_sign=direction, replicated_in_cohort=bool(consistent), opposite_in_cohort=bool(opposite)))
        if "EMT" in scores and len(fpres) >= 5 and G.shape[1] >= 30:
            est, lo, hi = C.boot_ci(lambda a, b: stats.spearmanr(a, b)[0], [scores["EMT"].to_numpy(), mfib.to_numpy()], rng)
            rows.append(dict(hypothesis="R6", cohort=acc, eligible=True, n_pos=int(G.shape[1]), effect=est, lo=lo, hi=hi, expected_sign=1,
                             replicated_in_cohort=bool(est >= 0.7 and lo >= 0.5), opposite_in_cohort=bool(hi < 0)))
        else:
            rows.append(dict(hypothesis="R6", cohort=acc, eligible=False, reason="EMT score or fibroblast genes unavailable"))
        pd.DataFrame(scores).assign(MCP_fibroblasts=mfib).to_csv(C.RES / f"e3_scores_{acc}.csv", index_label="sample_id")
    res = pd.DataFrame(rows)
    res.to_csv(C.RES / "e3_replication_tests.csv", index=False)
    pd.DataFrame(covg).to_csv(C.RES / "e3_gene_coverage.csv", index=False)
    verdict = {}
    ref = {"R1": -0.296, "R2": 0.188, "R3": -0.455, "R4": 0.482, "R5": "eps2 0.167 (dMMR higher)", "R6": 0.950}
    for hid in ["R1", "R2", "R3", "R4", "R5", "R6"]:
        h = res[(res.hypothesis == hid) & (res.eligible == True)]  # noqa: E712
        n_el, n_rep, n_opp = len(h), int(h.replicated_in_cohort.sum()) if len(h) else 0, int(h.opposite_in_cohort.sum()) if len(h) else 0
        v = "NOT TESTABLE" if n_el == 0 else "REPLICATED" if (n_rep >= 2 and n_opp == 0) else "PARTIALLY REPLICATED" if (n_rep == 1 and n_opp == 0) else "NOT REPLICATED"
        verdict[hid] = {"eligible_cohorts": n_el, "consistent": n_rep, "opposite": n_opp, "verdict": v, "GSE39582_effect": ref[hid]}
    X.write_json(C.REP / "e3_summary.json", {"extension": "E3", "label": "post-hoc extension, outcome-free, independent cohorts never pooled",
                                             "cohorts": cohort_info, "verdicts": verdict, "bootstrap": {"B": 2000, "seed": C.SEED}})
    print(json.dumps(cohort_info, indent=1))
    print(pd.DataFrame(covg).to_string())
    print(res.round(3).to_string())
    print(json.dumps(verdict, indent=1))


if __name__ == "__main__":
    main()
