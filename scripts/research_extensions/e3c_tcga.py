"""E3 addendum: replication of R1–R6 in TCGA-COAD (docs/RESEARCH_EXTENSION_ADDENDUM_TCGA.md). Outcome-free; never pooled."""
import json

import numpy as np
import pandas as pd
from scipy import stats

import rcommon as C
from rcommon import X
from e3b_positive_control_qc import qc

T = X.ROOT / "data" / "external" / "tcga_coad_xena"
HYP = [("R1", "TP53", "p53", -1), ("R2", "TP53", "E2F", +1), ("R3", "BRAF", "Wnt", -1), ("R4", "BRAF", "p53", +1), ("R5", "MSI", "p53", +1)]


def load():
    G = pd.read_csv(T / "HiSeqV2.gz", sep="\t", index_col=0)
    G = G.groupby(level=0).median()
    clin = pd.read_csv(T / "COAD_clinicalMatrix", sep="\t", index_col=0, dtype=str, keep_default_na=False)
    prim = [s for s in G.columns if s.endswith("-01") and s in clin.index and clin.loc[s, "sample_type"] == "Primary Tumor"]
    G = G[prim]
    mut = pd.read_csv(T / "COAD_mc3_gene_level.txt.gz", sep="\t", index_col=0)
    L = {}
    for g, name in (("TP53", "TP53"), ("BRAF", "BRAF")):
        s = pd.Series(np.nan, index=prim, dtype=object)
        common = [x for x in prim if x in mut.columns]
        s[common] = (mut.loc[g, common] > 0).map({True: True, False: False})
        L[name] = s
    msi = clin.loc[prim, "MSI_updated_Oct62011"].map({"MSI-H": True, "MSS": False})
    L["MSI"] = msi
    return G, L


def main():
    ctx = X.context()
    G, L = load()
    rng = np.random.default_rng(C.SEED)
    rs, med = qc(G)
    Z = C.zrows(G)
    scores, cov = {}, []
    for pid, lab in C.SETS:
        gs = ctx.panel["gene_sets"][pid]
        present = [g for g in gs["members_resolved"] if g in Z.index]
        ok = len(present) / int(gs["n_source"]) >= 0.8 and len(present) >= 20
        cov.append(dict(pathway=lab, n_present=len(present), coverage_rule_met=ok))
        if ok:
            scores[lab] = C.score(Z, present)
    mcp = pd.read_csv(X.EXT / "MCPcounter_genes.txt", sep="\t")
    fib = [g for g in mcp.loc[mcp["Cell population"] == "Fibroblasts", "HUGO symbols"] if g in G.index]
    mfib = G.loc[fib].mean(axis=0)
    rows = []
    for hid, ann, sc, direction in HYP:
        lab = L[ann]
        pos = scores[sc][(lab == True).values].to_numpy()   # noqa: E712
        neg = scores[sc][(lab == False).values].to_numpy()  # noqa: E712
        if len(pos) < 10 or len(neg) < 10:
            rows.append(dict(hypothesis=hid, cohort="TCGA-COAD", eligible=False, n_pos=len(pos), n_neg=len(neg), reason="groups < 10"))
            continue
        v, g = np.r_[pos, neg], np.r_[np.ones(len(pos)), np.zeros(len(neg))]
        est, lo, hi = C.boot_ci(lambda vv, gg: C.rank_biserial(vv[gg == 1], vv[gg == 0])[0], [v, g], rng, strat=g)
        rows.append(dict(hypothesis=hid, cohort="TCGA-COAD", eligible=True, n_pos=len(pos), n_neg=len(neg), effect=est, lo=lo, hi=hi, expected_sign=direction,
                         replicated_in_cohort=bool(lo > 0 if direction > 0 else hi < 0), opposite_in_cohort=bool(hi < 0 if direction > 0 else lo > 0)))
    est, lo, hi = C.boot_ci(lambda a, b: stats.spearmanr(a, b)[0], [scores["EMT"].to_numpy(), mfib.to_numpy()], rng)
    rows.append(dict(hypothesis="R6", cohort="TCGA-COAD", eligible=True, n_pos=int(G.shape[1]), effect=est, lo=lo, hi=hi, expected_sign=1,
                     replicated_in_cohort=bool(est >= 0.7 and lo >= 0.5), opposite_in_cohort=bool(hi < 0)))
    emt_dis = [g for g in ctx.panel["gene_sets"]["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]["members_resolved"] if g in Z.index and g not in fib]
    dis = C.boot_ci(lambda a, b: stats.spearmanr(a, b)[0], [C.score(Z, emt_dis).to_numpy(), mfib.to_numpy()], rng)
    tc = pd.DataFrame(rows)
    tc.to_csv(C.RES / "e3_tcga_replication_tests.csv", index=False)
    # combined verdicts: QC-passing original cohorts + TCGA
    s3 = json.loads((C.REP / "e3_summary.json").read_text())
    orig = pd.read_csv(C.RES / "e3_replication_tests.csv")
    orig = orig[~orig.cohort.isin(s3.get("qc_failed_cohorts", []))]
    allr = pd.concat([orig, tc], ignore_index=True)
    verdicts = {}
    for hid in ["R1", "R2", "R3", "R4", "R5", "R6"]:
        h = allr[(allr.hypothesis == hid) & (allr.eligible == True)]  # noqa: E712
        n_el = len(h)
        n_rep = int(h.replicated_in_cohort.astype(bool).sum()) if n_el else 0
        n_opp = int(h.opposite_in_cohort.astype(bool).sum()) if n_el else 0
        v = "NOT TESTABLE" if n_el == 0 else "REPLICATED" if (n_rep >= 2 and n_opp == 0) else "PARTIALLY REPLICATED" if (n_rep == 1 and n_opp == 0) else "NOT REPLICATED"
        verdicts[hid] = {"eligible_cohorts": n_el, "consistent": n_rep, "opposite": n_opp, "verdict": v, "cohorts": list(h.cohort)}
    X.write_json(C.REP / "e3_tcga_summary.json", {"label": "registered addendum; outcome-free; TCGA all stages; never pooled",
                                                   "n_primary_tumours": int(G.shape[1]), "qc": {"median_r": med, "pairs": rs, "passed": med >= 0.5},
                                                   "coverage": cov, "groups": {k: {"pos": int((v == True).sum()), "neg": int((v == False).sum())} for k, v in L.items()},  # noqa: E712
                                                   "fibroblast_markers_present": len(fib), "R6_gene_disjoint": dis,
                                                   "verdicts_qc_passing_plus_tcga": verdicts})
    print("n", G.shape, "QC", round(med, 3))
    print(tc.round(3).to_string())
    print("gene-disjoint R6", dis)
    print(json.dumps(verdicts, indent=1))


if __name__ == "__main__":
    main()
