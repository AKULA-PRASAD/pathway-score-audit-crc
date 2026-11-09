"""Post-review challenge tests (declared deviation) for E3 / TCGA addendum. Outcome-free.

C5  TCGA hypermutation sensitivity for R1–R4: hypermutated = > 300 non-silent mutated genes in MC3 (conventional
    threshold; cohort median ~113). Effects excluding hypermutated tumours, and OLS with log(mutated genes + 1).
C6  Cohort-specific random-gene-set baseline for R6: 200 random sets of the EMT set size (genes present), scored the same
    way, Spearman rho with the fibroblast-marker mean, in GSE39582, GSE41258, GSE143985 and TCGA-COAD.
"""
import numpy as np
import pandas as pd
from scipy import stats

import rcommon as C
from rcommon import X
from e3_replication import COHORTS, gene_matrix, labels, platform_map, read_matrix
from e3c_tcga import T, load


def rb_ci(pos, neg, rng):
    v, g = np.r_[pos, neg], np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    return C.boot_ci(lambda vv, gg: C.rank_biserial(vv[gg == 1], vv[gg == 0])[0], [v, g], rng, strat=g)


def ols(y, x, z):
    A = np.column_stack([np.ones(len(y)), x, z])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ b
    s2 = r @ r / (len(y) - A.shape[1])
    se = np.sqrt(np.diag(s2 * np.linalg.inv(A.T @ A)))
    return float(b[1]), float(b[1] - 1.96 * se[1]), float(b[1] + 1.96 * se[1])


def random_baseline(G, emt_genes, fib, rng, n=200):
    Z = C.zrows(G)
    pool = list(Z.index)
    size = len([g for g in emt_genes if g in Z.index])
    f = G.loc[[g for g in fib if g in G.index]].mean(axis=0).to_numpy()
    vals = [stats.spearmanr(C.score(Z, list(rng.choice(pool, size, replace=False))).to_numpy(), f)[0] for _ in range(n)]
    return {"median": float(np.median(vals)), "p2_5": float(np.percentile(vals, 2.5)), "p97_5": float(np.percentile(vals, 97.5)), "set_size": size}


def main():
    rng = np.random.default_rng(C.SEED)
    ctx = X.context()
    G, L = load()
    Z = C.zrows(G)
    sc = {lab: C.score(Z, [g for g in ctx.panel["gene_sets"][pid]["members_resolved"] if g in Z.index]) for pid, lab in C.SETS}
    mut = pd.read_csv(T / "COAD_mc3_gene_level.txt.gz", sep="\t", index_col=0)
    burden = mut.sum(axis=0).reindex(G.columns)
    hyper = burden > 300
    c5 = {"rule": "hypermutated = > 300 non-silent mutated genes (MC3)", "n_with_mc3": int(burden.notna().sum()), "n_hypermutated": int(hyper.sum()),
          "median_burden": float(burden.median())}
    for hid, ann, s in (("R1", "TP53", "p53"), ("R2", "TP53", "E2F"), ("R3", "BRAF", "Wnt"), ("R4", "BRAF", "p53")):
        lab = L[ann]
        ok = lab.notna() & burden.notna()
        nh = ok & ~hyper
        pos, neg = sc[s][(nh & (lab == True)).values].to_numpy(), sc[s][(nh & (lab == False)).values].to_numpy()  # noqa: E712
        e_excl = rb_ci(pos, neg, rng)
        y = sc[s][ok.values].to_numpy()
        x = lab[ok].astype(bool).to_numpy(float)
        z = np.log1p(burden[ok].to_numpy(float))
        c5[hid] = {"excluding_hypermutated": {"n_pos": len(pos), "n_neg": len(neg), "rank_biserial": e_excl},
                   "ols_adjusted_log_burden": ols(y, x, z),
                   "n_pos_hypermutated": int((ok & hyper & (lab == True)).sum()), "n_pos_total": int((ok & (lab == True)).sum())}  # noqa: E712
    # C6
    mcp = pd.read_csv(X.EXT / "MCPcounter_genes.txt", sep="\t")
    fib = mcp.loc[mcp["Cell population"] == "Fibroblasts", "HUGO symbols"].tolist()
    emt = ctx.panel["gene_sets"]["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]["members_resolved"]
    emt_dis = [g for g in emt if g not in fib]
    c6 = {"GSE39582": random_baseline(ctx.X, emt_dis, fib, rng), "TCGA-COAD": random_baseline(G, emt_dis, fib, rng)}
    for acc in ("GSE41258", "GSE143985"):
        a, expr = read_matrix(acc)
        Gc = gene_matrix(expr, platform_map(COHORTS[acc]["gpl"]), COHORTS[acc]["linear"])
        ids, _ = labels(acc, a)
        c6[acc] = random_baseline(Gc[[i for i in ids if i in Gc.columns]], emt_dis, fib, rng)
    X.write_json(C.REP / "tcga_challenge_tests.json", {"label": "post-review challenge tests (declared deviation); outcome-free",
                                                        "C5_tcga_hypermutation": c5, "C6_random_baseline_R6": c6})
    import json
    print(json.dumps({"C5": c5, "C6": c6}, indent=1, default=float))


if __name__ == "__main__":
    main()
