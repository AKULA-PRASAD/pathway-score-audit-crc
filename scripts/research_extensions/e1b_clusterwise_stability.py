"""E1 deviation 1: clusterwise Jaccard stability (Hennig, Comput Stat Data Anal 2007;52:258-271).

Added after the registered ARI criterion proved degenerate for single-module sets. Does NOT replace the registered
decision (H-E1 evaluability uses ARI as registered). For each module, the mean over 200 patient bootstraps of the best
Jaccard overlap with any bootstrap cluster (same clustering rule, before the size filter). Hennig's guide: >= 0.75 stable,
0.6-0.75 pattern present, < 0.5 dissolved.
"""
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

import rcommon as C
from rcommon import X
from e1_gene_architecture import CUT


def raw_clusters(Zs):
    R = np.corrcoef(Zs.to_numpy())
    D = np.clip(1 - R, 0, 2)
    np.fill_diagonal(D, 0)
    return pd.Series(fcluster(linkage(squareform(D, checks=False), "average"), t=CUT, criterion="distance"), index=Zs.index)


def main():
    ctx = X.context()
    Z = C.zrows(ctx.X)
    genes = pd.read_csv(C.RES / "e1_gene_contributions.csv")
    rng = np.random.default_rng(C.SEED)
    rows = []
    for pid, lab in C.SETS:
        g = genes[genes.pathway == lab]
        Zs = Z.loc[g.gene]
        mods = {m: set(g.gene[g.module == m]) for m in sorted(g.module.unique()) if m > 0}
        jac = {m: [] for m in mods}
        for _ in range(200):
            ix = rng.integers(0, Zs.shape[1], Zs.shape[1])
            Zb = Zs.iloc[:, ix]
            Zb = Zb.sub(Zb.mean(axis=1), axis=0).div(Zb.std(axis=1, ddof=1).replace(0, np.nan), axis=0).fillna(0)
            cl = raw_clusters(Zb)
            groups = [set(cl.index[cl == k]) for k in cl.unique()]
            for m, s in mods.items():
                jac[m].append(max(len(s & b) / len(s | b) for b in groups))
        for m, v in jac.items():
            mj = float(np.mean(v))
            rows.append(dict(module=f"{lab}_M{m}", n_genes=len(mods[m]), mean_jaccard=mj,
                             class_hennig="stable" if mj >= 0.75 else "pattern" if mj >= 0.6 else "weak" if mj >= 0.5 else "dissolved"))
    out = pd.DataFrame(rows)
    out.to_csv(C.RES / "e1_clusterwise_jaccard.csv", index=False)
    print(out.round(3).to_string())


if __name__ == "__main__":
    main()
