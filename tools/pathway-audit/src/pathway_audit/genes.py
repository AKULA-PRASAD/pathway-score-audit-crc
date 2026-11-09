"""Gene-level structure of a score: item–rest correlations, leave-one-gene-out, co-expression modules and their stability."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def _mean_std(z):
    s = z.mean(axis=0)
    return (s - s.mean()) / s.std(ddof=1)


def gene_contributions(z: pd.DataFrame, genes) -> pd.DataFrame:
    g = [x for x in genes if x in z.index]
    full = _mean_std(z.loc[g])
    rows = []
    for gene in g:
        rest = _mean_std(z.loc[[x for x in g if x != gene]])
        rows.append({"gene": gene, "item_rest_r": float(np.corrcoef(z.loc[gene], rest)[0, 1]),
                     "logo_r_with_full": float(np.corrcoef(rest, full)[0, 1])})
    out = pd.DataFrame(rows).sort_values("item_rest_r", ascending=False).reset_index(drop=True)
    r = np.corrcoef(z.loc[g].to_numpy())
    n = len(g)
    out.attrs["mean_inter_gene_r"] = float((r.sum() - n) / (n * (n - 1))) if n > 1 else float("nan")
    return out


def _clusters(zs: pd.DataFrame, cut: float) -> pd.Series:
    r = np.corrcoef(zs.to_numpy())
    d = np.clip(1 - r, 0, 2)
    np.fill_diagonal(d, 0)
    return pd.Series(fcluster(linkage(squareform(d, checks=False), "average"), t=cut, criterion="distance"), index=zs.index)


def modules(z: pd.DataFrame, genes, cut: float = 0.7, min_size: int = 10, n_boot: int = 200, seed: int = 0) -> pd.DataFrame:
    """Average-linkage modules (distance 1 − r, cut), with clusterwise Jaccard stability (Hennig 2007) over patient bootstraps."""
    g = [x for x in genes if x in z.index]
    zs = z.loc[g]
    lab = _clusters(zs, cut)
    sizes = lab.value_counts()
    mods = {i + 1: set(lab.index[lab == k]) for i, k in enumerate(sizes[sizes >= min_size].index)}
    rng = np.random.default_rng(seed)
    jac = {m: [] for m in mods}
    for _ in range(n_boot):
        ix = rng.integers(0, zs.shape[1], zs.shape[1])
        zb = zs.iloc[:, ix]
        zb = zb.sub(zb.mean(axis=1), axis=0).div(zb.std(axis=1, ddof=1).replace(0, np.nan), axis=0).fillna(0)
        cl = _clusters(zb, cut)
        groups = [set(cl.index[cl == k]) for k in cl.unique()]
        for m, s in mods.items():
            jac[m].append(max(len(s & b) / len(s | b) for b in groups))
    rows = []
    for m, s in mods.items():
        j = float(np.mean(jac[m])) if jac[m] else float("nan")
        rows.append({"module": m, "n_genes": len(s), "genes": sorted(s), "mean_jaccard": j,
                     "stability": "stable" if j >= 0.75 else "pattern" if j >= 0.6 else "weak" if j >= 0.5 else "dissolved"})
    return pd.DataFrame(rows)
