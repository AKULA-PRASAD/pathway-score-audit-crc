"""Associations with uncertainty: genotype effects with/without adjustment, partial correlations against a random-gene-set
null, and gene-disjoint marker tests. All bootstrap intervals are percentile intervals."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def boot_ci(fn, arrays, rng, B=1000, strat=None):
    n = len(arrays[0])
    est = fn(*arrays)
    vals = np.empty(B)
    idx = np.arange(n)
    for i in range(B):
        if strat is None:
            ix = rng.integers(0, n, n)
        else:
            ix = np.concatenate([rng.choice(idx[strat == k], (strat == k).sum(), replace=True) for k in np.unique(strat)])
        vals[i] = fn(*[a[ix] for a in arrays])
    return float(est), float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def rank_biserial(a, b):
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u / (len(a) * len(b)) - 1, p


def _ols_coef(y, x, adj):
    cols = [np.ones(len(y)), x] + ([adj] if adj is not None and adj.size else [])
    A = np.column_stack(cols)
    return np.linalg.lstsq(A, y, rcond=None)[0][1]


def genotype_effect(score: pd.Series, label: pd.Series, adjust: pd.DataFrame | None, rng, B=1000, min_group=10) -> dict:
    """Difference in score (SD units) between label==True and False, unadjusted and adjusted (OLS), bootstrap CIs.

    Verdict: 'robust' (adjusted CI excludes 0, same sign as unadjusted), 'explained by covariates' (unadjusted CI excludes
    0, adjusted includes 0 and |adjusted| < 0.5 |unadjusted|), 'not associated', or 'inconclusive'.
    """
    keep = label.notna()
    if adjust is not None:
        keep &= adjust.notna().all(axis=1)
    y = score[keep].to_numpy(float)
    x = label[keep].astype(bool).to_numpy(float)
    n1, n0 = int(x.sum()), int((1 - x).sum())
    out = {"n_pos": n1, "n_neg": n0}
    if n1 < min_group or n0 < min_group:
        out["verdict"] = f"not testable (groups {n1}/{n0} < {min_group})"
        return out
    out["rank_biserial"], out["p_mannwhitney"] = rank_biserial(y[x == 1], y[x == 0])
    out["unadjusted"] = boot_ci(lambda yy, xx: _ols_coef(yy, xx, None), [y, x], rng, B)
    if adjust is not None:
        A = adjust[keep].to_numpy(float)
        out["adjusted"] = boot_ci(lambda yy, xx, aa: _ols_coef(yy, xx, aa), [y, x, A], rng, B)
        u, a = out["unadjusted"], out["adjusted"]
        u_sig, a_sig = (u[1] > 0 or u[2] < 0), (a[1] > 0 or a[2] < 0)
        if a_sig and np.sign(a[0]) == np.sign(u[0]):
            out["verdict"] = "robust to adjustment"
        elif u_sig and not a_sig and abs(a[0]) < 0.5 * abs(u[0]):
            out["verdict"] = "explained by adjustment covariates"
        elif not u_sig and not a_sig:
            out["verdict"] = "not associated"
        else:
            out["verdict"] = "inconclusive"
    else:
        u = out["unadjusted"]
        out["verdict"] = "associated" if (u[1] > 0 or u[2] < 0) else "not associated"
    return out


def _rank_resid(y, Z):
    yr = stats.rankdata(y)
    A = np.column_stack([np.ones(len(y))] + [stats.rankdata(z) for z in Z.T])
    return yr - A @ np.linalg.lstsq(A, yr, rcond=None)[0]


def partial_spearman(a, b, Z):
    return float(np.corrcoef(_rank_resid(a, Z), _rank_resid(b, Z))[0, 1])


def random_set_null(z: pd.DataFrame, target: pd.Series, size: int, fn, rng, n=500, pool=None):
    """Distribution of fn(target, random-gene-set score) over n random sets of `size` genes."""
    pool = list(z.index if pool is None else pool)
    vals = []
    for _ in range(n):
        g = rng.choice(pool, size, replace=False)
        s = z.loc[g].mean(axis=0)
        vals.append(fn(target.to_numpy(), ((s - s.mean()) / s.std(ddof=1)).to_numpy()))
    return np.array(vals)


def gene_disjoint_marker_rho(z: pd.DataFrame, expr: pd.DataFrame, set_genes, marker_genes, rng, B=1000):
    """Spearman rho between the score built WITHOUT any marker genes and the mean of the marker genes."""
    sg = [g for g in set_genes if g in z.index and g not in set(marker_genes)]
    mk = [g for g in marker_genes if g in expr.index]
    if len(sg) < 10 or len(mk) < 3:
        return None
    s = z.loc[sg].mean(axis=0)
    m = expr.loc[mk].mean(axis=0)
    est = boot_ci(lambda a, b: stats.spearmanr(a, b)[0], [s.to_numpy(), m.to_numpy()], rng, B)
    return {"rho": est, "n_set_genes_used": len(sg), "n_markers": len(mk), "n_shared_removed": len(set(set_genes) & set(marker_genes))}
