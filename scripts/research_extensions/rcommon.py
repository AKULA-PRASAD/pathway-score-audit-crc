"""Shared helpers for registered research extensions E1–E3 (docs/RESEARCH_EXTENSION_REGISTER.md). Outcome-free; read-only on frozen data."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "exploratory"))
import xcommon as X  # noqa: E402

SEED = 20261108
REP = ROOT / "reports" / "research_extensions"
RES = ROOT / "results" / "research_extensions"
FIG = ROOT / "figures" / "research_extensions"
EXT = ROOT / "data" / "external" / "research_extensions"
SETS = X.PATHWAYS          # (id, label)
CANONICAL_P53 = ["CDKN1A", "MDM2", "TP53I3", "SESN1", "ZMAT3", "FDXR", "RRM2B", "TNFRSF10B", "BAX", "DDB2", "XPC"]


def zrows(M):
    sd = M.std(axis=1, ddof=1)
    M = M.loc[sd > 0]
    return M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1, ddof=1), axis=0)


def score(Z, genes):
    g = [x for x in genes if x in Z.index]
    s = Z.loc[g].mean(axis=0)
    return (s - s.mean()) / s.std(ddof=1)


def rank_biserial(a, b):
    from scipy import stats
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 2 * U / (len(a) * len(b)) - 1, p


def boot_ci(fn, arrays, rng, B=2000, strat=None):
    """Percentile bootstrap CI for fn(*arrays); resamples rows jointly (stratified by `strat` labels if given)."""
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
    return est, float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def dirs():
    for d in (REP, RES, FIG):
        d.mkdir(parents=True, exist_ok=True)
