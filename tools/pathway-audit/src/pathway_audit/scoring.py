"""Gene-set scoring: mean of per-gene z-scores across samples (the most common bulk pathway score)."""
from __future__ import annotations

import pandas as pd


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    """Per-gene z-scores across samples (ddof=1). Zero-variance genes are dropped."""
    sd = expr.std(axis=1, ddof=1)
    e = expr.loc[sd > 0]
    return e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1, ddof=1), axis=0)


def score_gene_set(z: pd.DataFrame, genes, min_genes: int = 10, min_coverage: float = 0.8, standardise: bool = True):
    """Return (score Series or None, coverage dict). Coverage rule: >= min_coverage of genes present and >= min_genes."""
    genes = list(dict.fromkeys(genes))
    present = [g for g in genes if g in z.index]
    cov = {"n_genes": len(genes), "n_present": len(present), "coverage": len(present) / max(len(genes), 1)}
    cov["coverage_ok"] = bool(cov["coverage"] >= min_coverage and len(present) >= min_genes)
    if not cov["coverage_ok"]:
        return None, cov
    s = z.loc[present].mean(axis=0)
    if standardise:
        s = (s - s.mean()) / s.std(ddof=1)
    return s, cov
