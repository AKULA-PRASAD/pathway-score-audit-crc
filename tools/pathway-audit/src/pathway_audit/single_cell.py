"""Which annotated cell types supply a gene set's transcripts? (single-cell source check)

Given single-cell UMI counts and cell-type labels, each gene set's share of transcripts per cell type is
compared with the same share for random gene sets matched on size and mean-expression decile. An
enrichment above 1 means the cell type supplies more of the set's transcripts than expression-matched
random genes would; below 1, less. With a group column (for example patient), shares are computed per
group and summarised by the median, so that one large sample cannot dominate.

This is the method of the registered single-cell extension of the colorectal case study
(docs/SINGLE_CELL_EXTENSION_PLAN.md in the project repository). It measures mRNA only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp


def _matched_random(target_idx, deciles, rng, n_random):
    """Sparse indicator (n_random x n_genes) of random sets matched to the target on mean-expression decile."""
    target = set(target_idx.tolist())
    pools = {d: np.array([i for i in np.nonzero(deciles == d)[0] if i not in target]) for d in range(10)}
    need = pd.Series(deciles[target_idx]).value_counts()
    rows, cols = [], []
    for r in range(n_random):
        for d, k in need.items():
            pick = rng.choice(pools[int(d)], int(min(k, len(pools[int(d)]))), replace=False)
            rows.extend([r] * len(pick))
            cols.extend(pick.tolist())
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_random, len(deciles)))


def cell_type_sources(counts, genes, cell_types, gene_sets, groups=None, min_detected=0.01, min_cells=20,
                      n_random=1000, seed=0):
    """Enrichment of each gene set's transcripts in each cell type versus expression-matched random gene sets.

    counts: cells x genes UMI matrix (scipy sparse or dense); genes: gene symbols (len = n_genes);
    cell_types: label per cell; gene_sets: {name: [genes]}; groups: optional label per cell (e.g. patient).
    Genes detected in fewer than `min_detected` of cells are excluded; groups need >= `min_cells` cells per cell type
    to contribute that cell type. Returns a long DataFrame (set, cell_type, n_genes, n_groups, share, random_share,
    enrichment), with medians across groups when `groups` is given.
    """
    X = sp.csr_matrix(counts, dtype=np.float64)
    genes = np.asarray(genes).astype(str)
    cell_types = np.asarray(cell_types).astype(str)
    if X.shape != (len(cell_types), len(genes)):
        raise ValueError("counts must be cells x genes, matching cell_types and genes")
    groups = np.asarray(groups).astype(str) if groups is not None else np.array(["all"] * X.shape[0])
    detected = np.asarray((X > 0).sum(axis=0)).ravel() / X.shape[0]
    keep = detected >= min_detected
    X, genes = X[:, keep], genes[keep]
    total = np.asarray(X.sum(axis=0)).ravel()
    order = pd.Series(total).rank(method="first").to_numpy()
    deciles = pd.qcut(order, 10, labels=False)
    types = sorted(set(cell_types))
    pos = pd.Series(np.arange(len(genes)), index=genes)
    pos = pos[~pos.index.duplicated()]
    rng = np.random.default_rng(seed)
    rows = []
    for name, members in gene_sets.items():
        tidx = np.array(sorted(pos.reindex([g for g in members if g in pos.index]).to_numpy()), dtype=int)
        if len(tidx) < 5:
            continue
        R = _matched_random(tidx, deciles, rng, n_random)
        per = []
        for g in sorted(set(groups)):
            m = groups == g
            sums, ok = [], []
            for t in types:
                sel = m & (cell_types == t)
                ok.append(sel.sum() >= min_cells)
                sums.append(np.asarray(X[sel].sum(axis=0)).ravel() if sel.any() else np.zeros(len(genes)))
            S = np.vstack(sums).T  # genes x types
            ok = np.array(ok)
            if ok.sum() < 2:
                continue
            S = S[:, ok]
            obs = S[tidx].sum(axis=0)
            obs_share = obs / obs.sum()
            rnd = R @ S
            rnd_share = (rnd / rnd.sum(axis=1, keepdims=True)).mean(axis=0)
            for t, sh, rs in zip(np.array(types)[ok], obs_share, rnd_share, strict=True):
                per.append({"group": g, "cell_type": t, "share": sh, "random_share": rs, "enrichment": sh / rs if rs > 0 else np.nan})
        if not per:
            continue
        df = pd.DataFrame(per)
        agg = df.groupby("cell_type").agg(n_groups=("group", "nunique"), share=("share", "median"),
                                          random_share=("random_share", "median"), enrichment=("enrichment", "median")).reset_index()
        agg.insert(0, "n_genes", len(tidx))
        agg.insert(0, "set", name)
        rows.append(agg)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def verdicts(table, high=2.0, low=0.5):
    """Plain-language summary: which cell types supply each set's transcripts (enrichment >= high) or are depleted (<= low)."""
    out = {}
    for name, t in table.groupby("set"):
        hi = t.loc[t["enrichment"] >= high].sort_values("enrichment", ascending=False)
        lo = t.loc[t["enrichment"] <= low]
        parts = []
        if len(hi):
            parts.append("mainly supplied by " + ", ".join(f"{r.cell_type} ({r.enrichment:.1f}x)" for r in hi.itertuples()))
        else:
            parts.append(f"no cell type supplies >= {high:g}x its expression-matched share")
        if len(lo):
            parts.append("depleted in " + ", ".join(f"{r.cell_type} ({r.enrichment:.2f}x)" for r in lo.itertuples()))
        out[name] = "; ".join(parts)
    return out
