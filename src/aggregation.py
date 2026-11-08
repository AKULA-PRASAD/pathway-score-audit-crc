"""Probe-to-gene aggregation strategies.

Every strategy here is OUTCOME-BLIND by construction: no function in this
module accepts, reads or can access recurrence status, survival time, or any
clinical variable. The only inputs are the probe-level expression matrix and
the annotation-derived probe-to-gene map.

Strategies
----------
representative_iqr
    Per gene, keep the single probe with the largest interquartile range
    across samples. IQR is a robust dispersion measure: it prefers the probe
    that actually varies between tumours over one sitting in background noise,
    without being dragged by outliers the way variance is.

representative_mean
    Per gene, keep the single probe with the largest mean expression.
    Prefers the best-hybridising probe; says nothing about whether it varies.

median
    Per gene, the median across its probes. Robust to one discordant probe,
    but for a gene whose probes disagree it returns a value that may
    correspond to no probe's biology.

mean
    Per gene, the arithmetic mean across its probes. Included as a comparator
    because it is the common default, NOT because it is recommended: it is the
    least robust to a discordant or cross-hybridising probe.

Ties
----
Representative strategies break ties deterministically on the lexicographically
smallest probe id, so a rebuild is reproducible. Tie counts are reported.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

STRATEGIES = ("representative_iqr", "representative_mean", "median", "mean")


@dataclass
class AggregationResult:
    matrix: pd.DataFrame                  # genes x samples
    strategy: str
    n_genes: int
    n_probes_used: int
    representative_probe: pd.Series | None = None
    diagnostics: dict = field(default_factory=dict)


def _probe_stats(expr: pd.DataFrame) -> pd.DataFrame:
    """Per-probe mean and IQR across samples. No outcome involvement."""
    values = expr.to_numpy()
    q1, q3 = np.percentile(values, [25, 75], axis=1)
    return pd.DataFrame(
        {"mean": values.mean(axis=1), "iqr": q3 - q1},
        index=expr.index,
    )


def _pick_representative(
    expr: pd.DataFrame, gene_map: pd.Series, criterion: str,
    selection_columns: list[str] | None = None,
) -> tuple[pd.Series, int]:
    """Choose one probe per gene by `criterion`; deterministic tie-break.

    `selection_columns` restricts WHICH SAMPLES are used to score probes. The
    emitted matrix always keeps every sample; only the selection is scoped.
    This is what makes a "freeze on discovery, apply everywhere" representation
    possible. It is outcome-blind: a sample subset is not an outcome.
    """
    scoring = expr if selection_columns is None else expr[selection_columns]
    stats = _probe_stats(scoring.loc[gene_map.index])
    df = pd.DataFrame({
        "probe": gene_map.index,
        "gene": gene_map.values,
        "score": stats[criterion].reindex(gene_map.index).values,
    })
    # descending score, then ascending probe id -> deterministic winner
    df = df.sort_values(["gene", "score", "probe"],
                        ascending=[True, False, True], kind="mergesort")
    winners = df.drop_duplicates("gene", keep="first")

    # how often was the winning score tied within its gene?
    top = df.groupby("gene")["score"].transform("max")
    tied = df[df["score"] == top].groupby("gene").size()
    n_ties = int((tied > 1).sum())

    return pd.Series(winners["probe"].values,
                     index=winners["gene"].values).sort_index(), n_ties


def aggregate(
    expr: pd.DataFrame,
    gene_map: pd.Series,
    strategy: str,
    selection_columns: list[str] | None = None,
) -> AggregationResult:
    """Collapse a probe x sample matrix to genes x samples.

    Parameters
    ----------
    expr
        Probe-level matrix (probes x samples), values exactly as deposited.
    gene_map
        probe_id -> gene_symbol for UNIQUELY mapped probes only.
    strategy
        One of STRATEGIES.
    selection_columns
        Samples used to SCORE probes for representative strategies. Defaults
        to every sample. The returned matrix always contains every sample
        regardless. Ignored by median/mean, which are per-sample operations
        with no selection step.
    """
    if strategy not in STRATEGIES:
        raise ValueError(f"unknown strategy {strategy!r}; expected {STRATEGIES}")

    missing = gene_map.index.difference(expr.index)
    if len(missing):
        raise ValueError(
            f"{len(missing)} mapped probes absent from the expression matrix "
            f"(first few: {list(missing[:5])})")

    sub = expr.loc[gene_map.index]

    if strategy.startswith("representative_"):
        criterion = "iqr" if strategy.endswith("iqr") else "mean"
        rep, n_ties = _pick_representative(expr, gene_map, criterion,
                                           selection_columns)
        matrix = expr.loc[rep.values].copy()
        matrix.index = pd.Index(rep.index, name="gene_symbol")
        matrix = matrix.sort_index()
        return AggregationResult(
            matrix=matrix, strategy=strategy, n_genes=len(matrix),
            n_probes_used=len(rep), representative_probe=rep,
            diagnostics={
                "criterion": criterion,
                "n_genes_with_tied_score": n_ties,
                "selection_samples": (len(expr.columns)
                                      if selection_columns is None
                                      else len(selection_columns)),
                "selection_scope": ("all samples" if selection_columns is None
                                    else "subset"),
            },
        )

    grouped = sub.groupby(gene_map.values, sort=True)
    matrix = grouped.median() if strategy == "median" else grouped.mean()
    matrix.index = pd.Index(matrix.index, name="gene_symbol")
    return AggregationResult(
        matrix=matrix, strategy=strategy, n_genes=len(matrix),
        n_probes_used=len(gene_map),
        diagnostics={"note": f"{strategy} across all probes of each gene"},
    )
