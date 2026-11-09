"""Positive-control co-expression QC.

Pairs of genes that are strongly co-expressed in essentially every bulk solid-tumour dataset (collagen triplet,
proliferation, T-cell receptor complex, intestinal epithelium). A dataset in which these pairs are uncorrelated has a
processing or labelling defect (for example probe values misaligned to probe IDs) and gene-level interpretation is unsafe.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

POSITIVE_CONTROL_PAIRS = [("COL1A1", "COL1A2"), ("COL1A1", "COL3A1"), ("COL1A2", "COL3A1"), ("SPARC", "COL5A1"),
                          ("MKI67", "TOP2A"), ("BUB1", "CCNB1"), ("CD2", "CD3E"), ("CDX2", "VIL1")]
PASS_MEDIAN_R = 0.5


def coexpression_qc(expr: pd.DataFrame, pairs=POSITIVE_CONTROL_PAIRS, threshold: float = PASS_MEDIAN_R) -> dict:
    rs = {}
    for a, b in pairs:
        if a in expr.index and b in expr.index:
            rs[f"{a}-{b}"] = float(np.corrcoef(expr.loc[a], expr.loc[b])[0, 1])
    med = float(np.median(list(rs.values()))) if rs else float("nan")
    return {"pairs": rs, "n_pairs": len(rs), "median_r": med, "threshold": threshold,
            "passed": bool(len(rs) >= 4 and med >= threshold),
            "note": "fewer than 4 control pairs present: QC not evaluable" if len(rs) < 4 else ""}
