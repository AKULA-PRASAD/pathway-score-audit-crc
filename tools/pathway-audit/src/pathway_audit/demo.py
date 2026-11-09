"""Built-in demonstration on synthetic data with known ground truth (no downloads).

Planted truth:
* STROMAL_PROGRAMME: driven by a tissue 'stroma' factor → expected verdict "composition-dominated".
* TARGET_PROGRAMME: lowered directly by a mutation → genotype effect "robust to adjustment".
* CONFOUNDED_PROGRAMME: driven by an MSI-like factor that is correlated with the mutation → "explained by adjustment covariates".
* RANDOM_SET: unrelated genes → "incoherent".
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def synthetic_cohort(n: int = 300, seed: int = 7):
    rng = np.random.default_rng(seed)
    stroma = rng.normal(size=n)
    immune = rng.normal(size=n)
    msi = (rng.random(n) < 0.2)
    mutation = np.where(msi, rng.random(n) < 0.7, rng.random(n) < 0.15)
    hit = rng.random(n) < 0.4
    rows = {}
    for g in ["COL1A1", "COL1A2", "COL3A1", "SPARC", "COL5A1"] + [f"STR{i}" for i in range(60)]:
        rows[g] = 8 + stroma + 0.4 * rng.normal(size=n)
    prolif = rng.normal(size=n)
    for g in ["MKI67", "TOP2A", "BUB1", "CCNB1"]:
        rows[g] = 7 + prolif + 0.4 * rng.normal(size=n)
    for g in ["CD2", "CD3E"] + [f"IMM{i}" for i in range(30)]:
        rows[g] = 5 + immune + 0.5 * rng.normal(size=n)
    epi = rng.normal(size=n)
    for g in ["CDX2", "VIL1"]:
        rows[g] = 9 + epi + 0.4 * rng.normal(size=n)
    for i in range(40):
        rows[f"TGT{i}"] = 6 - 1.0 * hit + 0.7 * rng.normal(size=n)
    for i in range(40):
        rows[f"CNF{i}"] = 6 + 1.2 * msi + 0.7 * rng.normal(size=n)
    for i in range(400):
        rows[f"G{i}"] = 5 + rng.normal(size=n)
    expr = pd.DataFrame(rows).T
    expr.columns = [f"S{i:03d}" for i in range(n)]
    ann = pd.DataFrame({"target_hit": np.where(hit, "M", "WT"), "mutation": np.where(mutation, "M", "WT"),
                        "msi": np.where(msi, "MSI", "MSS")}, index=expr.columns)
    comp = pd.DataFrame({"StromalScore": stroma + 0.3 * rng.normal(size=n), "ImmuneScore": immune + 0.3 * rng.normal(size=n)}, index=expr.columns)
    sets = {"STROMAL_PROGRAMME": [f"STR{i}" for i in range(40)] + ["COL1A1", "COL3A1"],
            "TARGET_PROGRAMME": [f"TGT{i}" for i in range(40)],
            "CONFOUNDED_PROGRAMME": [f"CNF{i}" for i in range(40)],
            "RANDOM_SET": [f"G{i}" for i in range(40)]}
    return expr, ann, comp, sets


def run_demo(out, seed: int = 7):
    from . import report
    from .audit import audit
    from .report_html import write_html
    expr, ann, comp, sets = synthetic_cohort(seed=seed)
    res = audit(expr, sets, ann, binary={"Target mutation": ("target_hit", "M", "WT"), "Driver mutation": ("mutation", "M", "WT")},
                adjust_for=["msi=MSI"], composition=comp, fibroblast_markers=["COL1A1", "COL1A2", "COL3A1", "SPARC", "COL5A1"],
                n_boot=300, n_null=100, module_boot=30, seed=seed)
    report.write(res, out, "pathway-audit demo (synthetic data with planted truth)")
    return res, write_html(res, out, "pathway-audit demo (synthetic data with planted truth)")
