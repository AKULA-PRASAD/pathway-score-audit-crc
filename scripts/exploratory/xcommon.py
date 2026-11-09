"""Shared, read-only helpers for the registered exploratory analyses (docs/EXPLORATORY_ANALYSIS_PLAN.md).

Reuses frozen project code without modifying it. Nothing here writes to a frozen location.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.pathway_scoring import score_panel  # noqa: E402
from src.phase9_corrected.data import load_context  # noqa: E402
from src.phase9_execution import build_clinical, rcs_knots, AGE_COL  # noqa: E402

OUT = ROOT / "reports" / "exploratory"
EXT = ROOT / "data" / "external" / "exploratory"
ANALYSIS = ROOT / "data" / "processed" / "analysis"
SEED = 20261108
PATHWAYS = [("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT"), ("HALLMARK_E2F_TARGETS", "E2F"),
            ("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "Wnt"), ("HALLMARK_P53_PATHWAY", "p53")]
LABEL = dict(PATHWAYS)


def context():
    """Frozen A-19 analysis cohort (452 patients) and its inputs."""
    return load_context(ROOT)


def cohort_scores(ctx, ids=None):
    """Frozen four-set scores, gene-level standardised on `ids` (default: all 452), then score-standardised on `ids`."""
    ids = list(ctx.ids if ids is None else ids)
    raw, diag = score_panel(ctx.X, ctx.panel, train_ids=ids)
    z = (raw - raw.loc[ids].mean()) / raw.loc[ids].std(ddof=1)
    z = z.rename(columns=LABEL)[[lab for _, lab in PATHWAYS]]
    return z.loc[ctx.ids], diag


def clinical_design(ctx, ids=None):
    ids = list(ctx.ids if ids is None else ids)
    age = pd.to_numeric(ctx.cov.loc[ids, AGE_COL]).to_numpy(float)
    knots, linear = rcs_knots(age)
    return build_clinical(ctx.cov, ids, knots, linear)


def endpoint(ctx):
    ep = pd.read_csv(ANALYSIS / "endpoint.csv").set_index("sample_id").loc[ctx.ids]
    return ep


def outdir(name):
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, obj):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        raise TypeError(type(o))
    Path(path).write_text(json.dumps(obj, indent=2, default=conv, ensure_ascii=False))


def bh(p):
    """Benjamini–Hochberg adjusted q-values."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = n - rank
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def soft_characteristics():
    """All `!Sample_characteristics_ch1` fields from the frozen GEO SOFT file (read-only), one row per GSM."""
    import gzip
    rows, cur = [], None
    with gzip.open(ROOT / "data" / "raw" / "gse39582" / "GSE39582_family.soft.gz", "rt", errors="replace") as f:
        for line in f:
            if line.startswith("^SAMPLE"):
                cur = {"sample_id": line.split("=", 1)[1].strip()}
                rows.append(cur)
            elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
                kv = line.split("=", 1)[1].strip()
                if ":" in kv:
                    k, v = kv.split(":", 1)
                    cur[k.strip()] = v.strip()
    return pd.DataFrame(rows).set_index("sample_id")
