"""Shared Phase 10 utilities: configuration, seeds, folds, transformers, metrics, hashing.

Phase 9 frozen code is imported, never modified. Every data-dependent transformer here is
fitted only on the rows it is given, so placing it inside a scikit-learn Pipeline that is tuned
by GridSearchCV gives strict nesting: each inner-training fold refits it.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "phase10" / "phase10.yaml"
REPORTS = ROOT / "reports" / "phase10"
TCGA_DIR = ROOT / "data" / "external" / "tcga_coad_xena"
GSE_EXPR = ROOT / "data" / "processed" / "analysis" / "expression_median_S0.parquet"

TRACK_CODE = {"mmr": 1, "recurrence": 2}


def load_config() -> dict:
    cfg = yaml.safe_load(CONFIG.read_text())
    return cfg


def grid(spec):
    """Expand {logspace: [a, b, n]} / {geomspace: [a, b, n]} or return a list unchanged."""
    if isinstance(spec, dict):
        (kind, (a, b, n)), = spec.items()
        return list(getattr(np, kind)(a, b, int(n)))
    return list(spec)


def seed(*parts: int) -> int:
    """Deterministic 32-bit seed from integer parts (numpy SeedSequence)."""
    return int(np.random.SeedSequence([int(p) for p in parts]).generate_state(1)[0])


# ---------------------------------------------------------------- folds
def make_folds(ids, strat, track: str, repeats: int, k: int, base: int) -> pd.DataFrame:
    """Patient-level stratified folds for every repeat; returns long table id/repeat/fold."""
    ids = list(ids)
    strat = np.asarray(strat)
    rows = []
    for r in range(repeats):
        skf = StratifiedKFold(k, shuffle=True, random_state=seed(base, TRACK_CODE[track], r))
        for f, (_, te) in enumerate(skf.split(np.zeros(len(ids)), strat)):
            rows += [(ids[i], r, f) for i in te]
    out = pd.DataFrame(rows, columns=["sample_id", "repeat", "fold"])
    return out.sort_values(["repeat", "sample_id"]).reset_index(drop=True)


def inner_cv(base: int, track: str, r: int, f: int, k: int) -> StratifiedKFold:
    """Inner splitter; identical for every model given the same training rows and labels."""
    return StratifiedKFold(k, shuffle=True, random_state=seed(base, TRACK_CODE[track], r, f, 99))


def frame_sha(df: pd.DataFrame) -> str:
    return hashlib.sha256(df.to_csv(index=False).encode()).hexdigest()


# ---------------------------------------------------------------- transformers
class VarianceFilter(BaseEstimator, TransformerMixin):
    """Keep the top-k columns by variance (ddof=1) of the TRAINING rows passed to fit.

    Outcome-free (B-02 rule). Ties are broken by column order (stable sort).
    """

    def __init__(self, k: int = 2000):
        self.k = k

    def fit(self, X, y=None):
        X = pd.DataFrame(X)
        var = X.var(axis=0, ddof=1).to_numpy()
        order = np.argsort(-var, kind="stable")[: min(self.k, X.shape[1])]
        self.columns_ = list(X.columns[np.sort(order)])
        self.n_fit_rows_ = X.shape[0]
        return self

    def transform(self, X):
        return pd.DataFrame(X)[self.columns_]

    def get_feature_names_out(self, input_features=None):
        return np.asarray(self.columns_, dtype=object)


class ColumnSelect(BaseEstimator, TransformerMixin):
    """Select fixed, pre-registered columns (no data-dependent fitting)."""

    def __init__(self, columns=()):
        self.columns = columns

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return pd.DataFrame(X)[list(self.columns)]

    def get_feature_names_out(self, input_features=None):
        return np.asarray(list(self.columns), dtype=object)


# ---------------------------------------------------------------- metrics
def auc_pair_counts(y, score):
    """(concordant + 0.5*tied, positive-negative pairs) for one fold."""
    y = np.asarray(y).astype(int)
    s = np.asarray(score, float)
    pos, neg = s[y == 1], s[y == 0]
    diff = pos[:, None] - neg[None, :]
    num = float((diff > 0).sum() + 0.5 * (diff == 0).sum())
    return num, float(len(pos) * len(neg))


def harrell_pair_counts(t, e, risk):
    """(concordant + 0.5*tied-risk, comparable pairs), lifelines 0.30.3 rules.

    A pair (i, j) is comparable when i had the event and T_j > T_i, or when j was censored at
    T_j == T_i. Two events at the same time are not comparable. Higher risk should mean
    shorter time; tied risk scores count 0.5.
    """
    t = np.asarray(t, float)
    e = np.asarray(e).astype(bool)
    r = np.asarray(risk, float)
    later = t[None, :] > t[:, None]
    same_time_censored = (t[None, :] == t[:, None]) & (~e)[None, :]
    comp = e[:, None] & (later | same_time_censored)
    d = r[:, None] - r[None, :]
    num = float(((d > 0) & comp).sum() + 0.5 * ((d == 0) & comp).sum())
    return num, float(comp.sum())


def pooled(counts):
    num = sum(c[0] for c in counts)
    den = sum(c[1] for c in counts)
    return num / den if den else float("nan")


def _logit(p, clip):
    p = np.clip(np.asarray(p, float), clip, 1 - clip)
    return np.log(p / (1 - p))


def calibration(y, p, clip=1e-6):
    """Calibration-in-the-large (intercept with logit(p) offset) and slope, by Newton steps.

    One prediction per patient is required (apply within a single repeat).
    """
    y = np.asarray(y, float)
    o = _logit(p, clip)
    a = 0.0
    for _ in range(100):                      # intercept with offset
        m = 1 / (1 + np.exp(-(a + o)))
        g, h = (y - m).sum(), (m * (1 - m)).sum()
        step = g / h
        a += step
        if abs(step) < 1e-10:
            break
    beta = np.zeros(2)
    Z = np.column_stack([np.ones_like(o), o])
    for _ in range(100):                      # logistic slope on the logit
        m = 1 / (1 + np.exp(-(Z @ beta)))
        g = Z.T @ (y - m)
        H = (Z * (m * (1 - m))[:, None]).T @ Z
        step = np.linalg.solve(H, g)
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return float(a), float(beta[1])


def summarise(values) -> dict:
    v = np.asarray([x for x in values if np.isfinite(x)], float)
    n = len(v)
    sd = float(v.std(ddof=1)) if n > 1 else float("nan")
    return {"mean": float(v.mean()) if n else float("nan"), "sd_across_repeats": sd,
            "mc_se": sd / np.sqrt(n) if n > 1 else float("nan"),
            "min": float(v.min()) if n else float("nan"), "max": float(v.max()) if n else float("nan"),
            "n_repeats": n}


# ---------------------------------------------------------------- provenance
def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_head() -> str:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def environment_record() -> dict:
    import lifelines
    import sklearn
    rec = {"python": sys.version.split()[0], "platform": platform.platform(),
           "numpy": np.__version__, "pandas": pd.__version__, "scikit-learn": sklearn.__version__,
           "lifelines": lifelines.__version__}
    try:
        import sksurv
        rec["scikit-survival"] = sksurv.__version__
    except ImportError:
        rec["scikit-survival"] = None
    return rec


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=False, default=float) + "\n")


def phase9_integrity() -> dict:
    """Recompute the Phase 9 freeze-record and archive-manifest hashes (read-only)."""
    rec = json.loads((ROOT / "config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json").read_text())
    bad_code = [p for p, h in rec["code"].items() if sha256_file(ROOT / p) != h]
    manifest = ROOT / "reports/phase9/ARCHIVE_MANIFEST_SHA256.txt"
    bad_arch, n_arch = [], 0
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        h, p = line.split(None, 1)
        p = p.strip().lstrip("*")
        n_arch += 1
        full = (manifest.parent / p) if not (ROOT / p).exists() else ROOT / p
        if not full.exists() or sha256_file(full) != h:
            bad_arch.append(p)
    return {"frozen_code_files": len(rec["code"]), "frozen_code_mismatches": bad_code,
            "archive_files": n_arch, "archive_mismatches": bad_arch}
