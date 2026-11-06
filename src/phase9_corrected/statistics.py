"""Interval statistics, EP-1 (both options), W5 (both actions) and ladder assessment.

Conventions follow documentation/STATISTICAL_CONVENTIONS.md (type-7 quantiles, 95 %,
Silverman 1986 eq. 3.31, Gaussian KDE at the endpoint, between-run SD ddof = 1).
Owner choices (EP-1 option, W5 action, ladder/stopping rule) are REQUIRED arguments;
None raises UnresolvedDecisionError.
"""
from __future__ import annotations

import math
from typing import Mapping, Optional

import numpy as np
from scipy.stats import norm

from .config import (BEYOND_LADDER_RULES, EP1_OPTIONS, PLAN_EP1_FRACTION, W5_ACTIONS,
                     W5_CORRECTIVE_RULES, W5_EXCLUSION_RULES, UnresolvedDecisionError)

ALPHA = 0.05
BW_MULT = (0.5, 1.0, 2.0)                                       # plan line 106-107
TRACE_GRID = (50, 100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000)   # [PROP] D5


class OwnerReviewRequired(RuntimeError):
    """A rule returned stop_for_owner_review: no B_out change may be derived (Decision 6B)."""


def quantile(x, p) -> float:
    return float(np.quantile(np.asarray(x, dtype=np.float64), p, method="linear"))


def percentile_ci(x, alpha=ALPHA):
    return quantile(x, alpha / 2), quantile(x, 1 - alpha / 2)


def basic_ci(x, theta, alpha=ALPHA):
    lo, hi = percentile_ci(x, alpha)
    return 2 * theta - hi, 2 * theta - lo


def z0(x, theta) -> float:
    x = np.asarray(x, dtype=np.float64)
    p = float((x < theta).sum()) / x.size                      # strict '<' (convention 4)
    if p in (0.0, 1.0):
        return math.copysign(math.inf, p - 0.5)
    return float(norm.ppf(p))


def z0_counts(x, theta) -> dict:
    """Convention 4: ties with theta are counted and reported, never split."""
    x = np.asarray(x, dtype=np.float64)
    return {"n_less": int((x < theta).sum()), "n_tie": int((x == theta).sum()), "B": int(x.size)}


def silverman_bandwidth(x) -> float:
    """Silverman (1986) eq. 3.31: 0.9 * min(s, IQR/1.34) * n^(-1/5)."""
    x = np.asarray(x, dtype=np.float64)
    s = float(np.std(x, ddof=1))
    iqr = quantile(x, 0.75) - quantile(x, 0.25)
    h = 0.9 * min(s, iqr / 1.34) * x.size ** (-0.2)
    if not h > 0:
        raise ValueError("Silverman bandwidth is zero (degenerate sample: SD or IQR is 0)")
    return h


def gaussian_kde_at(x, t, h) -> float:
    x = np.asarray(x, dtype=np.float64)
    u = (t - x) / h
    return float(np.exp(-0.5 * u * u).sum() / (x.size * h * math.sqrt(2 * math.pi)))


def kernel_quantile_se(x, p, q, h) -> float:
    """Plan section 2.2 (line 105): sqrt(p(1-p)/B) / f_hat(q_p)."""
    return math.sqrt(p * (1 - p) / len(x)) / gaussian_kde_at(x, q, h)


def endpoint_mcse(x, alpha=ALPHA) -> dict:
    """Per-run kernel MC-SE of both percentile endpoints at 0.5h, h, 2h."""
    lo, hi = percentile_ci(x, alpha)
    h = silverman_bandwidth(x)
    out = {"h": h, "width": hi - lo, "lower": {"q": lo, "p": alpha / 2},
           "upper": {"q": hi, "p": 1 - alpha / 2}}
    for e in ("lower", "upper"):
        out[e]["SE"] = {m: kernel_quantile_se(x, out[e]["p"], out[e]["q"], h * m) for m in BW_MULT}
    return out


def ep1_evaluate(run_values: Mapping[int, np.ndarray], option: Optional[str],
                 fraction=PLAN_EP1_FRACTION) -> dict:
    """EP-1: MC-SE(endpoint) <= 0.05 x width, every run, both endpoints (plan line 84).
    option A: ratio uses the Silverman bandwidth h only; B: max over {0.5h, h, 2h}.
    The between-run SD is NOT part of EP-1 under either option (it is reported separately)."""
    if option not in EP1_OPTIONS:
        raise UnresolvedDecisionError("EP-1 option (A or B) not supplied")
    rows = []
    for r, v in sorted(run_values.items()):
        t = endpoint_mcse(v)
        for e in ("lower", "upper"):
            ses = t[e]["SE"]
            se = ses[1.0] if option == "A" else max(ses.values())
            rho = se / t["width"]
            rows.append({"run": r, "endpoint": e, "SE_by_bw": {str(k): s for k, s in ses.items()},
                         "governing_SE": se, "width": t["width"], "ratio": rho,
                         "pass": rho <= fraction})
    rho_max = max(x["ratio"] for x in rows)
    return {"option": option, "fraction": fraction, "rows": rows, "rho_max": rho_max,
            "pass": all(x["pass"] for x in rows), "B_out": {r: len(v) for r, v in run_values.items()}}


def ladder_assessment(B_cur: int, ep1: dict, ladder: Optional[dict]) -> dict:
    """Next size under the DRAFT sizing rule B_req = B_cur * (rho_max / 0.05)^2, rounded up
    to the smallest rung >= B_req and > B_cur ([PLAN 118-119] at B_cur = 500; [PROP] above).
    The beyond-ladder behaviour is an owner decision. Returns a decision; executes nothing."""
    if not ladder or ladder.get("beyond_ladder_rule") not in BEYOND_LADDER_RULES:
        raise UnresolvedDecisionError("ladder and beyond_ladder_rule not supplied")
    if ep1["pass"]:
        return {"action": "stop_ep1_pass", "B_next": None}
    rungs = ladder["rungs"]
    if B_cur >= rungs[-1]:
        return {"action": "stop_label_imprecise", "B_next": None, "basis": "H3, plan line 456"}
    B_req = B_cur * (ep1["rho_max"] / ep1["fraction"]) ** 2
    cand = [L for L in rungs if L >= B_req and L > B_cur]
    if cand:
        return {"action": "top_up", "B_req": B_req, "B_next": cand[0]}
    if ladder["beyond_ladder_rule"] == "top_up_to_max_rung_then_reassess":
        return {"action": "top_up", "B_req": B_req, "B_next": rungs[-1],
                "note": "B_req exceeds the top rung"}
    return {"action": "stop_label_imprecise", "B_req": B_req, "B_next": None,
            "note": "B_req exceeds the top rung"}


def w5_grid_points(B_out: int, grid=TRACE_GRID) -> list:
    return [B for B in grid if 0.6 * B_out <= B <= B_out]


def w5_evaluate(values_by_k: Mapping[int, Optional[float]], B_out: int,
                width: Optional[float] = None, fraction=PLAN_EP1_FRACTION, grid=TRACE_GRID,
                exclusion_rule: Optional[str] = None,
                h1_triggered: Optional[bool] = None) -> dict:
    """[PROP D5] spread = max - min of the type-7 percentile endpoint on the first B
    replicates in k order, B in grid ∩ [0.6 B_out, B_out]; breach if spread > 0.05 * width.

    exclusion_rule=None (default): unchanged behaviour -- `width` is required and any
    excluded replicate (value None) is refused, the trace population being undefined.
    exclusion_rule='usable_within_k_positions_v1' (owner Decision 6B, 2026-11-06): see
    `_w5_usable_within_k_positions`."""
    ks = sorted(values_by_k)
    if ks != list(range(B_out)):
        raise ValueError("values_by_k must cover k = 0..B_out-1")
    if exclusion_rule is not None:
        if exclusion_rule not in W5_EXCLUSION_RULES:
            raise UnresolvedDecisionError(f"unknown W5 exclusion rule {exclusion_rule!r}")
        return _w5_usable_within_k_positions(values_by_k, B_out, width, fraction, grid,
                                             h1_triggered)
    if width is None:
        raise ValueError("width is required when no exclusion rule is supplied")
    if any(values_by_k[k] is None for k in ks):
        raise UnresolvedDecisionError("W5 trace with excluded outer replicates is not defined")
    x = np.array([values_by_k[k] for k in ks], dtype=np.float64)
    pts = w5_grid_points(B_out, grid)
    if len(pts) < 2:
        raise ValueError(f"fewer than two trace points for B_out = {B_out}")
    tr = {B: percentile_ci(x[:B]) for B in pts}
    spread = {"lower": max(v[0] for v in tr.values()) - min(v[0] for v in tr.values()),
              "upper": max(v[1] for v in tr.values()) - min(v[1] for v in tr.values())}
    bound = fraction * width
    return {"points": pts, "trace": {str(B): list(v) for B, v in tr.items()}, "spread": spread,
            "bound": bound, "breach": max(spread.values()) > bound}


def _w5_usable_within_k_positions(values_by_k, B_out, width, fraction, grid,
                                  h1_triggered) -> dict:
    """Owner Decision 6B rule 'usable_within_k_positions_v1':
      1. at grid point B, endpoints use only the usable values at original positions k < B
         (excluded = None); positions are kept, nothing is back-filled;
      2. n_usable is reported at every grid point;
      3. bound = fraction x percentile width of the usable replicates at full B_out (the
         primary interval); a caller-supplied width must equal it exactly;
      4. UNEVALUABLE (breach None) if H1 is triggered or any grid point -- or the full
         population -- has fewer than two usable replicates. No substitute is computed."""
    if not isinstance(h1_triggered, bool):
        raise ValueError("h1_triggered (bool) is required under usable_within_k_positions_v1")
    pts = w5_grid_points(B_out, grid)
    if len(pts) < 2:
        raise ValueError(f"fewer than two trace points for B_out = {B_out}")
    usable = {B: [values_by_k[k] for k in range(B) if values_by_k[k] is not None] for B in pts}
    full = [values_by_k[k] for k in range(B_out) if values_by_k[k] is not None]
    n_usable = {str(B): len(v) for B, v in usable.items()}
    out = {"exclusion_rule": "usable_within_k_positions_v1", "points": pts,
           "n_usable": n_usable, "n_usable_full": len(full), "h1_triggered": h1_triggered}
    reasons = []
    if h1_triggered:
        reasons.append("H1 triggered (> 5 % of outer replicates excluded, plan 6.5)")
    short = [B for B in pts if len(usable[B]) < 2]
    if short or len(full) < 2:
        reasons.append(f"fewer than two usable replicates at grid point(s) {short or [B_out]}")
    if reasons:
        return {**out, "status": "UNEVALUABLE", "breach": None, "reasons": reasons}
    lo, hi = percentile_ci(full)
    w_usable = hi - lo
    if width is not None and width != w_usable:
        raise ValueError(f"caller-supplied width {width!r} conflicts with the usable-population "
                         f"width {w_usable!r} required by usable_within_k_positions_v1")
    tr = {B: percentile_ci(np.array(usable[B], dtype=np.float64)) for B in pts}
    spread = {"lower": max(v[0] for v in tr.values()) - min(v[0] for v in tr.values()),
              "upper": max(v[1] for v in tr.values()) - min(v[1] for v in tr.values())}
    bound = fraction * w_usable
    return {**out, "status": "EVALUATED", "width": w_usable,
            "trace": {str(B): list(v) for B, v in tr.items()}, "spread": spread,
            "bound": bound, "breach": max(spread.values()) > bound}


def w5_decision(w5: dict, action: Optional[str], corrective_rule: Optional[str],
                B_cur: int, ladder: Optional[dict]) -> dict:
    if action not in W5_ACTIONS:
        raise UnresolvedDecisionError("W5 action (W5-A or W5-B) not supplied")
    if w5.get("status") == "UNEVALUABLE":                     # Decision 6B, both actions
        return {"action": "stop_for_owner_review", "B_next": None}
    if not w5["breach"]:
        return {"action": "none"}
    if action == "W5-B":
        return {"action": "investigate_and_report", "B_next": None}
    if corrective_rule not in W5_CORRECTIVE_RULES:
        raise UnresolvedDecisionError("W5-A corrective rule not supplied")
    if not ladder:
        raise UnresolvedDecisionError("ladder not supplied")
    nxt = [L for L in ladder["rungs"] if L > B_cur]
    if not nxt:
        return {"action": "report_persisting_breach", "B_next": None}
    return {"action": "top_up", "B_next": nxt[0]}


def combined_next_size(ep1_dec: dict, w5_dec: dict) -> Optional[int]:
    """Under W5-A the target is the larger of the EP-1 and W5 targets ([PROP] D5).
    If either decision is stop_for_owner_review, no size is derived (Decision 6B)."""
    for d in (ep1_dec, w5_dec):
        if d.get("action") == "stop_for_owner_review":
            raise OwnerReviewRequired("W5/EP-1 sizing requires owner review; no top-up may proceed")
    c = [d.get("B_next") for d in (ep1_dec, w5_dec) if d.get("B_next")]
    return max(c) if c else None


def run_summary(values, theta_hat: float) -> dict:
    v = np.asarray(values, dtype=np.float64)
    lo, hi = percentile_ci(v)
    blo, bhi = basic_ci(v, theta_hat)
    w = hi - lo
    m = 0.5 * (lo + hi)
    return {"n": int(v.size), "percentile": [lo, hi], "width": w, "basic": [blo, bhi],
            "z0": z0(v, theta_hat), "A7": 2 * abs(m - theta_hat) / w,
            "A8": abs(z0(v, theta_hat)), "z0_counts": z0_counts(v, theta_hat),
            "mcse": endpoint_mcse(v)}


def between_run(endpoints: Mapping[int, tuple]) -> dict:
    out = {}
    for i, e in enumerate(("lower", "upper")):
        vals = [endpoints[r][i] for r in sorted(endpoints)]
        out[e] = {"values": vals, "sd_ddof1": float(np.std(vals, ddof=1)),
                  "range": float(max(vals) - min(vals)), "df": len(vals) - 1}
    return out
