"""Simulation machinery for external-validation precision assessment.

ENTIRELY SYNTHETIC. No project data is read, no predictive model is fitted to
real data, and nothing here touches the Phase 3 bundle. The simulated linear
predictor is drawn from an assumed distribution -- it is NOT derived from
GSE39582 expression.

Framework follows Riley et al. (Stat Med 2022;41:1280-1295), which targets a
confidence-interval width for the CALIBRATION SLOPE, that being the measure
that typically needs the largest sample size. Assumptions required by that
framework -- the linear-predictor distribution and the event/censoring
distributions -- are explicit parameters here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Scenario:
    """Explicit assumptions. Every one is a choice, not a fact."""
    sigma_lp: float        # SD of the linear predictor (model strength)
    event_fraction: float  # proportion experiencing the event by end of follow-up
    name: str = ""


def simulate(n: int, sigma_lp: float, event_fraction: float,
             rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """One synthetic validation cohort.

    Exponential baseline hazard, proportional hazards in a normal linear
    predictor, administrative censoring chosen to hit `event_fraction`.
    Returns (time, event, lp).
    """
    lp = rng.normal(0.0, sigma_lp, n)
    # Exponential survival with rate exp(lp); baseline rate folded into tau.
    u = rng.uniform(size=n)
    t_event = -np.log(u) / np.exp(lp)
    # Administrative censoring time chosen so that ~event_fraction are observed.
    tau = np.quantile(t_event, event_fraction)
    time = np.minimum(t_event, tau)
    event = (t_event <= tau).astype(float)
    return time, event, lp


def cox_univariate(time: np.ndarray, event: np.ndarray, x: np.ndarray,
                   tol: float = 1e-10, max_iter: int = 100
                   ) -> tuple[float, float]:
    """Univariate Cox partial-likelihood fit (Breslow ties). Returns (beta, se).

    Fitting the outcome on the linear predictor gives the CALIBRATION SLOPE:
    1.0 under perfect calibration.
    """
    order = np.argsort(time)
    t, d, z = time[order], event[order], x[order]
    beta = 0.0
    for _ in range(max_iter):
        w = np.exp(beta * z)
        # risk sets: reverse cumulative sums over time-sorted data
        s0 = np.cumsum(w[::-1])[::-1]
        s1 = np.cumsum((w * z)[::-1])[::-1]
        s2 = np.cumsum((w * z * z)[::-1])[::-1]
        ev = d == 1
        if not ev.any():
            return np.nan, np.nan
        mu = s1[ev] / s0[ev]
        grad = np.sum(z[ev] - mu)
        info = np.sum(s2[ev] / s0[ev] - mu ** 2)
        if info <= 0 or not np.isfinite(info):
            return np.nan, np.nan
        step = grad / info
        beta += step
        if abs(step) < tol:
            break
    return beta, float(np.sqrt(1.0 / info))


def harrell_c(time: np.ndarray, event: np.ndarray, score: np.ndarray,
              rng: np.random.Generator | None = None,
              max_pairs: int = 2_000_000) -> float:
    """Harrell's concordance. Higher score = higher risk = shorter survival."""
    n = len(time)
    if n * (n - 1) // 2 > max_pairs and rng is not None:
        i = rng.integers(0, n, max_pairs)
        j = rng.integers(0, n, max_pairs)
        keep = i != j
        i, j = i[keep], j[keep]
    else:
        i, j = np.triu_indices(n, k=1)
    ti, tj = time[i], time[j]
    di, dj = event[i], event[j]
    si, sj = score[i], score[j]
    # comparable if the earlier time is an event
    comp = ((ti < tj) & (di == 1)) | ((tj < ti) & (dj == 1))
    if not comp.any():
        return np.nan
    i, j = i[comp], j[comp]
    ti, tj, si, sj = time[i], time[j], score[i], score[j]
    earlier_i = ti < tj
    conc = np.where(earlier_i, si > sj, sj > si)
    ties = si == sj
    return float((conc.sum() + 0.5 * ties.sum()) / len(i))


def n_for_events(target_events: int, event_fraction: float) -> int:
    """Sample size whose expected event count is `target_events`."""
    return int(round(target_events / event_fraction))
