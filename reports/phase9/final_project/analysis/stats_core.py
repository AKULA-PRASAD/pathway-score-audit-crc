"""Pure statistical functions for the Phase 9 section 9.1 analysis (plan v2.1).

No I/O. Each function cites the plan section that defines it. Where the plan leaves a
convention unstated, the convention adopted is named here and in
documentation/STATISTICAL_CONVENTIONS.md; none of them is a methodological choice that the
plan makes differently.
"""
from __future__ import annotations

import math
from statistics import NormalDist

import numpy as np

_N = NormalDist()
QUANTILE_METHOD = "linear"      # Hyndman-Fan type 7 = numpy default; plan 3.4 uses type 7 for knots,
                                # Phase 9 primary p2.5/p97.5 used numpy default (type 7)
ALPHA = 0.05                    # 95 % two-sided; see STATISTICAL_CONVENTIONS.md section 1


def quantile(x, p):
    return float(np.quantile(np.asarray(x, dtype=np.float64), p, method=QUANTILE_METHOD))


# ---------------------------------------------------------------- section 7 intervals
def percentile_ci(x, alpha=ALPHA):
    """Plan section 7 / 3.1: interval = percentile{ dC_corrected(k) } -- PRIMARY."""
    return quantile(x, alpha / 2), quantile(x, 1 - alpha / 2)


def basic_ci(x, theta, alpha=ALPHA):
    """Plan section 7: basic (reverse-percentile) -- DIAGNOSTIC.
    [2*theta - q(1-alpha/2), 2*theta - q(alpha/2)]."""
    lo_q, hi_q = percentile_ci(x, alpha)
    return 2 * theta - hi_q, 2 * theta - lo_q


def z0(x, theta):
    """Plan section 7: z0 = Phi^-1( #{dC_corr(k) < theta} / B_out ).  Strict '<' as written;
    ties with theta are counted and reported separately, never split.
    Returns dict; z0 is +/-inf (flagged) if the proportion is 0 or 1."""
    x = np.asarray(x, dtype=np.float64)
    n_less = int((x < theta).sum()); n_tie = int((x == theta).sum()); B = int(x.size)
    p = n_less / B
    if p <= 0.0 or p >= 1.0:
        z, boundary = (-math.inf if p <= 0 else math.inf), True
    else:
        z, boundary = _N.inv_cdf(p), False
    se = (math.sqrt(p * (1 - p) / B) / _N.pdf(z)) if not boundary else math.nan  # plan 7.1 SE(z0)
    return {"z0": z, "p_less": p, "n_less": n_less, "n_tie": n_tie, "B": B,
            "boundary": boundary, "SE_z0_plan_7_1": se}


def bc_z0_ci(x, theta, alpha=ALPHA):
    """Bias-corrected (z0-only, acceleration a = 0) percentile interval -- DIAGNOSTIC.
    Named in docs/PHASE9_DELTAC_INTERVAL_AUDIT.md line 98 ('bias-corrected (z0-only) variant');
    plan v2.1 section 7 lists z0 as a bias diagnostic and EXCLUDES BCa (no acceleration).
    alpha1 = Phi(2 z0 + z_{alpha/2}),  alpha2 = Phi(2 z0 + z_{1-alpha/2})."""
    zz = z0(x, theta)
    if zz["boundary"]:
        return {"lower": math.nan, "upper": math.nan, "alpha1": math.nan, "alpha2": math.nan,
                "boundary": True}
    a1 = _N.cdf(2 * zz["z0"] + _N.inv_cdf(alpha / 2))
    a2 = _N.cdf(2 * zz["z0"] + _N.inv_cdf(1 - alpha / 2))
    return {"lower": quantile(x, a1), "upper": quantile(x, a2), "alpha1": a1, "alpha2": a2,
            "boundary": False}


# ---------------------------------------------------------------- section 2.2 kernel MC-SE
def silverman_bandwidth(x):
    """Silverman's rule of thumb (Silverman 1986, eq. 3.31):
    h = 0.9 * min(s, IQR/1.34) * n^(-1/5), s = sample SD (ddof=1), IQR by type-7 quantiles."""
    x = np.asarray(x, dtype=np.float64); n = x.size
    s = float(np.std(x, ddof=1))
    iqr = quantile(x, 0.75) - quantile(x, 0.25)
    return 0.9 * min(s, iqr / 1.34) * n ** (-0.2)


def scipy_silverman_bandwidth(x):
    """Cross-check variant only: scipy.stats.gaussian_kde(bw_method='silverman') for d = 1,
    h = s * (n * 3/4)^(-1/5)  (~1.06 s n^-1/5). NOT the primary convention."""
    x = np.asarray(x, dtype=np.float64); n = x.size
    return float(np.std(x, ddof=1)) * (n * 3 / 4) ** (-0.2)


def gaussian_kde_at(x, t, h):
    """Gaussian kernel density estimate f_hat(t) = (1/(n h)) sum phi((t - x_i)/h)."""
    x = np.asarray(x, dtype=np.float64)
    u = (t - x) / h
    return float(np.exp(-0.5 * u * u).sum() / (x.size * h * math.sqrt(2 * math.pi)))


def kernel_quantile_se(x, p, q, h):
    """Plan section 2.2: SE(q_p) = sqrt(p(1-p)/B_out) / f_hat(q_p)."""
    f = gaussian_kde_at(x, q, h)
    return {"SE": math.sqrt(p * (1 - p) / len(x)) / f, "f_hat": f, "h": h}


# ---------------------------------------------------------------- section 6.4 bounding
def bounding_record(rec, B_in):
    """Plan section 6.4: if outer replicate k has failed (non-paired) inner replicates, impute the
    missing PAIRED optimism differences d_b = opt_M2,b - opt_M1,b at the min and at the max observed
    within that outer replicate; return [dC_corr^min(k), dC_corr^max(k)].
    With no missing b the bound collapses to the stored value (vacuous case)."""
    o1 = np.asarray(rec["inner_C"]["M1"], dtype=np.float64)
    o2 = np.asarray(rec["inner_C"]["M2"], dtype=np.float64)
    d = (o2[:, 0] - o2[:, 1]) - (o1[:, 0] - o1[:, 1])
    n_obs = d.size; n_miss = int(B_in - n_obs)
    if n_miss == 0:
        v = rec["dC_corrected"]
        return {"n_missing": 0, "dC_min": v, "dC_max": v, "width": 0.0}
    dopt_lo = (d.sum() + n_miss * d.min()) / B_in      # smallest dOpt -> largest dC_corr
    dopt_hi = (d.sum() + n_miss * d.max()) / B_in
    return {"n_missing": n_miss, "dC_min": rec["dC_app"] - dopt_hi,
            "dC_max": rec["dC_app"] - dopt_lo, "width": (dopt_hi - dopt_lo)}


# ---------------------------------------------------------------- section 1.2 / W3 inner noise
def inner_noise(records):
    """Plan section 1.2: Var(observed outer) = Var(sampling) + sigma_Delta^2 / B_in.
    Per outer k, the inner MC variance of dOpt(k) is var(d_b, ddof=1)/n_b, d_b the paired
    optimism differences. Returns the mean noise variance and the implied inflation of the
    observed variance and SD relative to the sampling component."""
    nv, vd = [], []
    for r in records:
        o1 = np.asarray(r["inner_C"]["M1"]); o2 = np.asarray(r["inner_C"]["M2"])
        d = (o2[:, 0] - o2[:, 1]) - (o1[:, 0] - o1[:, 1])
        vd.append(float(np.var(d, ddof=1)))
        nv.append(vd[-1] / d.size)
    noise = float(np.mean(nv))
    obs = float(np.var([r["dC_corrected"] for r in records], ddof=1))
    samp = obs - noise
    return {"mean_inner_noise_var": noise, "inner_MC_SE_dOpt": math.sqrt(noise),
            "observed_var": obs, "implied_sampling_var": samp,
            "variance_inflation": noise / samp if samp > 0 else math.inf,
            "SD_inflation": math.sqrt(obs / samp) - 1 if samp > 0 else math.inf,
            "sigma_Delta_per_inner": math.sqrt(float(np.mean(vd)))}


# ---------------------------------------------------------------- MC-1 running trace (W5)
MC1_GRID = (50, 100, 200, 300, 400, 500)   # Phase 6B MC-1 grid, as used by Phase 9 primary
MC1_FINAL = (300, 400, 500)                # 'final ~40 %' of B_out = 500


def running_trace(x, grid=MC1_GRID, alpha=ALPHA):
    """Percentile endpoints using the first B replicates in k order (k is the canonical,
    completion-order-independent replicate index)."""
    x = np.asarray(x, dtype=np.float64)
    return {int(B): percentile_ci(x[:B], alpha) for B in grid}
