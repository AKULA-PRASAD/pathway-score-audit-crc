"""Unit tests for analysis/stats_core.py -- synthetic inputs with hand-derived answers.

Run: .venv/bin/python -B -m pytest -p no:cacheprovider reports/phase9/final_project/validation/tests -q
"""
import math, sys
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analysis"))
import stats_core as sc  # noqa: E402

N = NormalDist()


# ---------------- quantiles / percentile interval
def test_quantile_type7_hand_value():
    # type 7: h = (n-1)p = 0.75 -> 1 + 0.75*(2-1) = 1.75
    assert sc.quantile([1, 2, 3, 4], 0.25) == pytest.approx(1.75, abs=0)


def test_percentile_ci_exact_on_grid():
    x = np.arange(101.0)            # (n-1)p = 2.5 and 97.5 -> exact
    assert sc.percentile_ci(x) == (2.5, 97.5)


def test_percentile_ci_order_invariant():
    rng = np.random.default_rng(1); x = rng.normal(size=500)
    assert sc.percentile_ci(x) == sc.percentile_ci(rng.permutation(x))


# ---------------- basic interval: direction of the reflection
def test_basic_ci_reflection_direction():
    x = np.arange(101.0)            # q.025 = 2.5, q.975 = 97.5
    lo, hi = sc.basic_ci(x, theta=10.0)
    assert (lo, hi) == (2 * 10 - 97.5, 2 * 10 - 2.5)
    assert lo < hi


def test_basic_equals_percentile_when_symmetric_about_theta():
    x = np.concatenate([-np.arange(1, 251.0), np.arange(1, 251.0)])  # symmetric about 0
    p = sc.percentile_ci(x); b = sc.basic_ci(x, 0.0)
    assert b[0] == pytest.approx(p[0]) and b[1] == pytest.approx(p[1])


def test_basic_moves_opposite_to_skew():
    # right-skewed bootstrap distribution: basic interval shifts LEFT of percentile
    rng = np.random.default_rng(2); x = rng.exponential(size=5000)
    theta = float(np.median(x))
    p = sc.percentile_ci(x); b = sc.basic_ci(x, theta)
    assert b[1] < p[1] and b[0] < p[0]


# ---------------- z0
def test_z0_centre_ties_and_boundary():
    z = sc.z0([1, 2, 3, 4], 2.5)
    assert z["z0"] == pytest.approx(0.0) and z["n_less"] == 2 and z["n_tie"] == 0
    t = sc.z0([1, 2, 3, 4], 2.0)        # strict '<': one below, one tie
    assert t["n_less"] == 1 and t["n_tie"] == 1 and t["z0"] == pytest.approx(N.inv_cdf(0.25))
    b = sc.z0([1, 2, 3, 4], 0.0)
    assert b["boundary"] and b["z0"] == -math.inf


def test_z0_se_formula_matches_plan_value():
    # plan 7.1: SE(z0) = sqrt(p(1-p)/B)/phi(z0) = 0.056 at B_out = 500, z0 = 0
    x = np.arange(500.0); z = sc.z0(x, 249.5)
    assert z["p_less"] == 0.5
    assert z["SE_z0_plan_7_1"] == pytest.approx(0.0560, abs=5e-4)


def test_bc_interval_reduces_to_percentile_when_z0_zero():
    x = np.arange(500.0)
    bc = sc.bc_z0_ci(x, 249.5); p = sc.percentile_ci(x)
    assert bc["alpha1"] == pytest.approx(0.025) and bc["alpha2"] == pytest.approx(0.975)
    assert bc["lower"] == pytest.approx(p[0]) and bc["upper"] == pytest.approx(p[1])


def test_bc_interval_shifts_down_when_z0_negative():
    x = np.arange(500.0); theta = 150.0      # p_less = 0.3 -> z0 < 0
    bc = sc.bc_z0_ci(x, theta)
    assert bc["alpha1"] < 0.025 and bc["alpha2"] < 0.975
    z = sc.z0(x, theta)["z0"]
    assert bc["alpha1"] == pytest.approx(N.cdf(2 * z + N.inv_cdf(0.025)))


# ---------------- bandwidth, KDE, kernel SE
def test_silverman_rule_hand_value():
    rng = np.random.default_rng(3); x = rng.normal(size=500)
    s = np.std(x, ddof=1); iqr = np.quantile(x, .75) - np.quantile(x, .25)
    assert sc.silverman_bandwidth(x) == pytest.approx(0.9 * min(s, iqr / 1.34) * 500 ** -0.2)


def test_kde_matches_scipy_at_same_bandwidth():
    rng = np.random.default_rng(4); x = rng.normal(size=500); h = sc.silverman_bandwidth(x)
    k = stats.gaussian_kde(x, bw_method=h / np.std(x, ddof=1))
    for t in (-2.0, 0.0, 1.7):
        assert sc.gaussian_kde_at(x, t, h) == pytest.approx(float(k(t)[0]), rel=1e-10)


def test_scipy_silverman_variant_matches_scipy():
    rng = np.random.default_rng(5); x = rng.normal(size=500)
    k = stats.gaussian_kde(x, bw_method="silverman")
    assert sc.scipy_silverman_bandwidth(x) == pytest.approx(k.factor * np.std(x, ddof=1), rel=1e-12)


def test_kernel_quantile_se_large_sample_theory():
    rng = np.random.default_rng(6); x = rng.normal(size=200_000); p = 0.025
    q = sc.quantile(x, p); se = sc.kernel_quantile_se(x, p, q, sc.silverman_bandwidth(x))["SE"]
    theory = math.sqrt(p * (1 - p) / x.size) / N.pdf(N.inv_cdf(p))
    assert se == pytest.approx(theory, rel=0.03)


def test_kernel_quantile_se_calibration_by_simulation():
    """At B = 500 the kernel SE of the 2.5 % quantile should approximate the empirical
    SD of that quantile over repeated samples (normal parent)."""
    rng = np.random.default_rng(7); reps = 2000; qs, ses = [], []
    for _ in range(reps):
        x = rng.normal(size=500); q = sc.quantile(x, 0.025)
        qs.append(q); ses.append(sc.kernel_quantile_se(x, 0.025, q, sc.silverman_bandwidth(x))["SE"])
    emp = float(np.std(qs, ddof=1)); med = float(np.median(ses))
    assert med == pytest.approx(emp, rel=0.15)


# ---------------- bounding analysis (6.4)
def _rec(n_pairs, seed=0):
    rng = np.random.default_rng(seed)
    m1 = np.c_[rng.uniform(.6, .7, n_pairs), rng.uniform(.6, .65, n_pairs)]
    m2 = np.c_[rng.uniform(.6, .75, n_pairs), rng.uniform(.6, .65, n_pairs)]
    d = (m2[:, 0] - m2[:, 1]) - (m1[:, 0] - m1[:, 1])
    dC_app = 0.05
    return {"inner_C": {"M1": m1.tolist(), "M2": m2.tolist()},
            "dC_app": dC_app, "dOpt": float(d.mean()), "dC_corrected": dC_app - float(d.mean())}, d


def test_bounding_vacuous_when_complete():
    r, _ = _rec(100)
    b = sc.bounding_record(r, 100)
    assert b["n_missing"] == 0 and b["width"] == 0.0 and b["dC_min"] == b["dC_max"] == r["dC_corrected"]


def test_bounding_hand_calculation_with_missing():
    r, d = _rec(97, seed=1)
    b = sc.bounding_record(r, 100)
    assert b["n_missing"] == 3
    assert b["dC_max"] == pytest.approx(0.05 - (d.sum() + 3 * d.min()) / 100)
    assert b["dC_min"] == pytest.approx(0.05 - (d.sum() + 3 * d.max()) / 100)
    assert b["dC_min"] <= b["dC_max"]
    assert b["width"] == pytest.approx(3 * (d.max() - d.min()) / 100)


# ---------------- inner noise (1.2 / W3)
def test_inner_noise_recovers_known_variance():
    rng = np.random.default_rng(8); recs = []; sd_d = 0.01
    for k in range(400):
        d = rng.normal(0, sd_d, 100)
        m1 = np.c_[np.full(100, .6), np.full(100, .6)]
        m2 = np.c_[.6 + d, np.full(100, .6)]
        recs.append({"inner_C": {"M1": m1.tolist(), "M2": m2.tolist()},
                     "dC_corrected": rng.normal(0, 0.006) - d.mean()})
    out = sc.inner_noise(recs)
    assert out["sigma_Delta_per_inner"] == pytest.approx(sd_d, rel=0.02)
    assert out["inner_MC_SE_dOpt"] == pytest.approx(sd_d / 10, rel=0.02)


# ---------------- running trace
def test_running_trace_uses_first_B_in_k_order():
    x = np.arange(500.0)
    tr = sc.running_trace(x)
    assert tr[50] == sc.percentile_ci(x[:50]) and tr[500] == sc.percentile_ci(x)
