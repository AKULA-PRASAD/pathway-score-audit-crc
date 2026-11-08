"""Phase 9 execution code and data transformations."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.phase9_execution import (AGE_COL, CLINICAL_BINARIES, ConvergenceFailure,  # noqa: E402
                                  REFERENCE_AGE, build_clinical, c_index,
                                  calibration_slope, fit_cox, linear_predictor,
                                  pipeline_fit, rcs_basis, rcs_knots)

ANALYSIS = ROOT / "data" / "processed" / "analysis"
PANELS = ROOT / "config" / "pathways"
EXEC = ROOT / "reports" / "phase9" / "execution"
M1_VARS = [AGE_COL, "sex", "tnm.stage", "tumor.location", "chemotherapy.adjuvant"]


def driver():
    spec = importlib.util.spec_from_file_location(
        "p9exec", ROOT / "scripts" / "38_phase9_execute.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def data():
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    cov = pd.read_csv(ANALYSIS / "covariates.csv", keep_default_na=False,
                      dtype=str).set_index("sample_id")
    ids = [s for s in s3 if all(cov.loc[s, v] != "" for v in M1_VARS)]
    ep = pd.read_csv(ANALYSIS / "endpoint.csv").set_index("sample_id").reindex(ids)
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[ids].astype(np.float64)
    p4 = json.loads((PANELS / "m2_panel_4set.json").read_text())
    return (ids, cov, X, ep["rfs_time_full"].to_numpy(float),
            ep["rfs_event_full"].to_numpy(float), p4)


# ---------- age spline ----------

def test_E01_rcs_basis_matches_harrell_formula():
    k = np.array([40.0, 60.0, 80.0])
    x = np.linspace(20, 100, 50)
    B = rcs_basis(x, k)
    t1, t2, t3 = k
    cu = lambda z: np.maximum(z, 0) ** 3
    want = (cu(x - t1) - cu(x - t2) * (t3 - t1) / (t3 - t2)
            + cu(x - t3) * (t2 - t1) / (t3 - t2)) / (t3 - t1) ** 2
    assert np.allclose(B[:, 0], x)
    assert np.allclose(B[:, 1], want)
    assert B.shape[1] == 2                      # 3 knots -> 2 df


def test_E02_rcs_is_linear_beyond_outer_knots():
    """The defining property of a RESTRICTED cubic spline."""
    k = np.array([40.0, 60.0, 80.0])
    for seg in (np.linspace(10, 39, 12), np.linspace(81, 110, 12)):
        B = rcs_basis(seg, k)
        lp = B @ np.array([0.3, -0.7])
        d2 = np.diff(lp, 2)
        assert np.abs(d2).max() < 1e-8, "spline is not linear in the tail"


def test_E03_knots_are_training_fold_percentiles(data):
    ids, cov, *_ = data
    tr = ids[:300]
    age = pd.to_numeric(cov.loc[tr, AGE_COL]).to_numpy(float)
    k, fallback = rcs_knots(age)
    assert np.allclose(k, np.percentile(age, [10, 50, 90]))   # type 7
    assert fallback is False
    k2, _ = rcs_knots(pd.to_numeric(cov.loc[ids[300:], AGE_COL]).to_numpy(float))
    assert not np.allclose(k, k2), "knots must be fold-specific"


def test_E04_duplicate_knots_fall_back_to_linear_and_are_recorded(data):
    _, cov, *_ = data
    k, fallback = rcs_knots(np.full(50, 70.0))
    assert fallback is True and len(k) < 3


# ---------- clinical design ----------

def test_E05_reference_levels_and_df(data):
    ids, cov, *_ = data
    k, _ = rcs_knots(pd.to_numeric(cov.loc[ids, AGE_COL]).to_numpy(float))
    d = build_clinical(cov, ids, k, linear_age=False)
    assert list(d.columns) == ["age_rcs1", "age_rcs2", "sex=Male", "tnm.stage=3.0",
                               "tumor.location=proximal", "chemotherapy.adjuvant=Y"]
    assert d.shape[1] == 6                       # 6 df per the protocol
    for c, (pos, _ref) in CLINICAL_BINARIES.items():
        col = d[f"{c}={pos}"]
        assert set(np.unique(col)) <= {0.0, 1.0}
        assert col.sum() == (cov.loc[ids, c] == pos).sum()
    assert REFERENCE_AGE == 69


def test_E06_clinical_block_identical_across_models(data):
    """P6B-A1: the clinical block must be the same in M1, M2p, M2s."""
    ids, cov, X, t, e, p4 = data
    p5 = json.loads((PANELS / "m2_panel_5set.json").read_text())
    grid = np.geomspace(1e-3, 10, 4)
    outs = []
    for panel in (None, p4, p5):
        _, d, _, _ = pipeline_fit(X, cov, panel, ids, ids, t, e, grid,
                                  np.random.default_rng(0))
        outs.append(d)
    clin = list(outs[0].columns)
    for d in outs[1:]:
        assert list(d.columns)[:6] == clin
        assert np.array_equal(d[clin].to_numpy(), outs[0][clin].to_numpy())
    assert outs[1].shape[1] == 10 and outs[2].shape[1] == 11


# ---------- penalty scope (P6B-A1) ----------

def test_E07_clinical_block_is_unpenalised(data):
    ids, cov, X, t, e, p4 = data
    grid = np.geomspace(1e-3, 10, 4)
    _, d, _, dg = pipeline_fit(X, cov, p4, ids, ids, t, e, grid,
                               np.random.default_rng(0))
    pen = np.array(dg["penalizer_vector"])
    assert pen.shape[0] == d.shape[1]
    assert np.all(pen[:6] == 0.0), "clinical block is penalised"
    assert np.all(pen[6:] == dg["lambda"]) and dg["lambda"] > 0


def test_E08_m1_has_no_penalised_term(data):
    ids, cov, X, t, e, _ = data
    _, d, _, dg = pipeline_fit(X, cov, None, ids, ids, t, e,
                               np.geomspace(1e-3, 10, 4), np.random.default_rng(0))
    assert dg["lambda"] == 0.0 and d.shape[1] == 6


def test_E09_penalty_shrinks_only_the_molecular_block(data):
    """A larger lambda must move molecular coefficients and leave clinical freer."""
    ids, cov, X, t, e, p4 = data
    _, d, _, _ = pipeline_fit(X, cov, p4, ids, ids, t, e,
                              np.array([1e-6]), np.random.default_rng(0))
    small = fit_cox(d, t, e, np.r_[np.zeros(6), np.full(4, 1e-6)])
    big = fit_cox(d, t, e, np.r_[np.zeros(6), np.full(4, 50.0)])
    mol = [c for c in d.columns if c.startswith("HALLMARK")]
    assert np.abs(big.params_[mol]).sum() < np.abs(small.params_[mol]).sum()
    assert np.abs(big.params_[mol]).max() < 0.01, "molecular block not shrunk"


# ---------- leakage ----------

def test_E10_pipeline_uses_training_fold_only(data):
    """Poison held-out rows: the fitted model must be bit-identical."""
    ids, cov, X, t, e, p4 = data
    tr, ho = ids[:300], ids[300:]
    grid = np.geomspace(1e-3, 10, 4)
    m1, d1, _, g1 = pipeline_fit(X, cov, p4, tr, tr, t[:300], e[:300], grid,
                                 np.random.default_rng(7))
    Xp = X.copy(); Xp.loc[:, ho] += 1000.0
    covp = cov.copy(); covp.loc[ho, AGE_COL] = "99"
    m2, d2, _, g2 = pipeline_fit(Xp, covp, p4, tr, tr, t[:300], e[:300], grid,
                                 np.random.default_rng(7))
    assert g1["age_knots"] == g2["age_knots"]
    assert np.array_equal(d1.to_numpy(), d2.to_numpy())
    assert np.allclose(m1.params_.to_numpy(), m2.params_.to_numpy(), atol=0, rtol=0)


def test_E11_heldout_predictions_use_training_parameters(data):
    ids, cov, X, t, e, p4 = data
    tr, ho = ids[:300], ids[300:]
    m, d_tr, d_ho, _ = pipeline_fit(X, cov, p4, tr, ho, t[:300], e[:300],
                                    np.geomspace(1e-3, 10, 4), np.random.default_rng(1))
    assert list(d_ho.index) == ho and list(d_tr.index) == tr
    assert list(d_ho.columns) == list(d_tr.columns)
    mol = [c for c in d_tr.columns if c.startswith("HALLMARK")]
    # training scores are standardised to unit variance; held-out are NOT
    assert abs(d_tr[mol].to_numpy().std(ddof=1) - 1.0) < 0.05
    assert abs(d_ho[mol[0]].std(ddof=1) - 1.0) > 1e-6


# ---------- resampling ----------

def test_E12_bootstrap_is_event_stratified():
    d = driver()
    rng = np.random.default_rng(0)
    e = np.r_[np.ones(136), np.zeros(316)]
    for _ in range(20):
        bi = d.boot_indices(e, rng)
        assert len(bi) == len(e)
        assert int(e[bi].sum()) == 136, "event count not preserved by stratification"


def test_E13_risk_set_comparisons_bounded_and_model_free():
    d = driver()
    rng = np.random.default_rng(2)
    for _ in range(30):
        n = int(rng.integers(40, 200))
        t = rng.integers(1, 60, n).astype(float)
        e = rng.binomial(1, 0.3, n).astype(float)
        s = rng.integers(0, 9, n)
        n_un, _ = d.risk_set_comparisons(t, e)
        n_st, per = d.risk_set_comparisons(t, e, s)
        assert 0 <= n_st <= n_un
        assert sum(v["events"] for v in per.values()) == int(e.sum())


def test_E14_zero_event_stratum_contributes_nothing():
    d = driver()
    t = np.array([5., 6., 7., 8.]); e = np.array([1., 0., 0., 0.])
    s = np.array(["a", "a", "b", "b"])
    n, per = d.risk_set_comparisons(t, e, s)
    assert "b" not in per or per["b"]["events"] == 0
    assert n == 1                       # the single event has one comparator in 'a'


# ---------- metrics ----------

def test_E15_c_index_orientation(data):
    _, _, _, t, e, _ = data
    rng = np.random.default_rng(0)
    lp = -t + rng.normal(0, 0.01, len(t))      # higher lp -> shorter time
    assert c_index(t, e, lp) > 0.9
    assert c_index(t, e, -lp) < 0.1


def test_E16_calibration_slope_is_one_on_its_own_fit(data):
    ids, cov, X, t, e, p4 = data
    m, d, _, _ = pipeline_fit(X, cov, p4, ids, ids, t, e,
                              np.geomspace(1e-3, 10, 4), np.random.default_rng(0))
    s = calibration_slope(t, e, linear_predictor(m, d))
    assert 0.8 < s < 1.3


# ---------- results integrity ----------

@pytest.fixture(scope="module")
def res():
    f = EXEC / "phase9_execution_results.json"
    if not f.exists():
        pytest.skip("execution results not yet written")
    return json.loads(f.read_text())


def test_E17_population_and_events(res):
    p = res["population"]
    assert p["n"] == 452 and p["events"] == 136 and p["censored"] == 316
    assert res["endpoint"]["event_col"] == "rfs_event_full"
    assert res["endpoint"]["time_col"] == "rfs_time_full"


def test_E18_m3_not_fitted(res):
    assert res["M3_fitted"] is False
    assert set(res["authorised_analyses"]) == {"M1", "M2_primary", "M2_sensitivity"}


def test_E19_no_prohibited_quantity(res):
    """The declaration list NAMES them; nothing else may mention them."""
    decl = res["prohibited_quantities_not_computed"]
    assert len(decl) == 5
    rest = {k: v for k, v in res.items()
            if k != "prohibited_quantities_not_computed"}
    blob = json.dumps(rest).lower()
    for banned in ("five_year", "5_year", "five-year", "time_dependent_auc",
                   "time-dependent auc", "brier", "cif_at", "cif at",
                   "cloglog", "absolute risk", "absolute_risk"):
        assert banned not in blob, banned
    # and the declaration must cover every prohibited class
    joined = " ".join(decl).lower()
    for must in ("absolute risk", "cif", "time-dependent auc", "brier", "cloglog"):
        assert must in joined, must


def test_E20_nonestimable_recorded_not_redrawn(res):
    ne = res["non_estimable_replicates"]
    assert ne["redrawn"] is False
    assert ne["threshold"] == 0.05
    assert set(ne["rate"]) == {"M1", "M2_primary", "M2_sensitivity"}


def test_E21_five_seeds_and_mc_reported(res):
    assert res["seeds"] == [20261102, 20261103, 20261104, 20261105, 20261106]
    for m, mc in res["monte_carlo"].items():
        assert len(mc["per_seed"]) == 5
        assert "monte_carlo_SE" in mc and "sampling_SE_approx" in mc
        assert "MC3_pass" in mc and "running_trace_MC1" in mc


def test_E22_apparent_and_corrected_both_present(res):
    for m in ("M1", "M2_primary", "M2_sensitivity"):
        assert "C" in res["apparent"][m]
        assert "corrected_C_mean" in res["monte_carlo"][m]
        # optimism must be positive-ish: apparent should exceed corrected
        assert res["apparent"][m]["C"] >= res["monte_carlo"][m]["corrected_C_mean"]


def test_E23_batch_nine_strata_no_pooling(res):
    b = res["batch_sensitivity"]
    assert "nine individual strata, NO pooling" in b["strata"]
    per = b["risk_set_comparisons"]["per_stratum"]
    assert len({k for k in per}) == 9
    assert sum(v["n"] for v in per.values()) == 452
    r = b["risk_set_comparisons"]
    assert 0 <= r["retention"] <= 1
    assert "no threshold applied" in r["interpretation"]
    assert len(b["leave_one_batch_out"]) == 9


def test_E24_input_hashes_recorded_and_current(res):
    import hashlib
    for rel, h in res["input_hashes"].items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == h, rel
