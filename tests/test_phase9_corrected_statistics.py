"""EP-1, W5, ladder and interval statistics on FIXED synthetic inputs, cross-checked
against the independent final-project implementation (stats_core.py, read-only)."""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from src.phase9_corrected import statistics as st
from src.phase9_corrected.config import UnresolvedDecisionError

ROOT = Path(__file__).resolve().parents[1]
LADDER = {"rungs": [750, 1000, 1500, 2000, 3000], "beyond_ladder_rule": "stop_and_label_imprecise"}


@pytest.fixture(scope="module")
def sc():
    p = ROOT / "reports/phase9/final_project/analysis/stats_core.py"
    spec = importlib.util.spec_from_file_location("fp_stats_core_ro", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def sample(seed, n=500, loc=0.0, sd=0.014):
    return np.random.default_rng(seed).normal(loc, sd, n)


def test_cross_check_against_final_project_stats_core(sc):
    x = sample(1)
    assert st.percentile_ci(x) == sc.percentile_ci(x)
    assert st.silverman_bandwidth(x) == sc.silverman_bandwidth(x)
    lo, hi = st.percentile_ci(x)
    h = st.silverman_bandwidth(x)
    for m in st.BW_MULT:
        assert math.isclose(st.kernel_quantile_se(x, 0.025, lo, h * m),
                            sc.kernel_quantile_se(x, 0.025, lo, h * m)["SE"], rel_tol=1e-14)
    assert st.basic_ci(x, 0.004) == sc.basic_ci(x, 0.004)
    assert math.isclose(st.z0(x, 0.004), sc.z0(x, 0.004)["z0"], rel_tol=1e-12)


def test_silverman_formula_by_hand():
    x = np.array([0.0, 1.0, 2.0, 3.0, 10.0])
    s = np.std(x, ddof=1)
    iqr = np.quantile(x, 0.75) - np.quantile(x, 0.25)
    assert st.silverman_bandwidth(x) == 0.9 * min(s, iqr / 1.34) * 5 ** (-0.2)


def test_ep1_requires_option():
    with pytest.raises(UnresolvedDecisionError):
        st.ep1_evaluate({0: sample(0)}, None)


def test_ep1_option_B_is_never_more_lenient_than_A():
    runs = {r: sample(10 + r) for r in range(5)}
    a, b = st.ep1_evaluate(runs, "A"), st.ep1_evaluate(runs, "B")
    for ra, rb in zip(a["rows"], b["rows"]):
        assert rb["ratio"] >= ra["ratio"]
        assert ra["governing_SE"] == ra["SE_by_bw"]["1.0"]
        assert rb["governing_SE"] == max(rb["SE_by_bw"].values())
    assert b["rho_max"] >= a["rho_max"]


def test_ep1_separates_options_on_constructed_data():
    """Heavy-tailed run: 0.5h estimate fails while h passes -> A passes, B fails."""
    found = None
    for seed in range(200):
        x = np.random.default_rng(seed).standard_t(3, 500) * 0.01
        a, b = st.ep1_evaluate({0: x}, "A"), st.ep1_evaluate({0: x}, "B")
        if a["pass"] and not b["pass"]:
            found = seed
            break
    assert found is not None


def test_ladder_rule_matches_plan_table_at_500():
    # plan lines 121-125: 6 % -> 750, 8 % -> 1500, 10 % -> 2000
    for ratio, rung in ((0.06, 750), (0.08, 1500), (0.10, 2000)):
        d = st.ladder_assessment(500, {"pass": False, "rho_max": ratio, "fraction": 0.05}, LADDER)
        assert d["B_next"] == rung
    assert st.ladder_assessment(500, {"pass": True, "rho_max": 0.01, "fraction": 0.05},
                                LADDER)["action"] == "stop_ep1_pass"


def test_ladder_beyond_top_rung_requires_owner_rule():
    ep = {"pass": False, "rho_max": 0.2, "fraction": 0.05}          # B_req = 8000
    with pytest.raises(UnresolvedDecisionError):
        st.ladder_assessment(500, ep, None)
    assert st.ladder_assessment(500, ep, LADDER)["action"] == "stop_label_imprecise"
    alt = dict(LADDER, beyond_ladder_rule="top_up_to_max_rung_then_reassess")
    assert st.ladder_assessment(500, ep, alt)["B_next"] == 3000
    assert st.ladder_assessment(3000, ep, alt)["action"] == "stop_label_imprecise"


def test_w5_grid_points():
    assert st.w5_grid_points(500) == [300, 400, 500]
    assert st.w5_grid_points(750) == [500, 750]
    assert st.w5_grid_points(3000) == [2000, 3000]


def test_w5_evaluate_matches_final_project_definition_at_500(sc):
    x = sample(4)
    w = st.w5_evaluate({k: v for k, v in enumerate(x)}, 500, st.percentile_ci(x)[1] - st.percentile_ci(x)[0])
    tr = sc.running_trace(x)
    fin = [tr[B] for B in (300, 400, 500)]
    assert w["spread"]["lower"] == max(t[0] for t in fin) - min(t[0] for t in fin)
    assert w["spread"]["upper"] == max(t[1] for t in fin) - min(t[1] for t in fin)


def test_w5_refuses_excluded_replicates_and_incomplete_k():
    x = {k: 0.0 for k in range(500)}
    x[7] = None
    with pytest.raises(UnresolvedDecisionError):
        st.w5_evaluate(x, 500, 0.05)
    with pytest.raises(ValueError):
        st.w5_evaluate({k: 0.0 for k in range(499)}, 500, 0.05)


def test_w5_actions():
    breach = {"breach": True}
    with pytest.raises(UnresolvedDecisionError):
        st.w5_decision(breach, None, None, 500, LADDER)
    assert st.w5_decision(breach, "W5-B", None, 500, LADDER)["action"] == "investigate_and_report"
    with pytest.raises(UnresolvedDecisionError):
        st.w5_decision(breach, "W5-A", None, 500, LADDER)
    d = st.w5_decision(breach, "W5-A", "next_rung_above_current", 500, LADDER)
    assert d == {"action": "top_up", "B_next": 750}
    assert st.combined_next_size({"B_next": 1500}, d) == 1500
    assert st.w5_decision({"breach": False}, "W5-A", "next_rung_above_current", 500, LADDER)["action"] == "none"


def test_run_summary_A7_identity_and_between_run():
    x = sample(5, loc=-0.005)
    s = st.run_summary(x, 0.0043)
    lo, hi = s["percentile"]
    blo, bhi = s["basic"]
    assert math.isclose(abs(lo - blo) / s["width"], s["A7"], rel_tol=1e-12)
    assert math.isclose(abs(hi - bhi) / s["width"], s["A7"], rel_tol=1e-12)
    br = st.between_run({0: (1.0, 2.0), 1: (1.2, 2.1), 2: (0.9, 1.9)})
    assert br["lower"]["df"] == 2 and math.isclose(br["lower"]["sd_ddof1"], np.std([1.0, 1.2, 0.9], ddof=1))


def test_z0_conventions():
    assert st.z0(np.array([1.0, 2.0, 3.0]), 0.0) == -math.inf
    assert st.z0(np.array([1.0, 2.0, 3.0, 4.0]), 2.5) == 0.0
    assert st.z0(np.array([1.0, 2.0, 2.0, 4.0]), 2.0) < 0          # ties not counted ("<")
