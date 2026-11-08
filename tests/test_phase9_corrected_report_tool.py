"""Final-interval reporting tool (scripts/46_phase9_corrected_report.py).

`analyse` is checked on synthetic records against independent NumPy computations; the
command-line gate is checked for its refusals. No study data or checkpoints are read and
nothing is written inside the repository.
"""
import importlib.util
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "46_phase9_corrected_report.py"
CONFIG = ROOT / "config" / "analysis" / "PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json"


def load_tool():
    spec = importlib.util.spec_from_file_location("p9c_report_tool_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fake_records(run_id, L, rng, shift=0.0):
    recs = {}
    for k in range(L):
        inner = []
        for b in range(100):
            c1, c2 = rng.uniform(0.55, 0.7, 2)
            inner.append({"b": b, "paired_complete": True,
                          "fits": {"M1": {"C_fit": c1 + 0.02, "C_eval": c1},
                                   "M2": {"C_fit": c2 + 0.03, "C_eval": c2}}})
        dapp = float(rng.normal(0.02, 0.01))
        dopt = float(np.mean([(r["fits"]["M2"]["C_fit"] - r["fits"]["M2"]["C_eval"])
                              - (r["fits"]["M1"]["C_fit"] - r["fits"]["M1"]["C_eval"]) for r in inner]))
        recs[k] = {"k": k, "run_id": run_id, "status": "ok", "B_in_attempted": 100,
                   "B_in_paired_complete": 100, "inner": inner, "failures": [],
                   "dC_app": dapp, "dC_app_as_run": dapp, "dOpt": dopt,
                   "dC_corrected": float(rng.normal(0.012 + shift, 0.016)),
                   "C_corrected": {"M1": 0.60, "M2": 0.61}, "n_used": {"M1": 100, "M2": 100}}
    return recs


@pytest.fixture(scope="module")
def result():
    from src.phase9_corrected.config import load_config
    tool = load_tool()
    cfg = load_config(str(CONFIG), ROOT)
    rng = np.random.default_rng(7)
    recs = {r: fake_records(r, 500, rng) for r in cfg.run_ids}
    return tool.analyse(cfg, recs, 500, 0.0043), recs


def test_percentile_intervals_match_numpy(result):
    res, recs = result
    for r, x in res["runs"].items():
        v = np.array([recs[int(r)][k]["dC_corrected"] for k in range(500)])
        lo, hi = np.quantile(v, [0.025, 0.975], method="linear")
        assert x["percentile_primary"]["lower"] == pytest.approx(lo, rel=0, abs=1e-15)
        assert x["percentile_primary"]["upper"] == pytest.approx(hi, rel=0, abs=1e-15)
        assert x["n_usable"] == 500 and x["H1"] == "not triggered"


def test_between_run_and_pooled_match_numpy(result):
    res, recs = result
    lows = [res["runs"][str(r)]["percentile_primary"]["lower"] for r in sorted(recs)]
    assert res["between_run"]["lower"]["sd_ddof1"] == pytest.approx(np.std(lows, ddof=1), abs=1e-15)
    allv = np.concatenate([[recs[r][k]["dC_corrected"] for k in range(500)] for r in sorted(recs)])
    assert res["pooled_descriptive"]["n"] == 2500
    assert res["pooled_descriptive"]["percentile"][0] == pytest.approx(np.quantile(allv, 0.025), abs=1e-15)
    assert "DESCRIPTIVE ONLY" in res["pooled_descriptive"]["note"]


def test_kernel_se_matches_independent_formula(result):
    res, recs = result
    v = np.array([recs[0][k]["dC_corrected"] for k in range(500)])
    s = np.std(v, ddof=1)
    iqr = np.quantile(v, 0.75) - np.quantile(v, 0.25)
    h = 0.9 * min(s, iqr / 1.34) * 500 ** -0.2
    q = np.quantile(v, 0.025)
    f = np.mean(np.exp(-0.5 * ((q - v) / h) ** 2)) / (h * math.sqrt(2 * math.pi))
    se = math.sqrt(0.025 * 0.975 / 500) / f
    got = res["runs"]["0"]["kernel_mcse"]["lower"]["SE"]
    got = got.get("1.0", got.get(1.0))
    assert got == pytest.approx(se, rel=1e-12)


def test_diagnostics_and_flags_present(result):
    res, _ = result
    x = res["runs"]["0"]
    assert x["H2"]["vacuous"] is True and x["H2"]["H2"] == "not triggered"
    assert x["inner_fit_failures"] == 0 and x["H6_unequal_paired_counts"] == 0
    assert math.isclose(x["A7"], 2 * abs(0.5 * (x["percentile_primary"]["lower"]
                                                + x["percentile_primary"]["upper"]) - 0.0043)
                        / x["percentile_primary"]["width"], rel_tol=1e-12)
    comp = x["A7_components"]
    mid = 0.5 * (x["percentile_primary"]["lower"] + x["percentile_primary"]["upper"])
    assert comp["tail_asymmetry_mid_minus_median"] + comp["median_offset_median_minus_theta"] == \
        pytest.approx(mid - 0.0043, abs=1e-15)
    assert res["conditional_statement"].startswith("Conditional on the observed 136")
    assert res["sizing_decision_at_level"]["action"] in (
        "top_up", "stop_ep1_pass", "stop_label_imprecise", "report_persisting_breach", "stop_for_owner_review")


def run_cli(args, flags=()):
    return subprocess.run([sys.executable, *flags, str(SCRIPT), *args], capture_output=True,
                          text=True, cwd=ROOT, timeout=120)


def test_cli_usage_and_flag_refusals(tmp_path):
    pfx = tmp_path / "pyc"
    pfx.mkdir()
    r = run_cli([], ("-I", "-B", "-X", f"pycache_prefix={pfx}"))
    assert r.returncode == 2 and "usage" in r.stderr
    r = run_cli([str(CONFIG), str(tmp_path / "out")], ("-B",))
    assert r.returncode == 4 and "pycache_prefix" in r.stderr
    assert not (tmp_path / "out").exists()
