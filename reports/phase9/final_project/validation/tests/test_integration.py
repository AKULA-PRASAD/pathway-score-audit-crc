"""Integration / end-to-end tests on the REAL Phase 9 outputs (read-only).

Run after analysis/s0_verify_starting_state.py, analysis/s1_interval_analysis.py and
validation/independent_check.py.
"""
import csv, hashlib, json, math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[5]
FP = ROOT / "reports/phase9/final_project"
EXPECTED_CKPT = {0: "885bbd1bd891b570", 1: "1bb27c56612a5f9d", 2: "0f12dce12bfce95e",
                 3: "4115bd4cd014d345", 4: "2202d559b0418392"}
PROV = "PROVISIONAL"


@pytest.fixture(scope="module")
def res():
    return json.loads((FP / "results/s1_interval_results.json").read_text())


def test_checkpoints_unchanged():
    for r, h in EXPECTED_CKPT.items():
        p = ROOT / f"reports/phase9/dc_interval/run{r}_outer.jsonl"
        assert hashlib.sha256(p.read_bytes()).hexdigest().startswith(h)


def test_results_were_computed_on_the_verified_checkpoints(res):
    assert {int(k): v for k, v in res["checkpoint_hashes_16"].items()} == EXPECTED_CKPT


def test_starting_state_verification_passed():
    s0 = json.loads((FP / "results/s0_starting_state.json").read_text())
    assert s0["all_ok"] is True and s0["driver_check_frozen"] == "PASS"


def test_independent_check_passed():
    ic = json.loads((FP / "validation/independent_check_results.json").read_text())
    assert ic["all_ok"] and ic["n_checks"] >= 150


def test_every_M2_dependent_block_is_labelled_provisional(res):
    assert PROV in res["label"] and PROV in res["EP1"]["label"]
    assert all(PROV in x["label"] for x in res["runs"])


def test_interval_structure(res):
    for x in res["runs"]:
        p, b = x["percentile_primary"], x["basic_diagnostic"]
        assert p["lower"] < p["upper"] and b["lower"] < b["upper"]
        assert math.isclose(p["width"], p["upper"] - p["lower"])
        # basic is the reflection of percentile about theta_hat
        assert math.isclose(b["lower"] + p["upper"], 2 * res["theta_hat"]["value"], abs_tol=1e-15)
        assert math.isclose(b["upper"] + p["lower"], 2 * res["theta_hat"]["value"], abs_tol=1e-15)
        assert x["n_usable"] == 500 and x["z0_diagnostic"]["B"] == 500


def test_ep1_logic_consistent_with_grid(res):
    e = res["EP1"]
    assert e["n_evaluations"] == len(e["grid"]) == 40
    assert e["n_pass"] == sum(g["pass"] for g in e["grid"])
    for g in e["grid"]:
        assert g["pass"] == (g["MCSE"] <= 0.05 * g["width"])
    # status may only read PASS if every evaluation passes; never relaxed
    if not e["all_pass"]:
        assert "PASS under every" not in e["reported_status"]
        assert e["plan_2_4_ladder_if_strict"]["executed"] is False


def test_h3_not_claimed_without_ladder(res):
    assert res["hard_failures"]["H3"].startswith("NOT TRIGGERED")


def test_hard_failure_flags(res):
    assert res["hard_failures"]["H1"] is False and res["hard_failures"]["H2"] is False
    assert res["hard_failures"]["H6"] is False


def test_tables_match_results(res):
    with open(FP / "tables/T3_between_run.csv") as fh:
        rows = {r["quantity"]: r for r in csv.DictReader(fh)}
    for k in ("percentile_lower", "percentile_upper"):
        assert float(rows[k]["sd_ddof1"]) == res["between_run"][k]["sd_ddof1"]
