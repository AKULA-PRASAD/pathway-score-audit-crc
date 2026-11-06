"""Corrected replicate engine on SYNTHETIC data (controlled integration tests).

Includes an equivalence test against the FROZEN driver's own outer_replicate: the frozen
module is loaded in memory (file untouched), its data and fitting functions are swapped
for the synthetic stand-ins, and the new engine with the historical streams must give a
bit-identical result. This tests orchestration, pairing and arithmetic, not science."""
import json
from pathlib import Path

import numpy as np
import pytest

from src.phase9_corrected.config import PLAN_BASE, canonical_json
from src.phase9_corrected.data import LAM_GRID
from src.phase9_corrected.engine import compute_outer_replicate
from src.phase9_corrected.fitting import RealPipeline
from src.phase9_corrected.rng import HistoricalAsRunStreams, PlanStreams
from src.phase9_corrected.verification import (load_frozen_driver, v1_self_consistency,
                                               v2_record_assignment, v3_m1_invariance,
                                               v4_as_run_apparent, v5_spot_check)
from tests.p9c_support import (FakeFailure, FakePipeline, synthetic_context,
                               synthetic_real_pipeline_context)

ROOT = Path(__file__).resolve().parents[1]
T_TEST = 7
B_IN = 8


@pytest.fixture(scope="module")
def ctx():
    return synthetic_context()


def frozen_with_fakes(ctx, fake):
    drv = load_frozen_driver(ROOT)
    drv._G.update(ids=ctx.ids, cov=ctx.cov, t=ctx.t, e=ctx.e, X=ctx.X, panel=ctx.panel)
    drv.pipeline_fit = fake.pipeline_fit
    drv.linear_predictor = fake.linear_predictor
    drv.ConvergenceFailure = FakeFailure
    drv.NonEstimableReplicate = FakeFailure
    return drv


def run(ctx, streams, k, fake=None, B_in=B_IN, variant=False):
    return json.loads(canonical_json(compute_outer_replicate(
        ctx, streams, fake or FakePipeline(), streams.run_id, k, B_in, LAM_GRID,
        as_run_variant=variant)))


@pytest.mark.parametrize("k", [0, 3, 11])
def test_engine_with_historical_streams_reproduces_frozen_driver(ctx, k):
    fake = FakePipeline()
    drv = frozen_with_fakes(ctx, fake)
    old = json.loads(json.dumps(drv.outer_replicate((1, k, B_IN))))
    new = run(ctx, HistoricalAsRunStreams(PLAN_BASE, 1, B_IN), k, fake)
    assert new["outer"]["legacy_multiplicity_sha16"] == old["multiplicity_sha"]
    assert float(new["dC_corrected"]).hex() == float(old["dC_corrected"]).hex()
    assert new["dC_app"] == old["dC_app"] and new["dOpt"] == old["dOpt"]
    assert new["opt_mean"] == {"M1": old["opt_M1_mean"], "M2": old["opt_M2_mean"]}
    assert new["C_app"] == old["C_app"]
    assert {m: new["apparent"][m]["lam"] for m in ("M1", "M2")} == old["lambda_apparent"]
    for m in ("M1", "M2"):
        got = [np.round([r["fits"][m]["C_fit"], r["fits"][m]["C_eval"]], 8).tolist()
               for r in new["inner"]]
        assert got == old["inner_C"][m]
    hist = {k: old}
    assert v3_m1_invariance([new], hist)["status"] == "PASS"
    assert v4_as_run_apparent([new], hist, "as_run", None)["status"] == "PASS"


def test_corrected_streams_change_only_M2_cv_dependent_quantities(ctx):
    fake = FakePipeline()
    drv = frozen_with_fakes(ctx, fake)
    k = 5
    old = json.loads(json.dumps(drv.outer_replicate((2, k, B_IN))))
    new = run(ctx, PlanStreams(PLAN_BASE, 2, "independent", T_TEST), k, fake)
    hist = {k: old}
    assert new["outer"]["legacy_multiplicity_sha16"] == old["multiplicity_sha"]   # V1 outer
    assert v3_m1_invariance([new], hist)["status"] == "PASS"                     # V3
    # inner CV stream changed => at least one M2 inner lambda differs from the as-run path
    asrun = run(ctx, HistoricalAsRunStreams(PLAN_BASE, 2, B_IN), k, fake)
    assert any(a["fits"]["M2"]["lam"] != b["fits"]["M2"]["lam"]
               for a, b in zip(new["inner"], asrun["inner"]))
    assert all(a["resample"] == b["resample"] for a, b in zip(new["inner"], asrun["inner"]))
    # V4 must NOT be applied to the independent primary
    assert v4_as_run_apparent([new], hist, "independent", False)["status"] == "NOT_APPLICABLE"
    assert v4_as_run_apparent([new], hist, "independent", None)["status"] == "BLOCKED"


def test_as_run_variant_matches_history_while_primary_differs(ctx):
    fake = FakePipeline()
    drv = frozen_with_fakes(ctx, fake)
    k = 2
    old = json.loads(json.dumps(drv.outer_replicate((0, k, B_IN))))
    new = run(ctx, PlanStreams(PLAN_BASE, 0, "independent", T_TEST), k, fake, variant=True)
    assert v4_as_run_apparent([new], {k: old}, "independent", True)["status"] == "PASS"
    assert new["streams"]["apparent_cv"]["state_sha256"] != new["streams"]["outer"]["state_sha256"]


def test_cv_consumption_only_by_M2_and_stream_assignment(ctx):
    s = PlanStreams(PLAN_BASE, 0, "independent", T_TEST)
    recs = [run(ctx, s, k) for k in range(3)]
    for r in recs:
        assert r["apparent"]["M1"]["cv_rng_consumed"] is False
        assert r["apparent"]["M2"]["cv_rng_consumed"] is True
    assert v2_record_assignment(recs, s)["status"] == "PASS"
    # a record produced under the wrong policy is detected
    wrong = run(ctx, HistoricalAsRunStreams(PLAN_BASE, 0, B_IN), 1)
    assert v2_record_assignment([wrong], s)["status"] == "FAIL"


def test_v1_self_consistency_and_tamper_detection(ctx):
    s = PlanStreams(PLAN_BASE, 0, "independent", T_TEST)
    recs = [run(ctx, s, k) for k in range(3)]
    assert v1_self_consistency(recs, ctx.e, s)["status"] == "PASS"
    t = json.loads(json.dumps(recs))
    t[1]["inner"][4]["resample"]["index_sha256"] = "0" * 64
    t[2]["outer"]["multiplicity"][0] += 1
    t[0]["inner"][0]["fits"]["M2"]["resample_index_sha256"] = "f" * 64
    res = v1_self_consistency(t, ctx.e, s)
    assert res["status"] == "FAIL" and res["n_mismatches"] == 3


def test_pairing_preserving_exclusion_and_failure_log(ctx):
    calls = {"n": 0}

    def fail_when(train_ids, panel):            # M2 call 1 is apparent; fail M2 inner call 6
        if panel is not None:
            calls["n"] += 1
            if calls["n"] == 7:
                return "synthetic separation"
        return None
    s = PlanStreams(PLAN_BASE, 0, "independent", T_TEST)
    rec = run(ctx, s, 1, FakePipeline(fail_when), B_in=30)
    n_fail = sum(not r["paired_complete"] for r in rec["inner"])
    assert n_fail == 1 and not rec["inner"][5]["paired_complete"]
    assert rec["n_used"]["M1"] == rec["n_used"]["M2"] == 30 - n_fail     # H6
    inner_fail = [f for f in rec["failures"] if f["stage"] == "inner"]
    assert len(inner_fail) == n_fail
    for f in inner_fail:
        assert f["model"] == "M2" and f["reason"] == "synthetic separation"
        assert set(f) >= {"run_id", "k", "b", "n_distinct", "events", "binary_cells",
                          "failure_type", "lambda_at_failure"}
        row = rec["inner"][f["b"]]
        assert row["fits"]["M1"]["ok"] and not row["fits"]["M2"]["ok"]   # M1 kept for audit
    # opt means are over the paired set only
    pairs = [r for r in rec["inner"] if r["paired_complete"]]
    o1 = np.mean([r["fits"]["M1"]["C_fit"] - r["fits"]["M1"]["C_eval"] for r in pairs])
    assert np.isclose(rec["opt_mean"]["M1"], o1, rtol=0, atol=1e-15)
    assert rec["status"] == "ok" and rec["B_in_paired_complete"] == 29


def test_inner_shortfall_below_95_percent(ctx):
    calls = {"n": 0}

    def fail_when(train_ids, panel):
        if panel is not None:
            calls["n"] += 1
            if calls["n"] in (3, 4):                # two of 20 inner M2 fits -> 90 % paired
                return "synthetic separation"
        return None
    rec = run(ctx, PlanStreams(PLAN_BASE, 0, "independent", T_TEST), 1, FakePipeline(fail_when), B_in=20)
    assert rec["B_in_paired_complete"] == 18 and rec["status"] == "unreliable_inner_shortfall"
    assert "dC_corrected" not in rec and len(rec["failures"]) == 2


def test_apparent_failure_marks_unreliable(ctx):
    s = PlanStreams(PLAN_BASE, 0, "independent", T_TEST)
    always = FakePipeline(lambda ids, panel: "apparent boom" if panel is not None else None)
    rec = run(ctx, s, 0, always)
    assert rec["status"] == "unreliable_apparent_failed" and "inner" not in rec
    assert rec["apparent"]["M1"]["ok"] and not rec["apparent"]["M2"]["ok"]
    assert [f["stage"] for f in rec["failures"]] == ["apparent"]


def test_determinism_and_v5(ctx):
    s = PlanStreams(PLAN_BASE, 3, "independent", T_TEST)
    a = {k: run(ctx, s, k) for k in range(4)}
    b = {k: run(ctx, s, k) for k in range(4)}
    assert canonical_json(a) == canonical_json(b)
    assert v5_spot_check(a, lambda k: run(ctx, s, k), selected_k=range(4))["status"] == "PASS"
    other = PlanStreams(PLAN_BASE, 3, "independent", T_TEST + 1)
    assert v5_spot_check(a, lambda k: run(ctx, other, k), selected_k=range(4))["status"] == "FAIL"
    # different runs are different
    c = run(ctx, PlanStreams(PLAN_BASE, 4, "independent", T_TEST), 0)
    assert c["outer"]["multiplicity_sha256"] != a[0]["outer"]["multiplicity_sha256"]


def test_real_frozen_pipeline_m1_m2_and_lambda_selection_on_synthetic_data():
    """The REAL frozen pipeline (src/phase9_execution.py) on a small synthetic cohort:
    M1 selects lambda = 0 without drawing; M2 runs the ridge CV, draws from its CV
    generator and selects lambda from the frozen grid."""
    ctx = synthetic_real_pipeline_context()
    s = PlanStreams(PLAN_BASE, 0, "independent", T_TEST)
    rec = json.loads(canonical_json(compute_outer_replicate(ctx, s, RealPipeline(), 0, 0, 2,
                                                            LAM_GRID)))
    assert rec["status"] in ("ok", "unreliable_inner_shortfall", "unreliable_apparent_failed")
    a = rec["apparent"]
    assert a["M1"]["ok"] and a["M1"]["lam"] == 0.0 and a["M1"]["cv_rng_consumed"] is False
    assert a["M2"]["ok"] and a["M2"]["lam"] in [float(x) for x in LAM_GRID]
    assert a["M2"]["cv_rng_consumed"] is True and a["M2"]["inner_cv_c"] is not None
    for row in rec["inner"]:
        assert row["fits"]["M1"]["cv_rng_consumed"] is False
        if row["fits"]["M2"]["ok"]:
            assert row["fits"]["M2"]["cv_rng_consumed"] is True
    rec2 = json.loads(canonical_json(compute_outer_replicate(ctx, s, RealPipeline(), 0, 0, 2,
                                                             LAM_GRID)))
    assert rec == rec2                                               # bit-identical rerun
