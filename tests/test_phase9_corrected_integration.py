"""Controlled end-to-end integration on SYNTHETIC data (not the scientific analysis).

Five synthetic runs x 12 outer x 6 inner replicates through the real orchestration loop
(execute_replicates + CheckpointStore), with the fast stand-in pipeline. Checks exact
reproducibility across repeated runs, worker counts and an interrupted/resumed run, then
the verification interfaces and the statistical summaries on the outputs."""
import json
from pathlib import Path

import numpy as np
import pytest

from src.phase9_corrected import statistics as st
from src.phase9_corrected.checkpoint import CheckpointStore
from src.phase9_corrected.config import PLAN_BASE
from src.phase9_corrected.engine import compute_outer_replicate
from src.phase9_corrected.data import LAM_GRID
from src.phase9_corrected.orchestration import StopCondition, execute_replicates
from src.phase9_corrected.rng import PlanStreams, enumerate_stream_states
from src.phase9_corrected.verification import (v1_self_consistency, v2_record_assignment,
                                               v5_spot_check)
from tests.p9c_support import FakePipeline, synthetic_context

T_TEST, B_OUT, B_IN, RUNS = 7, 12, 6, range(5)


def execute(base_dir, run_id, workers, stop_after=None, ctx=None):
    ctx = ctx or synthetic_context()
    s = PlanStreams(PLAN_BASE, run_id, "independent", T_TEST)
    calls = {"n": 0}

    def on_record(rec, store):
        calls["n"] += 1
        if stop_after and calls["n"] >= stop_after:
            raise StopCondition("simulated interruption")
    with CheckpointStore(base_dir / "out", run_id, "c" * 64, "m" * 64, root=base_dir) as store:
        try:
            execute_replicates(store, ctx, s, FakePipeline(), B_IN, B_OUT, workers, False, on_record)
        except StopCondition:
            pass
        return store.path, dict(store.done)


@pytest.fixture(scope="module")
def five_runs(tmp_path_factory):
    a = tmp_path_factory.mktemp("a")
    return a, {r: execute(a, r, 1) for r in RUNS}


def test_repeat_and_worker_count_give_byte_identical_checkpoints(five_runs, tmp_path):
    a, res = five_runs
    for r in (0, 3):
        p1 = execute(tmp_path / f"w1_{r}", r, 1)[0]
        p3 = execute(tmp_path / f"w3_{r}", r, 3)[0]
        assert p1.read_bytes() == res[r][0].read_bytes() == p3.read_bytes()


def test_interrupted_run_resumes_to_identical_bytes(five_runs, tmp_path):
    a, res = five_runs
    p, done = execute(tmp_path, 2, 3, stop_after=5)
    assert 5 <= len(done) < B_OUT
    p, done = execute(tmp_path, 2, 3)
    assert sorted(done) == list(range(B_OUT))
    assert p.read_bytes() == res[2][0].read_bytes()


def test_records_complete_and_runs_distinct(five_runs):
    _, res = five_runs
    for r in RUNS:
        recs = res[r][1]
        assert sorted(recs) == list(range(B_OUT))
        for k, rec in recs.items():
            assert rec["run_id"] == r and rec["k"] == k and rec["status"] == "ok"
            assert len(rec["inner"]) == B_IN and rec["B_in_paired_complete"] == B_IN
            assert len(rec["outer"]["multiplicity"]) == 60
            assert {row["b"] for row in rec["inner"]} == set(range(B_IN))
    sha = {rec["outer"]["multiplicity_sha256"] for r in RUNS for rec in res[r][1].values()}
    assert len(sha) == len(RUNS) * B_OUT          # every outer resample distinct


def test_gate_interfaces_on_integration_output(five_runs):
    _, res = five_runs
    ctx = synthetic_context()
    for r in RUNS:
        s = PlanStreams(PLAN_BASE, r, "independent", T_TEST)
        recs = [res[r][1][k] for k in range(B_OUT)]
        assert v1_self_consistency(recs, ctx.e, s)["status"] == "PASS"
        assert v2_record_assignment(recs, s)["status"] == "PASS"
    s0 = PlanStreams(PLAN_BASE, 0, "independent", T_TEST)
    v5 = v5_spot_check(res[0][1], lambda k: compute_outer_replicate(
        ctx, s0, FakePipeline(), 0, k, B_IN, LAM_GRID), selected_k=range(0, B_OUT, 3))
    assert v5["status"] == "PASS"
    dom = enumerate_stream_states(s0, list(RUNS), B_OUT, B_IN)
    assert dom["pass"]


def test_statistical_summaries_on_integration_output(five_runs):
    _, res = five_runs
    vals = {r: np.array([res[r][1][k]["dC_corrected"] for k in range(B_OUT)]) for r in RUNS}
    for r, v in vals.items():
        # dC_corrected = dC_app - dOpt, and C_corrected difference equals dC_corrected
        for rec in res[r][1].values():
            assert rec["dC_corrected"] == rec["dC_app"] - rec["dOpt"]
            cc = rec["C_corrected"]["M2"] - rec["C_corrected"]["M1"]
            assert abs(cc - rec["dC_corrected"]) < 1e-15
        s = st.run_summary(v, 0.0043)
        assert s["percentile"][0] <= s["percentile"][1]
    for opt in ("A", "B"):
        ep = st.ep1_evaluate(vals, opt)
        assert len(ep["rows"]) == 10 and set(ep["B_out"].values()) == {B_OUT}
    br = st.between_run({r: st.percentile_ci(v) for r, v in vals.items()})
    assert br["lower"]["df"] == 4
