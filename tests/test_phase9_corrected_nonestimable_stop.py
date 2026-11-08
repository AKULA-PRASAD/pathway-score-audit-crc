"""Owner Decision 11 (2026-11-06): literal 5 % stop, dispatch refusal, three execution
paths, checkpoint continuity, completion eligibility (Stage B), Decision 7B sizing-file
verification by recomputation, Decision 9 code identity against a Decision 13 freeze,
global historical-preservation block, H6, durability.

Synthetic only: the fit (compute_outer_replicate), data loading, gate, H4, manifest and
(unless a test says otherwise) code identity are stubbed; freeze, sizing and window records
are synthetic and live in temporary directories. No scientific data, no real bootstrap."""
import hashlib
import importlib.util
import json
import py_compile
import shutil
import sys
import uuid
from pathlib import Path

import numpy as np
import pytest

from src.phase9_corrected import orchestration as orch
from src.phase9_corrected.checkpoint import CheckpointCorruption, CheckpointStore
from src.phase9_corrected.config import ConfigurationError, CorrectedConfig
from src.phase9_corrected.engine import ImplementationDefectH6
from src.phase9_corrected.integrity import HardFailureH4, historical_preservation_check
from tests.p9c_support import approved_test_config, make_synthetic_gate, synthetic_gate

ROOT = Path(__file__).resolve().parents[1]
DARWIN = sys.platform == "darwin"
REAL_VERIFY = orch.verify_code_identity
SYNTH_FREEZE = {"freeze_record_sha256": "f" * 64,
                "freeze": {"code": {"src/phase9_corrected/__init__.py": "0" * 64}}}
W5A = {"action": "W5-A", "corrective_rule": "next_rung_above_current",
       "definition": "adjudication_D5_as_proposed", "exclusion_rule": "usable_within_k_positions_v1"}
LADDER = {"rungs": [750, 1000, 1500, 2000, 3000],
          "beyond_ladder_rule": "top_up_to_max_rung_then_reassess"}

# ---------------------------------------------------------------- fakes
B = {"bad": {}, "raise_at": None, "exc": RuntimeError, "calls": 0, "drift": 0.0, "side": None}


def fake_record(ctx, streams, pipeline, run_id, k, B_in, lam_grid, as_run_variant=False):
    B["calls"] += 1
    if B["side"]:
        B["side"](run_id, k)
    if B["raise_at"] is not None and k == B["raise_at"]:
        raise B["exc"](f"simulated failure at k={k}")
    v = float(np.random.default_rng([run_id, k]).normal(0.0, 0.01) + B["drift"] * k / 500)
    return {"run_id": int(run_id), "k": int(k), "status": B["bad"].get(k, "ok"),
            "B_in_attempted": 0, "inner": [], "dC_corrected": v, "dC_app": v}


class StepClock:
    def __init__(self, step=0.001):
        self.step, self.n = step, 0

    def __call__(self):
        v = self.n * self.step
        self.n += 1
        return v


class Env:
    def __init__(self, out_rel):
        self.out_rel, self.out = out_rel, ROOT / out_rel
        self.h4_fail_stage, self.pres_fail_next, self.pres_result = None, False, None
        self.cfg_over = {}

    def cfg(self, n_workers=1, **over):
        kw = {"n_workers": n_workers, "runtime_budget_hours_per_run": 4.0, **self.cfg_over, **over}
        return CorrectedConfig(approved_test_config(ROOT, self.out_rel, **kw), ROOT)

    def run(self, run_id=0, target=500, clock=None, n_workers=1, **kw):
        return orch.run_corrected(self.cfg(n_workers), run_id, target,
                                  clock=("TEST_CLOCK", clock or StepClock()), **kw)

    def ckpt(self, r=0):
        return self.out / "runs" / f"run{r}.jsonl"

    def msha(self):
        return hashlib.sha256((self.out / "manifest/input_manifest.json").read_bytes()).hexdigest()

    def elig(self, r=0, target=500, **kw):
        return orch.evaluate_eligibility(self.cfg(), r, target, self.msha(), **kw)

    def halts(self, r, m, kind):
        return self.out / "halts" / f"run{r}.inv{m:04d}.{kind}.json"

    def rewrite(self, r, m, kind, fn):
        p = self.halts(r, m, kind)
        o = json.loads(p.read_text())
        fn(o)
        p.write_text(json.dumps(o))


@pytest.fixture
def env(monkeypatch):
    e = Env(f"scratch_p9c_ne_{uuid.uuid4().hex[:10]}")
    B.update(bad={}, raise_at=None, exc=RuntimeError, calls=0, drift=0.0, side=None)

    def h4(root, stage, raise_on_fail=True):
        if stage == e.h4_fail_stage:
            raise HardFailureH4(f"HARD FAILURE H4 ({stage}) -- simulated")
        return {"status": "PASS", "all_ok": True}

    def preservation(root):
        if e.pres_fail_next:
            e.pres_fail_next = False
            return e.pres_result or {"all_ok": False, "status": "FAIL", "files": {
                "reports/phase9/dc_interval/run0_outer.jsonl":
                    {"expected_prefix": "885bbd1bd891b570", "sha256": "ab" * 32, "ok": False}}}
        return {"all_ok": True, "status": "PASS"}

    monkeypatch.setattr(orch, "h4_check", h4)
    monkeypatch.setattr(orch, "historical_preservation_check", preservation)
    monkeypatch.setattr(orch, "build_manifest", lambda config, root: {
        "schema": "test-manifest", "config": config.sha256(),
        "code": {"src/phase9_corrected/__init__.py": "0" * 64}})
    monkeypatch.setattr(orch, "run_pre_execution_gate", synthetic_gate)
    monkeypatch.setattr(orch, "load_context", lambda root: None)
    monkeypatch.setattr(orch, "RealPipeline", lambda: None)
    monkeypatch.setattr(orch, "compute_outer_replicate", fake_record)
    monkeypatch.setattr(orch, "verify_code_identity",
                        lambda config, manifest=None, extra_files=(): SYNTH_FREEZE)
    yield e
    shutil.rmtree(e.out, ignore_errors=True)
    assert not e.out.exists()


def halt_at_55(env, **kw):
    with pytest.raises(orch.RuntimeBudgetHalt):
        env.run(clock=StepClock(600), **kw)


def ack_window(env, r, m):
    """A synthetic Stage A PASS result bound to the window Stage B reports."""
    w = [x for x in env.elig(r, max_invocation=None)["windows"] if x["invocation"] == m][0]
    res = {"schema": "p9c-window-1", "run_id": r, "invocation": m,
           "region_sha256": w["region_sha256"], "k_list": w["k_list"],
           "code_identity_verified": True, "status": "PASS",
           "per_k": {str(k): {"record_identical": True, "dC_hex_equal": True} for k in w["k_list"]}}
    p = orch.window_result_path(env.out, r, m)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(res))
    return res


def stray_tmp(env, r=0):
    (env.out / "halts" / f"run{r}.inv0099.end.json.tmp").write_text("{")


# ================================================================ 1. the literal rule
@pytest.mark.parametrize("n,b,stop", [(1, 1, True), (10, 1, True), (20, 1, False), (20, 2, True),
                                      (500, 25, False), (500, 26, True)])
def test_strict_integer_predicate(n, b, stop):
    recs = [{"k": k, "status": "ok" if k >= b else "unreliable_apparent_failed"} for k in range(n)]
    s = orch.nonestimable_status(recs)
    assert s["result"] is stop and s["lhs"] == 20 * b and s["rhs"] == n
    assert s["rule"] == "20*n_nonestimable > n_completed"


def test_prefix_scan_finds_earliest_violation_even_if_later_diluted():
    recs = [{"k": 0, "status": "unreliable_inner_shortfall"}] + \
           [{"k": k, "status": "ok"} for k in range(1, 100)]
    assert orch.nonestimable_status(recs)["result"] is False         # diluted final rate
    assert orch.nonestimable_prefix_violation(recs) == {"position": 1, "k": 0, "n_nonestimable": 1}


def test_unknown_status_fails_closed():
    assert orch.nonestimable_status([{"k": 0, "status": "a_future_status"}])["result"] is True


def test_engine_status_set_is_pinned():
    src = (ROOT / "src/phase9_corrected/engine.py").read_text()
    import re
    assert sorted(set(re.findall(r'"status"\] = "(\w+)"', src)) | {"ok"}) == \
        ["ok", "unreliable_apparent_failed", "unreliable_inner_shortfall"]


@pytest.mark.parametrize("bad_k,stop_k", [({0}, 0), ({1}, 1), ({9}, 9)])
def test_one_early_failure_stops_at_that_record(env, bad_k, stop_k):
    B["bad"] = {k: "unreliable_apparent_failed" for k in bad_k}
    with pytest.raises(orch.NonEstimableStop):
        env.run()
    end = json.loads(env.halts(0, 1, "end").read_text())
    assert end["outcome"] == "stop_nonestimable" and end["nonestimable"]["triggering_k"] == stop_k
    assert end["stop_reasons"] == ["nonestimable_rate"]


def test_failure_at_record_20_continues(env):
    B["bad"] = {19: "unreliable_apparent_failed"}
    assert env.run()["n_records"] == 500


def test_cumulative_crossing_across_resume(env):
    B["bad"] = {29: "unreliable_apparent_failed"}             # 1 of 55: fine
    halt_at_55(env)
    B["bad"][58] = "unreliable_apparent_failed"                # 2 by n = 59 (40 > 59? no)
    B["bad"][59] = "unreliable_apparent_failed"                # 3 at n = 60: 60 > 60? no
    B["bad"][60] = "unreliable_apparent_failed"                # 4 at n = 61: 80 > 61 -> stop
    with pytest.raises(orch.NonEstimableStop):
        env.run(resume_authorization="R1")
    end = json.loads(env.halts(0, 2, "end").read_text())
    assert end["nonestimable"]["triggering_k"] == 60 and end["nonestimable"]["n_completed"] == 61


def test_repeated_invocations_after_stop_are_refused_and_bytes_unchanged(env):
    B["bad"] = {0: "unreliable_apparent_failed"}
    with pytest.raises(orch.NonEstimableStop):
        env.run()
    before, calls = env.ckpt().read_bytes(), B["calls"]
    for _ in range(15):
        with pytest.raises(orch.DispatchRefused, match="non-estimable"):
            env.run()
    assert env.ckpt().read_bytes() == before and B["calls"] == calls
    assert not env.halts(0, 2, "start").exists()


def test_topup_after_stop_refused(env):
    B["bad"] = {0: "unreliable_apparent_failed"}
    with pytest.raises(orch.NonEstimableStop):
        env.run()
    with pytest.raises(orch.DispatchRefused, match="non-estimable"):
        env.run(target=750, topup_authorization="T")


@pytest.mark.parametrize("workers", [1, 3])
def test_stop_point_worker_invariant_and_nothing_appended_after(env, workers):
    B["bad"] = {30: "unreliable_apparent_failed", 33: "unreliable_apparent_failed"}
    with pytest.raises(orch.NonEstimableStop):
        env.run(n_workers=workers)
    recs = CheckpointStore.load_records(env.ckpt())
    assert sorted(recs) == list(range(34))                     # in-flight results discarded
    end = json.loads(env.halts(0, 1, "end").read_text())
    assert end["nonestimable"]["triggering_k"] == 33


def test_stop_bytes_identical_across_worker_counts(env):
    B["bad"] = {30: "unreliable_apparent_failed", 33: "unreliable_apparent_failed"}
    out = {}
    for w in (1, 3):
        with pytest.raises(orch.NonEstimableStop):
            env.run(n_workers=w)
        out[w] = env.ckpt().read_bytes()
        shutil.rmtree(env.out)
    assert out[1] == out[3]


# ================================================================ 2. simultaneous stops
class TriggerClock:
    """0 until the 6th call (start + evaluations at done = 50..54), then huge."""

    def __init__(self):
        self.n = 0

    def __call__(self):
        self.n += 1
        return 0.0 if self.n <= 6 else 1e9


def simultaneous(env):
    B["bad"] = {29: "unreliable_apparent_failed", 39: "unreliable_apparent_failed",
                54: "unreliable_inner_shortfall"}           # 3rd at n = 55: 60 > 55
    with pytest.raises(orch.RuntimeBudgetHalt) as ex:
        env.run(clock=TriggerClock())
    return ex.value


def test_simultaneous_stop_records_both_reasons(env):
    ex = simultaneous(env)
    assert ex.reasons == ["nonestimable_rate", "runtime_projection"]
    ends = list((env.out / "halts").glob("run0.inv*.end.json"))
    assert len(ends) == 1
    end = json.loads(ends[0].read_text())
    assert end["outcome"] == "runtime_halt"
    assert end["stop_reasons"] == ["nonestimable_rate", "runtime_projection"]
    ne = end["nonestimable"]
    assert (ne["n_nonestimable"], ne["n_completed"], ne["lhs"], ne["rhs"], ne["result"]) == (3, 55, 60, 55, True)
    assert ne["triggering_k"] == 54 and ne["triggering_position"] == 55
    assert ne["k_nonestimable"] == [29, 39, 54]
    assert end["details"]["reason"] == "runtime_projection"
    assert end["checkpoint_size_bytes"] == env.ckpt().stat().st_size


def test_simultaneous_stop_cannot_be_resumed_by_runtime_authorisation(env):
    simultaneous(env)
    for ref in (None, "R1"):
        with pytest.raises(orch.DispatchRefused, match="non-estimable"):
            env.run(resume_authorization=ref)
    assert not env.halts(0, 1, "auth").exists()
    assert orch.invocation_state(env.out, 0)["state"] == "HALTED"
    e = env.elig()
    assert e["class"] == "INELIGIBLE" and e["permanent"] is True
    codes = {r["code"] for r in e["reasons"]}
    assert {"nonestimable_prefix_violation", "history_nonestimable_stop"} <= codes


def test_label_independent_prefix_scan(env):
    simultaneous(env)
    env.rewrite(0, 1, "end", lambda o: (o.pop("stop_reasons"), o.pop("nonestimable")))
    e = env.elig()
    assert e["class"] == "INELIGIBLE"
    assert "nonestimable_prefix_violation" in {r["code"] for r in e["reasons"]}


def test_diluted_forged_history_is_ineligible(env):
    B["bad"] = {0: "unreliable_apparent_failed"}
    with pytest.raises(orch.NonEstimableStop):
        env.run()
    store_recs = CheckpointStore.load_records(env.ckpt())
    with open(env.ckpt(), "ab") as fh:                         # forged, well-formed appends
        for k in range(1, 40):
            r = {**store_recs[0], "k": k, "status": "ok"}
            fh.write(orch._canon({"kind": "record", "run_id": 0, "k": k,
                                  "record_sha256": orch._sha(orch.canonical_json(r)),
                                  "record": r}) + b"\n")
    e = env.elig()
    codes = {r["code"] for r in e["reasons"]}
    assert e["class"] == "INELIGIBLE" and "nonestimable_prefix_violation" in codes
    assert "chain_break" in codes


# ================================================================ 3. eligibility classes (13)
def g_scenario(env):
    """Runtime halt, authorisation committed, crash before the next start (CLEAN, 55 rec)."""
    halt_at_55(env)
    real = orch.log_event

    def log(out, ev):
        if ev.get("event") == "resume_authorized":
            raise RuntimeError("crash after authorisation commit")
        return real(out, ev)
    orch.log_event = log
    try:
        with pytest.raises(RuntimeError):
            env.run(resume_authorization="R1")
    finally:
        orch.log_event = real


def complete_with_window(env):
    B["raise_at"] = 120
    with pytest.raises(RuntimeError):
        env.run()
    B["raise_at"] = None
    env.run(resume_authorization="ACK1")


def record_preservation_failure(env):
    env.pres_fail_next = True
    with pytest.raises(SystemExit):
        env.run(run_id=1)


def check(env, e, cls, permanent, dispatch):
    assert (e["class"], e["permanent"], e["dispatch"]) == (cls, permanent, dispatch), e["reasons"]
    assert e["permanent"] == (e["class"] == "INELIGIBLE")


def test_c01_incomplete_clean(env):
    g_scenario(env)
    check(env, env.elig(), "INCOMPLETE", False, "permitted_without_authorization")
    assert env.run()["n_records"] == 500                       # actual dispatch permitted


def test_c02_incomplete_halted(env):
    halt_at_55(env)
    check(env, env.elig(), "INCOMPLETE", False, "permitted_with_authorization")
    with pytest.raises(orch.InvocationStateError):
        env.run()
    assert env.run(resume_authorization="R")["n_records"] == 500


def test_c03_incomplete_unterminated(env):
    B["raise_at"] = 100
    with pytest.raises(RuntimeError):
        env.run()
    B["raise_at"] = None
    check(env, env.elig(), "INCOMPLETE", False, "permitted_with_authorization")
    with pytest.raises(orch.InvocationStateError):
        env.run()
    assert env.run(resume_authorization="A")["n_records"] == 500


def test_c04_incomplete_ambiguous(env):
    halt_at_55(env)
    stray_tmp(env)
    check(env, env.elig(), "UNVERIFIABLE", False, "blocked_pending_owner")
    with pytest.raises(orch.InvocationStateError, match="AMBIGUOUS"):
        env.run(resume_authorization="R")


def test_c05_complete_ambiguous(env):
    env.run()
    stray_tmp(env)
    check(env, env.elig(), "UNVERIFIABLE", False, "blocked_pending_owner")
    with pytest.raises(orch.InvocationStateError):
        env.run()


def test_c06_complete_window_unrecomputed_then_recomputed(env):
    complete_with_window(env)
    e = env.elig()
    check(env, e, "UNVERIFIABLE", False, "not_applicable")
    assert "window_not_recomputed" in {r["code"] for r in e["reasons"]}
    assert env.run()["n_records"] == 500                       # nothing blocked (no-op resume)
    ack_window(env, 0, 1)
    e2 = env.elig()
    check(env, e2, "ELIGIBLE", False, "not_applicable")
    assert "window_recomputed_bit_identical" in {d["code"] for d in e2["disclosures"]}


def test_c06b_window_fail_and_unbound_results(env):
    complete_with_window(env)
    res = ack_window(env, 0, 1)
    res["per_k"][str(res["k_list"][0])]["record_identical"] = False
    res["status"] = "FAIL"
    orch.window_result_path(env.out, 0, 1).write_text(json.dumps(res))
    check(env, env.elig(), "INELIGIBLE", True, "blocked_permanently")
    res["status"], res["region_sha256"] = "PASS", "0" * 64
    orch.window_result_path(env.out, 0, 1).write_text(json.dumps(res))
    e = env.elig()
    assert e["class"] == "UNVERIFIABLE" and "window_result_not_bound" in {r["code"] for r in e["reasons"]}


def test_c07_incomplete_preservation_pending(env):
    halt_at_55(env)
    record_preservation_failure(env)
    check(env, env.elig(), "UNVERIFIABLE", False, "blocked_pending_owner")
    with pytest.raises(orch.DispatchRefused, match="preservation"):
        env.run(resume_authorization="R")


def test_c08_h4_plus_interrupted(env):
    B["raise_at"] = 100
    with pytest.raises(RuntimeError):
        env.run()
    B["raise_at"] = None
    env.h4_fail_stage = "during run 0 k=150"
    with pytest.raises(HardFailureH4):
        env.run(resume_authorization="A")
    env.h4_fail_stage = None
    check(env, env.elig(), "INELIGIBLE", True, "blocked_permanently")
    with pytest.raises(orch.DispatchRefused):
        env.run()


def test_c09_h6_plus_missing_authorisation(env):
    halt_at_55(env)
    B.update(raise_at=80, exc=ImplementationDefectH6)
    with pytest.raises(ImplementationDefectH6):
        env.run(resume_authorization="R1")
    B.update(raise_at=None, exc=RuntimeError)
    assert json.loads(env.halts(0, 2, "end").read_text())["outcome"] == "hard_failure_h6"
    env.halts(0, 1, "auth").unlink()                           # authorisation missing
    e = env.elig()
    check(env, e, "INELIGIBLE", True, "blocked_permanently")
    codes = {r["code"] for r in e["reasons"]}
    assert {"hard_failure_h6", "invocation_state_ambiguous"} <= codes


def test_c10_nonestimable_plus_broken_chain(env):
    B["bad"] = {0: "unreliable_apparent_failed"}
    with pytest.raises(orch.NonEstimableStop):
        env.run()
    with open(env.ckpt(), "ab") as fh:
        fh.write(env.ckpt().read_bytes().split(b"\n")[1] + b"\n")   # duplicate line appended
    e = env.elig()
    check(env, e, "INELIGIBLE", True, "blocked_permanently")
    assert {"history_nonestimable_stop", "duplicate_k"} <= {r["code"] for r in e["reasons"]}


def test_c11_nonestimable_plus_ambiguous(env):
    B["bad"] = {0: "unreliable_apparent_failed"}
    with pytest.raises(orch.NonEstimableStop):
        env.run()
    stray_tmp(env)
    e = env.elig()
    check(env, e, "INELIGIBLE", True, "blocked_permanently")
    assert {"history_nonestimable_stop", "invocation_state_ambiguous"} <= {r["code"] for r in e["reasons"]}


def test_c12_complete_preservation_pending(env):
    env.run()
    record_preservation_failure(env)
    check(env, env.elig(), "UNVERIFIABLE", False, "blocked_pending_owner")
    with pytest.raises(orch.DispatchRefused, match="preservation"):
        env.run()


def test_c13_all_pass(env):
    env.run()
    e = env.elig()
    check(env, e, "ELIGIBLE", False, "not_applicable")
    assert e["reasons"] == [] and orch._HEX64.match(e["eligibility_result_sha256"])
    assert env.elig()["eligibility_result_sha256"] == e["eligibility_result_sha256"]   # pure


def test_eligibility_target_mismatch_and_incomplete_target(env):
    env.run()
    assert env.elig(target=750)["class"] == "INCOMPLETE"
    assert env.elig(target=400)["class"] == "INELIGIBLE"


def test_stage_b_performs_no_fitting(env):
    env.run()
    calls = B["calls"]
    env.elig()
    assert B["calls"] == calls


# ================================================================ 4. checkpoint continuity
def test_valid_multi_invocation_chain_eligible(env):
    halt_at_55(env)
    halt_at_55(env, resume_authorization="R1")
    env.run(resume_authorization="R2")
    check(env, env.elig(), "ELIGIBLE", False, "not_applicable")


def test_external_append_between_invocations(env):
    halt_at_55(env)
    recs = CheckpointStore.load_records(env.ckpt())
    r = {**recs[0], "k": 55}
    with open(env.ckpt(), "ab") as fh:
        fh.write(orch._canon({"kind": "record", "run_id": 0, "k": 55,
                              "record_sha256": orch._sha(orch.canonical_json(r)), "record": r}) + b"\n")
    with pytest.raises(orch.DispatchRefused, match="has changed"):
        env.run(resume_authorization="R1")
    assert not env.halts(0, 1, "auth").exists()                # authorisation not consumed
    assert "chain_break" in {x["code"] for x in env.elig()["reasons"]}


def test_external_append_during_invocation(env):
    def side(run_id, k):
        if k == 10:
            r = {"run_id": 0, "k": 999, "status": "ok", "B_in_attempted": 0, "inner": [],
                 "dC_corrected": 0.0, "dC_app": 0.0}
            with open(env.ckpt(), "ab") as fh:
                fh.write(orch._canon({"kind": "record", "run_id": 0, "k": 999,
                                      "record_sha256": orch._sha(orch.canonical_json(r)),
                                      "record": r}) + b"\n")
    B["side"] = side
    env.run()
    e = env.elig()
    assert e["class"] == "INELIGIBLE"
    # the sealed end boundary (size from the end record) no longer equals the boundary
    # implied by its own record count: detected as a chain break
    assert {"chain_break", "region_mismatch"} & {x["code"] for x in e["reasons"]}


@pytest.mark.parametrize("damage,code", [
    ("drop_record", "chain_break"), ("wrong_run_id", "checkpoint_region_invalid"),
    ("bad_hash", "checkpoint_region_invalid"), ("end_size", "chain_break"),
    ("end_hash", "chain_break"), ("rewrite", "chain_break")])
def test_checkpoint_and_boundary_damage(env, damage, code):
    env.run()
    raw = env.ckpt().read_bytes()
    lines = raw.split(b"\n")[:-1]
    if damage == "drop_record":
        env.ckpt().write_bytes(b"\n".join(lines[:200] + lines[201:]) + b"\n")
    elif damage in ("wrong_run_id", "bad_hash"):
        o = json.loads(lines[5])
        if damage == "wrong_run_id":
            o["run_id"] = 3
        else:
            o["record"]["dC_corrected"] = 1.0
        lines[5] = orch._canon(o)
        env.ckpt().write_bytes(b"\n".join(lines) + b"\n")
    elif damage == "end_size":
        env.rewrite(0, 1, "end", lambda o: o.update(checkpoint_size_bytes=o["checkpoint_size_bytes"] - 1))
    elif damage == "end_hash":
        env.rewrite(0, 1, "end", lambda o: o.update(checkpoint_sha256="0" * 64))
    else:
        o = json.loads(lines[7])
        o["record"]["dC_corrected"] = 1.0
        o["record_sha256"] = orch._sha(orch.canonical_json(o["record"]))
        lines[7] = orch._canon(o)
        env.ckpt().write_bytes(b"\n".join(lines) + b"\n")
    e = env.elig()
    assert e["class"] == "INELIGIBLE" and code in {x["code"] for x in e["reasons"]}, e["reasons"]


def test_partial_tail_blocks_and_is_unverifiable(env):
    halt_at_55(env)
    with open(env.ckpt(), "ab") as fh:
        fh.write(b'{"kind":"rec')
    e = env.elig()
    check(env, e, "UNVERIFIABLE", False, "blocked_pending_owner")
    assert "partial_tail" in {x["code"] for x in e["reasons"]}
    with pytest.raises(CheckpointCorruption):
        env.run(resume_authorization="R1")


def test_unterminated_acknowledged_and_boundary_mismatch(env):
    g_scenario(env)                                            # auth committed, no next start
    recs = CheckpointStore.load_records(env.ckpt())
    r = {**recs[0], "k": 55}
    with open(env.ckpt(), "ab") as fh:
        fh.write(orch._canon({"kind": "record", "run_id": 0, "k": 55,
                              "record_sha256": orch._sha(orch.canonical_json(r)), "record": r}) + b"\n")
    with pytest.raises(orch.DispatchRefused, match="has changed"):
        env.run()


def test_auth_and_end_records_carry_sizes(env):
    halt_at_55(env)
    env.run(resume_authorization="R1")
    a = json.loads(env.halts(0, 1, "auth").read_text())
    e1 = json.loads(env.halts(0, 1, "end").read_text())
    assert a["checkpoint_size_bytes"] == e1["checkpoint_size_bytes"]
    assert json.loads(env.halts(0, 2, "start").read_text())["checkpoint_size_bytes"] == a["checkpoint_size_bytes"]


def test_orphan_records_refused(env):
    env.run()                                                  # creates manifest; run 0 done
    # orphan: run 1 checkpoint written outside any invocation
    msha = env.msha()
    s1 = CheckpointStore(env.out, 1, env.cfg().sha256(), msha, ROOT)
    with s1:
        s1.append(fake_record(None, None, None, 1, 0, 0, None))
    with pytest.raises(orch.DispatchRefused, match="orphan_records"):
        env.run(run_id=1)
    assert env.elig(1)["class"] == "INELIGIBLE"


# ================================================================ 5. paths A / B / C and sizing
def complete_five(env):
    for r in range(5):
        env.run(run_id=r)


def write_sizing(env, B_cur):
    art = orch.build_sizing_decision(env.cfg(), B_cur, env.msha(), SYNTH_FREEZE)
    orch.durable_create(orch._sizing_path(env.out, B_cur), art)
    return art


@pytest.fixture
def topup_env(env):
    env.cfg_over = {"w5": W5A, "ladder": LADDER}
    B["drift"] = 0.05                                          # forces a W5-A breach
    complete_five(env)
    art = write_sizing(env, 500)
    assert art["decision"]["action"] == "top_up", art["decision"]["action"]
    env.next_target = art["decision"]["combined_B_next"]
    return env


def test_fresh_path_rules(env):
    with pytest.raises(orch.DispatchRefused, match="fresh run"):
        env.run(target=750, topup_authorization="T")
    with pytest.raises(orch.DispatchRefused, match="fresh run"):
        env.run(topup_authorization="T")
    assert env.run()["n_records"] == 500
    assert json.loads(env.halts(0, 1, "start").read_text())["path"] == "fresh"


def test_valid_topup_and_progression(topup_env):
    env, T = topup_env, topup_env.next_target
    assert T == 750
    res = env.run(run_id=0, target=T, topup_authorization="OWNER-TOPUP-750")
    assert res["n_records"] == T
    s = json.loads(env.halts(0, 2, "start").read_text())
    assert s["path"] == "topup" and orch._HEX64.match(s["sizing_decision_sha256"])
    for r in range(1, 5):
        env.run(run_id=r, target=T, topup_authorization="OWNER-TOPUP-750")
    for r in range(5):
        check(env, env.elig(r, T), "ELIGIBLE", False, "not_applicable")
    assert env.elig(0, 500)["class"] == "INELIGIBLE"           # evaluated below its level
    art2 = write_sizing(env, T)
    assert art2["decision"]["action"] == "top_up" and art2["decision"]["combined_B_next"] == 1000
    env.run(run_id=0, target=1000, topup_authorization="OWNER-TOPUP-1000")
    check(env, env.elig(0, 1000), "ELIGIBLE", False, "not_applicable")
    # repeated top-up at the same level with a different reference / target decrease
    with pytest.raises(orch.DispatchRefused, match="original top-up reference"):
        env.run(run_id=0, target=1000, topup_authorization="SOMETHING-ELSE")
    assert env.run(run_id=0, target=1000, topup_authorization="OWNER-TOPUP-1000")["n_records"] == 1000
    with pytest.raises(orch.DispatchRefused, match="below the current level"):
        env.run(run_id=0, target=750, topup_authorization="OWNER-TOPUP-750")


def test_topup_refusals(topup_env):
    env, T = topup_env, topup_env.next_target
    with pytest.raises(orch.DispatchRefused, match="differs from the independent recomputation|not 1000"):
        env.run(run_id=0, target=1000, topup_authorization="X")        # not the 7B target
    with pytest.raises(ConfigurationError, match="ladder rung"):
        env.run(run_id=0, target=600, topup_authorization="X")
    with pytest.raises(ConfigurationError, match="authorisation"):
        env.run(run_id=0, target=T)
    with pytest.raises(ConfigurationError, match="blank"):
        env.run(run_id=0, target=T, topup_authorization="   ")


def tamper(env, fn):
    p = orch._sizing_path(env.out, 500)
    art = json.loads(p.read_text())
    fn(art)
    p.unlink()
    orch.durable_create(p, art)                                # internally consistent forgery


@pytest.mark.parametrize("what,match", [
    ("decision", "independent recomputation"), ("eligibility", "not reproduced"),
    ("package", "code identity"), ("freeze", "code identity"), ("missing", "missing")])
def test_sizing_file_cannot_self_certify(topup_env, what, match):
    env, T = topup_env, topup_env.next_target
    if what == "decision":
        tamper(env, lambda a: a["decision"]["ep1"].update(rho_max=0.001))
    elif what == "eligibility":
        tamper(env, lambda a: a["runs"]["2"].update(eligibility_result_sha256="0" * 64))
    elif what == "package":
        tamper(env, lambda a: a["code_identity"]["package"].update({"src/x.py": "1" * 64}))
    elif what == "freeze":
        tamper(env, lambda a: a["code_identity"].update(freeze_record_sha256="e" * 64))
    else:
        orch._sizing_path(env.out, 500).unlink()
    with pytest.raises(orch.DispatchRefused, match=match):
        env.run(run_id=0, target=T, topup_authorization="OWNER")
    assert len(CheckpointStore.load_records(env.ckpt(0))) == 500


def test_sizing_refuses_checkpoint_drift(topup_env):
    env, T = topup_env, topup_env.next_target
    lines = env.ckpt(2).read_bytes()
    with open(env.ckpt(2), "ab") as fh:                        # another run drifts
        fh.write(lines.split(b"\n")[1] + b"\n")
    with pytest.raises(orch.DispatchRefused, match="run 2"):
        env.run(run_id=0, target=T, topup_authorization="OWNER")


def test_different_target_across_runs_refused(topup_env):
    env, T = topup_env, topup_env.next_target
    env.run(run_id=0, target=T, topup_authorization="OWNER")
    with pytest.raises(orch.DispatchRefused):
        env.run(run_id=1, target=1000, topup_authorization="OWNER")


def test_topup_from_incomplete_level_refused(env):
    env.cfg_over = {"w5": W5A, "ladder": LADDER}
    g_scenario(env)                                            # CLEAN with 55 of 500
    with pytest.raises(orch.DispatchRefused, match="not complete"):
        env.run(target=750, topup_authorization="T")


def test_forged_level_increase_without_sizing_is_ineligible(topup_env):
    env, T = topup_env, topup_env.next_target
    env.run(run_id=0, target=T, topup_authorization="OWNER")
    env.rewrite(0, 2, "start", lambda o: o.update(sizing_decision_sha256=None))
    e = env.elig(0, T)
    assert e["class"] == "INELIGIBLE" and "target_progression_violation" in {x["code"] for x in e["reasons"]}


# ================================================================ 6. preservation
def test_preservation_evidence_from_actual_check():
    ev = orch.preservation_evidence(historical_preservation_check(ROOT), "post")
    assert ev["malformed"] is False and ev["n_failing"] == 0 and ev["n_checked"] == len(ev["files"])
    assert [f["path"] for f in ev["files"]] == sorted(f["path"] for f in ev["files"])


def test_preservation_evidence_multiple_and_missing_files():
    res = {"all_ok": False, "files": {
        "b.json": {"expected_prefix": "11", "sha256": "22" * 32, "ok": False},
        "a.json": {"expected_prefix": "33", "sha256": None, "ok": False},
        "c.json": {"expected_prefix": "44", "sha256": "44" * 32, "ok": True}}}
    ev = orch.preservation_evidence(res, "post")
    assert [f["path"] for f in ev["files"]] == ["a.json", "b.json", "c.json"]
    assert ev["n_failing"] == 2 and ev["files"][0]["exists"] is False
    assert orch._evidence_malformed(ev) is False
    assert orch.preservation_evidence({"all_ok": False}, "post")["malformed"] is True


def test_post_run_preservation_failure_records_structured_evidence(env):
    env.pres_result = {"all_ok": False, "status": "FAIL", "files": {
        "x/a.json": {"expected_prefix": "aa", "sha256": "bb" * 32, "ok": False},
        "x/b.json": {"expected_prefix": "cc", "sha256": None, "ok": False}}}
    calls = {"n": 0}
    real = orch.historical_preservation_check

    def second_call_fails(root):
        calls["n"] += 1
        if calls["n"] == 2:
            env.pres_fail_next = True
        return real(root)
    orch.historical_preservation_check = second_call_fails
    try:
        with pytest.raises(orch.HistoricalPreservationFailure):
            env.run()
    finally:
        orch.historical_preservation_check = real
    d = json.loads(env.halts(0, 1, "end").read_text())["details"]
    assert d["malformed"] is False and d["n_failing"] == 2
    assert [f["path"] for f in d["files"]] == ["x/a.json", "x/b.json"]
    with pytest.raises(orch.DispatchRefused, match="preservation"):
        env.run()


def test_malformed_preservation_evidence_still_blocks(env):
    env.run()
    env.pres_result = {"all_ok": False}                        # no file list
    env.pres_fail_next = True
    with pytest.raises(SystemExit):
        env.run(run_id=1)
    e = env.elig()
    assert {"preservation_failure_pending", "preservation_evidence_malformed"} <= \
        {x["code"] for x in e["reasons"]}
    with pytest.raises(orch.DispatchRefused):
        env.run()


def test_global_preservation_block_path_a(env):
    env.pres_fail_next = True
    with pytest.raises(SystemExit):
        env.run()                                              # fresh out dir: recorded
    assert (env.out / "preservation" / "failure.0001.json").exists()
    with pytest.raises(orch.DispatchRefused, match="preservation"):
        env.run()                                              # restored files: still blocked
    with pytest.raises(orch.DispatchRefused, match="preservation"):
        env.run(run_id=3)


def test_global_preservation_block_path_b(env):
    halt_at_55(env)
    record_preservation_failure(env)
    with pytest.raises(orch.DispatchRefused, match="preservation"):
        env.run(resume_authorization="R1")


def test_global_preservation_block_path_c(topup_env):
    env, T = topup_env, topup_env.next_target
    env.pres_fail_next = True
    with pytest.raises(SystemExit):
        env.run(run_id=4, target=T, topup_authorization="OWNER")
    with pytest.raises(orch.DispatchRefused, match="preservation"):
        env.run(run_id=0, target=T, topup_authorization="OWNER")


# ================================================================ 7. H4 / H6 terminal behaviour
def test_h6_terminal_and_refusal(env):
    B.update(raise_at=10, exc=ImplementationDefectH6)
    with pytest.raises(ImplementationDefectH6):
        env.run()
    B.update(raise_at=None, exc=RuntimeError)
    end = json.loads(env.halts(0, 1, "end").read_text())
    assert end["outcome"] == "hard_failure_h6" and "hard_failure_h6" in orch.OUTCOMES
    assert orch.invocation_state(env.out, 0)["state"] == "CLEAN"
    for ref in (None, "R"):
        with pytest.raises(orch.DispatchRefused):
            env.run(resume_authorization=ref)
    check(env, env.elig(), "INELIGIBLE", True, "blocked_permanently")


@pytest.mark.parametrize("stage", ["during run 0 k=20", "post"])
def test_h4_during_or_after_execution_is_permanent(env, stage):
    env.h4_fail_stage = stage
    with pytest.raises(HardFailureH4):
        env.run()
    env.h4_fail_stage = None
    with pytest.raises(orch.DispatchRefused):
        env.run()
    check(env, env.elig(), "INELIGIBLE", True, "blocked_permanently")


def test_h4_and_gate_before_dispatch_leave_no_history(env, monkeypatch):
    env.h4_fail_stage = "pre"
    with pytest.raises(HardFailureH4):
        env.run()
    env.h4_fail_stage = None
    monkeypatch.setattr(orch, "run_pre_execution_gate", make_synthetic_gate(fail_check="V2_numpy"))
    with pytest.raises(SystemExit, match="gate"):
        env.run()
    assert not (env.out / "halts").exists()
    monkeypatch.setattr(orch, "run_pre_execution_gate", synthetic_gate)
    assert env.run()["n_records"] == 500                       # not disqualifying


# ================================================================ 8. code identity / freeze chain
FROZEN_FILES = [f"src/phase9_corrected/{p.name}" for p in (ROOT / "src/phase9_corrected").glob("*.py")] + \
    ["scripts/41_phase9_corrected_driver.py", "src/phase9_execution.py",
     "src/pathway_scoring.py", "src/score_standardisation.py",
     "scripts/42_phase9_corrected_launch.sh",                  # DEC-E1: launcher is frozen
     # Decision 13 (B1): the expanded mandatory inventory
     "scripts/43_phase9_corrected_verify.py", "scripts/44_phase9_freeze_prepare.sh",
     "scripts/44_phase9_freeze_capture.py", "scripts/40_phase9_dc_interval.py",
     "src/__init__.py", "src/paths.py", "src/reproducibility.py"]


@pytest.fixture(autouse=True)
def _e1_runtime_flags_stub(monkeypatch, tmp_path_factory):
    """DEC-E1 flag-reader stub: pytest is not launched with -I -B -X pycache_prefix, so the
    real interpreter flags are replaced by E1-compliant values (fresh empty external dir).
    The real flag checks are exercised in tests/test_phase9_corrected_d11_remediation.py."""
    prefix = tmp_path_factory.mktemp("p9c_pyc_prefix")
    monkeypatch.setattr(orch, "_runtime_flags", lambda: {
        "isolated": True, "dont_write_bytecode": True, "pycache_prefix": str(prefix)})


def make_root(tmp_path, policy="checked_hash_pyc", env_edit=None, freeze_edit=None):
    root = tmp_path / "repo"
    for rel in FROZEN_FILES:
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, root / rel)
    envd = orch._jsonable(orch.current_environment())
    if env_edit:
        env_edit(envd)
    freeze = {"schema": orch.FREEZE_SCHEMA,
              "code": {rel: hashlib.sha256((root / rel).read_bytes()).hexdigest() for rel in FROZEN_FILES},
              "environment": envd, "bytecode_policy": policy}
    if freeze_edit:
        freeze_edit(freeze)
    fp = root / orch.FREEZE_RECORD_RELPATH
    fp.parent.mkdir(parents=True, exist_ok=True)
    fp.write_bytes(json.dumps(freeze, sort_keys=True).encode())
    return root, hashlib.sha256(fp.read_bytes()).hexdigest(), freeze


def cfg_for(root, anchor, **d13):
    """DEC-B2: `freeze_record_sha256` is the anchor; `sha256` is the driver's own hash."""
    raw = approved_test_config(root, "scratch_unused")
    drv = root / "scripts/41_phase9_corrected_driver.py"
    driver_sha = hashlib.sha256(drv.read_bytes()).hexdigest() if drv.exists() else "0" * 64
    raw["owner_decisions_outside_schema"] = {"D13_driver_freeze": {
        "approved": True, "date": "2026-11-06", "reference": "TEST FIXTURE",
        "driver_path": "scripts/41_phase9_corrected_driver.py",
        "freeze_record_sha256": anchor, "sha256": driver_sha, **d13}}
    return CorrectedConfig(raw, root)


def test_freeze_chain_passes_on_matching_synthetic_freeze(tmp_path):
    root, anchor, freeze = make_root(tmp_path)
    info = REAL_VERIFY(cfg_for(root, anchor))
    assert info["freeze_record_sha256"] == anchor
    assert REAL_VERIFY(cfg_for(root, anchor), manifest={"code": {FROZEN_FILES[0]: freeze["code"][FROZEN_FILES[0]]}})


@pytest.mark.parametrize("case,match", [
    ("unapproved", "not approved"), ("badhex", "64 hex"), ("missing", "missing or unreadable"),
    ("altered", "does not match the approved"), ("schema", "schema invalid"),
    ("extra_pkg", "package file set"), ("modified", "differs from the frozen hash"),
    ("env", "environment distributions"), ("manifest", "manifest code hash"),
    ("manifest_empty", "no code hashes"), ("policy", "bytecode_policy")])
def test_freeze_chain_fails_closed(tmp_path, case, match):
    edit = {"env": lambda e: e["distributions"].update({"numpy": ["0.0.1"]}),
            "schema": None, "policy": None}.get(case)
    fedit = {"schema": lambda f: f.update(schema="other"),
             "policy": lambda f: f.update(bytecode_policy="whatever")}.get(case)
    root, anchor, freeze = make_root(tmp_path, env_edit=edit if case == "env" else None, freeze_edit=fedit)
    kw, manifest = {}, None
    if case == "unapproved":
        kw = {"approved": False}
    elif case == "badhex":
        kw = {"sha256": "XYZ"}
    elif case == "missing":
        (root / orch.FREEZE_RECORD_RELPATH).unlink()
    elif case == "altered":
        (root / orch.FREEZE_RECORD_RELPATH).write_text("{}")
    elif case == "extra_pkg":
        (root / "src/phase9_corrected/extra_module.py").write_text("x = 1\n")
    elif case == "modified":
        with open(root / "src/phase9_corrected/statistics.py", "a") as fh:
            fh.write("\n# modified\n")
    elif case == "manifest":
        manifest = {"code": {FROZEN_FILES[0]: "0" * 64}}
    elif case == "manifest_empty":
        manifest = {"code": {}}
    with pytest.raises(orch.CodeIdentityRefusal, match=match):
        REAL_VERIFY(cfg_for(root, anchor, **kw), manifest=manifest)


def test_thread_environment_is_pinned(tmp_path, monkeypatch):
    root, anchor, _ = make_root(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "3")
    with pytest.raises(orch.CodeIdentityRefusal, match="thread_env"):
        REAL_VERIFY(cfg_for(root, anchor))


def test_bytecode_policies(tmp_path, monkeypatch):
    root, anchor, _ = make_root(tmp_path, policy="no_pyc")
    monkeypatch.setattr(orch, "_dont_write_bytecode", lambda: False)
    with pytest.raises(orch.CodeIdentityRefusal, match="PYTHONDONTWRITEBYTECODE"):
        REAL_VERIFY(cfg_for(root, anchor))
    monkeypatch.setattr(orch, "_dont_write_bytecode", lambda: True)
    REAL_VERIFY(cfg_for(root, anchor))
    src = root / "src/phase9_corrected/config.py"
    py_compile.compile(str(src), cfile=importlib.util.cache_from_source(str(src)),
                       invalidation_mode=py_compile.PycInvalidationMode.TIMESTAMP)
    with pytest.raises(orch.CodeIdentityRefusal, match="stale bytecode"):
        REAL_VERIFY(cfg_for(root, anchor))


def test_checked_hash_bytecode_policy(tmp_path):
    root, anchor, _ = make_root(tmp_path, policy="checked_hash_pyc")
    src = root / "src/phase9_corrected/config.py"
    pyc = importlib.util.cache_from_source(str(src))
    py_compile.compile(str(src), cfile=pyc, invalidation_mode=py_compile.PycInvalidationMode.TIMESTAMP)
    with pytest.raises(orch.CodeIdentityRefusal, match="not checked-hash"):
        REAL_VERIFY(cfg_for(root, anchor))
    py_compile.compile(str(src), cfile=pyc, invalidation_mode=py_compile.PycInvalidationMode.CHECKED_HASH)
    REAL_VERIFY(cfg_for(root, anchor))


def test_optional_distribution_record_check(tmp_path):
    import importlib.metadata as md
    dist = md.distribution("numpy")
    rec = next(f for f in dist.files if f.name == "RECORD")
    good = hashlib.sha256(Path(dist.locate_file(rec)).read_bytes()).hexdigest()
    root, anchor, _ = make_root(tmp_path, freeze_edit=lambda f: f.update(distribution_records={"numpy": good}))
    REAL_VERIFY(cfg_for(root, anchor))
    root2, anchor2, _ = make_root(tmp_path / "b", freeze_edit=lambda f: f.update(distribution_records={"numpy": "0" * 64}))
    with pytest.raises(orch.CodeIdentityRefusal, match="RECORD"):
        REAL_VERIFY(cfg_for(root2, anchor2))


def test_real_freeze_is_the_approved_anchor_and_verification_fails_closed():
    """Updated after Decision 13 (approved 2026-11-06). Before approval the real freeze record
    was absent and a configuration anchored to any hash failed as 'missing or unreadable'. Now
    the approved record is installed: it must match the approved anchor, and a configuration
    anchored to any other hash must still be refused."""
    cfg = json.loads((ROOT / "config/analysis/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json").read_text())
    anchor = cfg["owner_decisions_outside_schema"]["D13_driver_freeze"]["freeze_record_sha256"]
    assert hashlib.sha256((ROOT / orch.FREEZE_RECORD_RELPATH).read_bytes()).hexdigest() == anchor
    with pytest.raises(orch.CodeIdentityRefusal, match="does not match the approved Decision 13 hash"):
        REAL_VERIFY(cfg_for(ROOT, "a" * 64))
    with pytest.raises(orch.CodeIdentityRefusal, match="not approved"):
        REAL_VERIFY(CorrectedConfig(approved_test_config(ROOT, "scratch_unused"), ROOT))


def test_code_identity_enforced_on_all_three_paths(topup_env, monkeypatch):
    env, T = topup_env, topup_env.next_target
    fresh = Env(f"scratch_p9c_ne_{uuid.uuid4().hex[:10]}")
    monkeypatch.setattr(orch, "verify_code_identity", REAL_VERIFY)
    try:
        with pytest.raises(orch.CodeIdentityRefusal):              # Path A
            fresh.run()
        assert not fresh.out.exists()                              # nothing written
        with pytest.raises(orch.CodeIdentityRefusal):              # Path B (no-op resume)
            env.run(run_id=0)
        with pytest.raises(orch.CodeIdentityRefusal):              # Path C
            env.run(run_id=0, target=T, topup_authorization="OWNER")
    finally:
        shutil.rmtree(fresh.out, ignore_errors=True)


# ================================================================ 9. durability / diagnostics / driver
def test_directory_levels_created_and_synced_in_order(tmp_path, monkeypatch):
    seen = []
    real = orch._sync_path
    monkeypatch.setattr(orch, "_sync_path", lambda p: (seen.append(Path(p)), real(p)))
    target = tmp_path / "a" / "b" / "c" / "x.json"
    orch.durable_create(target, {"x": 1})
    assert seen[:3] == [tmp_path, tmp_path / "a", tmp_path / "a" / "b"]
    assert seen[3:] == [tmp_path / "a" / "b" / "c", tmp_path / "a" / "b" / "c"]
    assert json.loads(target.read_text()) == {"x": 1}


def test_diagnostic_report_is_accounting_only(env):
    B["bad"] = {0: "unreliable_apparent_failed"}
    with pytest.raises(orch.NonEstimableStop):
        env.run()
    name, text = orch.diagnostic_report(env.elig())
    assert name.startswith("diagnostics/NOT_ELIGIBLE__NOT_A_SCIENTIFIC_RESULT")
    assert text.splitlines()[0] == orch.NOT_ELIGIBLE_MARKER
    body = json.loads(text.split("\n", 1)[1])
    assert not {"interval", "ep1", "w5", "percentile", "H1", "H2"} & set(body)
    shutil.rmtree(env.out)
    B["bad"] = {}
    env.run()
    with pytest.raises(ValueError):                            # never for an ELIGIBLE run
        orch.diagnostic_report(env.elig())


@pytest.fixture(scope="module")
def driver():
    spec = importlib.util.spec_from_file_location(
        "p9c_driver_d11", ROOT / "scripts/41_phase9_corrected_driver.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.mark.parametrize("tail", [["--topup-authorization"], ["--topup-authorization", "  "],
                                  ["--topup-authorization", "--resume-authorization", "R"]])
def test_driver_rejects_missing_or_blank_topup_reference(driver, tmp_path, tail, monkeypatch):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(approved_test_config(ROOT, "scratch_unused")))
    called = []
    monkeypatch.setattr(orch, "run_corrected", lambda *a, **k: called.append(1))
    assert driver.main(["run", str(p), "0", "750", *tail]) == 2 and not called
