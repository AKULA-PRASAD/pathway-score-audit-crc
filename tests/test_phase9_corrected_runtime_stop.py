"""Owner Decision 10 (2026-11-06): runtime budget, warm-up, clock, invocation sentinel and
halt/unterminated authorisation. Fake clocks and synthetic records only: the scientific
fit (compute_outer_replicate), data loading, gate, H4 and manifest are stubbed, so no
scientific data is read and no real bootstrap replicate is computed."""
import importlib.util
import json
import math
import shutil
import sys
import time
import uuid
from pathlib import Path

import pytest

from src.phase9_corrected import orchestration as orch
from src.phase9_corrected.checkpoint import CheckpointStore
from src.phase9_corrected.config import ConfigurationError, CorrectedConfig
from src.phase9_corrected.integrity import HardFailureH4
from tests.p9c_support import approved_test_config, synthetic_gate

ROOT = Path(__file__).resolve().parents[1]
REAL_LINK, REAL_UNLINK = orch.os.link, orch.os.unlink
DARWIN = sys.platform == "darwin"

# ---------------------------------------------------------------- fakes
BEHAVIOUR = {"bad": {}, "raise_at": None, "exc": RuntimeError, "calls": 0}


def fake_record(ctx, streams, pipeline, run_id, k, B_in, lam_grid, as_run_variant=False):
    BEHAVIOUR["calls"] += 1
    if BEHAVIOUR["raise_at"] is not None and k == BEHAVIOUR["raise_at"]:
        raise BEHAVIOUR["exc"](f"simulated crash at k={k}")
    return {"run_id": int(run_id), "k": int(k), "status": BEHAVIOUR["bad"].get(k, "ok"),
            "value": k * 0.5}


class StepClock:
    """Returns 0, step, 2*step, ... on successive calls; counts calls."""

    def __init__(self, step):
        self.step, self.n = step, 0

    def __call__(self):
        v = self.n * self.step
        self.n += 1
        return v


class SeqClock:
    """Returns the given values in order (start value first)."""

    def __init__(self, values):
        self.values, self.n = list(values), 0

    def __call__(self):
        v = self.values[self.n]
        self.n += 1
        return v


class Env:
    def __init__(self, out_rel):
        self.out_rel = out_rel
        self.out = ROOT / out_rel
        self.h4_fail_stage = None
        self.preservation_calls = 0
        self.preservation_fail_on_call = None

    def cfg(self, n_workers=1, budget=4.0):
        return CorrectedConfig(approved_test_config(ROOT, self.out_rel, n_workers=n_workers,
                                                    runtime_budget_hours_per_run=budget), ROOT)

    def run(self, clock_step=0.001, n_workers=1, clock=None, **kw):
        clk = clock if clock is not None else StepClock(clock_step)
        return orch.run_corrected(self.cfg(n_workers), 0, 500, clock=("TEST_CLOCK", clk), **kw)

    @property
    def ckpt(self):
        return self.out / "runs" / "run0.jsonl"

    def state(self):
        return orch.invocation_state(self.out, 0)

    def records(self, kind):
        return {int(p.name.split(".inv")[1][:4]): json.loads(p.read_text())
                for p in (self.out / "halts").glob(f"run0.inv*.{kind}.json")}


@pytest.fixture
def env(monkeypatch):
    e = Env(f"scratch_p9c_rt_{uuid.uuid4().hex[:10]}")
    BEHAVIOUR.update(bad={}, raise_at=None, exc=RuntimeError, calls=0)

    def h4(root, stage, raise_on_fail=True):
        if stage == e.h4_fail_stage:
            raise HardFailureH4(f"HARD FAILURE H4 ({stage}) -- simulated")
        return {"status": "PASS", "all_ok": True}

    def preservation(root):
        e.preservation_calls += 1
        ok = e.preservation_calls != e.preservation_fail_on_call
        return {"all_ok": ok, "status": "PASS" if ok else "FAIL"}

    monkeypatch.setattr(orch, "h4_check", h4)
    monkeypatch.setattr(orch, "historical_preservation_check", preservation)
    monkeypatch.setattr(orch, "build_manifest", lambda config, root: {"schema": "test-manifest",
                                                                      "config": config.sha256()})
    monkeypatch.setattr(orch, "run_pre_execution_gate", synthetic_gate)
    monkeypatch.setattr(orch, "load_context", lambda root: None)
    monkeypatch.setattr(orch, "RealPipeline", lambda: None)
    monkeypatch.setattr(orch, "compute_outer_replicate", fake_record)
    monkeypatch.setattr(orch, "verify_code_identity",            # Decision 11: synthetic freeze
                        lambda config, manifest=None, extra_files=(): {
                            "freeze_record_sha256": "f" * 64, "freeze": {"code": {}}})
    yield e
    shutil.rmtree(e.out, ignore_errors=True)
    assert not e.out.exists()


def uninterrupted_bytes(env, n_workers=1):
    env.run(n_workers=n_workers)
    b = env.ckpt.read_bytes()
    shutil.rmtree(env.out)
    return b


def halt_first_invocation(env):
    with pytest.raises(orch.RuntimeBudgetHalt):
        env.run(clock=StepClock(600))


# ================================================================ 1-8: RuntimeMonitor (unit)
@pytest.mark.parametrize("n_workers,W", [(1, 50), (8, 50), (32, 64)])
def test_warmup_size(n_workers, W):
    assert orch.RuntimeMonitor(4.0, 0, 500, n_workers, StepClock(1)).W == W


class GuardedClock:
    def __init__(self):
        self.reads, self.armed = 0, False

    def __call__(self):
        if self.armed:
            self.reads += 1
        return 0.0 if not self.armed else 1.0


@pytest.mark.parametrize("n_workers", [1, 8, 32])
def test_projection_not_evaluated_before_warmup(n_workers):
    c = GuardedClock()
    m = orch.RuntimeMonitor(4.0, 0, 500, n_workers, c)
    m.start()
    c.armed = True
    for d in range(1, m.W):                            # first record ... W-1: clock never read
        assert m.observe(d) is None
    assert c.reads == 0
    assert m.observe(m.W) is None and c.reads == 1     # first evaluation exactly at W


def test_warmup_counts_only_this_invocation_after_resume():
    c = GuardedClock()
    m = orch.RuntimeMonitor(4.0, 300, 500, 8, c)       # cumulative checkpoint already 300
    m.start()
    c.armed = True
    assert m.n_new == 200 and m.W == 50
    for total in range(301, 350):                      # cumulative >= W, but done_now < W
        assert m.observe(total) is None
    assert c.reads == 0
    m.observe(350)
    assert c.reads == 1


@pytest.mark.parametrize("n_new,expected", [(250, 14400.0), (500, 28800.0), (2500, 144000.0)])
def test_threshold_scaling(n_new, expected):
    m = orch.RuntimeMonitor(4.0, 0, n_new, 8, StepClock(1))
    assert m.threshold_s.hex() == expected.hex()


def test_eight_worker_arrival_pattern_does_not_halt():
    t_rep = 216.0                                      # 3.6 min single-worker replicate
    times = [0.0] + [math.ceil(d / 8) * t_rep for d in range(50, 501)]
    m = orch.RuntimeMonitor(4.0, 0, 500, 8, SeqClock(times))
    m.start()
    out = [m.observe(d) for d in range(1, 501)]
    assert all(o is None for o in out)
    assert math.isclose(7 * t_rep / 50 * 500, 15120.0)  # projection at d = 50 = 4.2 h < 8 h


def test_old_first_record_projection_would_have_halted():
    """Documents the defect: projecting from record 1 with 8 workers overstates ~8x."""
    t_rep = 216.0
    assert t_rep / 1 * 500 > 2 * 4.0 * 3600            # 30 h > 8 h
    m = orch.RuntimeMonitor(4.0, 0, 500, 8, SeqClock([0.0, t_rep]))
    m.start()
    assert m.observe(1) is None                          # suppressed now


def test_projection_halt_above_threshold_and_details():
    m = orch.RuntimeMonitor(4.0, 0, 500, 1, SeqClock([0.0, 3000.0]))
    m.start()
    h = m.observe(50)                                  # 3000/50*500 = 30000 > 28800
    assert h["reason"] == "runtime_projection" and h["done_now"] == 50
    assert h["projection_s"] == 30000.0 and h["threshold_s"] == 28800.0


def test_simulated_sleep_clock_jump_halts():
    vals = [0.0] + [27.2 * d for d in range(50, 60)] + [27.2 * 60 + 30000.0]
    m = orch.RuntimeMonitor(4.0, 0, 500, 8, SeqClock(vals))
    m.start()
    assert all(m.observe(d) is None for d in range(1, 60))
    assert m.observe(60)["projection_s"] > 28800.0


def test_short_invocation_never_evaluates():
    c = GuardedClock()
    m = orch.RuntimeMonitor(4.0, 470, 500, 8, c)
    m.start()
    c.armed = True
    assert m.n_new == 30 and m.runtime_check == "not_evaluated (n_new < W)"
    assert all(m.observe(t) is None for t in range(471, 501)) and c.reads == 0
    assert orch.RuntimeMonitor(4.0, 0, 500, 8, c).runtime_check == "evaluated_after_warmup"


@pytest.mark.parametrize("bad", [0, -1, True, "4", None])
def test_budget_validation_monitor_and_config(bad):
    with pytest.raises(ConfigurationError):
        orch.RuntimeMonitor(bad, 0, 500, 1, StepClock(1))
    c = approved_test_config(ROOT, "scratch_p9c_test_out", runtime_budget_hours_per_run=bad)
    probs = CorrectedConfig(c, ROOT).validate_structure() + CorrectedConfig(c, ROOT).missing_owner_items()
    assert any("runtime_budget_hours_per_run" in p for p in probs)


def test_monitor_rejects_more_records_than_target():
    with pytest.raises(ConfigurationError):
        orch.RuntimeMonitor(4.0, 600, 500, 1, StepClock(1))


# ================================================================ 9: clock selection
@pytest.mark.skipif(not DARWIN, reason="macOS clock")
def test_darwin_clock_is_monotonic_raw():
    name, fn = orch.select_clock()
    assert name == "CLOCK_MONOTONIC_RAW"
    a, b = fn(), fn()
    assert isinstance(a, float) and b >= a


@pytest.mark.skipif(not DARWIN, reason="macOS clock")
def test_darwin_raw_clock_relation_to_uptime():
    """Relationship only: MONOTONIC_RAW >= UPTIME_RAW (UPTIME_RAW excludes sleep). This does
    NOT by itself prove that MONOTONIC_RAW includes sleep (Decision 10F)."""
    assert time.clock_gettime(time.CLOCK_MONOTONIC_RAW) >= time.clock_gettime(time.CLOCK_UPTIME_RAW)


def test_linux_clock_selection(monkeypatch):
    monkeypatch.setattr(time, "CLOCK_BOOTTIME", time.CLOCK_MONOTONIC, raising=False)
    assert orch.select_clock("linux")[0] == "CLOCK_BOOTTIME"
    monkeypatch.delattr(time, "CLOCK_BOOTTIME")
    with pytest.raises(ConfigurationError, match="not available"):
        orch.select_clock("linux")


@pytest.mark.parametrize("plat", ["win32", "cygwin", "freebsd14"])
def test_unsupported_platform_refuses(plat):
    with pytest.raises(ConfigurationError, match="unsupported platform"):
        orch.select_clock(plat)


@pytest.mark.skipif(not DARWIN, reason="macOS clock")
def test_default_clock_recorded_in_start_record_and_event(env):
    orch.run_corrected(env.cfg(), 0, 500)              # default: approved platform clock
    assert env.records("start")[1]["clock"] == "CLOCK_MONOTONIC_RAW"
    ev = [json.loads(l) for l in (env.out / "logs/events.jsonl").read_text().splitlines()]
    assert [e for e in ev if e["event"] == "run_start"][0]["clock"] == "CLOCK_MONOTONIC_RAW"


# ================================================================ 5, 19: normal completion
def test_normal_completion_records_and_clean_state(env):
    clk = StepClock(0.001)
    res = env.run(clock=clk)
    assert res["n_records"] == 500 and res["invocation"] == 1
    assert clk.n == 1 + 451 + 1 + 1                    # start, 451 evaluations, end, result
    st = env.state()
    assert st["state"] == "CLEAN" and st["last_outcome"] == "completed"
    s, e = env.records("start")[1], env.records("end")[1]
    assert s["n_new"] == 500 and s["W"] == 50 and s["threshold_s"] == 28800.0
    assert s["runtime_check"] == "evaluated_after_warmup" and s["clock"] == "TEST_CLOCK"
    assert e["outcome"] == "completed" and e["checkpoint_sha256"] == res["checkpoint_sha256"]
    assert not list((env.out / "halts").glob("*.tmp"))
    res2 = env.run()                                    # clean subsequent invocation: no work
    assert res2["invocation"] == 2 and res2["n_records"] == 500
    assert env.records("start")[2]["runtime_check"] == "not_evaluated (n_new < W)"


# ================================================================ 6, 12-14: halt and authorisation
def test_runtime_halt_checkpoint_integrity_and_state(env):
    halt_first_invocation(env)
    assert issubclass(orch.RuntimeBudgetHalt, orch.StopCondition)
    recs = CheckpointStore.load_records(env.ckpt)       # strict load: no partial tail
    assert sorted(recs) == list(range(55))              # halted at done_now = 55 (see StepClock)
    assert not (env.out / "runs" / "run0.lock").exists()
    end = env.records("end")[1]
    assert end["outcome"] == "runtime_halt" and end["details"]["reason"] == "runtime_projection"
    assert end["checkpoint_sha256"] == orch.hashlib.sha256(env.ckpt.read_bytes()).hexdigest()
    assert env.state()["state"] == "HALTED"


def test_unauthorised_resume_refused_and_nothing_changes(env):
    halt_first_invocation(env)
    before = env.ckpt.read_bytes()
    n_before = BEHAVIOUR["calls"]
    with pytest.raises(orch.InvocationStateError, match="runtime-budget halt"):
        env.run()
    assert env.ckpt.read_bytes() == before and BEHAVIOUR["calls"] == n_before
    assert sorted(env.records("start")) == [1]          # no new invocation was started


@pytest.mark.parametrize("resume_workers", [1, 3])
def test_authorised_resume_is_byte_identical(env, resume_workers):
    ref = uninterrupted_bytes(env)
    halt_first_invocation(env)
    env.run(n_workers=resume_workers, resume_authorization="OWNER-REF-1")
    assert env.ckpt.read_bytes() == ref
    a = env.records("auth")[1]
    assert a["authorization_kind"] == "runtime_halt_authorization"
    assert a["reference"] == "OWNER-REF-1" and a["invocation"] == 1
    assert env.state()["state"] == "CLEAN"


def test_warmup_after_resume_counts_only_new_invocation(env):
    halt_first_invocation(env)
    clk = StepClock(0.001)
    env.run(clock=clk, resume_authorization="R")
    s = env.records("start")[2]
    assert s["n_start"] == 55 and s["n_new"] == 445
    assert clk.n == 1 + (500 - (55 + 50) + 1) + 1 + 1   # evaluations start at done_now = 50


def test_checkpoint_mismatch_after_halt_refused(env):
    halt_first_invocation(env)
    lines = env.ckpt.read_bytes().split(b"\n")
    env.ckpt.write_bytes(b"\n".join(lines[:-2]) + b"\n")   # drop the last record (valid prefix)
    with pytest.raises(orch.InvocationStateError, match="checkpoint changed"):
        env.run(resume_authorization="R")
    assert env.state()["state"] == "HALTED" and not env.records("auth")


def test_authorisation_refused_when_nothing_pending(env):
    with pytest.raises(orch.InvocationStateError, match="no runtime halt"):
        env.run(resume_authorization="R")               # fresh run: nothing pending
    shutil.rmtree(env.out)
    halt_first_invocation(env)
    env.run(resume_authorization="R1")
    with pytest.raises(orch.InvocationStateError, match="no runtime halt"):
        env.run(resume_authorization="R1")              # duplicate submission


@pytest.mark.parametrize("blank", ["", "   "])
def test_blank_reference_refused(env, blank):
    with pytest.raises(orch.InvocationStateError, match="blank"):
        env.run(resume_authorization=blank)


def test_second_halt_needs_new_authorisation_and_refs_are_single_use(env):
    halt_first_invocation(env)
    with pytest.raises(orch.RuntimeBudgetHalt):
        env.run(clock=StepClock(600), resume_authorization="R1")
    st = env.state()
    assert st["state"] == "HALTED" and st["pending"] == 2
    with pytest.raises(orch.InvocationStateError, match="already used"):
        env.run(resume_authorization="R1")
    env.run(resume_authorization="R2")
    assert {m: a["reference"] for m, a in env.records("auth").items()} == {1: "R1", 2: "R2"}
    assert env.state()["state"] == "CLEAN"


# ================================================================ 11, 17, 18: UNTERMINATED
@pytest.mark.parametrize("exc", [RuntimeError, KeyboardInterrupt])
def test_crash_leaves_unterminated_and_needs_acknowledgement(env, exc):
    ref = uninterrupted_bytes(env)
    BEHAVIOUR.update(raise_at=120, exc=exc)
    with pytest.raises(exc):
        env.run()
    st = env.state()
    assert st["state"] == "UNTERMINATED" and st["pending"] == 1 and not env.records("end")
    BEHAVIOUR.update(raise_at=None)
    with pytest.raises(orch.InvocationStateError, match="cause unknown"):
        env.run()
    env.run(resume_authorization="ACK-1")
    assert env.records("auth")[1]["authorization_kind"] == "unterminated_acknowledgement"
    assert env.ckpt.read_bytes() == ref and env.state()["state"] == "CLEAN"


def test_crash_during_authorised_resume(env):
    ref = uninterrupted_bytes(env)
    halt_first_invocation(env)
    BEHAVIOUR.update(raise_at=120)
    with pytest.raises(RuntimeError):
        env.run(resume_authorization="R1")
    assert env.state()["state"] == "UNTERMINATED" and env.state()["pending"] == 2
    BEHAVIOUR.update(raise_at=None)
    with pytest.raises(orch.InvocationStateError):
        env.run(resume_authorization="R1")              # old reference cannot be reused
    env.run(resume_authorization="R2")
    assert env.ckpt.read_bytes() == ref


def test_unterminated_prefix_change_refused(env):
    halt_first_invocation(env)                          # inv1: 55 records
    BEHAVIOUR.update(raise_at=120)
    with pytest.raises(RuntimeError):
        env.run(resume_authorization="R1")              # inv2 starts with 55 records
    BEHAVIOUR.update(raise_at=None)
    lines = env.ckpt.read_bytes().split(b"\n")
    env.ckpt.write_bytes(b"\n".join(lines[:41]) + b"\n")  # header + 40 records
    with pytest.raises(orch.InvocationStateError, match="has changed"):
        env.run(resume_authorization="R2")


def test_crash_after_final_record_before_completion_record(env, monkeypatch):
    real = orch._write_end

    def crash_on_completed(out, run_id, m, outcome, *a, **k):
        if outcome == "completed":
            raise RuntimeError("simulated crash before the completion record")
        return real(out, run_id, m, outcome, *a, **k)
    monkeypatch.setattr(orch, "_write_end", crash_on_completed)
    with pytest.raises(RuntimeError):
        env.run()
    assert len(CheckpointStore.load_records(env.ckpt)) == 500
    assert env.state()["state"] == "UNTERMINATED"
    monkeypatch.setattr(orch, "_write_end", real)
    with pytest.raises(orch.InvocationStateError):
        env.run()
    res = env.run(resume_authorization="ACK")
    assert res["n_records"] == 500 and env.state()["state"] == "CLEAN"
    assert env.records("start")[2]["n_new"] == 0


def test_short_invocation_logged_and_not_evaluated(env):
    BEHAVIOUR.update(raise_at=470)
    with pytest.raises(RuntimeError):
        env.run()
    BEHAVIOUR.update(raise_at=None)
    clk = StepClock(10_000.0)                           # would halt instantly if evaluated
    env.run(clock=clk, resume_authorization="ACK")
    s = env.records("start")[2]
    assert s["n_new"] == 30 and s["runtime_check"] == "not_evaluated (n_new < W)"
    assert clk.n == 3                                   # start, end, result: no evaluation
    ev = [json.loads(l) for l in (env.out / "logs/events.jsonl").read_text().splitlines()]
    assert [e for e in ev if e["event"] == "run_start"][-1]["runtime_check"] == \
        "not_evaluated (n_new < W)"


# ================================================================ 20: non-runtime outcomes
def test_nonestimable_stop_recorded_distinctly_existing_handling(env):
    BEHAVIOUR["bad"] = {0: "unreliable_apparent_failed"}
    with pytest.raises(orch.NonEstimableStop) as ex:
        env.run()
    assert isinstance(ex.value, orch.StopCondition)
    assert env.records("end")[1]["outcome"] == "stop_nonestimable"
    st = env.state()
    assert st["state"] == "CLEAN" and st["last_outcome"] == "stop_nonestimable"
    with pytest.raises(orch.DispatchRefused, match="non-estimable"):
        env.run(resume_authorization="R")              # Decision 11B: refused before dispatch
    with pytest.raises(orch.DispatchRefused, match="non-estimable"):
        env.run()                                       # Decision 11B: permanent under Option A


def test_h4_failure_recorded_distinctly(env):
    env.h4_fail_stage = "during run 0 k=10"
    with pytest.raises(HardFailureH4):
        env.run()
    assert env.records("end")[1]["outcome"] == "hard_failure_h4"
    st = env.state()
    assert st["state"] == "CLEAN" and st["last_outcome"] == "hard_failure_h4"
    with pytest.raises(orch.InvocationStateError):
        env.run(resume_authorization="R")


def test_post_run_preservation_failure_recorded_distinctly(env):
    env.preservation_fail_on_call = 2                   # call 1 = pre (ok), call 2 = post
    with pytest.raises(orch.HistoricalPreservationFailure):
        env.run()
    assert env.records("end")[1]["outcome"] == "historical_preservation_failure"
    assert env.state()["state"] == "CLEAN"


# ================================================================ 10: durability / atomic writes
@pytest.mark.skipif(not DARWIN, reason="F_FULLFSYNC is macOS-specific")
def test_full_fsync_used(env, monkeypatch):
    calls = []
    real = orch.fcntl.fcntl
    monkeypatch.setattr(orch.fcntl, "fcntl",
                        lambda fd, op, *a: calls.append(op) or real(fd, op, *a))
    env.run()
    assert calls and all(op == orch.fcntl.F_FULLFSYNC for op in calls)


@pytest.mark.skipif(not DARWIN, reason="F_FULLFSYNC is macOS-specific")
def test_full_fsync_unavailable_fails_closed_before_dispatch(env, monkeypatch):
    monkeypatch.delattr(orch.fcntl, "F_FULLFSYNC")
    with pytest.raises(orch.DurabilityError, match="F_FULLFSYNC unavailable"):
        env.run()
    assert BEHAVIOUR["calls"] == 0                      # nothing dispatched
    assert CheckpointStore.load_records(env.ckpt) == {}


def test_full_sync_failure_fails_closed(env, monkeypatch):
    def boom(fd):
        raise orch.DurabilityError("simulated sync failure")
    monkeypatch.setattr(orch, "_full_sync", boom)
    with pytest.raises(orch.DurabilityError):
        env.run()
    assert BEHAVIOUR["calls"] == 0


def failing_link_for(suffix):
    def link(src, dst, *a, **k):
        if str(dst).endswith(suffix):
            raise OSError("simulated crash before commit")
        return REAL_LINK(src, dst, *a, **k)
    return link


def test_crash_before_start_commit_is_ambiguous(env, monkeypatch):
    monkeypatch.setattr(orch.os, "link", failing_link_for("start.json"))
    with pytest.raises(orch.DurabilityError):
        env.run()
    assert BEHAVIOUR["calls"] == 0
    monkeypatch.setattr(orch.os, "link", REAL_LINK)
    st = env.state()
    assert st["state"] == "AMBIGUOUS" and "incomplete record" in st["problems"][0]


def test_crash_while_writing_halt_record_is_ambiguous(env, monkeypatch):
    monkeypatch.setattr(orch.os, "link", failing_link_for("end.json"))
    with pytest.raises(orch.DurabilityError):
        env.run(clock=StepClock(600))
    assert env.state()["state"] == "AMBIGUOUS"
    for ref in (None, "R"):
        with pytest.raises(orch.InvocationStateError, match="AMBIGUOUS"):
            orch.run_corrected(env.cfg(), 0, 500, clock=("T", StepClock(0.001)),
                               resume_authorization=ref)


def test_crash_after_link_before_tmp_removal_is_ambiguous(env, monkeypatch):
    def unlink(p, *a, **k):
        if str(p).endswith("end.json.tmp"):
            raise OSError("simulated crash after commit")
        return REAL_UNLINK(p, *a, **k)
    monkeypatch.setattr(orch.os, "unlink", unlink)
    with pytest.raises(orch.DurabilityError):
        env.run(clock=StepClock(600))
    monkeypatch.setattr(orch.os, "unlink", REAL_UNLINK)
    assert env.records("end")[1]["outcome"] == "runtime_halt"
    assert env.state()["state"] == "AMBIGUOUS"


def test_crash_during_authorisation_creation_is_ambiguous(env, monkeypatch):
    halt_first_invocation(env)
    monkeypatch.setattr(orch.os, "link", failing_link_for("auth.json"))
    with pytest.raises(orch.DurabilityError):
        env.run(resume_authorization="R1")
    monkeypatch.setattr(orch.os, "link", REAL_LINK)
    assert env.state()["state"] == "AMBIGUOUS" and not env.records("auth")


def test_crash_after_authorisation_before_resume_event(env, monkeypatch):
    ref = uninterrupted_bytes(env)
    halt_first_invocation(env)
    real = orch.log_event

    def log(out, ev):
        if ev.get("event") == "resume_authorized":
            raise RuntimeError("simulated crash after authorisation commit")
        return real(out, ev)
    monkeypatch.setattr(orch, "log_event", log)
    with pytest.raises(RuntimeError):
        env.run(resume_authorization="R1")
    monkeypatch.setattr(orch, "log_event", real)
    assert env.state()["state"] == "CLEAN"              # the auth record is authoritative
    env.run()
    assert env.ckpt.read_bytes() == ref


def test_create_only_records_cannot_be_overwritten(tmp_path):
    p = tmp_path / "halts" / "run0.inv0001.start.json"
    orch.durable_create(p, {"a": 1})
    with pytest.raises(orch.DurabilityError, match="create-only"):
        orch.durable_create(p, {"a": 2})
    assert json.loads(p.read_text()) == {"a": 1}


@pytest.mark.parametrize("files,problem", [
    ({"run0.inv0001.auth.json": {"kind": "auth", "run_id": 0, "invocation": 1}}, "without a start"),
    ({"run0.inv0001.start.json": {"kind": "start", "run_id": 0, "invocation": 1},
      "run0.inv0001.end.json": {"kind": "end", "run_id": 0, "invocation": 1, "outcome": "runtime_halt"},
      "run0.inv0002.start.json": {"kind": "start", "run_id": 0, "invocation": 2}}, "more than one"),
    ({"run0.inv0001.start.json": {"kind": "start", "run_id": 0, "invocation": 1},
      "run0.inv0001.end.json": {"kind": "end", "run_id": 0, "invocation": 1, "outcome": "completed"},
      "run0.inv0001.auth.json": {"kind": "auth", "run_id": 0, "invocation": 1,
                                 "authorization_kind": "runtime_halt_authorization"}}, "does not match"),
    ({"run0.inv0001.start.json": {"kind": "end", "run_id": 0, "invocation": 1}}, "identity"),
])
def test_inconsistent_records_are_ambiguous(tmp_path, files, problem):
    d = tmp_path / "halts"
    d.mkdir()
    for n, o in files.items():
        (d / n).write_text(json.dumps(o))
    st = orch.invocation_state(tmp_path, 0)
    assert st["state"] == "AMBIGUOUS" and any(problem in p for p in st["problems"])


def test_other_runs_records_do_not_affect_state(tmp_path):
    d = tmp_path / "halts"
    d.mkdir()
    (d / "run1.inv0001.start.json").write_text(json.dumps({"kind": "start", "run_id": 1,
                                                           "invocation": 1}))
    assert orch.invocation_state(tmp_path, 0)["state"] == "CLEAN"
    assert orch.invocation_state(tmp_path, 1)["state"] == "UNTERMINATED"


# ================================================================ 21: determinism / workers
def test_worker_count_invariance_uninterrupted(env):
    a = uninterrupted_bytes(env, n_workers=1)
    b = uninterrupted_bytes(env, n_workers=3)
    assert a == b


def test_topup_still_requires_its_own_authorisation(env):
    env.run()
    with pytest.raises(ConfigurationError, match="top-up"):
        orch.run_corrected(env.cfg(), 0, 750, clock=("T", StepClock(0.001)))


# ================================================================ 23: driver
@pytest.fixture(scope="module")
def driver():
    spec = importlib.util.spec_from_file_location(
        "p9c_driver_under_test", ROOT / "scripts/41_phase9_corrected_driver.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture
def cfg_file(tmp_path):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(approved_test_config(ROOT, "scratch_p9c_drv_unused")))
    return str(p)


@pytest.mark.parametrize("tail", [["--resume-authorization"],
                                  ["--resume-authorization", "   "],
                                  ["--resume-authorization", "--topup-authorization", "X"]])
def test_driver_rejects_missing_or_blank_reference(driver, cfg_file, tail, monkeypatch):
    called = []
    monkeypatch.setattr(orch, "run_corrected", lambda *a, **k: called.append(1))
    assert driver.main(["run", cfg_file, "0", "500", *tail]) == 2
    assert not called


def test_driver_passes_reference_through(driver, cfg_file, monkeypatch):
    seen = {}

    def fake(cfg, run_id, target, **kw):
        seen.update(run_id=run_id, target=target, **kw)
        return {"ok": True}
    monkeypatch.setattr(orch, "run_corrected", fake)
    assert driver.main(["run", cfg_file, "0", "500", "--resume-authorization", "OWNER-7"]) == 0
    assert seen == {"run_id": 0, "target": 500, "topup_authorization": None,
                    "resume_authorization": "OWNER-7"}
    assert driver.main(["run", cfg_file, "0", "750", "--topup-authorization", "T-1"]) == 0
    assert seen["topup_authorization"] == "T-1" and seen["resume_authorization"] is None


def test_driver_reports_refusal(driver, cfg_file, monkeypatch):
    def refuse(*a, **k):
        raise orch.InvocationStateError("no runtime halt or unterminated invocation is pending")
    monkeypatch.setattr(orch, "run_corrected", refuse)
    assert driver.main(["run", cfg_file, "0", "500", "--resume-authorization", "R"]) == 3
