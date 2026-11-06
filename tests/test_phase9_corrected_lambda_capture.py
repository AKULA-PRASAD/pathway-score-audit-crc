"""Owner Decision 12 (Option C, 2026-11-06): capture of the ORIGINAL tune_lambda result via
a per-call private namespace. Uses the REAL frozen pipeline (src/phase9_execution.py) on
the small synthetic cohort of tests/p9c_support.py. Failures are injected deterministically
by test-only stand-ins; no scientific data, no real bootstrap, no gate."""
import hashlib
import json
import shutil
import uuid
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src import phase9_execution as px
from src.phase9_corrected.checkpoint import CheckpointStore
from src.phase9_corrected.config import PLAN_BASE, canonical_json
from src.phase9_corrected.data import LAM_GRID
from src.phase9_corrected.engine import compute_outer_replicate
from src.phase9_corrected.fitting import (LAMBDA_CAPTURED, LAMBDA_FIXED_ZERO, LAMBDA_NOT_APPLICABLE,
                                          LambdaCaptureDefect, RealPipeline, TuneCapture,
                                          fit_and_score)
from src.phase9_corrected.orchestration import execute_replicates
from src.phase9_corrected.rng import PlanStreams
from src.score_standardisation import NonEstimableReplicate
from tests.p9c_support import synthetic_real_pipeline_context

ROOT = Path(__file__).resolve().parents[1]
S = PlanStreams(PLAN_BASE, 0, "independent", 7)
ORIG_TUNE, ORIG_FIT, ORIG_FITCOX = px.tune_lambda, px.pipeline_fit, px.fit_cox
ORIG_CODE_BYTES = ORIG_FIT.__code__.co_code


class Uncaptured(RealPipeline):
    capturing_fit = None                         # the pre-Decision-12 behaviour


@pytest.fixture(scope="module")
def ctx():
    return synthetic_real_pipeline_context()


def gen():
    return S.generator(S.apparent_cv(0))


def m2_fit(pipeline, ctx, rng=None):
    return fit_and_score(pipeline, ctx, "M2", "apparent", ctx.ids, None, ctx.t, ctx.e,
                         None, None, rng or gen(), LAM_GRID)


def canon(o):
    return canonical_json(o)


# ---------------------------------------------------------------- failure injection (tests only)
def final_fit_fails(pred=lambda design: True):
    """fit_cox stand-in that raises ConvergenceFailure ONLY on pipeline_fit's final call
    (the only call passing `strata=`; tune_lambda's internal fits never do)."""
    def fc(design, t, e, penalizer, *a, **kw):
        if "strata" in kw and pred(design):
            raise px.ConvergenceFailure("injected final-fit failure")
        return ORIG_FITCOX(design, t, e, penalizer, *a, **kw)
    return fc


def is_m2(design):
    return any(str(c).startswith("P") for c in design.columns)


@pytest.fixture
def restore_px():
    yield
    px.fit_cox, px.tune_lambda = ORIG_FITCOX, ORIG_TUNE
    px.standardise_scores_training_fold = _ORIG_STD
    assert px.pipeline_fit is ORIG_FIT and px.pipeline_fit.__code__.co_code == ORIG_CODE_BYTES


_ORIG_STD = px.standardise_scores_training_fold


# ================================================================ A. capture correctness (1-7)
def test_successful_m2_capture_equals_diag_and_score(ctx):
    pl, cap = RealPipeline(), TuneCapture()
    f = pl.capturing_fit(cap)
    rng = gen()
    m, d_tr, d_ap, dg = f(ctx.X, ctx.cov, ctx.panel, ctx.ids, ctx.ids, ctx.t, ctx.e, LAM_GRID, rng)
    assert cap.calls == 1 and cap.raised is None
    assert cap.result == (dg["lambda"], dg["inner_cv_c"])                  # tests 1 and 2
    assert isinstance(cap.result[0], float) and cap.result[0] in [float(x) for x in LAM_GRID]
    assert cap.args[5] is rng and cap.args[4] is LAM_GRID and cap.kwargs == {}   # test 4
    assert cap.args[3] == dg["n_clinical_df"]
    pl.check_success_capture(cap, dg)


def test_exactly_one_original_tuning_call_per_m2_fit(ctx, restore_px):
    calls = []

    def counting(*a, **k):
        calls.append(1)
        return ORIG_TUNE(*a, **k)
    px.tune_lambda = counting
    m2_fit(RealPipeline(), ctx)
    assert len(calls) == 1
    calls.clear()
    m2_fit(Uncaptured(), ctx)
    assert len(calls) == 1                                                  # test 3 / 18


def test_module_namespace_and_code_object_never_modified(ctx, restore_px):
    seen = []

    def watching_fit_cox(*a, **k):                     # runs DURING the captured fit
        seen.append(px.tune_lambda is ORIG_TUNE and px.pipeline_fit is ORIG_FIT)
        return ORIG_FITCOX(*a, **k)
    before = dict(px.__dict__)
    px.fit_cox = watching_fit_cox
    out = m2_fit(RealPipeline(), ctx)
    px.fit_cox = ORIG_FITCOX
    assert out.ok and seen and all(seen)                                    # test 5 (during)
    assert px.tune_lambda is ORIG_TUNE and px.pipeline_fit is ORIG_FIT      # test 5 (after)
    assert {k: v for k, v in px.__dict__.items()} == before                 # namespace intact
    f = RealPipeline().capturing_fit(TuneCapture())
    assert f.__code__ is ORIG_FIT.__code__ and ORIG_FIT.__code__.co_code == ORIG_CODE_BYTES
    assert f.__defaults__ == ORIG_FIT.__defaults__ and f.__closure__ is ORIG_FIT.__closure__
    assert f.__kwdefaults__ == ORIG_FIT.__kwdefaults__                     # test 6
    assert f.__globals__ is not px.__dict__ and f.__globals__["tune_lambda"] is not ORIG_TUNE


def test_missing_repeated_or_inconsistent_capture_fails_closed(ctx, monkeypatch):
    pl = RealPipeline()
    tuned = {"lambda": 0.1, "inner_cv_c": 0.6}
    for calls, result in ((0, None), (2, (0.1, 0.6)), (1, (0.2, 0.6)), (1, (0.1, 0.5)), (1, "x")):
        c = TuneCapture()
        c.calls, c.result = calls, result
        with pytest.raises(LambdaCaptureDefect):
            pl.check_success_capture(c, tuned)
    c = TuneCapture()
    c.calls = 1
    with pytest.raises(LambdaCaptureDefect):                               # tuning without diag
        pl.check_success_capture(c, {"lambda": 0.0})
    c.result = (0.1, 0.6)
    with pytest.raises(LambdaCaptureDefect):                               # non-estimable AFTER tuning
        pl.failure_diagnostics(c, NonEstimableReplicate("x"))
    with pytest.raises(LambdaCaptureDefect):
        pl.failure_diagnostics(TuneCapture(), RuntimeError("unclassifiable"))
    # integration: a capturing_fit that bypasses the capture while tuning still occurs
    monkeypatch.setattr(RealPipeline, "capturing_fit", lambda self, cap: px.pipeline_fit)
    with pytest.raises(LambdaCaptureDefect):
        m2_fit(RealPipeline(), ctx)                                         # test 7


# ================================================================ B. failure scenarios (8-14)
def test_m2_final_fit_failure_records_original_lambda_and_score(ctx, restore_px):
    ok = m2_fit(RealPipeline(), ctx)                     # the same fit without the failure
    px.fit_cox = final_fit_fails()
    bad = m2_fit(RealPipeline(), ctx)
    assert ok.ok and not bad.ok and bad.extra["failure_type"] == "ConvergenceFailure"
    assert bad.extra["lambda_status"] == LAMBDA_CAPTURED
    assert bad.lambda_at_failure.hex() == ok.lam.hex()                      # test 8
    assert bad.extra["tuning_score_at_failure"].hex() == ok.inner_cv_c.hex()


def test_m2_failure_before_tuning_not_applicable(ctx, restore_px):
    def std_fails(scores, train_ids):
        raise NonEstimableReplicate("injected non-finite scaling")
    px.standardise_scores_training_fold = std_fails
    out = m2_fit(RealPipeline(), ctx)
    assert not out.ok and out.extra["lambda_status"] == LAMBDA_NOT_APPLICABLE
    assert out.lambda_at_failure is None and out.extra["tuning_score_at_failure"] is None   # 9
    assert out.cv_rng_consumed is False


def test_m1_failure_fixed_zero(ctx, restore_px):
    px.fit_cox = final_fit_fails()
    out = fit_and_score(RealPipeline(), ctx, "M1", "apparent", ctx.ids, None, ctx.t, ctx.e,
                        None, None, gen(), LAM_GRID)
    assert not out.ok and out.extra["lambda_status"] == LAMBDA_FIXED_ZERO
    assert out.lambda_at_failure == 0.0 and out.extra["tuning_score_at_failure"] is None   # 10


def test_m2_all_pathways_dropped_fixed_zero(ctx, restore_px):
    def all_dropped(scores, train_ids):
        return (pd.DataFrame(index=scores.index), {"n_dropped": scores.shape[1],
                                                   "dropped_pathways": list(scores.columns)})
    px.standardise_scores_training_fold = all_dropped
    ok = m2_fit(RealPipeline(), ctx)
    assert ok.ok and ok.lam == 0.0 and ok.cv_rng_consumed is False           # no tuning, no defect
    px.fit_cox = final_fit_fails()
    out = m2_fit(RealPipeline(), ctx)
    assert out.extra["lambda_status"] == LAMBDA_FIXED_ZERO and out.lambda_at_failure == 0.0   # 11


def run_rec(ctx, pipeline, k=0, B_in=3):
    return json.loads(canon(compute_outer_replicate(ctx, S, pipeline, 0, k, B_in, LAM_GRID)))


def test_engine_pairing_and_independent_failure_diagnostics(ctx, restore_px):
    base = run_rec(ctx, RealPipeline())
    assert base["status"] == "ok" and not base["failures"]
    order = []

    def fc(design, t, e, penalizer, *a, **kw):
        if "strata" in kw:
            order.append(is_m2(design))
            i = len(order) - 1                 # 0 app M1, 1 app M2, then (M1, M2) per b
            if i in (3, 4, 5):                 # b0: M2 fails; b1: both fail
                raise px.ConvergenceFailure("injected")
        return ORIG_FITCOX(design, t, e, penalizer, *a, **kw)
    px.fit_cox = fc
    rec = run_rec(ctx, RealPipeline())
    rows = rec["inner"]
    assert [r["paired_complete"] for r in rows] == [False, False, True]       # test 12 pairing
    assert rows[0]["fits"]["M1"]["ok"] and rows[0]["fits"]["M1"]["lam"] == 0.0
    assert rec["B_in_paired_complete"] == 1
    f = {(x["b"], x["model"]): x for x in rec["failures"]}
    assert set(f) == {(0, "M2"), (1, "M1"), (1, "M2")}                       # test 13
    for b in (0, 1):
        assert f[(b, "M2")]["lambda_status"] == LAMBDA_CAPTURED
        assert f[(b, "M2")]["lambda_at_failure"] == base["inner"][b]["fits"]["M2"]["lam"]
        assert f[(b, "M2")]["tuning_score_at_failure"] == base["inner"][b]["fits"]["M2"]["inner_cv_c"]
        assert rows[b]["fits"]["M2"]["lambda_at_failure"] == f[(b, "M2")]["lambda_at_failure"]
    assert f[(1, "M1")]["lambda_status"] == LAMBDA_FIXED_ZERO and f[(1, "M1")]["lambda_at_failure"] == 0.0
    # successful fits are identical to the failure-free run
    assert rows[2] == base["inner"][2] and rec["apparent"] == base["apparent"]


def test_exceptions_inside_capture_propagate_unchanged(ctx, restore_px):
    boom = RuntimeError("raised inside tune_lambda")

    def raising(*a, **k):
        raise boom
    px.tune_lambda = raising
    with pytest.raises(RuntimeError) as ei:
        m2_fit(RealPipeline(), ctx)
    assert ei.value is boom                                                # test 14: same object
    cf = px.ConvergenceFailure("a failure type raised from inside tuning (impossible today)")

    def raising_cf(*a, **k):
        raise cf
    px.tune_lambda = raising_cf
    with pytest.raises(LambdaCaptureDefect) as ei2:                        # fails closed, chained
        m2_fit(RealPipeline(), ctx)
    assert ei2.value.__cause__ is cf


# ================================================================ C. RNG and determinism (15-18)
def test_rng_state_consumption_and_outputs_identical_with_and_without_capture(ctx):
    r1, r2 = gen(), gen()
    a, b = m2_fit(RealPipeline(), ctx, r1), m2_fit(Uncaptured(), ctx, r2)
    assert r1.bit_generator.state == r2.bit_generator.state                 # test 15
    assert a.cv_rng_consumed is b.cv_rng_consumed is True                  # test 16
    assert canon(a.to_dict()) == canon(b.to_dict())                         # test 17
    m1a = fit_and_score(RealPipeline(), ctx, "M1", "apparent", ctx.ids, None, ctx.t, ctx.e,
                        None, None, gen(), LAM_GRID)
    assert m1a.cv_rng_consumed is False


def test_full_record_byte_identical_with_and_without_capture(ctx):
    assert canon(run_rec(ctx, RealPipeline())) == canon(run_rec(ctx, Uncaptured()))   # test 17


# ================================================================ D. workers and resume (19-22)
def det_pred(design):
    """Deterministic, content-based injection: a subset of M2 final fits fail."""
    h = hashlib.sha256("".join(map(str, design.index)).encode()).digest()[0]
    return is_m2(design) and h % 3 == 0


@pytest.fixture
def scratch():
    d = ROOT / f"scratch_p9c_d12_{uuid.uuid4().hex[:10]}"
    yield d
    shutil.rmtree(d, ignore_errors=True)
    assert not d.exists()


def run_store(ctx, out, pipeline, workers, target=4, interrupt_after=None, B_in=2):
    class Stop(Exception):
        pass
    def on_record(rec, st):
        if interrupt_after is not None and len(st.done) == interrupt_after:
            raise Stop()
    with CheckpointStore(out, 0, "c" * 64, "m" * 64, ROOT) as st:
        try:
            execute_replicates(st, ctx, S, pipeline, B_in, target, workers, False, on_record)
        except Stop:
            pass
    return (out / "runs" / "run0.jsonl").read_bytes()


def test_workers_resume_and_capture_are_byte_identical(ctx, scratch, restore_px):
    px.fit_cox = final_fit_fails(det_pred)
    results = {}
    for name, pl, w, intr in (("cap1", RealPipeline(), 1, None), ("cap3", RealPipeline(), 3, None),
                              ("unc1", Uncaptured(), 1, None)):
        results[name] = run_store(ctx, scratch / name, pl, w, interrupt_after=intr)
    # interrupted after 2 records at 3 workers, then resumed at 1 worker
    run_store(ctx, scratch / "res", RealPipeline(), 3, interrupt_after=2)
    results["resumed"] = run_store(ctx, scratch / "res", RealPipeline(), 1)
    assert results["cap1"] == results["cap3"] == results["resumed"]          # tests 19-22
    recs = CheckpointStore.load_records(scratch / "cap1" / "runs" / "run0.jsonl")
    fails = [f for r in recs.values() for f in r["failures"]]
    assert fails, "the deterministic injection must produce at least one failure"
    assert all(f["lambda_status"] == LAMBDA_CAPTURED and f["lambda_at_failure"] is not None
               for f in fails if f["model"] == "M2")
    # uncaptured differs ONLY in the diagnostic fields of failed fits
    unc = CheckpointStore.load_records(scratch / "unc1" / "runs" / "run0.jsonl")
    for k in recs:
        a, b = json.loads(json.dumps(recs[k])), json.loads(json.dumps(unc[k]))
        for rr in (a, b):
            for f in rr["failures"]:
                for key in ("lambda_at_failure", "lambda_status", "tuning_score_at_failure"):
                    f.pop(key, None)
            for row in rr.get("inner", []):
                for fd in row["fits"].values():
                    if not fd["ok"]:
                        for key in ("lambda_at_failure", "lambda_status", "tuning_score_at_failure"):
                            fd.pop(key, None)
            for fd in rr.get("apparent", {}).values():
                if not fd["ok"]:
                    for key in ("lambda_at_failure", "lambda_status", "tuning_score_at_failure"):
                        fd.pop(key, None)
        assert a == b


def test_successful_checkpoint_bytes_identical_with_and_without_capture(ctx, scratch):
    assert run_store(ctx, scratch / "a", RealPipeline(), 1, target=3) == \
        run_store(ctx, scratch / "b", Uncaptured(), 3, target=3)


# ================================================================ E. V5-style recomputation (23)
def test_v5_style_recomputation_identical(ctx, scratch, restore_px):
    px.fit_cox = final_fit_fails(det_pred)
    run_store(ctx, scratch / "v5", RealPipeline(), 1)
    stored = CheckpointStore.load_records(scratch / "v5" / "runs" / "run0.jsonl")
    for k in (0, 3):
        again = json.loads(canon(compute_outer_replicate(ctx, S, RealPipeline(), 0, k, 2, LAM_GRID)))
        assert again == stored[k]


# ================================================================ D12-T1: degenerate tuning
import math                                                                # noqa: E402

from src.phase9_corrected import orchestration as orch                     # noqa: E402
from src.phase9_corrected.config import CorrectedConfig                    # noqa: E402
from src.phase9_corrected.data import AnalysisContext                      # noqa: E402
from tests.p9c_support import approved_test_config, synthetic_gate         # noqa: E402


def degenerate_ctx(ctx):
    """Only 4 events: every tuning fold fails the 'training events >= 5' guard
    (src/phase9_execution.py:152), so tune_lambda returns (grid[0], -inf)."""
    e = np.zeros(len(ctx.e))
    e[[3, 17, 41, 66]] = 1.0
    return AnalysisContext(ctx.ids, ctx.cov, ctx.t, e, ctx.X, ctx.panel)


def test_d12_t1_degenerate_tuning_then_final_fit_failure(ctx, scratch, restore_px):
    dctx = degenerate_ctx(ctx)
    calls = []

    def counting(*a, **k):
        calls.append(1)
        return ORIG_TUNE(*a, **k)
    px.tune_lambda = counting
    px.fit_cox = final_fit_fails()
    out = m2_fit(RealPipeline(), dctx)
    assert len(calls) == 1                                         # no extra tuning call
    assert out.extra["lambda_status"] == LAMBDA_CAPTURED
    assert out.lambda_at_failure.hex() == float(LAM_GRID[0]).hex()  # original selection
    sc = out.extra["tuning_score_at_failure"]
    assert math.isinf(sc) and sc < 0
    assert b'"tuning_score_at_failure":-Infinity' in canon(out.to_dict())
    # whole replicate: apparent failure record carries the -inf score; checkpoint round trip
    rec = compute_outer_replicate(dctx, S, RealPipeline(), 0, 0, 2, LAM_GRID)
    assert rec["status"] == "unreliable_apparent_failed"
    f = [x for x in rec["failures"] if x["model"] == "M2"][0]
    assert f["lambda_status"] == LAMBDA_CAPTURED and f["tuning_score_at_failure"] == -math.inf
    b1 = canon(rec)
    assert b1 == canon(compute_outer_replicate(dctx, S, RealPipeline(), 0, 0, 2, LAM_GRID))
    with CheckpointStore(scratch, 0, "c" * 64, "m" * 64, ROOT) as st:
        st.append(rec)
    loaded = CheckpointStore.load_records(scratch / "runs" / "run0.jsonl")[0]
    assert loaded == json.loads(b1) and canon(loaded) == b1          # exact, hash-stable
    assert hashlib.sha256(canon(loaded)).hexdigest() == hashlib.sha256(b1).hexdigest()


# ================================================================ D12-T2: LambdaCaptureDefect -> UNTERMINATED
@pytest.fixture
def run_env(monkeypatch, scratch):
    out_rel = scratch.relative_to(ROOT).as_posix()
    hits = {"k": None}

    def fake_record(ctx_, streams, pipeline, run_id, k, B_in, lam_grid, as_run_variant=False):
        if k == 3:
            hits["k"] = k
            raise LambdaCaptureDefect("synthetic capture defect at k=3")
        return {"run_id": int(run_id), "k": int(k), "status": "ok", "failures": []}
    monkeypatch.setattr(orch, "h4_check", lambda root, stage, raise_on_fail=True: {"all_ok": True})
    monkeypatch.setattr(orch, "historical_preservation_check", lambda root: {"all_ok": True})
    monkeypatch.setattr(orch, "build_manifest", lambda c, r: {"schema": "t", "code": {"x": "0" * 64}})
    monkeypatch.setattr(orch, "run_pre_execution_gate", synthetic_gate)
    monkeypatch.setattr(orch, "load_context", lambda root: None)
    monkeypatch.setattr(orch, "RealPipeline", lambda: None)
    monkeypatch.setattr(orch, "compute_outer_replicate", fake_record)
    monkeypatch.setattr(orch, "verify_code_identity", lambda config, manifest=None, extra_files=():
                        {"freeze_record_sha256": "f" * 64, "freeze": {"code": {}}})
    cfg = CorrectedConfig(approved_test_config(ROOT, out_rel, n_workers=1,
                                               runtime_budget_hours_per_run=4.0), ROOT)

    class Clock:
        n = 0

        def __call__(self):
            Clock.n += 1
            return Clock.n * 0.001

    def run(**kw):
        return orch.run_corrected(cfg, 0, 500, clock=("TEST_CLOCK", Clock()), **kw)
    return run, scratch, hits


def test_d12_t2_capture_defect_leaves_invocation_unterminated(run_env):
    run, out, hits = run_env
    with pytest.raises(LambdaCaptureDefect):
        run()
    assert hits["k"] == 3
    recs = CheckpointStore.load_records(out / "runs" / "run0.jsonl")
    assert sorted(recs) == [0, 1, 2] and all(not r["failures"] for r in recs.values())  # not a model failure
    assert not (out / "halts" / "run0.inv0001.end.json").exists()
    st = orch.invocation_state(out, 0)
    assert st["state"] == "UNTERMINATED" and st["pending"] == 1
    with pytest.raises(orch.InvocationStateError, match="cause unknown"):
        run()                                                            # ordinary resume refused
    with pytest.raises(LambdaCaptureDefect):                             # acknowledged: recurs
        run(resume_authorization="OWNER-ACK-1")
    auth = json.loads((out / "halts" / "run0.inv0001.auth.json").read_text())
    assert auth["authorization_kind"] == "unterminated_acknowledgement"
    assert sorted(CheckpointStore.load_records(out / "runs" / "run0.jsonl")) == [0, 1, 2]
    assert orch.invocation_state(out, 0)["state"] == "UNTERMINATED"


def test_d12_t2_real_engine_propagates_defect_not_failure(ctx, monkeypatch):
    monkeypatch.setattr(RealPipeline, "capturing_fit", lambda self, cap: px.pipeline_fit)
    with pytest.raises(LambdaCaptureDefect):
        compute_outer_replicate(ctx, S, RealPipeline(), 0, 0, 2, LAM_GRID)


# ================================================================ D12-T3: as-run apparent variant
def test_d12_t3_as_run_variant_failure_diagnostics(ctx, restore_px):
    base = json.loads(canon(compute_outer_replicate(ctx, S, RealPipeline(), 0, 0, 2, LAM_GRID,
                                                    as_run_variant=True)))
    assert canon(base) == canon(compute_outer_replicate(ctx, S, Uncaptured(), 0, 0, 2, LAM_GRID,
                                                        as_run_variant=True))   # success unchanged
    assert base["streams"]["apparent_cv"]["state_sha256"] != \
        base["streams"]["apparent_cv_as_run"]["state_sha256"]                 # distinct streams
    order, tunes = [], []

    def fc(design, t, e, penalizer, *a, **kw):
        if "strata" in kw:
            order.append(is_m2(design))
            if len(order) in (2, 4):          # 1 app M1, 2 app M2, 3 as-run M1, 4 as-run M2
                raise px.ConvergenceFailure("injected apparent final-fit failure")
        return ORIG_FITCOX(design, t, e, penalizer, *a, **kw)

    def counting(*a, **k):
        tunes.append(1)
        return ORIG_TUNE(*a, **k)
    px.fit_cox, px.tune_lambda = fc, counting
    rec = json.loads(canon(compute_outer_replicate(ctx, S, RealPipeline(), 0, 0, 2, LAM_GRID,
                                                   as_run_variant=True)))
    assert rec["status"] == "unreliable_apparent_failed" and len(tunes) == 2   # one per stage
    f = {(x["stage"], x["model"]): x for x in rec["failures"]}
    assert set(f) == {("apparent", "M2"), ("apparent_as_run", "M2")}
    for stage, src in (("apparent", "apparent"), ("apparent_as_run", "apparent_as_run")):
        assert f[(stage, "M2")]["lambda_status"] == LAMBDA_CAPTURED
        assert f[(stage, "M2")]["lambda_at_failure"] == base[src]["M2"]["lam"]
        assert f[(stage, "M2")]["tuning_score_at_failure"] == base[src]["M2"]["inner_cv_c"]
        assert rec[src]["M2"]["cv_rng_consumed"] is True and rec[src]["M1"]["cv_rng_consumed"] is False
        assert rec[src]["M1"] == base[src]["M1"]                             # M1 untouched
    assert rec["streams"] == base["streams"]
