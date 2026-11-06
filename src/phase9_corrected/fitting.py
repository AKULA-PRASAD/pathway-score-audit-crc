"""Model fitting interface.

Reuses the frozen pipeline `src.phase9_execution` UNCHANGED (pipeline_fit, c_index,
linear_predictor). The random generator passed in is consumed only by the molecular
ridge CV (`tune_lambda`, src/phase9_execution.py:196-201): M1 has no molecular block
and draws nothing. Whether the generator was consumed is recorded for every fit
(evidence for V2 stream assignment and V3 M1 invariance).

A `pipeline` object is injectable so the orchestration can be tested on synthetic
data with a fast stand-in; the execution entry point always uses `RealPipeline`.

Owner Decision 12 (Option C, 2026-11-06) -- lambda_at_failure
--------------------------------------------------------------
For every single fit, `RealPipeline.capturing_fit(cap)` builds a function from the
ORIGINAL `pipeline_fit.__code__` (same defaults, keyword defaults and closure) bound to a
PRIVATE COPY of the `src.phase9_execution` namespace in which only the name
`tune_lambda` is rebound to a pass-through. The pass-through calls the original
`tune_lambda` exactly once with the original arguments and generator, records its
return value, and returns it unchanged; exceptions propagate unchanged. The module
namespace and the frozen code object are never modified; no global, thread-local or
shared capture state exists (the capture object lives for one fit only).
On failure the record gains `lambda_status` / `tuning_score_at_failure` and the original
selected lambda in `lambda_at_failure`. Successful-fit outcomes are unchanged. Any
inconsistency between the capture and the pipeline's own diagnostics raises
`LambdaCaptureDefect` (fail closed; not a model-failure type).
"""
from __future__ import annotations

import copy
import math
import types
from dataclasses import asdict, dataclass, field
from typing import Optional


class LambdaCaptureDefect(RuntimeError):
    """Decision 12: the capture is missing, repeated or inconsistent (implementation defect)."""


LAMBDA_CAPTURED = "captured"
LAMBDA_NOT_APPLICABLE = "not_applicable_before_tuning"
LAMBDA_FIXED_ZERO = "fixed_zero_no_molecular_block"


class TuneCapture:
    """Per-fit capture object (created by fit_and_score for exactly one fit)."""
    __slots__ = ("calls", "result", "raised", "args", "kwargs")

    def __init__(self):
        self.calls, self.result, self.raised, self.args, self.kwargs = 0, None, None, None, None


def _same_float(a, b) -> bool:
    return (isinstance(a, float) and isinstance(b, float)
            and (a == b or (math.isnan(a) and math.isnan(b))))


class RealPipeline:
    name = "src.phase9_execution"

    def __init__(self):
        from src import phase9_execution as px
        from src.score_standardisation import NonEstimableReplicate
        self._px = px
        self._convergence_failure = px.ConvergenceFailure
        self._nonestimable = NonEstimableReplicate
        self.pipeline_fit = px.pipeline_fit
        self.c_index = px.c_index
        self.linear_predictor = px.linear_predictor
        self.failure_types = (px.ConvergenceFailure, NonEstimableReplicate)

    def capturing_fit(self, cap: TuneCapture):
        """The frozen pipeline_fit code object, bound to a private namespace copy whose
        `tune_lambda` is a pass-through recording into `cap` (Decision 12, Option C)."""
        px = self._px
        original_fit, original_tune = px.pipeline_fit, px.tune_lambda
        g = dict(px.__dict__)                      # private copy; px.__dict__ is untouched

        def tune_lambda(*args, **kwargs):          # called by the frozen code only
            cap.calls += 1
            cap.args, cap.kwargs = args, kwargs
            try:
                res = original_tune(*args, **kwargs)
            except BaseException as ex:            # original exception object, re-raised
                cap.raised = ex
                raise
            cap.result = res
            return res                             # returned unchanged
        g["tune_lambda"] = tune_lambda
        f = types.FunctionType(original_fit.__code__, g, original_fit.__name__,
                               original_fit.__defaults__, original_fit.__closure__)
        f.__kwdefaults__ = (dict(original_fit.__kwdefaults__)
                            if original_fit.__kwdefaults__ is not None else None)
        f.__qualname__ = original_fit.__qualname__
        if f.__code__ is not original_fit.__code__:
            raise LambdaCaptureDefect("constructed function does not use the frozen code object")
        return f

    def check_success_capture(self, cap: TuneCapture, diag: dict):
        """A successful fit: the capture must agree exactly with the pipeline's own
        diagnostics (tuning happened <=> inner_cv_c present)."""
        if cap.calls > 1 or cap.raised is not None:
            raise LambdaCaptureDefect(f"tune_lambda called {cap.calls} times / raised")
        tuned = "inner_cv_c" in diag
        if tuned != (cap.calls == 1):
            raise LambdaCaptureDefect("tuning call and pipeline diagnostics disagree")
        if tuned:
            r = cap.result
            if not (isinstance(r, tuple) and len(r) == 2 and _same_float(r[0], diag["lambda"])
                    and _same_float(r[1], diag["inner_cv_c"])):
                raise LambdaCaptureDefect("captured tuning result differs from diag")

    def failure_diagnostics(self, cap: TuneCapture, ex: BaseException) -> dict:
        """Classify a model failure by what the frozen code demonstrably did:
        ConvergenceFailure escapes pipeline_fit only from the final fit_cox (tune_lambda
        catches its own, src/phase9_execution.py:157-158); NonEstimableReplicate arises
        only in score standardisation, before tuning."""
        if cap.calls > 1 or cap.raised is not None:
            raise LambdaCaptureDefect(f"tune_lambda called {cap.calls} times / raised "
                                      f"{type(cap.raised).__name__ if cap.raised else None}") from ex
        if cap.calls == 1:
            r = cap.result
            if isinstance(ex, self._nonestimable) or not (
                    isinstance(r, tuple) and len(r) == 2
                    and all(isinstance(x, float) for x in r)):
                raise LambdaCaptureDefect("tuning ran but the capture is unusable") from ex
            return {"lambda_at_failure": r[0], "lambda_status": LAMBDA_CAPTURED,
                    "tuning_score_at_failure": r[1]}
        if isinstance(ex, self._nonestimable):
            return {"lambda_at_failure": None, "lambda_status": LAMBDA_NOT_APPLICABLE,
                    "tuning_score_at_failure": None}
        if isinstance(ex, self._convergence_failure):
            return {"lambda_at_failure": 0.0, "lambda_status": LAMBDA_FIXED_ZERO,
                    "tuning_score_at_failure": None}
        raise LambdaCaptureDefect(f"unclassifiable failure {type(ex).__name__}") from ex


@dataclass
class FitOutcome:
    model: str
    stage: str                       # "apparent" | "apparent_as_run" | "inner"
    ok: bool
    C_fit: Optional[float] = None    # C on the fitting sample (apparent C, or C_boot)
    C_eval: Optional[float] = None   # inner only: C of the inner fit evaluated on D*_k
    lam: Optional[float] = None
    age_knots: Optional[list] = None
    age_linear_fallback: Optional[bool] = None
    inner_cv_c: Optional[float] = None
    cv_rng_consumed: Optional[bool] = None
    reason: Optional[str] = None
    lambda_at_failure: Optional[float] = None   # not observable through pipeline_fit
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        d.pop("extra")
        d.update(self.extra)
        return d


def _state(rng):
    return copy.deepcopy(rng.bit_generator.state)


def fit_and_score(pipeline, ctx, model: str, stage: str, train_ids, eval_ids,
                  t_tr, e_tr, t_ev, e_ev, rng, lam_grid) -> FitOutcome:
    """Fit on train_ids; C on the training sample and (if eval given) on eval_ids."""
    panel = ctx.models[model]
    capturing = getattr(pipeline, "capturing_fit", None)
    cap = TuneCapture() if capturing is not None else None     # one capture object per fit
    fit = capturing(cap) if cap is not None else pipeline.pipeline_fit
    before = _state(rng)
    try:
        m, d_tr, d_ev, dg = fit(ctx.X, ctx.cov, panel, train_ids,
                                eval_ids if eval_ids is not None else train_ids,
                                t_tr, e_tr, lam_grid, rng)
        c_fit = pipeline.c_index(t_tr, e_tr, pipeline.linear_predictor(m, d_tr))
        c_ev = (pipeline.c_index(t_ev, e_ev, pipeline.linear_predictor(m, d_ev))
                if eval_ids is not None else None)
        if cap is not None:
            pipeline.check_success_capture(cap, dg)
        return FitOutcome(model, stage, True, c_fit, c_ev, float(dg["lambda"]),
                          [float(x) for x in dg["age_knots"]],
                          bool(dg.get("age_linear_fallback", False)),
                          (float(dg["inner_cv_c"]) if "inner_cv_c" in dg else None),
                          _state(rng) != before)
    except pipeline.failure_types as ex:
        lam_diag = pipeline.failure_diagnostics(cap, ex) if cap is not None else {
            "lambda_at_failure": None, "lambda_status": None, "tuning_score_at_failure": None}
        return FitOutcome(model, stage, False, reason=getattr(ex, "reason", str(ex)),
                          cv_rng_consumed=_state(rng) != before,
                          lambda_at_failure=lam_diag["lambda_at_failure"],
                          extra={"failure_type": type(ex).__name__,
                                 "lambda_status": lam_diag["lambda_status"],
                                 "tuning_score_at_failure": lam_diag["tuning_score_at_failure"]})
