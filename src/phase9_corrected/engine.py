"""One outer replicate of the nested bootstrap (plan section 3.1, lines 138-161).

Execution path for replicate (run_id, k):
  1. outer resample   I*_k = strat_boot(e, G(outer(k)))                     [PLAN 142]
  2. apparent fits    for M in (M1, M2): fit on D*_k with G(apparent_cv(k))  [PLAN 145-148]
                      fresh generator per model                               [PROP D2]
     optional variant same, with G(outer(k)) (as-run apparent convention)    [PROP D2]
  3. inner loop b     J*_{k,b} = strat_boot(eO, G(inner(k,b)))               [PLAN 151, 363]
                      for M: fit on D*_k[J], CV with G(inner_cv(k,b)) fresh  [PLAN 364]
                      opt(M,k,b) = C(D*_k[J]) - C(D*_k)                       [PLAN 154]
  4. pairing          b enters only if both models succeed                    [PLAN 271]
  5. aggregation      dOpt = mean opt(M2) - mean opt(M1); dC_corr = dC_app - dOpt [PLAN 156-158]
Arithmetic is identical to the frozen driver (lines 162-166) so that results are
bit-comparable wherever the streams coincide.
"""
from __future__ import annotations

import numpy as np

from . import SCHEMA_VERSION
from .fitting import fit_and_score
from .resampling import (INNER_MULT_TAG, OUTER_MULT_TAG, binary_cell_counts, draw_inner,
                         draw_outer, tagged_sha256)
from .rng import seedseq_descriptor, state_fingerprint

CLINICAL_BINARIES = {"sex": ("Male", "Female"), "tnm.stage": ("3.0", "2.0"),
                     "tumor.location": ("proximal", "distal"),
                     "chemotherapy.adjuvant": ("Y", "N")}   # src/phase9_execution.py:33-35
MODEL_ORDER = ("M1", "M2")


class ImplementationDefectH6(RuntimeError):
    """Paired-complete counts differ between M1 and M2 (plan line 459)."""


def _stream(ss) -> dict:
    d = seedseq_descriptor(ss)
    d["state_sha256"] = state_fingerprint(ss)
    return d


def compute_outer_replicate(ctx, streams, pipeline, run_id: int, k: int, B_in: int,
                            lam_grid, as_run_variant: bool = False,
                            record_binary_cells: bool = True) -> dict:
    ids, t, e = ctx.ids, ctx.t, ctx.e
    n = len(ids)
    ss_out = streams.outer(k)
    outer = draw_outer(e, streams.generator(ss_out))
    I = outer.indices
    tO, eO = t[I], e[I]
    oids = [ids[i] for i in I]
    idpos = {s: i for i, s in enumerate(ids)}

    def consumed(sample_ids, tag):
        """Multiplicity hash recomputed from the id list actually passed to the fit
        (independent of the index arrays) -- evidence that each fit consumed its resample."""
        return tagged_sha256(tag, np.bincount([idpos[s] for s in sample_ids], minlength=n))

    rec = {"schema": SCHEMA_VERSION, "run_id": int(run_id), "k": int(k),
           "B_in_attempted": int(B_in), "stream_policy": streams.kind,
           "apparent_convention": getattr(streams, "convention", "as_run"),
           "T": getattr(streams, "T", None),
           "streams": {"outer": _stream(ss_out), "apparent_cv": _stream(streams.apparent_cv(k))},
           "outer": outer.evidence(), "failures": [], "status": "ok"}

    # ---- 2. apparent on D*_k, both models, SAME outer sample
    def apparent(stage, ss):
        res = {}
        for name in MODEL_ORDER:
            fo = fit_and_score(pipeline, ctx, name, stage, oids, None, tO, eO, None, None,
                               streams.generator(ss), lam_grid)
            fo.extra["consumed_multiplicity_sha256"] = consumed(oids, OUTER_MULT_TAG)
            res[name] = fo
            if not fo.ok:
                rec["failures"].append({"run_id": int(run_id), "k": int(k), "b": None,
                                        "stage": stage, "model": name, "reason": fo.reason,
                                        "failure_type": fo.extra.get("failure_type"),
                                        "n_distinct": outer.n_distinct, "events": outer.events,
                                        "lambda_at_failure": fo.lambda_at_failure,       # D12
                                        "lambda_status": fo.extra.get("lambda_status"),
                                        "tuning_score_at_failure": fo.extra.get("tuning_score_at_failure"),
                                        "binary_cells": (binary_cell_counts(ctx.cov, oids, CLINICAL_BINARIES)
                                                         if record_binary_cells else None)})
        return res

    app = apparent("apparent", streams.apparent_cv(k))
    rec["apparent"] = {m: fo.to_dict() for m, fo in app.items()}
    if as_run_variant:
        ss_ar = streams.apparent_cv_as_run(k)
        rec["streams"]["apparent_cv_as_run"] = _stream(ss_ar)
        var = apparent("apparent_as_run", ss_ar)
        rec["apparent_as_run"] = {m: fo.to_dict() for m, fo in var.items()}
        if all(fo.ok for fo in var.values()):
            rec["dC_app_as_run"] = var["M2"].C_fit - var["M1"].C_fit
    if not all(fo.ok for fo in app.values()):
        rec["status"] = "unreliable_apparent_failed"
        return rec
    rec["C_app"] = {m: app[m].C_fit for m in MODEL_ORDER}
    rec["dC_app"] = app["M2"].C_fit - app["M1"].C_fit

    # ---- 3. inner optimism, evaluated on D*_k
    inner_rows, used = [], {m: [] for m in MODEL_ORDER}
    for b in range(B_in):
        ss_in, ss_cv = streams.inner(k, b), streams.inner_cv(k, b)
        J = draw_inner(eO, I, n, streams.generator(ss_in))
        jids = [oids[i] for i in J.positions]
        tJ, eJ = tO[J.positions], eO[J.positions]
        row = {"b": int(b), "resample": J.evidence(),
               "streams": {"inner": _stream(ss_in), "inner_cv": _stream(ss_cv)}, "fits": {}}
        got = {}
        for name in MODEL_ORDER:
            fo = fit_and_score(pipeline, ctx, name, "inner", jids, oids, tJ, eJ, tO, eO,
                               streams.generator(ss_cv), lam_grid)
            fd = fo.to_dict()
            fd["resample_index_sha256"] = J.index_sha256     # pairing evidence, per fit
            fd["consumed_multiplicity_sha256"] = consumed(jids, INNER_MULT_TAG)
            row["fits"][name] = fd
            if fo.ok:
                got[name] = (fo.C_fit, fo.C_eval)
            else:
                rec["failures"].append({"run_id": int(run_id), "k": int(k), "b": int(b),
                                        "stage": "inner", "model": name, "reason": fo.reason,
                                        "failure_type": fo.extra.get("failure_type"),
                                        "n_distinct": J.n_distinct, "events": J.events,
                                        "lambda_at_failure": fo.lambda_at_failure,       # D12
                                        "lambda_status": fo.extra.get("lambda_status"),
                                        "tuning_score_at_failure": fo.extra.get("tuning_score_at_failure"),
                                        "binary_cells": (binary_cell_counts(ctx.cov, jids, CLINICAL_BINARIES)
                                                         if record_binary_cells else None)})
        row["paired_complete"] = len(got) == 2              # 4. pairing-preserving exclusion
        if row["paired_complete"]:
            for name in MODEL_ORDER:
                used[name].append([got[name][0], got[name][1]])
        inner_rows.append(row)
    rec["inner"] = inner_rows
    paired_ok = sum(r["paired_complete"] for r in inner_rows)
    rec["B_in_paired_complete"] = int(paired_ok)
    rec["n_used"] = {m: len(used[m]) for m in MODEL_ORDER}
    if len(used["M1"]) != len(used["M2"]):
        raise ImplementationDefectH6(f"run {run_id} k {k}: {rec['n_used']}")
    if paired_ok < 0.95 * B_in:
        rec["status"] = "unreliable_inner_shortfall"
        return rec

    # ---- 5. aggregation (driver lines 162-166, same arithmetic)
    o1 = np.array(used["M1"]); o2 = np.array(used["M2"])
    opt1 = o1[:, 0] - o1[:, 1]; opt2 = o2[:, 0] - o2[:, 1]
    rec["opt_mean"] = {"M1": float(opt1.mean()), "M2": float(opt2.mean())}
    rec["C_corrected"] = {"M1": float(rec["C_app"]["M1"] - opt1.mean()),
                          "M2": float(rec["C_app"]["M2"] - opt2.mean())}
    rec["dOpt"] = float(opt2.mean() - opt1.mean())
    rec["dC_corrected"] = float(rec["dC_app"] - rec["dOpt"])
    rec["dOpt_per_replicate_min"] = float((opt2 - opt1).min())
    rec["dOpt_per_replicate_max"] = float((opt2 - opt1).max())
    return rec
