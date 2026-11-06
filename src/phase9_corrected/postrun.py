"""Decision 13 P1 -- post-run verification runner (entry point scripts/43_phase9_corrected_verify.py,
started only through `scripts/42_phase9_corrected_launch.sh verify <owner-config>`).

BOUNDARY (Decision 13, B3). The runner never calls `run_corrected`, never dispatches,
resumes or tops up a bootstrap replicate, never appends to or modifies a checkpoint, the
manifest or a sizing-decision file, and never bypasses the launcher. Its only writes are
create-only, under the approved output_dir:
  * Stage A window results   <output_dir>/windows/run<r>.inv<m>.json   (p9c-window-1)
  * the verification report  <output_dir>/verification/L<level>.report.json / .txt
An existing file is never replaced: a window result that differs from the fresh
recomputation, or any prior report at the same level, is a refusal / FAIL for owner review.
Stage A, V1-inner and V5 perform the individual REAL-DATA refits / resample draws they
require; executing them is a scientific computation reserved for Decision 14.

Sequence (each step's evidence is recorded; nothing is inferred from a stored verdict):
  0. preconditions -- executable configuration; code identity (DEC-B2 / DEC-E1); stored
     manifest code == freeze (and every invocation's start record bound to this
     configuration, manifest and freeze); H4; historical preservation; no prior report.
  1. history       -- orchestration.assess_run for every run (Decision 11 checks).
  2. Stage A       -- every acknowledged unattested window is recomputed here, record by
                      record, and compared byte-for-byte (canonical JSON) and by float.hex.
  3. Stage B       -- orchestration.evaluate_eligibility at the campaign level, with the
                      Stage A results; every run must be ELIGIBLE.
  4. sizing        -- every level's decision is recomputed with the approved
                      `compute_sizing_decision` (EP-1, Decision 5B / section 6.5 MC-SE,
                      W5-A, Decision 6B, Decision 7B); stored sizing files must equal the
                      recomputation and their snapshots must reproduce; all runs must have
                      moved to the decided level; the final level must be terminal.
  5. record checks, per run: V1-self, V1-outer-historical, V1-inner, V2-record,
     D12-lambda (Decision 12 statuses), V3, V4; then V5 (run 0). V5 FAIL is hard failure H5.

Verdict semantics (Decision 13 P4). Statuses: PASS | FAIL | NOT_APPLICABLE | NOT_COMPARABLE,
each with counts and reasons. Raw BLOCKED / missing / malformed / empty evidence is FAIL,
never NOT_APPLICABLE. NOT_APPLICABLE is accepted only where APPROVED_NOT_APPLICABLE names
the exact configuration condition; anything else that a check reports as not applicable
is FAIL. NOT_COMPARABLE is never PASS. Overall: FAIL if any FAIL; OWNER_REVIEW_REQUIRED if
any NOT_COMPARABLE (including the owner-approved V4 classification, D2) or a sizing
decision requires owner review; PASS otherwise. Only PASS permits reporting.

Review remediation (Decision 13, second review cycle):
  D1  before EVERY output write (each Stage A window result, the JSON report, the text
      report) the runner re-checks loaded-module identity against the freeze and refuses
      importable shadow entries at the repository root; any failure refuses the write.
      LIMITATION: this detects a shadowing or unexpected repository-local module only
      AFTER it was imported -- it cannot authenticate code before that code runs; it
      prevents such a process from writing verification evidence.
  D2  (owner-approved: NOT_COMPARABLE) inside the historical comparison domain (k < 500),
      a V4 record with no historical row, a non-ok historical row, a failed corrected fit or
      a missing corrected as-run field is NOT_COMPARABLE (result `V4:not_comparable`, with
      runs, k positions and reasons); it is never a successful comparison and never a
      scientific FAIL. Actual comparison mismatches remain FAIL. If every in-domain record
      is skipped, V4 itself is NOT_COMPARABLE. Positions k >= 500 lie outside the
      comparison domain and are counted separately. An error in the accounting itself makes
      V4 NOT_COMPARABLE with the exception recorded (owner review), never PASS.
  D3  (gate remediation) every invocation's start record binds its own freshly executed
      gate evidence (path, SHA-256, PASS). The runner re-verifies each one with
      orchestration.start_gate_problems: the file exists, matches the recorded SHA-256 and
      is a complete six-check PASS for this configuration, manifest and NumPy version
      (`gate-evidence` result per invocation; any problem is FAIL, and assess_run makes the
      run non-ELIGIBLE). The files are not cryptographically authenticated: a coordinated
      rewrite of evidence and start record is outside what hashing alone can detect.
  D4  V5 keeps the approved `!=` comparison. Under it a TOP-LEVEL NaN field always differs
      (NaN != NaN), whereas NaN nested in a list or dict compares equal (json returns one
      shared NaN object and container comparison checks identity first) and +/-Infinity
      compares equal. The report classifies non-finite content of every V5 field
      accordingly; the result is unchanged.
  D5  run 0 is not assumed: V5's approved run must be selected, otherwise V5 is an explicit
      owner-review outcome; selected runs must exist and unselected runs must not. Every
      file in runs/ and halts/ is classified: complete artifacts, incomplete artifacts
      (checkpoint lock, *.tmp) and unrecognised run-like names. Any incomplete or
      unrecognised artifact, or an unselected run in any artifact, refuses verification.
  D6  an exception during V5 recomputation is FAIL with its type, message, run and stage;
      it is not H5 (H5 is the approved V5 FAIL condition only).

Non-finite values (DEC-A1): the JSON report is canonical JSON with NaN / Infinity /
-Infinity tokens (not RFC 8259); the text rendering spells each one out and the report
lists every non-finite value's path.
"""
from __future__ import annotations

import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from . import orchestration as orch
from . import verification as V
from .config import ConfigurationError, CorrectedConfig, canonical_json
from .data import LAM_GRID
from .fitting import LAMBDA_CAPTURED, LAMBDA_FIXED_ZERO, LAMBDA_NOT_APPLICABLE

REPORT_SCHEMA = "p9c-verify-2"
WINDOW_SCHEMA = "p9c-window-1"
VERIFICATION_DIR = "verification"
STATUSES = ("PASS", "FAIL", "NOT_APPLICABLE", "NOT_COMPARABLE")
VERDICTS = ("PASS", "FAIL", "OWNER_REVIEW_REQUIRED")
RECORD_CHECKS = ("V1-self", "V1-outer-historical", "V1-inner", "V2-record", "D12-lambda",
                 "V3", "V4")
V5_RUN = 0                                       # plan line 383: k = 0, 25, ..., 475 of run 0
TERMINAL_ACTIONS = ("stop_ep1_pass", "report_persisting_breach", "stop_label_imprecise")
LAMBDA_STATUSES = (LAMBDA_CAPTURED, LAMBDA_NOT_APPLICABLE, LAMBDA_FIXED_ZERO)
APPROVED_GATE_ID = orch.APPROVED_GATE_ID
GATE_CHECKS = orch.GATE_CHECKS
HISTORICAL_B_OUT = orch.HISTORICAL_B_OUT        # the historical runs hold k = 0..499
# importable top-level repository entries permitted on sys.path[0] (= the repository root)
ROOT_IMPORTABLE_ALLOWED = frozenset({"src", "tests"})
_IMPORTABLE_SUFFIXES = (".py", ".pyc", ".pyo", ".so", ".pyd", ".dylib")
D2_POLICY = ("owner-approved: V4 records inside the historical comparison domain (k < 500) "
             "that cannot be compared are NOT_COMPARABLE; k >= 500 is outside the domain")
_RUN_ARTIFACT = {"runs": re.compile(r"^run(\d+)\.(jsonl|lock|jsonl\.tmp)$"),
                 "halts": re.compile(r"^run(\d+)\.inv(\d{4})\.(start|end|auth)\.json(\.tmp)?$")}


def _v4_not_applicable(config: CorrectedConfig) -> Optional[str]:
    """Owner rule (Decision 13 P4, items 5-6): V4 may be NOT_APPLICABLE only when the
    approved configuration has as_run_variant = false under the independent convention."""
    if config.as_run_variant is False and config.convention == "independent":
        return "as_run_variant = false under the independent apparent-CV convention"
    return None


# check -> function(config) returning the exact approved condition, or None
APPROVED_NOT_APPLICABLE = {"V4": _v4_not_applicable}
# owner-approved NOT_COMPARABLE classifications (they still require owner review, P4)
APPROVED_NOT_COMPARABLE: dict = {"V4": "D2: V4 records inside the historical comparison domain "
                                       "that cannot be compared"}


class VerificationRefused(RuntimeError):
    """A precondition failed; nothing was computed and no report was written."""


class HardFailureH5(RuntimeError):
    """V5 failed: the spot-check recomputation is not identical (hard failure H5)."""


@dataclass
class VerificationEnv:
    """Injection points (tests substitute synthetic doubles; the defaults are the real,
    data-reading implementations, used only by the launched entry point)."""
    load_context: Optional[Callable] = None        # root -> AnalysisContext
    pipeline: Optional[Callable] = None            # () -> pipeline
    recompute: Optional[Callable] = None           # (ctx, streams, pipeline, run_id, k, B_in, lam_grid, as_run_variant) -> record
    load_historical: Optional[Callable] = None     # (root, run_id) -> {k: historical row}
    load_frozen_driver: Optional[Callable] = None  # root -> module
    streams: Optional[Callable] = None             # (config, run_id) -> PlanStreams

    @classmethod
    def real(cls) -> "VerificationEnv":
        from .data import load_context
        from .engine import compute_outer_replicate
        from .fitting import RealPipeline
        from .rng import PlanStreams
        return cls(load_context=load_context, pipeline=RealPipeline,
                   recompute=compute_outer_replicate, load_historical=V.load_historical_run,
                   load_frozen_driver=V.load_frozen_driver,
                   streams=lambda c, r: PlanStreams(c.base, r, c.convention, c.T))


# ----------------------------------------------------------------- results and rendering
def result(check, status, *, run_id=None, reasons=(), counts=None, evidence=None, raw_status=None):
    assert status in STATUSES, status
    return {"check": check, "run_id": run_id, "status": status, "raw_status": raw_status,
            "reasons": list(reasons), "counts": dict(counts or {}), "evidence": evidence}


def nonfinite_items(obj, path="$") -> list:
    out = []
    if isinstance(obj, float) and not math.isfinite(obj):
        out.append((path, obj))
    elif isinstance(obj, dict):
        for k in sorted(obj, key=str):
            out.extend(nonfinite_items(obj[k], f"{path}.{k}"))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            out.extend(nonfinite_items(v, f"{path}[{i}]"))
    return out


def nonfinite_paths(obj, path="$") -> list:
    return [p for p, _ in nonfinite_items(obj, path)]


def render_value(v) -> str:
    """DEC-A1 explicit rendering: never an empty cell, never a silent conversion."""
    if isinstance(v, float):
        if math.isnan(v):
            return "NaN"
        if math.isinf(v):
            return "+Infinity" if v > 0 else "-Infinity"
        return repr(v)
    if v is None:
        return "null"
    return json.dumps(v, sort_keys=True, allow_nan=True, default=str)


def render_text(report: dict) -> str:
    lines = [f"PHASE 9 CORRECTED -- POST-RUN VERIFICATION ({report['schema']})",
             f"verdict: {report['verdict']}    reporting permitted: {report['reporting_permitted']}",
             f"level B_out: {render_value(report['level'])}    hard failure H5: {report['hard_failure_h5']}",
             f"freeze_record_sha256: {report['freeze_record_sha256']}",
             f"manifest_sha256: {report['manifest_sha256']}", "",
             "status counts: " + ", ".join(f"{s}={report['status_counts'][s]}" for s in STATUSES), ""]
    for r in report["results"]:
        run = "" if r["run_id"] is None else f" run {r['run_id']}"
        lines.append(f"[{r['status']}] {r['check']}{run}"
                     + (f"  (raw {r['raw_status']})" if r.get("raw_status") else ""))
        for k in sorted(r["counts"]):
            lines.append(f"    {k}: {render_value(r['counts'][k])}")
        for why in r["reasons"]:
            lines.append(f"    reason: {render_value(why)}")
    nf = report["nonfinite_values"]
    lines += ["", f"non-finite values in this report (DEC-A1): {nf['count']}"]
    lines += [f"    {p} = {render_value(v)}"
              for p, v in nonfinite_items({k: v for k, v in report.items() if k != "nonfinite_values"})[:200]]
    for name, text in sorted((report.get("policies") or {}).items()):
        lines += ["", f"policy {name}: {text}"]
    ro = report.get("environment_record_only") or {}
    lines += ["", "record-only environment (never compared at launch):",
              f"    freeze: {render_value(ro.get('freeze'))}"]
    lines.append("")
    lines.append("Synthetic software tests do not establish scientific validity; this report "
                 "is valid only for the frozen package and data it names.")
    return "\n".join(lines) + "\n"


def _create_exclusive(path: Path, data: bytes):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        os.write(fd, data)
        os.fsync(fd)
    finally:
        os.close(fd)


def report_paths(out: Path, level) -> tuple:
    tag = f"L{int(level):05d}" if isinstance(level, int) else "Lunknown"
    d = Path(out) / VERIFICATION_DIR
    return d / f"{tag}.report.json", d / f"{tag}.report.txt"


# ----------------------------------------------------------------- normalisation (P4)
_EMPTY_COUNT_FIELDS = ("n_compared", "n_outer", "n_fits")


def normalise(raw: dict, config: CorrectedConfig, run_id=None) -> list:
    """Map a verification.py result onto the four P4 statuses. Never upgrades evidence."""
    check, st = raw.get("check"), raw.get("status")
    counts = {k: raw[k] for k in raw if k.startswith("n_") and isinstance(raw[k], (int, float))}
    evidence = {k: raw[k] for k in raw if k not in counts and k not in ("check", "status")}
    out = []
    if st == "PASS":
        empty = [f for f in _EMPTY_COUNT_FIELDS if f in raw and raw[f] == 0]
        if empty:
            out.append(result(check, "FAIL", run_id=run_id, raw_status=st, counts=counts,
                              evidence=evidence, reasons=[f"no evidence: {empty[0]} = 0"]))
        else:
            out.append(result(check, "PASS", run_id=run_id, raw_status=st, counts=counts,
                              evidence=evidence))
    elif st == "FAIL":
        out.append(result(check, "FAIL", run_id=run_id, raw_status=st, counts=counts,
                          evidence=evidence, reasons=["check failed"]))
    elif st == "NOT_APPLICABLE":
        rule = APPROVED_NOT_APPLICABLE.get(_base_check(check))
        cond = rule(config) if rule else None
        if cond:
            out.append(result(check, "NOT_APPLICABLE", run_id=run_id, raw_status=st, counts=counts,
                              evidence={**evidence, "approved_condition": cond},
                              reasons=[f"approved inapplicability: {cond}"]))
        else:
            out.append(result(check, "FAIL", run_id=run_id, raw_status=st, counts=counts,
                              evidence=evidence,
                              reasons=["inapplicability is not established by the approved "
                                       "specification (missing or empty evidence is never "
                                       f"NOT_APPLICABLE): {raw.get('reason', 'no comparable records')}"]))
    elif st == "NOT_COMPARABLE":
        out.append(result(check, "NOT_COMPARABLE", run_id=run_id, raw_status=st, counts=counts,
                          evidence=evidence, reasons=[raw.get("reason", "not comparable")]))
    else:                                          # BLOCKED, missing, malformed, unknown
        out.append(result(check or "unknown", "FAIL", run_id=run_id, raw_status=st, counts=counts,
                          evidence=evidence,
                          reasons=[f"evidence missing or malformed: {raw.get('reason', st)}"]))
    n_nc = raw.get("n_not_comparable") or 0
    if n_nc:
        out.append(result(f"{check}:not_comparable", "NOT_COMPARABLE", run_id=run_id,
                          counts={"n_not_comparable": n_nc},
                          reasons=[f"{k}: {why}" for k, why in (raw.get("not_comparable") or [])][:50]))
    return out


def _base_check(check: str) -> str:
    return (check or "").split(":")[0]


# ----------------------------------------------------------------- D1: identity guard
def root_shadow_entries(root: Path) -> list:
    """Importable entries at the repository root (sys.path[0] of every entry point) other
    than the permitted packages: a top-level module file or extension, or a directory that
    is a REGULAR package (has __init__). Namespace directories (no __init__) cannot shadow a
    regular module found later on sys.path and are not reported."""
    bad = []
    for p in sorted(Path(root).iterdir(), key=lambda x: x.name):
        name = p.name
        if p.is_dir():
            inits = [q.name for q in p.iterdir() if q.name.startswith("__init__.")
                     and (q.name.endswith(_IMPORTABLE_SUFFIXES) or ".cpython-" in q.name
                          or ".abi3." in q.name)]
            if inits and name not in ROOT_IMPORTABLE_ALLOWED:
                bad.append(name)
        elif name.endswith(_IMPORTABLE_SUFFIXES) or ".cpython-" in name:
            bad.append(name)
    return bad


def required_set(root: Path) -> set:
    root = Path(root)
    return ({f"{orch.PACKAGE_DIR}/{p.name}" for p in (root / orch.PACKAGE_DIR).glob("*.py")}
            | set(orch.REQUIRED_FILES))


def identity_guard(root: Path, code: dict, what: str):
    """D1: run immediately before every runner output write. Fails closed."""
    try:
        shadows = root_shadow_entries(root)
        if shadows:
            raise orch.CodeIdentityRefusal(f"importable entries at the repository root shadow "
                                           f"imports: {shadows}")
        orch.verify_loaded_modules(root, code, required_set(root))
    except Exception as ex:                                     # any error fails closed
        raise VerificationRefused(f"verification refused: identity re-check before writing "
                                  f"{what} failed: {type(ex).__name__}: {ex}") from ex


# ----------------------------------------------------------------- D2: V4 skip accounting
def v4_skip_accounting(records, historical: dict, config: CorrectedConfig) -> Optional[dict]:
    """Independently enumerate every V4 skip condition (verification.v4_as_run_apparent:
    no historical row, historical status not ok, corrected as-run field absent). None when
    V4 compares nothing by design (NOT_APPLICABLE configuration)."""
    if config.convention == "as_run":
        field = "dC_app"
    elif config.as_run_variant is True:
        field = "dC_app_as_run"
    else:
        return None
    stage = "apparent" if field == "dC_app" else "apparent_as_run"
    reasons = {"beyond_historical_range": [], "historical_row_missing": [],
               "historical_status_not_ok": [], "corrected_fit_failed": [],
               "corrected_field_missing": []}
    n_compared = 0
    for rec in records:
        k = rec["k"]
        h = historical.get(k)
        if k >= HISTORICAL_B_OUT and h is None:
            reasons["beyond_historical_range"].append(k)       # outside the comparison domain
        elif h is None:
            reasons["historical_row_missing"].append(k)
        elif h.get("status") != "ok":
            reasons["historical_status_not_ok"].append(k)
        elif field not in rec:
            fits = rec.get(stage) or {}
            failed = any(isinstance(fd, dict) and fd.get("ok") is False for fd in fits.values())
            reasons["corrected_fit_failed" if failed else "corrected_field_missing"].append(k)
        else:
            n_compared += 1
    policy_relevant = sum(len(reasons[x]) for x in ("historical_row_missing",
                                                    "historical_status_not_ok",
                                                    "corrected_fit_failed",
                                                    "corrected_field_missing"))
    return {"field": field, "n_records": len(records), "n_compared_expected": n_compared,
            "n_skipped": {x: len(v) for x, v in reasons.items()},
            "skipped_k": {x: v[:50] for x, v in reasons.items()},
            "n_policy_relevant_skips": policy_relevant}


# ----------------------------------------------------------------- D3: gate evidence
expected_v2_domain_states = orch.expected_v2_domain_states


def invocation_gate_results(config: CorrectedConfig, msha: str) -> list:
    """One `gate-evidence` result per recorded invocation, independently re-verified."""
    out, res = Path(config.output_dir), []
    for r in config.run_ids:
        invs = orch.invocation_state(out, r)["invocations"]
        for m in sorted(invs):
            st = invs[m].get("start")
            if not st:
                continue
            cat, problems = orch.start_gate_problems(out, st, config, msha)
            ge = st.get("gate_evidence") if isinstance(st.get("gate_evidence"), dict) else {}
            res.append(result("gate-evidence", "FAIL" if cat else "PASS", run_id=r,
                              counts={"invocation": m},
                              evidence={"path": ge.get("path"), "sha256": ge.get("sha256"),
                                        "category": cat},
                              reasons=problems or ["six-check PASS bound to the start record "
                                                   "(not cryptographically authenticated)"]))
    return res


# ----------------------------------------------------------------- D4: NaN disclosure for V5
def _nonfinite_kinds(v, container=None) -> set:
    """top_level_nan | nested_nan_in_list | nested_nan_in_dict | other_nonfinite (innermost
    container named)."""
    if isinstance(v, float) and not math.isfinite(v):
        if math.isnan(v):
            return {"top_level_nan" if container is None else f"nested_nan_in_{container}"}
        return {"other_nonfinite"}
    if isinstance(v, dict):
        return set().union(*(_nonfinite_kinds(x, "dict") for x in v.values())) if v else set()
    if isinstance(v, list):
        return set().union(*(_nonfinite_kinds(x, "list") for x in v)) if v else set()
    return set()


def v5_nonfinite_disclosure(rows, stored: dict, recompute_cache: dict) -> dict:
    """Classify the non-finite content of every V5 field under the approved `!=`:
    top_level_nan (always differs), nested_nan_in_list / nested_nan_in_dict (compares equal:
    shared json NaN object, identity-first container comparison) and other_nonfinite
    (+/-Infinity, compares equal). `nan_only_differences` lists differing fields whose
    canonical JSON is identical, i.e. that differ only because of a top-level NaN."""
    by_kind, nan_only = {}, []
    for row in rows:
        k = row.get("k")
        old, new = stored.get(k), recompute_cache.get(k)
        if old is None or new is None:
            continue
        for f in sorted(set(old) | set(new)):
            for kind in _nonfinite_kinds(old.get(f)) | _nonfinite_kinds(new.get(f)):
                by_kind.setdefault(kind, set()).add(f)
        diff = [f for f in row.get("differing_fields") or []
                if f in old and f in new and canonical_json(old[f]) == canonical_json(new[f])]
        if diff:
            nan_only.append({"k": k, "fields": diff})
    return {"fields_by_kind": {k: sorted(v) for k, v in sorted(by_kind.items())},
            "nan_only_differences": nan_only}


# ----------------------------------------------------------------- Decision 12 statuses
def d12_lambda_check(records, run_id=None) -> dict:
    """Every failed fit carries a Decision 12 Option C status consistent with its failure
    type and lambda/score fields; successful fits carry none; the replicate's `failures`
    list mirrors the fit records. Missing statuses are FAIL (capture missing)."""
    bad, counts = [], {s: 0 for s in LAMBDA_STATUSES}
    n_fits = n_failed = 0
    grid = [float(x) for x in LAM_GRID]
    for rec in records:
        k = rec.get("k")
        fails = {(f.get("stage"), f.get("model"), f.get("b")): f for f in rec.get("failures") or []}
        seen = set()
        fits = []
        for stage in ("apparent", "apparent_as_run"):
            for m, fd in (rec.get(stage) or {}).items():
                fits.append((stage, m, None, fd))
        for row in rec.get("inner") or []:
            for m, fd in (row.get("fits") or {}).items():
                fits.append(("inner", m, row.get("b"), fd))
        for stage, m, b, fd in fits:
            n_fits += 1
            where = (k, stage, m, b)
            if fd.get("ok") is True:
                if fd.get("lambda_status") is not None:
                    bad.append((*where, "successful fit carries a failure lambda status"))
                continue
            n_failed += 1
            st, lam, score = fd.get("lambda_status"), fd.get("lambda_at_failure"), fd.get("tuning_score_at_failure")
            ftype = fd.get("failure_type")
            if st not in LAMBDA_STATUSES:
                bad.append((*where, f"lambda_status {st!r} missing or unknown"))
                continue
            counts[st] += 1
            if st == LAMBDA_CAPTURED:
                if not (m == "M2" and isinstance(lam, float) and any(lam == g for g in grid)
                        and isinstance(score, float) and fd.get("cv_rng_consumed") is True):
                    bad.append((*where, "captured status inconsistent with its lambda/score/model"))
            elif st == LAMBDA_NOT_APPLICABLE:
                if not (lam is None and score is None and ftype == "NonEstimableReplicate"):
                    bad.append((*where, "not_applicable_before_tuning inconsistent"))
            elif st == LAMBDA_FIXED_ZERO:
                if not (lam == 0.0 and score is None and ftype == "ConvergenceFailure"):
                    bad.append((*where, "fixed_zero_no_molecular_block inconsistent"))
            f = fails.get((stage, m, b))
            seen.add((stage, m, b))
            if f is None or f.get("lambda_status") != st or not _same(f.get("lambda_at_failure"), lam) \
                    or not _same(f.get("tuning_score_at_failure"), score):
                bad.append((*where, "failures entry missing or differs from the fit record"))
        for key in fails:
            if key not in seen:
                bad.append((k, *key, "failures entry without a failed fit"))
    return {"check": "D12-lambda", "status": "PASS" if not bad else "FAIL",
            "n_records": len(records), "n_fits": n_fits, "n_failed_fits": n_failed,
            **{f"n_{s}": c for s, c in counts.items()}, "n_mismatches": len(bad),
            "mismatches": [list(x) for x in bad[:50]]}


def _same(a, b) -> bool:
    if isinstance(a, float) and isinstance(b, float):
        return a == b or (math.isnan(a) and math.isnan(b))
    return a == b and type(a) is type(b)


def _v4_results(raw: dict, recs, hist, config: CorrectedConfig, run_id) -> list:
    """D2 (owner-approved NOT_COMPARABLE). Comparison mismatches stay FAIL; in-domain records
    that cannot be compared are NOT_COMPARABLE; k >= 500 is disclosed separately."""
    try:
        acct = v4_skip_accounting(recs, hist, config)
    except Exception as ex:                                     # structured, never PASS
        nc = result("V4:skip-accounting", "NOT_COMPARABLE", run_id=run_id,
                    evidence={"exception_type": type(ex).__name__, "message": str(ex),
                              "stage": "V4 skip accounting"},
                    reasons=["V4 skip accounting failed; V4 evidence is incomplete (owner review)"])
        main = normalise(raw, config, run_id=run_id)            # a real FAIL is never hidden
        if main[0]["status"] in ("PASS", "NOT_APPLICABLE"):    # incomplete evidence: never PASS
            main[0] = result("V4", "NOT_COMPARABLE", run_id=run_id, raw_status=raw.get("status"),
                             counts=main[0]["counts"], evidence=main[0]["evidence"],
                             reasons=["V4 skip accounting failed; the comparison cannot be "
                                      "confirmed complete"])
        return main + [nc]
    if acct is None:                                            # NOT_APPLICABLE configuration
        return normalise(raw, config, run_id=run_id)
    if raw.get("status") in ("PASS", "FAIL") and raw.get("n_compared") != acct["n_compared_expected"]:
        return [result("V4", "FAIL", run_id=run_id, raw_status=raw.get("status"),
                       evidence={"skip_accounting": acct, "n_compared": raw.get("n_compared")},
                       reasons=["V4 comparison count disagrees with the independent skip accounting"])]
    n_rel = acct["n_policy_relevant_skips"]
    if raw.get("status") == "NOT_APPLICABLE" and acct["n_compared_expected"] == 0 and n_rel:
        main = [result("V4", "NOT_COMPARABLE", run_id=run_id, raw_status="NOT_APPLICABLE",
                       counts={"n_compared": 0, "n_not_comparable": n_rel},
                       evidence={"skip_accounting": acct},
                       reasons=["every in-domain V4 record is skipped: NOT_COMPARABLE (D2)"])]
    else:
        main = normalise({**raw, "skip_accounting": acct,
                          "n_beyond_historical_domain": acct["n_skipped"]["beyond_historical_range"]},
                         config, run_id=run_id)
    if n_rel:
        reasons = [f"{why}: k={ks}" for why, ks in acct["skipped_k"].items()
                   if ks and why != "beyond_historical_range"]
        main.append(result("V4:not_comparable", "NOT_COMPARABLE", run_id=run_id,
                           counts={"n_not_comparable": n_rel,
                                   **{f"n_{w}": c for w, c in acct["n_skipped"].items()
                                      if w != "beyond_historical_range"}},
                           evidence={"skipped_k": {w: v for w, v in acct["skipped_k"].items()
                                                   if w != "beyond_historical_range"}},
                           reasons=reasons))
    return main


# ----------------------------------------------------------------- the runner
def run_verification(config: CorrectedConfig, env: Optional[VerificationEnv] = None) -> dict:
    env = env or VerificationEnv.real()
    root, out = Path(config.root), Path(config.output_dir)

    def refuse(msg):
        raise VerificationRefused(f"verification refused: {msg}")

    # ---- 0. preconditions (no real data is read before all of these pass)
    try:
        config.validate_for_execution()
    except ConfigurationError as ex:
        refuse(str(ex))
    if not out.is_dir():
        refuse(f"output_dir {out} does not exist")
    freeze_info = orch.verify_code_identity(config)
    anchor, code = freeze_info["freeze_record_sha256"], freeze_info["freeze"]["code"]
    mpath = out / "manifest" / "input_manifest.json"
    try:
        mraw = mpath.read_bytes()
        manifest = json.loads(mraw)
    except (OSError, ValueError):
        refuse(f"manifest {mpath} missing or unreadable")
    msha = orch._sha(mraw)
    orch.verify_code_identity(config, manifest=manifest)        # manifest code == freeze
    orch.h4_check(root, "verify")                                 # raises HardFailureH4
    pres = orch.historical_preservation_check(root)
    if not pres.get("all_ok"):
        refuse("historical preservation check failed")

    # D5: the selected runs must exist; runs outside the selection must not
    on_disk, incomplete, unrecognised, unrelated = set(), [], [], []
    for d, rx in _RUN_ARTIFACT.items():
        if not (out / d).is_dir():
            continue
        for q in sorted((out / d).iterdir(), key=lambda x: x.name):
            m = rx.match(q.name)
            if m:
                on_disk.add(int(m.group(1)))
                if q.name.endswith((".lock", ".tmp")):
                    incomplete.append(f"{d}/{q.name}")
            elif q.name.startswith("run"):
                unrecognised.append(f"{d}/{q.name}")
            else:
                unrelated.append(f"{d}/{q.name}")
    if unrecognised:
        refuse(f"unrecognised run artifacts; campaign completeness cannot be established: {unrecognised}")
    if incomplete:
        refuse("incomplete run artifacts present (checkpoint lock: an active writer or a stale "
               f"lock; *.tmp: an interrupted write); owner inspection required: {incomplete}")
    missing = [r for r in config.run_ids if r not in on_disk]
    extra = sorted(on_disk - set(config.run_ids))
    if missing:
        refuse(f"selected runs have no campaign records: {missing}")
    if extra:
        refuse(f"the campaign contains runs outside the selected run_ids: {extra}")

    results = []
    add = results.append
    if unrelated:
        add(result("campaign-artifacts", "PASS", counts={"n_unrelated_files": len(unrelated)},
                   evidence={"unrelated_files": unrelated[:50]},
                   reasons=["files in runs/ or halts/ that are not run artifacts; not counted as runs"]))
    results.extend(invocation_gate_results(config, msha))      # D3: every invocation
    # identity binding of every invocation (the first one included) to config/manifest/freeze
    for r in config.run_ids:
        invs = orch.invocation_state(out, r)["invocations"]
        bad = [m for m in sorted(invs) if "start" in invs[m] and (
            invs[m]["start"].get("freeze_record_sha256") != anchor
            or invs[m]["start"].get("manifest_sha256") != msha
            or invs[m]["start"].get("config_sha256") != config.sha256())]
        if not invs:
            add(result("identity-binding", "FAIL", run_id=r, reasons=["run has no invocations"]))
        elif bad:
            add(result("identity-binding", "FAIL", run_id=r, counts={"n_invocations": len(invs)},
                       reasons=[f"invocation {m} is not bound to this configuration, manifest "
                                f"and freeze" for m in bad]))
        else:
            add(result("identity-binding", "PASS", run_id=r, counts={"n_invocations": len(invs)}))

    # ---- 1. history
    assess = {r: orch.assess_run(config, r, msha) for r in config.run_ids}
    levels = {r: assess[r]["levels"].get("L_cur") for r in config.run_ids}
    level = levels[config.run_ids[0]] if len(set(levels.values())) == 1 else None
    jpath, tpath = report_paths(out, level)
    if jpath.exists() or tpath.exists():
        refuse(f"a verification report for this level already exists ({jpath.name}); "
               "it is never replaced -- owner review required")
    if level is None:
        add(result("campaign-level", "FAIL", reasons=[f"runs are at different levels: {levels}"]))
    else:
        add(result("campaign-level", "PASS", counts={"level": level}))

    # ---- 2. Stage A (recomputed here; stored PASS files are never trusted)
    ctx_box = {}

    def ctx():
        if "ctx" not in ctx_box:
            ctx_box["ctx"] = env.load_context(root)
            ctx_box["pipeline"] = env.pipeline()
        return ctx_box["ctx"], ctx_box["pipeline"]

    def recompute(r, k):
        c, p = ctx()
        return env.recompute(c, env.streams(config, r), p, r, k, config.B_in, LAM_GRID,
                             bool(config.as_run_variant))
    stage_a = {}
    for r in config.run_ids:
        a = assess[r]
        by_k = {x["k"]: x for x in a["records"]}
        stage_a[r] = {}
        for w in a["windows"]:
            m = w["invocation"]
            per_k, problems = {}, []
            for k in w["k_list"]:
                old = by_k.get(k)
                if old is None:
                    problems.append(f"record k={k} absent from the checkpoint")
                    per_k[str(k)] = {"record_identical": False, "dC_hex_equal": False}
                    continue
                new = json.loads(canonical_json(recompute(r, k)))
                per_k[str(k)] = {"record_identical": canonical_json(old) == canonical_json(new),
                                 "dC_hex_equal": V._hex(old.get("dC_corrected")) == V._hex(new.get("dC_corrected"))}
            ok = bool(w["k_list"]) and not problems and all(
                v["record_identical"] and v["dC_hex_equal"] for v in per_k.values())
            res = {"schema": WINDOW_SCHEMA, "run_id": r, "invocation": m,
                   "region_sha256": w["region_sha256"], "k_list": w["k_list"],
                   "code_identity_verified": True, "status": "PASS" if ok else "FAIL",
                   "per_k": per_k}
            p = orch.window_result_path(out, r, m)
            counts = {"n_records": len(w["k_list"]),
                      "n_identical": sum(v["record_identical"] and v["dC_hex_equal"] for v in per_k.values())}
            if p.exists():
                try:
                    existing = json.loads(p.read_text())
                except (OSError, ValueError):
                    existing = None
                if existing is None or canonical_json(existing) != canonical_json(res):
                    add(result("StageA-window", "FAIL", run_id=r, counts=counts,
                               evidence={"invocation": m},
                               reasons=["an existing window result differs from the fresh "
                                        "recomputation; it is not replaced (owner review)"]))
                    stage_a[r][m] = existing if isinstance(existing, dict) else {"status": "UNREADABLE"}
                    continue
            else:
                identity_guard(root, code, f"Stage A window result run {r} inv {m}")   # D1
                orch.durable_create(p, res)
            stage_a[r][m] = res
            add(result("StageA-window", "PASS" if ok else "FAIL", run_id=r, counts=counts,
                       evidence={"invocation": m, "region_sha256": w["region_sha256"]},
                       reasons=problems or ([] if ok else ["recomputed records are not identical"])))

    # ---- 3. Stage B
    eligible = {}
    for r in config.run_ids:
        if level is None:
            eligible[r] = False
            continue
        el = orch.evaluate_eligibility(config, r, level, msha, window_results=stage_a[r],
                                       freeze_record_sha256=anchor)
        eligible[r] = el["class"] == "ELIGIBLE"
        add(result("StageB-eligibility", "PASS" if eligible[r] else "FAIL", run_id=r,
                   counts={"n_reasons": len(el["reasons"])},
                   evidence={"class": el["class"], "dispatch": el["dispatch"],
                             "eligibility_result_sha256": el["eligibility_result_sha256"]},
                   reasons=[] if eligible[r] else [f"{x['category']}: {x['code']}" for x in el["reasons"]]))

    # ---- 4. sizing (approved logic, recomputed; cross-run compliance)
    owner_review = False
    checkpoints = {}
    for r in config.run_ids:
        try:
            checkpoints[r] = (out / "runs" / f"run{r}.jsonl").read_bytes()
        except OSError:
            checkpoints[r] = None
    chain = _level_chain(config, assess)
    if level is not None and chain is not None and all(eligible.values()):
        for i, L in enumerate(chain):
            nxt = chain[i + 1] if i + 1 < len(chain) else None
            add_res, review = _verify_level(config, L, nxt, msha, freeze_info, checkpoints)
            results.extend(add_res)
            owner_review = owner_review or review
    else:
        add(result("sizing", "FAIL", reasons=["sizing not verifiable: campaign level, level "
                                              "chain or eligibility not established"]))

    # ---- 5. record-level checks (real data only now)
    hard_h5 = False
    prereq_ok = not any(x["status"] == "FAIL" for x in results)
    if not prereq_ok:
        for r in config.run_ids:
            for chk in RECORD_CHECKS:
                add(result(chk, "NOT_COMPARABLE", run_id=r,
                           reasons=["not executed: a prerequisite check FAILED"]))
        add(result("V5", "NOT_COMPARABLE", run_id=V5_RUN if V5_RUN in config.run_ids else None,
                   reasons=["not executed: a prerequisite check FAILED"]))
    else:
        c, _ = ctx()
        e = c.e
        frozen = env.load_frozen_driver(root)
        identity_guard(root, code, "record checks (historical driver loaded)")      # B2 / D1
        for r in config.run_ids:
            recs = sorted((x for x in assess[r]["records"] if x["k"] < level), key=lambda x: x["k"])
            streams = env.streams(config, r)
            hist = env.load_historical(root, r)
            raws = [("V1-self", lambda: V.v1_self_consistency(recs, e, streams)),
                    ("V1-outer-historical", lambda: V.v1_outer_vs_historical(recs, hist)),
                    ("V1-inner", lambda: V.v1_inner_vs_frozen_path(recs, e, frozen, config.base, r,
                                                                    config.B_in)),
                    ("V2-record", lambda: V.v2_record_assignment(recs, streams)),
                    ("D12-lambda", lambda: d12_lambda_check(recs, r)),
                    ("V3", lambda: V.v3_m1_invariance(recs, hist)),
                    ("V4", lambda: V.v4_as_run_apparent(recs, hist, config.convention,
                                                        config.as_run_variant))]
            for name, fn in raws:
                try:
                    raw = fn()
                except Exception as ex:                  # malformed evidence fails closed
                    raw = {"check": name, "status": "MALFORMED",
                           "reason": f"{type(ex).__name__}: {ex}"}
                raw["check"] = name
                if name == "V4":
                    results.extend(_v4_results(raw, recs, hist, config, r))
                    continue
                results.extend(normalise(raw, config, run_id=r))
        if V5_RUN not in config.run_ids:                       # D5: no substitute run
            add(result("V5", "NOT_COMPARABLE", run_id=None,
                       reasons=[f"the approved V5 specification uses run {V5_RUN} (plan line 383), "
                                f"which is not in run_ids {list(config.run_ids)}: owner decision "
                                "required"]))
        else:
            stored0 = {x["k"]: x for x in assess[V5_RUN]["records"]}
            cache = {}

            def rc5(k):
                cache[k] = json.loads(canonical_json(recompute(V5_RUN, k)))
                return cache[k]
            exc = None
            try:
                raw5 = V.v5_spot_check(stored0, rc5)             # approved `!=` comparison (D4)
            except Exception as ex:                              # D6: FAIL, not H5
                exc = {"exception_type": type(ex).__name__, "message": str(ex),
                       "run_id": V5_RUN, "stage": "V5 recomputation"}
                raw5 = {"check": "V5", "status": "MALFORMED",
                        "reason": f"V5 recomputation raised {type(ex).__name__}: {ex}"}
            raw5["check"] = "V5"
            res5 = normalise(raw5, config, run_id=V5_RUN)
            if exc:
                res5[0]["evidence"] = {**(res5[0]["evidence"] or {}), **exc}
            disc = v5_nonfinite_disclosure(raw5.get("rows") or [], stored0, cache)
            if disc["fields_by_kind"]:
                res5[0]["evidence"] = {**(res5[0]["evidence"] or {}), "nonfinite": disc}
            if raw5.get("status") == "FAIL" and disc["nan_only_differences"]:
                res5[0]["reasons"].append(
                    "some fields differ under the approved != comparison only because they hold "
                    "a top-level NaN (NaN != NaN; their canonical JSON is identical); nested NaN "
                    "and +/-Infinity compare equal; the result is unchanged")
            results.extend(res5)
            hard_h5 = raw5.get("status") == "FAIL"

    # ---- verdict and create-only output
    sc = {s: sum(x["status"] == s for x in results) for s in STATUSES}
    if sc["FAIL"]:
        verdict = "FAIL"
    elif sc["NOT_COMPARABLE"] or owner_review or any(
            x["status"] == "NOT_APPLICABLE" and _base_check(x["check"]) not in APPROVED_NOT_APPLICABLE
            for x in results):
        verdict = "OWNER_REVIEW_REQUIRED"
    else:
        verdict = "PASS"
    report = {"schema": REPORT_SCHEMA, "level": level, "verdict": verdict,
              "reporting_permitted": verdict == "PASS", "hard_failure_h5": hard_h5,
              "freeze_record_sha256": anchor, "manifest_sha256": msha,
              "config_sha256": config.sha256(), "status_counts": sc, "results": results,
              "approved_not_applicable": {k: f.__doc__.strip().splitlines()[0]
                                          for k, f in APPROVED_NOT_APPLICABLE.items()},
              "approved_not_comparable": dict(APPROVED_NOT_COMPARABLE),
              "policies": {"D2_V4_skipped_records": D2_POLICY},
              "environment_record_only": {"freeze": freeze_info["freeze"].get("environment_record_only")},
              "boundary": "verification only: no dispatch, resume, top-up, checkpoint, manifest "
                          "or sizing-file write"}
    nf = nonfinite_paths(report)
    report["nonfinite_values"] = {"count": len(nf), "paths": nf[:200]}
    data = canonical_json(report) + b"\n"
    identity_guard(root, code, f"report {jpath.name}")                              # D1
    _create_exclusive(jpath, data)
    identity_guard(root, code, f"report {tpath.name}")                              # D1
    _create_exclusive(tpath, render_text(report).encode())
    if hard_h5:
        raise HardFailureH5(f"HARD FAILURE H5: V5 spot-check recomputation differs ({jpath})")
    return report


def _level_chain(config: CorrectedConfig, assess: dict) -> Optional[list]:
    """The sequence of levels every run went through (B_out_initial, then each top-up
    target), from the start records. All runs must share it."""
    chains = set()
    for r, a in assess.items():
        st = orch.invocation_state(Path(config.output_dir), r)
        seq = []
        for m in sorted(st["invocations"]):
            t = (st["invocations"][m].get("start") or {}).get("target_B_out")
            if not seq or seq[-1] != t:
                seq.append(t)
        chains.add(tuple(seq))
    if len(chains) != 1:
        return None
    seq = list(chains.pop())
    return seq if seq and seq[0] == config.B_out_initial else None


def _verify_level(config, L, nxt, msha, freeze_info, checkpoints) -> tuple:
    """Recompute the sizing decision at level L from the checkpoint prefixes (k < L) with
    the approved compute_sizing_decision. If the campaign moved on to `nxt`, the stored
    sizing file must exist, reproduce exactly and decide top_up -> nxt; at the final level
    the recomputed decision must be terminal."""
    out = Path(config.output_dir)
    res, review = [], False
    sp = orch._sizing_path(out, L)
    records = {}
    for r in config.run_ids:
        raw = checkpoints.get(r)
        if raw is None:
            return [result("sizing", "FAIL", counts={"level": L},
                           reasons=[f"run {r} checkpoint missing"])], False
        try:
            by_k = orch._records_by_k(raw)
        except ValueError as ex:
            return [result("sizing", "FAIL", counts={"level": L},
                           reasons=[f"run {r} checkpoint malformed: {ex}"])], False
        if any(k not in by_k for k in range(L)):
            return [result("sizing", "FAIL", counts={"level": L},
                           reasons=[f"run {r} lacks records below {L}"])], False
        records[r] = {k: by_k[k] for k in range(L)}
    try:
        dec = orch.compute_sizing_decision(config, L, records)
    except Exception as ex:                                 # fail closed
        return [result("sizing", "FAIL", counts={"level": L},
                       reasons=[f"recomputation failed: {type(ex).__name__}: {ex}"])], False
    ev = {"level": L, "action": dec["action"], "combined_B_next": dec["combined_B_next"],
          "rho_max": (dec.get("ep1") or {}).get("rho_max")}
    if nxt is None:
        if sp.exists():
            res.extend(_verify_sizing_file(config, L, sp, dec, msha, freeze_info, checkpoints))
        if dec["action"] == "top_up":
            res.append(result("sizing", "FAIL", counts={"level": L}, evidence=ev,
                              reasons=[f"the recomputed decision requires a top-up to "
                                       f"{dec['combined_B_next']}: the campaign is not complete"]))
        elif dec["action"] == "stop_for_owner_review":
            review = True
            res.append(result("sizing", "PASS", counts={"level": L}, evidence=ev,
                              reasons=["recomputed decision: stop_for_owner_review (owner review "
                                       "required before reporting)"]))
        elif dec["action"] in TERMINAL_ACTIONS:
            res.append(result("sizing", "PASS", counts={"level": L}, evidence=ev,
                              reasons=[f"terminal decision {dec['action']}"]))
        else:
            res.append(result("sizing", "FAIL", counts={"level": L}, evidence=ev,
                              reasons=[f"unrecognised decision action {dec['action']!r}"]))
        return res, review
    if not sp.exists():
        return [result("sizing", "FAIL", counts={"level": L}, evidence=ev,
                       reasons=[f"the campaign moved to {nxt} without a sizing file at {L}"])], False
    res.extend(_verify_sizing_file(config, L, sp, dec, msha, freeze_info, checkpoints))
    if dec["action"] != "top_up" or dec["combined_B_next"] != nxt:
        res.append(result("sizing", "FAIL", counts={"level": L}, evidence=ev,
                          reasons=[f"recomputed decision {dec['action']} -> {dec['combined_B_next']} "
                                   f"does not authorise the move to {nxt}"]))
    else:
        res.append(result("sizing", "PASS", counts={"level": L}, evidence=ev,
                          reasons=[f"top_up -> {nxt} reproduced for every run"]))
    return res, review


def _verify_sizing_file(config, L, sp, dec, msha, freeze_info, checkpoints) -> list:
    out = Path(config.output_dir)
    bad = []
    try:
        sraw = sp.read_bytes()
        art = json.loads(sraw)
    except (OSError, ValueError):
        return [result("sizing-file", "FAIL", counts={"level": L}, reasons=["unreadable"])]
    code = freeze_info["freeze"]["code"]
    if art.get("kind") != "sizing_decision" or art.get("schema") != "p9c-sizing-1" \
            or art.get("B_cur") != L or art.get("run_ids") != list(config.run_ids):
        bad.append("schema, level or run_ids invalid")
    if art.get("config_sha256") != config.sha256() or art.get("manifest_sha256") != msha:
        bad.append("configuration or manifest identity differs")
    ci = art.get("code_identity") or {}
    if ci.get("freeze_record_sha256") != freeze_info["freeze_record_sha256"] or \
            ci.get("package") != {k: code[k] for k in sorted(code)}:
        bad.append("code identity does not resolve to the verified freeze")
    if orch._canon(art.get("decision")) != orch._canon(dec):
        bad.append("stored decision differs from the recomputation")
    for r in config.run_ids:
        e = (art.get("runs") or {}).get(str(r)) or {}
        raw = checkpoints.get(r) or b""
        size = e.get("checkpoint_size_bytes")
        if not isinstance(size, int) or len(raw) < size or orch._sha(raw[:size]) != e.get("checkpoint_sha256"):
            bad.append(f"run {r} snapshot is not a prefix of the checkpoint")
            continue
        el = orch.evaluate_eligibility(config, r, L, msha, max_invocation=e.get("latest_invocation"),
                                       checkpoint_raw=raw[:size],
                                       freeze_record_sha256=freeze_info["freeze_record_sha256"])
        if el["class"] != "ELIGIBLE" or el["eligibility_result_sha256"] != e.get("eligibility_result_sha256"):
            bad.append(f"run {r} snapshot eligibility is not reproduced ({el['class']})")
    return [result("sizing-file", "FAIL" if bad else "PASS", counts={"level": L},
                   evidence={"sizing_sha256": orch._sha(sraw)}, reasons=bad)]


# ----------------------------------------------------------------- entry point
EXIT = {"PASS": 0, "FAIL": 1, "OWNER_REVIEW_REQUIRED": 3}


def main(argv, root: Path) -> int:
    from .config import load_config
    from .integrity import HardFailureH4
    if len(argv) != 1:
        print("usage: 43_phase9_corrected_verify.py <owner-config.json>")
        return 2
    try:
        cfg = load_config(argv[0], root)
        rep = run_verification(cfg)
    except HardFailureH5 as ex:
        print(str(ex))
        return 5
    except HardFailureH4 as ex:
        print(str(ex))
        return 4
    except (VerificationRefused, orch.CodeIdentityRefusal, ConfigurationError, OSError, ValueError) as ex:
        print(f"REFUSED: {ex}")
        return 2
    print(json.dumps({"verdict": rep["verdict"], "level": rep["level"],
                      "status_counts": rep["status_counts"]}))
    return EXIT[rep["verdict"]]
