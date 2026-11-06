"""Execution orchestration for one corrected run (fail-fast; writes only to output_dir).

Order: validate config -> H4 pre -> historical preservation -> output-dir newness ->
manifest (write-once) -> gate -> checkpoint (lock, strict load) -> invocation state check
(+ explicit authorisation if pending) -> durable invocation-start record -> dispatch
missing k -> append each completed record IN k ORDER -> H4 after every record -> plan 9.5
stop conditions -> H4 post -> durable invocation-end record.

Records are collected with the ORDERED Pool.imap, so the checkpoint byte stream is
independent of worker count and completion order (resumption appends the missing k in
ascending order, reproducing the uninterrupted file exactly).

Owner Decision 10 (2026-11-06) -- runtime budget and halt control
-----------------------------------------------------------------
* Budget: `runtime_budget_hours_per_run` = expected wall-clock hours per run per 500 NEW
  replicates at the configured worker count. Per invocation (one `run` command):
      threshold_s = 2 x budget_hours x 3600 x n_new / 500
      W           = max(50, 2 x n_workers)
      projection  = elapsed_s / done_now x n_new, evaluated only when done_now >= W
  where done_now counts records appended IN THIS INVOCATION (never the cumulative
  checkpoint count). Halt if projection > threshold. Pool throughput is measured
  directly (done_now / elapsed); nothing is divided by the worker count.
* No elapsed guard and no hang detection (10C, 10D). When n_new < W the projection is
  never evaluated; this is an accepted limitation recorded as
  runtime_check = "not_evaluated (n_new < W)" (10E).
* Clock (10F): macOS CLOCK_MONOTONIC_RAW, Linux CLOCK_BOOTTIME, otherwise refuse.
* Invocation sentinel (10G): <output_dir>/halts/run<r>.inv<m>.start.json is written
  durably BEFORE any replicate is dispatched; run<r>.inv<m>.end.json records the outcome.
  A start without an end is UNTERMINATED (cause unknown). A runtime halt (end outcome
  runtime_halt) or an UNTERMINATED invocation must be resolved with an explicit,
  single-use `resume_authorization` (run<r>.inv<m>.auth.json) before the run continues.
* Non-runtime outcomes (5 % non-estimable stop, H4, historical preservation) are recorded
  with their own outcome code and keep their existing handling: no authorisation is
  required by this mechanism (10H).
* Durability: every sentinel/halt/authorisation record is create-only (tmp -> fsync ->
  link -> fsync dir -> unlink tmp -> fsync dir), using F_FULLFSYNC on macOS and fsync on
  Linux. Any failure raises DurabilityError (fail closed); leftover *.tmp files make the
  state AMBIGUOUS, which is refused until inspected.

Owner Decision 11 (2026-11-06) -- 5 % stop, dispatch control, completion eligibility
------------------------------------------------------------------------------------
* 5 % rule (literal, Option A): stop when 20 * n_nonestimable > n_completed, cumulative
  over all records, evaluated after every append; any status != "ok" counts. A violation
  at ANY prefix of the checkpoint is permanent (prefix scan), whatever the end-record label.
* Both stop reasons are evaluated on every record; runtime keeps outcome precedence
  (authorisation semantics) and every end record carries `stop_reasons`, a full
  `nonestimable` snapshot and `checkpoint_size_bytes`.
* Code identity (Decision 9 condition 4) is verified on EVERY path before any gate,
  manifest, data load or dispatch, against the approved Decision 13 freeze record
  (fail closed while it does not exist). Path C additionally verifies the Decision 7B
  sizing file by independent recomputation.
* Before dispatch, the run's full history is assessed with the same predicates as the
  completion-eligibility evaluation (chain continuity, region contents, 5 % prefix scan,
  H4/H6, global preservation block, identities, authorisation references, target
  progression); a blocked assessment refuses BEFORE any authorisation is consumed.
* Paths: A fresh (500, no records, no top-up reference), B resume at the current level
  (same level reference above 500), C top-up (previous level complete, strictly higher
  rung, new reference, verified sizing file).
* Eligibility classes (pure Stage B): INELIGIBLE > UNVERIFIABLE > INCOMPLETE > ELIGIBLE.

Owner decisions 2026-11-06 (DEC-A1..F1, D11-D1/D2)
-------------------------------------------------
* D11-D1 / DEC-C1: window results are assessed in `assess_run`, so a window FAIL blocks
  dispatch exactly as it disqualifies eligibility; pending windows permit only continuation
  of the already-authorised level and block eligibility, new top-ups and reporting.
* D11-D2: Path C assesses the complete current history of every run, independent of
  checkpoint growth (`_later_invocations_problem`).
* DEC-E1: scientific execution is launched only through scripts/42_phase9_corrected_launch.sh
  (external SHA-256 check of the freeze inventory before Python starts, then
  `python -I -B -X pycache_prefix=<fresh external dir>`); `verify_code_identity` checks the
  resulting flags. MANUAL ANCHOR STEP before each execution campaign: (1) the SHA-256 of
  config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json equals the owner-approved
  `freeze_record_sha256`; (2) the SHA-256 of the launcher equals its entry in that record.
  The launcher cannot authenticate itself.
* Decision 13 gate remediation: the six-check pre-execution gate is executed FRESHLY for
  every invocation attempt, into its own new evidence file
  <output_dir>/gate/run<r>.inv<mmmm>.attempt<aa>.json (gate.py never finds a cache there, so
  it never reuses one). Dispatch requires the written evidence to pass
  `gate_evidence_problems` (all six checks PASS, correct key and domain), not merely an
  `all_pass` flag. A failed attempt keeps its evidence and writes no start record. The
  start record binds the evidence path, SHA-256 and result; `assess_run` re-verifies it for
  every invocation (missing -> blocking, altered / invalid -> disqualifier).
* DEC-F1 (accepted limitation): a crash between H4 detection and the durable end record
  leaves the invocation UNTERMINATED with the H4 event lost if the frozen inputs are later
  restored. PROCEDURAL SAFEGUARD: before acknowledging ANY UNTERMINATED invocation, the owner
  reviews the process output and logs for an H4 message. Decision 10 semantics are unchanged.
"""
from __future__ import annotations

import fcntl
import functools
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from . import SCHEMA_VERSION
from .checkpoint import CheckpointStore
from .config import ConfigurationError, CorrectedConfig, canonical_json
from .data import LAM_GRID, load_context
from .engine import ImplementationDefectH6, compute_outer_replicate
from .fitting import RealPipeline
from .gate import run_pre_execution_gate
from .integrity import (PROTECTED_HISTORICAL, HardFailureH4, build_manifest, h4_check,
                        historical_preservation_check, write_manifest_once)
from .rng import PHASE9_PRIMARY_SEEDS, PlanStreams

_W = {}


class StopCondition(SystemExit):
    """Plan section 9.5 stop condition: halt and report."""


class NonEstimableStop(StopCondition):
    """Plan 9.5: outer non-estimable rate > 5 % (outcome code stop_nonestimable)."""

    def __init__(self, msg, nonestimable: Optional[dict] = None, reasons=None):
        super().__init__(msg)
        self.nonestimable = nonestimable
        self.reasons = list(reasons or ["nonestimable_rate"])


class RuntimeBudgetHalt(StopCondition):
    """Plan 9.5: projected runtime > twice the budget (Decision 10). When the 5 % rule
    fires on the same record, both reasons are carried (Decision 11)."""

    def __init__(self, msg, details: dict, nonestimable: Optional[dict] = None, reasons=None):
        super().__init__(msg)
        self.details = details
        self.nonestimable = nonestimable
        self.reasons = list(reasons or ["runtime_projection"])


class HistoricalPreservationFailure(SystemExit):
    """A protected historical artifact changed during the run."""

    def __init__(self, msg, evidence: Optional[dict] = None):
        super().__init__(msg)
        self.evidence = evidence


class DurabilityError(RuntimeError):
    """A durable record could not be committed; nothing further may be dispatched."""


class InvocationStateError(ConfigurationError):
    """The run's invocation state forbids this invocation (pending or ambiguous)."""


class DispatchRefused(InvocationStateError):
    """Decision 11: the run's history or the requested path forbids any further dispatch."""


class CodeIdentityRefusal(InvocationStateError):
    """Decision 9 condition 4: code / environment do not match the approved freeze."""


def _work(k):
    return compute_outer_replicate(_W["ctx"], _W["streams"], _W["pipeline"], _W["run_id"], k,
                                   _W["B_in"], LAM_GRID, as_run_variant=_W["variant"])


def execute_replicates(store: CheckpointStore, ctx, streams, pipeline, B_in: int,
                       target_B_out: int, n_workers: int, as_run_variant: bool,
                       on_record=None) -> list:
    """Compute every k in [0, target) not yet in `store`, append in ascending k.
    `on_record(rec, store)` runs after each append (H4, stop conditions)."""
    todo = [k for k in range(target_B_out) if k not in store.done]
    _W.update(ctx=ctx, streams=streams, B_in=B_in, variant=bool(as_run_variant),
              pipeline=pipeline, run_id=store.run_id)
    if not todo:
        return todo
    if n_workers == 1:
        it = map(_work, todo)
        for rec in it:
            store.append(rec)
            if on_record:
                on_record(rec, store)
        return todo
    import multiprocessing as mp
    with mp.get_context("fork").Pool(n_workers) as pool:     # workers inherit _W (fork)
        for rec in pool.imap(_work, todo, chunksize=1):
            store.append(rec)
            if on_record:
                try:
                    on_record(rec, store)
                except BaseException:
                    pool.terminate()
                    raise
    return todo


def allowed_targets(config: CorrectedConfig) -> list:
    return [config.B_out_initial] + list((config.get("ladder") or {}).get("rungs", []))


def log_event(output_dir: Path, event: dict):
    """Informational event log (not authoritative; the halts/ records are)."""
    p = Path(output_dir) / "logs" / "events.jsonl"
    p.parent.mkdir(parents=True, exist_ok=True)
    event = {"utc": datetime.now(timezone.utc).isoformat(), **event}
    with open(p, "a") as fh:
        fh.write(json.dumps(event) + "\n")


def assert_output_dir_new_or_ours(out: Path):
    """The output directory must be new, or a previous run of THIS pipeline (it holds our
    manifest). Any other existing content is refused."""
    if not out.exists() or (out / "manifest" / "input_manifest.json").exists():
        return
    names = {p.name for p in out.iterdir()}
    if names and not names <= {PRESERVATION_DIR}:     # only pre-dispatch preservation records
        raise ConfigurationError(f"output_dir {out} exists, is not empty and holds no p9c manifest")


# ======================================================================= Decision 10F
def select_clock(platform: Optional[str] = None) -> tuple:
    """Approved sleep-inclusive, non-adjustable clocks. Returns (name, callable)."""
    plat = sys.platform if platform is None else platform
    if plat == "darwin":
        name = "CLOCK_MONOTONIC_RAW"
    elif plat.startswith("linux"):
        name = "CLOCK_BOOTTIME"
    else:
        raise ConfigurationError(f"unsupported platform {plat!r}: no approved sleep-inclusive "
                                 "clock (Decision 10F); refusing to start")
    cid = getattr(time, name, None)
    if cid is None:
        raise ConfigurationError(f"{name} is not available in this Python build; refusing to start")
    return name, functools.partial(time.clock_gettime, cid)


# ======================================================================= Decision 10A/B/E
RUNTIME_REF_REPLICATES = 500
WARMUP_MIN = 50
RUNTIME_NOT_EVALUATED = "not_evaluated (n_new < W)"
RUNTIME_EVALUATED = "evaluated_after_warmup"


class RuntimeMonitor:
    """Plan 9.5 projected-runtime rule, per invocation (Decision 10)."""

    def __init__(self, budget_hours, n_start: int, target_B_out: int, n_workers: int,
                 clock: Callable[[], float]):
        if isinstance(budget_hours, bool) or not isinstance(budget_hours, (int, float)) \
                or not budget_hours > 0:
            raise ConfigurationError("runtime_budget_hours_per_run must be a positive number")
        self.budget_hours = float(budget_hours)
        self.n_start = int(n_start)
        self.n_new = int(target_B_out) - self.n_start
        if self.n_new < 0:
            raise ConfigurationError("checkpoint already holds more records than the target")
        self.n_workers = int(n_workers)
        self.W = max(WARMUP_MIN, 2 * self.n_workers)
        self.threshold_s = 2 * self.budget_hours * 3600 * self.n_new / RUNTIME_REF_REPLICATES
        self.clock = clock
        self.t0 = None

    @property
    def runtime_check(self) -> str:
        return RUNTIME_EVALUATED if self.n_new >= self.W else RUNTIME_NOT_EVALUATED

    def describe(self) -> dict:
        return {"n_start": self.n_start, "n_new": self.n_new, "W": self.W,
                "n_workers": self.n_workers, "budget_hours_per_500": self.budget_hours,
                "threshold_s": self.threshold_s, "runtime_check": self.runtime_check}

    def start(self):
        self.t0 = self.clock()

    def elapsed(self) -> float:
        return self.clock() - self.t0

    def observe(self, n_records_total: int) -> Optional[dict]:
        """Call after each append with the store's TOTAL record count; only records
        appended in this invocation count. Returns halt details, or None."""
        done_now = int(n_records_total) - self.n_start
        if done_now < self.W:
            return None                                # warm-up: the clock is not read
        el = self.elapsed()
        proj = el / done_now * self.n_new
        if proj > self.threshold_s:
            return {"reason": "runtime_projection", "done_now": done_now, "elapsed_s": el,
                    "projection_s": proj, **self.describe()}
        return None


# ======================================================================= Decision 10G: durability
def _full_sync(fd: int):
    """F_FULLFSYNC on macOS (forces the drive cache), fsync on Linux. Fail closed."""
    if sys.platform == "darwin":
        op = getattr(fcntl, "F_FULLFSYNC", None)
        if op is None:
            raise DurabilityError("F_FULLFSYNC unavailable: refusing to downgrade durability")
        try:
            fcntl.fcntl(fd, op)
        except OSError as ex:
            raise DurabilityError(f"F_FULLFSYNC failed: {ex}") from ex
    elif sys.platform.startswith("linux"):
        try:
            os.fsync(fd)
        except OSError as ex:
            raise DurabilityError(f"fsync failed: {ex}") from ex
    else:
        raise DurabilityError(f"no approved durability mechanism on {sys.platform!r}")


def _sync_path(path: Path):
    try:
        fd = os.open(path, os.O_RDONLY)
    except OSError as ex:
        raise DurabilityError(f"cannot open {path} for sync: {ex}") from ex
    try:
        _full_sync(fd)
    finally:
        os.close(fd)


def durable_create(path: Path, obj: dict):
    """Create-only durable record: tmp(O_EXCL) -> write -> sync -> link (fails if the final
    name exists) -> sync dir -> unlink tmp -> sync dir. Any failure raises DurabilityError
    and may leave the *.tmp behind, which makes the invocation state AMBIGUOUS."""
    path = Path(path)
    d = path.parent
    missing = []
    q = d
    while not q.exists():
        missing.append(q)
        q = q.parent
    for level in reversed(missing):              # one level at a time, parent synced each time
        try:
            level.mkdir()
        except OSError as ex:
            raise DurabilityError(f"cannot create {level}: {ex}") from ex
        _sync_path(level.parent)
    tmp = path.with_name(path.name + ".tmp")
    data = canonical_json(obj) + b"\n"
    try:
        fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except OSError as ex:
        raise DurabilityError(f"cannot create {tmp}: {ex}") from ex
    try:
        view = memoryview(data)
        while view:
            n = os.write(fd, view)
            view = view[n:]
        _full_sync(fd)
    except OSError as ex:
        raise DurabilityError(f"write to {tmp} failed: {ex}") from ex
    finally:
        os.close(fd)
    try:
        os.link(tmp, path)
    except OSError as ex:
        raise DurabilityError(f"cannot commit {path} (create-only): {ex}") from ex
    _sync_path(d)
    try:
        os.unlink(tmp)
    except OSError as ex:
        raise DurabilityError(f"cannot remove {tmp} after commit: {ex}") from ex
    _sync_path(d)


# ======================================================================= Decision 10G/H: state
HALTS_DIR = "halts"
OUTCOMES = ("completed", "runtime_halt", "stop_nonestimable", "hard_failure_h4",
            "historical_preservation_failure", "hard_failure_h6")      # h6: Decision 11
AUTH_FOR_HALT = "runtime_halt_authorization"
AUTH_FOR_UNTERMINATED = "unterminated_acknowledgement"
_REC_NAME = re.compile(r"^run(\d+)\.inv(\d{4})\.(start|end|auth)\.json(\.tmp)?$")


def _inv_path(out: Path, run_id: int, m: int, kind: str) -> Path:
    return Path(out) / HALTS_DIR / f"run{run_id}.inv{m:04d}.{kind}.json"


def invocation_state(out: Path, run_id: int, max_invocation: Optional[int] = None) -> dict:
    """Derive the run's state from its halts/ records only.
    CLEAN        nothing pending (or every halt/unterminated invocation resolved)
    HALTED       latest invocation ended in runtime_halt, not yet authorised
    UNTERMINATED latest invocation has a start record and no end record, not acknowledged
    AMBIGUOUS    leftover *.tmp, unparseable/inconsistent record, orphan end/auth,
                 more than one pending, or a pending invocation that is not the latest"""
    d = Path(out) / HALTS_DIR
    inv, problems, seen_m = {}, [], set()
    if d.exists():
        for f in sorted(d.iterdir()):
            mt = _REC_NAME.match(f.name)
            if not mt or int(mt.group(1)) != run_id:
                if f.name.startswith(f"run{run_id}.inv"):
                    problems.append(f"unrecognised record name {f.name}")
                continue
            m, kind = int(mt.group(2)), mt.group(3)
            if max_invocation is not None and m > max_invocation:
                continue                                 # snapshot view (Decision 11 Path C)
            seen_m.add(m)
            if mt.group(4):
                problems.append(f"incomplete record {f.name} (interrupted durable write)")
                continue
            try:
                obj = json.loads(f.read_text())
            except (OSError, json.JSONDecodeError) as ex:
                problems.append(f"unreadable record {f.name}: {ex}")
                continue
            if obj.get("kind") != kind or obj.get("run_id") != run_id or obj.get("invocation") != m:
                problems.append(f"record {f.name} identity mismatch")
                continue
            inv.setdefault(m, {})[kind] = obj
    pending = []
    for m in sorted(inv):
        v = inv[m]
        end, auth = v.get("end"), v.get("auth")
        if "start" not in v:
            problems.append(f"invocation {m}: end/auth record without a start record")
            continue
        if end is not None and end.get("outcome") not in OUTCOMES:
            problems.append(f"invocation {m}: unknown outcome {end.get('outcome')!r}")
        if auth is not None:
            want = (AUTH_FOR_HALT if end is not None and end.get("outcome") == "runtime_halt"
                    else AUTH_FOR_UNTERMINATED if end is None else None)
            if auth.get("authorization_kind") != want:
                problems.append(f"invocation {m}: authorisation does not match its outcome")
        if auth is None and (end is None or end.get("outcome") == "runtime_halt"):
            pending.append(m)
    if len(pending) > 1:
        problems.append(f"more than one unresolved invocation: {pending}")
    if pending and inv and pending[-1] != max(inv):
        problems.append(f"unresolved invocation {pending[-1]} is not the latest")
    latest = max(inv) if inv else None
    if problems:
        state = "AMBIGUOUS"
    elif pending:
        end = inv[pending[0]].get("end")
        state = "HALTED" if end is not None else "UNTERMINATED"
    else:
        state = "CLEAN"
    last_outcome = None
    if latest is not None and "start" in inv[latest]:
        e = inv[latest].get("end")
        last_outcome = e["outcome"] if e else ("unterminated_acknowledged"
                                               if inv[latest].get("auth") else "unterminated")
    known = [v["end"]["elapsed_s"] for v in inv.values()
             if v.get("end") and isinstance(v["end"].get("elapsed_s"), (int, float))]
    return {"run_id": run_id, "state": state, "pending": pending[0] if pending else None,
            "problems": problems, "invocations": inv, "last_outcome": last_outcome,
            "next_invocation": (max(seen_m) + 1) if seen_m else 1,
            "cumulative_elapsed_s_known": float(sum(known)),
            "n_invocations_elapsed_unknown": sum(1 for v in inv.values() if not v.get("end"))}


def _all_authorization_refs(out: Path) -> set:
    d = Path(out) / HALTS_DIR
    refs = set()
    if d.exists():
        for f in d.glob("run*.inv*.auth.json"):
            try:
                refs.add(json.loads(f.read_text())["reference"])
            except (OSError, ValueError, KeyError):
                pass                                     # an unreadable auth is AMBIGUOUS anyway
    return refs


def _prefix_sha256(path: Path, n_bytes: int) -> Optional[str]:
    raw = Path(path).read_bytes()
    if len(raw) < n_bytes:
        return None
    return hashlib.sha256(raw[:n_bytes]).hexdigest()


def resolve_authorization(out: Path, run_id: int, state: dict, reference: Optional[str],
                          store: CheckpointStore, config_sha256: str, manifest_sha256: str):
    """Enforce the state model before a new invocation; record a single-use authorisation
    when one is pending and validly supplied. Returns the authorisation record or None."""
    if reference is not None:
        reference = str(reference).strip()
        if not reference:
            raise InvocationStateError("--resume-authorization reference is blank")
    st = state["state"]
    if st == "AMBIGUOUS":
        raise InvocationStateError(f"run {run_id}: invocation state AMBIGUOUS -- owner inspection "
                                   f"required: {state['problems']}")
    if st == "CLEAN":
        if reference is not None:
            raise InvocationStateError(f"run {run_id}: no runtime halt or unterminated invocation "
                                       "is pending; --resume-authorization refused")
        return None
    m = state["pending"]
    rec = state["invocations"][m]
    if reference is None:
        what = ("ended in a runtime-budget halt" if st == "HALTED"
                else "has no terminal record (interrupted, cause unknown)")
        raise InvocationStateError(f"run {run_id}: invocation {m} {what}; resuming requires "
                                   "--resume-authorization REF")
    if reference in _all_authorization_refs(out):
        raise InvocationStateError(f"authorisation reference {reference!r} was already used")
    start = rec["start"]
    if start.get("config_sha256") != config_sha256 or start.get("manifest_sha256") != manifest_sha256:
        raise InvocationStateError(f"run {run_id}: configuration or manifest differs from "
                                   f"invocation {m}; authorisation refused")
    cur = store.file_sha256()
    if st == "HALTED":
        if cur != rec["end"].get("checkpoint_sha256"):
            raise InvocationStateError(f"run {run_id}: checkpoint changed since the runtime halt "
                                       f"of invocation {m}; authorisation refused")
        kind = AUTH_FOR_HALT
    else:
        if _prefix_sha256(store.path, start["checkpoint_size_bytes"]) != start["checkpoint_sha256"]:
            raise InvocationStateError(f"run {run_id}: checkpoint content present when invocation "
                                       f"{m} started has changed; acknowledgement refused")
        kind = AUTH_FOR_UNTERMINATED
    auth = {"kind": "auth", "run_id": run_id, "invocation": m, "authorization_kind": kind,
            "reference": reference, "utc": datetime.now(timezone.utc).isoformat(),
            "checkpoint_sha256": cur, "checkpoint_n_records": len(store.done),
            "checkpoint_size_bytes": store.path.stat().st_size,
            "config_sha256": config_sha256, "manifest_sha256": manifest_sha256}
    durable_create(_inv_path(out, run_id, m, "auth"), auth)
    log_event(out, {"event": "resume_authorized", "run_id": run_id, "invocation": m,
                    "authorization_kind": kind, "reference": reference})
    return auth


def _write_end(out, run_id, m, outcome, monitor, store, prior_elapsed, details=None,
               stop_reasons=None, nonestimable=None):
    el = monitor.elapsed() if monitor.t0 is not None else None
    rec = {"kind": "end", "run_id": run_id, "invocation": m, "outcome": outcome,
           "utc": datetime.now(timezone.utc).isoformat(), "elapsed_s": el,
           "cumulative_elapsed_s_known": prior_elapsed + (el or 0.0),
           "n_records": len(store.done), "done_now": len(store.done) - monitor.n_start,
           "checkpoint_sha256": store.file_sha256(),
           "checkpoint_size_bytes": store.path.stat().st_size,
           "stop_reasons": list(stop_reasons or []),
           "nonestimable": (nonestimable if nonestimable is not None
                            else nonestimable_status(store.done.values())),
           "details": details}
    durable_create(_inv_path(out, run_id, m, "end"), rec)


# ======================================================================= Decision 11: helpers
NONESTIMABLE_RULE = "20*n_nonestimable > n_completed"
PRESERVATION_DIR = "preservation"
SIZING_DIR = "sizing"
WINDOWS_DIR = "windows"
DIAGNOSTICS_DIR = "diagnostics"
GATE_DIR = "gate"
APPROVED_GATE_ID = "V1-V5_H4_as_implemented_p9c-1"
GATE_CHECKS = ("H4", "historical_preservation", "V2_numpy", "V2_structural", "V2_domain",
               "V1_outer_historical")
HISTORICAL_B_OUT = 500                              # historical runs hold k = 0..499
NOT_ELIGIBLE_MARKER = "NOT ELIGIBLE — NOT A SCIENTIFIC RESULT"
CODE_IDENTITY_NOTICE = "CODE IDENTITY NOT VERIFIED — NO ANALYSIS PERFORMED"
_HEX64 = re.compile(r"^[0-9a-f]{64}$")


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _np_default(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o).__name__}")


def _canon(obj) -> bytes:
    """canonical_json (config.py) with numpy scalars converted; used for every hash below."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=True,
                      default=_np_default).encode()


def _jsonable(obj):
    return json.loads(_canon(obj))


def _nonblank(x) -> bool:
    return isinstance(x, str) and bool(x.strip())


def expected_v2_domain_states(run_ids, k_max: int, b_max: int, convention: str) -> int:
    """States rng.enumerate_stream_states adds over the declared domain: per run and k,
    outer (+ apparent_cv under the independent convention) and, per b, inner, inner_cv and
    inner_reserve; then the Phase 9 primary seeds."""
    per_k = (2 if convention == "independent" else 1) + 3 * int(b_max)
    return len(run_ids) * int(k_max) * per_k + len(PHASE9_PRIMARY_SEEDS)


def gate_evidence_problems(raw: Optional[bytes], config: CorrectedConfig, msha: str) -> list:
    """Decision 13: every reason why gate evidence bytes do NOT establish that all six
    approved pre-execution checks passed for this configuration, manifest and NumPy
    version. An empty list means the evidence is a complete, consistent six-check PASS.
    Never raises: any error while examining malformed evidence (wrong types, excessive
    nesting, ...) is itself reported as a problem (fail closed, structured)."""
    try:
        return _gate_evidence_problems(raw, config, msha)
    except Exception as ex:                          # incl. RecursionError, AttributeError
        return [f"gate evidence malformed: {type(ex).__name__}"]


def _gate_evidence_problems(raw: Optional[bytes], config: CorrectedConfig, msha: str) -> list:
    if config.get("verification_gate") != APPROVED_GATE_ID:
        return [f"configured verification_gate is not {APPROVED_GATE_ID}"]
    if raw is None:
        return ["gate evidence file missing or unreadable"]
    try:
        g = json.loads(raw)
    except ValueError as ex:
        return [f"gate evidence is not valid JSON ({ex.__class__.__name__})"]
    if not isinstance(g, dict):
        return ["gate evidence is not a JSON object"]
    bad = []
    if set(g) != {"key", "all_pass", *GATE_CHECKS}:
        bad.append(f"gate evidence does not hold exactly the approved checks: {sorted(g)}")
    key = {"config_sha256": config.sha256(), "manifest_sha256": msha, "numpy": np.__version__}
    if g.get("key") != key:
        bad.append("gate key does not match this configuration, manifest and NumPy version")
    for c in ("H4", "historical_preservation"):
        if g.get(c) != "PASS":
            bad.append(f"{c} is not PASS")
    for c in ("V2_numpy", "V2_structural", "V2_domain", "V1_outer_historical"):
        if not isinstance(g.get(c), dict) or g[c].get("status") != "PASS":
            bad.append(f"{c} is not PASS")
    if bad:
        return bad
    nd, sd = g["V2_numpy"].get("detail"), g["V2_structural"].get("detail")
    if not isinstance(nd, dict) or nd.get("all_pass") is not True:
        bad.append("V2_numpy detail does not pass")
    if not isinstance(sd, dict) or sd.get("holds") is not True:
        bad.append("V2_structural detail does not hold")
    dom, v1 = g["V2_domain"], g["V1_outer_historical"]
    rungs = (config.get("ladder") or {}).get("rungs") or []
    k_max = max([config.B_out_initial] + list(rungs))
    exp_domain = {"run_ids": list(config.run_ids), "k_max": k_max, "b_max": config.B_in_max,
                  "convention": config.convention, "T": config.T}
    if dom.get("domain") != exp_domain:
        bad.append(f"V2_domain domain {dom.get('domain')} differs from the approved {exp_domain}")
    try:
        n_exp = expected_v2_domain_states(config.run_ids, k_max, config.B_in_max, config.convention)
    except (TypeError, ValueError):
        n_exp = None
    if dom.get("pass") is not True or dom.get("collisions") != [] or n_exp is None \
            or dom.get("n_states") != n_exp or dom.get("n_distinct") != n_exp:
        bad.append(f"V2_domain counts inconsistent (expected {n_exp} distinct states, no collisions)")
    if v1.get("n_compared") != HISTORICAL_B_OUT * len(config.run_ids) or v1.get("mismatches") != []:
        bad.append("V1_outer_historical does not cover every historical replicate of every run")
    if g.get("all_pass") is not True:
        bad.append("all_pass is not true")
    return bad


def _gate_attempt_path(out: Path, run_id: int, m: int, a: int) -> Path:
    return Path(out) / GATE_DIR / f"run{run_id}.inv{m:04d}.attempt{a:02d}.json"


_GATE_ATTEMPT_NAME = re.compile(r"^gate/run(\d+)\.inv(\d{4})\.attempt(\d{2})\.json$")


def start_gate_problems(out: Path, start: dict, config: CorrectedConfig, msha: str) -> tuple:
    """Re-verify the gate evidence bound into a start record. Returns (category, problems):
    category None when valid, "blocking" when the evidence cannot be found or read,
    "disqualifier" when it is altered (SHA-256 mismatch), belongs to another invocation or
    does not establish a PASS. Never raises."""
    try:
        return _start_gate_problems(out, start, config, msha)
    except Exception as ex:
        return "disqualifier", [f"gate evidence reference malformed: {type(ex).__name__}"]


def _start_gate_problems(out: Path, start: dict, config: CorrectedConfig, msha: str) -> tuple:
    ge = start.get("gate_evidence")
    if not isinstance(ge, dict) or not isinstance(ge.get("path"), str) or not ge.get("sha256"):
        return "blocking", ["start record holds no gate evidence reference"]
    rel = ge["path"]
    mm = _GATE_ATTEMPT_NAME.match(rel)
    if not mm:
        return "disqualifier", [f"gate evidence path not permitted: {rel!r}"]
    if (int(mm.group(1)), int(mm.group(2))) != (start.get("run_id"), start.get("invocation")) \
            or int(mm.group(3)) != ge.get("attempt"):
        return "disqualifier", [f"gate evidence {rel} does not belong to run {start.get('run_id')} "
                                f"invocation {start.get('invocation')} attempt {ge.get('attempt')}"]
    try:
        raw = (Path(out) / rel).read_bytes()
    except (OSError, ValueError):
        return "blocking", [f"gate evidence {rel} missing or unreadable"]
    if _sha(raw) != ge["sha256"]:
        return "disqualifier", [f"gate evidence {rel} differs from the SHA-256 in the start record"]
    bad = gate_evidence_problems(raw, config, msha)
    if ge.get("all_pass") is not True:
        bad = ["start record does not record a gate PASS"] + bad
    return ("disqualifier", bad) if bad else (None, [])


def run_fresh_gate(config: CorrectedConfig, out: Path, run_id: int, m: int, msha: str) -> dict:
    """Execute the complete gate for this invocation attempt into a NEW evidence file (so
    gate.py cannot reuse any cache), then establish the PASS from the written evidence.
    A failed attempt keeps its evidence and raises SystemExit; nothing is dispatched."""
    a = 1
    while _gate_attempt_path(out, run_id, m, a).exists():
        a += 1
    p = _gate_attempt_path(out, run_id, m, a)
    p.parent.mkdir(parents=True, exist_ok=True)
    run_pre_execution_gate(config, msha, out_path=p)
    try:
        raw = p.read_bytes()
    except OSError:
        raw = None
    bad = gate_evidence_problems(raw, config, msha)
    rel = p.relative_to(out).as_posix()
    if bad:
        log_event(out, {"event": "gate_failed", "run_id": run_id, "invocation": m,
                        "evidence": rel, "problems": bad[:10]})
        raise SystemExit(f"pre-execution gate FAILED for run {run_id} invocation {m} -- nothing "
                         f"dispatched; evidence preserved at {rel}: {bad[:3]}")
    return {"path": rel, "sha256": _sha(raw), "all_pass": True, "attempt": a}


def nonestimable_status(records_in_order) -> dict:
    """Decision 11A snapshot over records in append order (== ascending k)."""
    recs = list(records_in_order)
    bad = [r["k"] for r in recs if r.get("status") != "ok"]
    by = {}
    for r in recs:
        if r.get("status") != "ok":
            by[str(r.get("status"))] = by.get(str(r.get("status")), 0) + 1
    n, b = len(recs), len(bad)
    return {"rule": NONESTIMABLE_RULE, "n_nonestimable": b, "n_completed": n,
            "lhs": 20 * b, "rhs": n, "result": 20 * b > n,
            "triggering_k": recs[-1]["k"] if recs else None, "triggering_position": n,
            "k_nonestimable": bad, "by_status": dict(sorted(by.items()))}


def nonestimable_prefix_violation(records_in_order) -> Optional[dict]:
    """First prefix n with 20*b_n > n, or None. Permanent under Option A."""
    b = 0
    for i, r in enumerate(records_in_order, 1):
        if r.get("status") != "ok":
            b += 1
        if 20 * b > i:
            return {"position": i, "k": r["k"], "n_nonestimable": b}
    return None


# ----------------------------------------------------------------- preservation evidence
def preservation_evidence(res, stage: str) -> dict:
    """Structured evidence from the ACTUAL historical_preservation_check() result."""
    files = res.get("files") if isinstance(res, dict) else None
    well_formed = isinstance(files, dict) and bool(files) and all(
        isinstance(v, dict) and {"expected_prefix", "sha256", "ok"} <= set(v) for v in files.values())
    if not well_formed:
        return {"check": "historical_preservation", "stage": stage, "malformed": True,
                "raw": repr(res)[:2000]}
    rows = sorted(({"path": rel, "expected_prefix": v["expected_prefix"],
                    "observed_sha256": v["sha256"], "exists": v["sha256"] is not None,
                    "ok": bool(v["ok"])} for rel, v in files.items()), key=lambda r: r["path"])
    return {"check": "historical_preservation", "stage": stage, "malformed": False,
            "registry_sha256": _sha(_canon(PROTECTED_HISTORICAL)),
            "n_checked": len(rows), "n_failing": sum(not r["ok"] for r in rows), "files": rows}


def _evidence_malformed(ev) -> bool:
    return not (isinstance(ev, dict) and ev.get("malformed") is False
                and isinstance(ev.get("files"), list) and ev.get("n_failing", 0) >= 1)


def _record_preservation_failure(out: Path, evidence: dict):
    """Decision 11 item 17: durable, global record of a pre-dispatch preservation failure."""
    d = Path(out) / PRESERVATION_DIR
    seqs = [int(m.group(1)) for f in (d.iterdir() if d.exists() else [])
            for m in [re.match(r"^failure\.(\d{4})\.json(\.tmp)?$", f.name)] if m]
    seq = (max(seqs) + 1) if seqs else 1
    durable_create(d / f"failure.{seq:04d}.json",
                   {"kind": "preservation_failure", "stage": "pre_dispatch", "sequence": seq,
                    "utc": datetime.now(timezone.utc).isoformat(), "evidence": evidence})


def preservation_failures(out: Path, run_ids) -> list:
    """Every recorded historical-preservation failure, across ALL runs (global block)."""
    out = Path(out)
    found = []
    d = out / PRESERVATION_DIR
    if d.exists():
        for f in sorted(d.iterdir()):
            entry = {"source": f"{PRESERVATION_DIR}/{f.name}", "sha256": None, "malformed": True}
            if not f.name.endswith(".tmp"):
                try:
                    raw = f.read_bytes()
                    entry["sha256"] = _sha(raw)
                    entry["malformed"] = _evidence_malformed(json.loads(raw).get("evidence"))
                except (OSError, ValueError, AttributeError):
                    pass
            found.append(entry)
    for r in run_ids:
        for f in sorted((out / HALTS_DIR).glob(f"run{r}.inv*.end.json")):
            try:
                raw = f.read_bytes()
                e = json.loads(raw)
            except (OSError, ValueError):
                continue                                # unreadable end -> AMBIGUOUS elsewhere
            if e.get("outcome") == "historical_preservation_failure":
                found.append({"source": f"{HALTS_DIR}/{f.name}", "sha256": _sha(raw),
                              "malformed": _evidence_malformed(e.get("details"))})
    return found


def _authorization_reference_counts(out: Path) -> dict:
    counts = {}
    for f in sorted((Path(out) / HALTS_DIR).glob("run*.inv*.auth.json")):
        try:
            ref = json.loads(f.read_text())["reference"]
        except (OSError, ValueError, KeyError):
            continue
        counts[ref] = counts.get(ref, 0) + 1
    return counts


# ----------------------------------------------------------------- checkpoint parsing
def _parse_checkpoint(raw: bytes):
    """-> (offsets, header, records, problems). offsets[n] = byte length of header + n
    records. Lines are canonical JSON without embedded newlines (config.canonical_json)."""
    lines = raw.split(b"\n")[:-1]
    offsets, pos = [], 0
    for ln in lines:
        pos += len(ln) + 1
        offsets.append(pos)
    problems, records, header = [], [], None
    try:
        header = json.loads(lines[0]) if lines else None
    except ValueError:
        problems.append("header line is not valid JSON")
    for i, ln in enumerate(lines[1:], 1):
        try:
            o = json.loads(ln)
        except ValueError:
            problems.append(f"line {i + 1}: not valid JSON")
            records.append(None)
            continue
        r = o.get("record") if isinstance(o, dict) else None
        if (o.get("kind") != "record" or not isinstance(r, dict)
                or _sha(canonical_json(r)) != o.get("record_sha256")
                or o.get("k") != r.get("k") or o.get("run_id") != r.get("run_id")):
            problems.append(f"line {i + 1}: record envelope or hash invalid")
        records.append(r)
    return offsets, header, records, problems


class _Reasons:
    ORDER = ("disqualifier", "blocking", "pending_verification", "incomplete", "disclosure")

    def __init__(self):
        self.items = []

    def add(self, category, code, **detail):
        assert category in self.ORDER
        self.items.append({"category": category, "code": code, "detail": _jsonable(detail)})

    def has(self, category):
        return any(i["category"] == category for i in self.items)

    def codes(self, category):
        return [i["code"] for i in self.items if i["category"] == category]


def _sizing_path(out: Path, B_cur: int) -> Path:
    return Path(out) / SIZING_DIR / f"level{int(B_cur):05d}.json"


def _assess_windows(R: "_Reasons", run_id: int, windows: list, wres: dict):
    """DEC-C1 / D11-D1. Every acknowledged unattested window must have a BOUND, COMPLETE
    PASS result before scientific eligibility. No result, an unbound result or an incomplete
    result is `pending_verification` (never PASS); a bound FAIL (or a PASS contradicted by
    its own per-k evidence) is a permanent disqualifier that blocks ALL dispatch."""
    for w in windows:
        r = wres.get(w["invocation"])
        bound = (isinstance(r, dict) and r.get("run_id") == run_id
                 and r.get("invocation") == w["invocation"]
                 and r.get("region_sha256") == w["region_sha256"] and r.get("k_list") == w["k_list"])
        if r is None:
            R.add("pending_verification", "window_not_recomputed", invocation=w["invocation"])
        elif not bound:
            R.add("pending_verification", "window_result_not_bound", invocation=w["invocation"])
        else:
            per_k = r.get("per_k")
            per_k = per_k if isinstance(per_k, dict) else {}       # malformed -> incomplete
            vals = list(per_k.values())
            # contradiction = an entry that explicitly reports a difference; missing or
            # malformed entries are incomplete evidence (pending), never a disqualifier
            contradicted = any(isinstance(v, dict) and (v.get("record_identical") is False
                                                        or v.get("dC_hex_equal") is False)
                               for v in vals)
            all_ok = (set(per_k) == {str(k) for k in w["k_list"]}
                      and all(isinstance(v, dict) and v.get("record_identical") is True
                              and v.get("dC_hex_equal") is True for v in vals))
            if r.get("status") == "FAIL" or (contradicted and r.get("status") == "PASS"):
                R.add("disqualifier", "window_recomputation_failed", invocation=w["invocation"])
            elif r.get("status") == "PASS" and all_ok and r.get("code_identity_verified") is True:
                R.add("disclosure", "window_recomputed_bit_identical", invocation=w["invocation"],
                      n_records=len(w["k_list"]))
            else:
                R.add("pending_verification", "window_result_incomplete", invocation=w["invocation"])


def assess_run(config: CorrectedConfig, run_id: int, manifest_sha256: str,
               max_invocation: Optional[int] = None,
               checkpoint_raw: Optional[bytes] = None,
               window_results: Optional[dict] = None) -> dict:
    """Read-only history assessment shared by the pre-dispatch refusal and by Stage B.
    Decision 11 items 5, 6, 8, 9, 10, the 5 % prefix scan and (D11-D1) window results:
    a window FAIL disqualifies before any dispatch; pending windows never block
    continuation (DEC-C1)."""
    out, R = Path(config.output_dir), _Reasons()
    st = invocation_state(out, run_id, max_invocation=max_invocation)
    invs = st["invocations"]
    ms = sorted(m for m in invs if "start" in invs[m])
    info = {"state": st["state"], "pending": st["pending"], "invocations": ms, "windows": [],
            "levels": {}, "records": [], "checkpoint_sha256": None, "checkpoint_size_bytes": None,
            "window_results": {}}
    if st["state"] == "AMBIGUOUS":
        R.add("blocking", "invocation_state_ambiguous", problems=st["problems"])
    for p in preservation_failures(out, config.run_ids):
        R.add("blocking", "preservation_failure_pending", source=p["source"])
        if p["malformed"]:
            R.add("blocking", "preservation_evidence_malformed", source=p["source"])
    for ref, c in sorted(_authorization_reference_counts(out).items()):
        if c > 1:
            R.add("disqualifier", "duplicate_authorization_reference", reference=ref, count=c)
    for m in ms:
        e = invs[m].get("end") or {}
        if e.get("outcome") == "hard_failure_h4":
            R.add("disqualifier", "hard_failure_h4_during_execution", invocation=m)
        if e.get("outcome") == "hard_failure_h6":
            R.add("disqualifier", "hard_failure_h6", invocation=m)
        if e.get("outcome") == "stop_nonestimable" or "nonestimable_rate" in (e.get("stop_reasons") or []):
            R.add("disqualifier", "history_nonestimable_stop", invocation=m)
    cs, ms_ = config.sha256(), manifest_sha256
    for m in ms:                                     # Decision 13: gate evidence per invocation
        cat, gp = start_gate_problems(out, invs[m]["start"], config, ms_)
        if cat:
            R.add(cat, "gate_evidence_unverifiable" if cat == "blocking" else "gate_evidence_invalid",
                  invocation=m, problems=gp[:5])
    for m in ms:
        for kind in ("start", "auth"):
            rec = invs[m].get(kind)
            if rec and (rec.get("config_sha256") != cs or rec.get("manifest_sha256") != ms_):
                R.add("disqualifier", "identity_mismatch", invocation=m, record=kind)
    ck = out / "runs" / f"run{run_id}.jsonl"
    raw = checkpoint_raw
    if raw is None and ck.exists():
        raw = ck.read_bytes()
    if raw is None:
        if ms:
            R.add("blocking", "checkpoint_missing")
        else:
            R.add("incomplete", "not_started")
        return {"reasons": R, **info}
    info["checkpoint_sha256"], info["checkpoint_size_bytes"] = _sha(raw), len(raw)
    if raw and not raw.endswith(b"\n"):
        R.add("blocking", "partial_tail", bytes_after_last_newline=len(raw) - (raw.rfind(b"\n") + 1))
        return {"reasons": R, **info}
    offsets, header, recs, problems = _parse_checkpoint(raw)
    for pr in problems:
        R.add("disqualifier", "checkpoint_region_invalid", problem=pr)
    exp_header = {"kind": "header", "schema": SCHEMA_VERSION, "run_id": run_id,
                  "config_sha256": cs, "manifest_sha256": ms_}
    if header != exp_header:
        R.add("disqualifier", "checkpoint_header_mismatch")
    if problems or header != exp_header:
        return {"reasons": R, **info}
    ks = [r["k"] for r in recs]
    if len(set(ks)) != len(ks):
        R.add("disqualifier", "duplicate_k")
    info["records"] = recs
    v = nonestimable_prefix_violation(recs)
    if v:
        R.add("disqualifier", "nonestimable_prefix_violation", **v)
    if not ms:
        if recs:
            R.add("disqualifier", "orphan_records", n_records=len(recs))
        else:
            R.add("incomplete", "not_started")
        return {"reasons": R, **info}

    def boundary(n):
        return offsets[n] if 0 <= n < len(offsets) else None

    def chain_break(m, what, **d):
        R.add("disqualifier", "chain_break", invocation=m, what=what, **d)

    prev = None                                     # (size, sha, n) of the last sealed boundary
    L_cur, completed, refs = None, set(), {}
    for idx, m in enumerate(ms):
        v_ = invs[m]
        s, end, auth = v_["start"], v_.get("end"), v_.get("auth")
        n0, size0, sha0 = s.get("n_start"), s.get("checkpoint_size_bytes"), s.get("checkpoint_sha256")
        if idx == 0 and (n0 != 0 or size0 != boundary(0)):
            chain_break(m, "first invocation did not start from the header-only checkpoint")
        if idx > 0 and (prev is None or (size0, sha0) != prev[:2]):
            chain_break(m, "start boundary != previous sealed boundary")
        if boundary(n0 if isinstance(n0, int) else -1) != size0 or size0 is None \
                or _sha(raw[:size0]) != sha0:
            chain_break(m, "start boundary not reproducible from the checkpoint")
            prev = None
            break
        close = end or auth
        if end and auth and end.get("checkpoint_sha256") != auth.get("checkpoint_sha256"):
            chain_break(m, "end and authorisation boundaries differ")
        if close is not None:
            n1 = end["n_records"] if end else close.get("checkpoint_n_records")
            size1 = close.get("checkpoint_size_bytes") or boundary(n1 if isinstance(n1, int) else -1)
            if not isinstance(n1, int) or size1 is None or boundary(n1) != size1 \
                    or _sha(raw[:size1]) != close.get("checkpoint_sha256"):
                chain_break(m, "closing boundary not reproducible from the checkpoint")
                prev = None
                break
        else:
            n1, size1 = len(recs), len(raw)          # open region (latest, unterminated)
        region = recs[n0:n1]
        expect_n = end["done_now"] if end else n1 - n0
        target = s.get("target_B_out")
        done_before = set(ks[:n0])
        todo = [k for k in range(target if isinstance(target, int) else 0) if k not in done_before]
        if expect_n != len(region) or [r["k"] for r in region] != todo[:len(region)] \
                or any(r.get("run_id") != run_id for r in region):
            R.add("disqualifier", "region_mismatch", invocation=m)
        if end is None and auth is not None:
            info["windows"].append({"invocation": m, "region_start": size0, "region_end": size1,
                                    "region_sha256": _sha(raw[size0:size1]),
                                    "k_list": [r["k"] for r in region]})
        prev = (size1, _sha(raw[:size1]), n1) if close is not None else None
        # ---- target progression (Paths A / B / C), retrospectively
        ref = s.get("topup_authorization")
        ref = ref.strip() if isinstance(ref, str) else None
        ok = True
        if L_cur is None:
            ok = target == config.B_out_initial and n0 == 0 and not ref
        elif target == L_cur:
            ok = (ref == refs.get(L_cur)) if L_cur > config.B_out_initial else not ref
        elif isinstance(target, int) and target > L_cur:
            rungs = (config.get("ladder") or {}).get("rungs") or []
            ok = (L_cur in completed and n0 == L_cur and target in rungs and bool(ref)
                  and ref not in refs.values())
            if ok:
                sp, ssha = _sizing_path(out, L_cur), s.get("sizing_decision_sha256")
                try:
                    sraw = sp.read_bytes()
                    art = json.loads(sraw)
                    rr = art["runs"][str(run_id)]
                    ok = (_sha(sraw) == ssha and art.get("B_cur") == L_cur
                          and art["decision"]["combined_B_next"] == target
                          and rr["checkpoint_size_bytes"] == size0
                          and rr["checkpoint_sha256"] == sha0)
                except (OSError, ValueError, KeyError, TypeError):
                    ok = False
        else:
            ok = False
        if not ok:
            R.add("disqualifier", "target_progression_violation", invocation=m, target=target)
        if target != L_cur:
            refs[target] = ref
            L_cur = target
        if end and end.get("outcome") == "completed" and end.get("n_records") == target:
            completed.add(target)
    else:
        if prev is not None and prev[0] != len(raw):
            chain_break(ms[-1], "bytes appended after the last sealed boundary")
    info["levels"] = {"L_cur": L_cur, "completed": sorted(completed),
                      "refs": {str(k): v for k, v in refs.items()}}
    wres = (window_results if window_results is not None
            else load_window_results(out, run_id, [w["invocation"] for w in info["windows"]]))
    info["window_results"] = wres
    _assess_windows(R, run_id, info["windows"], wres)
    if st["state"] == "HALTED":
        R.add("incomplete", "pending_runtime_halt", invocation=st["pending"])
    elif st["state"] == "UNTERMINATED":
        R.add("incomplete", "pending_unterminated", invocation=st["pending"])
    return {"reasons": R, **info}


def _dispatch(R: _Reasons, state: str, complete: bool) -> str:
    if R.has("disqualifier"):
        return "blocked_permanently"
    if R.has("blocking"):
        return "blocked_pending_owner"
    if complete:
        return "not_applicable"
    return ("permitted_with_authorization" if state in ("HALTED", "UNTERMINATED")
            else "permitted_without_authorization")


# ----------------------------------------------------------------- Stage A results (inputs)
def window_result_path(out: Path, run_id: int, m: int) -> Path:
    return Path(out) / WINDOWS_DIR / f"run{run_id}.inv{m:04d}.json"


def load_window_results(out: Path, run_id: int, invocations) -> dict:
    res = {}
    for m in invocations:
        p = window_result_path(out, run_id, m)
        if p.exists():
            try:
                res[m] = json.loads(p.read_text())
            except (OSError, ValueError):
                res[m] = {"status": "UNREADABLE"}
    return res


# ----------------------------------------------------------------- Stage B (pure)
def evaluate_eligibility(config: CorrectedConfig, run_id: int, approved_target: int,
                         manifest_sha256: str, max_invocation: Optional[int] = None,
                         checkpoint_raw: Optional[bytes] = None,
                         window_results: Optional[dict] = None,
                         freeze_record_sha256: Optional[str] = None) -> dict:
    """Decision 11D Stage B: pure, deterministic, read-only. Performs no fitting.
    Classes: INELIGIBLE > UNVERIFIABLE > INCOMPLETE > ELIGIBLE (all reasons retained).

    Decision 11 item 4 as amended by DEC-C1: "UNVERIFIABLE implies dispatch
    blocked_pending_owner, except when every UNVERIFIABLE reason is a pending window
    verification (no result, an unbound result or an incomplete result). In that case
    continuation of the already-authorised level is permitted under the Decision 10 rules,
    while scientific eligibility, any new sizing-driven top-up and reporting remain blocked
    until every required window has a bound PASS result. A window FAIL is a permanent
    disqualifier and blocks all dispatch." Pending windows are category
    `pending_verification`, which is never blocking; the `dispatch` value therefore follows
    the same predicates as run_corrected's pre-dispatch refusal."""
    a = assess_run(config, run_id, manifest_sha256, max_invocation, checkpoint_raw,
                   window_results=window_results)
    R, out = a["reasons"], Path(config.output_dir)
    st = invocation_state(out, run_id, max_invocation=max_invocation)
    invs, ms = st["invocations"], a["invocations"]
    complete = False
    if ms and not R.has("blocking"):
        last = ms[-1]
        end, s = invs[last].get("end"), invs[last]["start"]
        if not end or end.get("outcome") != "completed":
            if not R.has("incomplete"):
                R.add("incomplete", "latest_not_completed", invocation=last)
        else:
            t = s.get("target_B_out")
            if t > approved_target:
                R.add("disqualifier", "target_exceeds_approved", target=t, approved=approved_target)
            elif t < approved_target:
                R.add("incomplete", "target_not_reached", target=t, approved=approved_target)
            elif [r["k"] for r in a["records"]] != list(range(approved_target)):
                R.add("disqualifier", "k_positions_incorrect")
            else:
                complete = True
    wres = a["window_results"]
    if R.has("disqualifier"):
        cls = "INELIGIBLE"
    elif R.has("blocking") or (R.has("pending_verification") and not R.has("incomplete")):
        cls = "UNVERIFIABLE"
    elif R.has("incomplete"):
        cls = "INCOMPLETE"
    else:
        cls = "ELIGIBLE"
    dispatch = _dispatch(R, st["state"], complete)
    halts_files = []
    for f in sorted((out / HALTS_DIR).glob(f"run{run_id}.inv*")):
        mt = _REC_NAME.match(f.name)
        if mt and (max_invocation is None or int(mt.group(2)) <= max_invocation):
            halts_files.append([f.name, _sha(f.read_bytes())])
    evidence = {
        "checkpoint_sha256": a["checkpoint_sha256"], "checkpoint_size_bytes": a["checkpoint_size_bytes"],
        "n_records": len(a["records"]), "invocation_records": halts_files,
        "window_results": sorted([str(m), _sha(_canon(r))] for m, r in wres.items()),
        "preservation_records": sorted([p["source"], p["sha256"]] for p in
                                       preservation_failures(out, config.run_ids)),
        "config_sha256": config.sha256(), "manifest_sha256": manifest_sha256,
        "freeze_record_sha256": freeze_record_sha256}
    key = lambda i: (i["code"], _canon(i["detail"]))
    obj = {"schema": "p9c-eligibility-1", "run_id": run_id, "evaluated_target": approved_target,
           "class": cls, "permanent": cls == "INELIGIBLE", "dispatch": dispatch,
           "reasons": sorted((i for i in R.items if i["category"] != "disclosure"), key=key),
           "disclosures": sorted((i for i in R.items if i["category"] == "disclosure"), key=key),
           "evidence": evidence}
    obj = _jsonable(obj)
    return {**obj, "windows": a["windows"], "eligibility_result_sha256": _sha(_canon(obj))}


def diagnostic_report(result: dict) -> tuple:
    """Item 14: accounting-only diagnostic for a run that is not ELIGIBLE. Returns
    (relative filename, text); writing it is the runner's job. Contains no interval,
    EP-1, W5 or H1/H2 values."""
    if result.get("class") == "ELIGIBLE":
        raise ValueError("diagnostic reports are only for runs that are not ELIGIBLE")
    name = (f"{DIAGNOSTICS_DIR}/NOT_ELIGIBLE__NOT_A_SCIENTIFIC_RESULT__run{result['run_id']}"
            f"__{result['class']}.json")
    body = {"marker": NOT_ELIGIBLE_MARKER, "run_id": result["run_id"], "class": result["class"],
            "permanent": result["permanent"], "dispatch": result["dispatch"],
            "reasons": result["reasons"], "disclosures": result["disclosures"],
            "evidence": result["evidence"]}
    return name, NOT_ELIGIBLE_MARKER + "\n" + json.dumps(body, indent=1, sort_keys=True) + "\n"


# ----------------------------------------------------------------- Decision 9 / 13 code identity
FREEZE_RECORD_RELPATH = "config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json"
FREEZE_SCHEMA = "p9c-freeze-2"          # Decision 13: + executable, launch_flags, record-only
FROZEN_PIPELINE_MODULES = ("src/phase9_execution.py", "src/pathway_scoring.py",
                           "src/score_standardisation.py")
DRIVER_RELPATH = "scripts/41_phase9_corrected_driver.py"
LAUNCHER_RELPATH = "scripts/42_phase9_corrected_launch.sh"       # DEC-E1: part of the freeze
VERIFY_RELPATH = "scripts/43_phase9_corrected_verify.py"         # Decision 13 P1: post-run runner
FREEZE_PREPARE_RELPATH = "scripts/44_phase9_freeze_prepare.sh"   # Decision 13 P2 (provenance)
FREEZE_CAPTURE_RELPATH = "scripts/44_phase9_freeze_capture.py"   # Decision 13 P2 (provenance)
HISTORICAL_DRIVER_RELPATH = "scripts/40_phase9_dc_interval.py"   # loaded by V1-inner
HISTORICAL_DRIVER_ALIAS = "p9_dc_driver_frozen_ro"               # its module name when loaded
PACKAGE_DIR = "src/phase9_corrected"
# Decision 13: every repository file the approved execution and verification paths need,
# besides the package directory (whose exact file set is checked separately). Mandatory:
# each must be in the freeze inventory and match its hash, whether or not it is imported.
REQUIRED_FILES = (DRIVER_RELPATH, LAUNCHER_RELPATH, VERIFY_RELPATH, FREEZE_PREPARE_RELPATH,
                  FREEZE_CAPTURE_RELPATH, HISTORICAL_DRIVER_RELPATH, "src/__init__.py",
                  "src/paths.py", "src/reproducibility.py", *FROZEN_PIPELINE_MODULES)
THREAD_ENV_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                   "MKL_NUM_THREADS")
# Decision 13 P3: BLOCKING environment fields (compared with the freeze at every launch).
# Record-only fields (macOS build, jq version) live in the freeze's `environment_record_only`
# and are never compared here. Neither set is complete environment identity: shared
# libraries loaded by the interpreter (libpython, BLAS/LAPACK binaries) are not hashed.
ENV_KEYS = ("python_implementation", "python_version", "sys_platform", "machine", "platform",
            "prefix", "distributions", "numpy_build", "thread_env", "executable", "launch_flags")
# DEC-E1: the Decision 13 freeze must use "fresh_compile_prefix" (launch with
# -I -B -X pycache_prefix=<fresh empty external directory>, so in-tree caches are never
# read). The earlier policies remain recognised; the E1 runtime launch flags below are
# required under EVERY policy.
BYTECODE_POLICIES = ("no_pyc", "checked_hash_pyc", "fresh_compile_prefix")
E1_BYTECODE_POLICY = "fresh_compile_prefix"


def current_environment() -> dict:
    dists = {}
    for d in importlib.metadata.distributions():
        name = (d.metadata["Name"] or "").lower().replace("_", "-")
        dists.setdefault(name, set()).add(d.version)
    nb = np.show_config(mode="dicts")["Build Dependencies"]
    return {"python_implementation": sys.implementation.name,
            "python_version": ".".join(map(str, sys.version_info[:3])),
            "sys_platform": sys.platform, "machine": platform.machine(),
            # Decision 13 P3: OS family and architecture only. platform.platform() embeds the
            # macOS point version, which is RECORD-ONLY (environment_record_only), not blocking.
            "platform": {"system": platform.system(), "machine": platform.machine(),
                         "architecture": platform.architecture()[0]},
            "prefix": sys.prefix,
            "distributions": {k: sorted(v) for k, v in sorted(dists.items())},
            "numpy_build": {lib: {f: nb.get(lib, {}).get(f) for f in ("name", "version", "detection method")}
                            for lib in ("blas", "lapack")},
            "thread_env": {k: os.environ.get(k) for k in THREAD_ENV_VARS},
            "executable": _executable_identity(),
            "launch_flags": _launch_flag_identity()}


def _executable_identity() -> dict:
    """Decision 13 P3: resolved interpreter path and the SHA-256 of that file."""
    p = Path(sys.executable).resolve()
    try:
        h = hashlib.sha256(p.read_bytes()).hexdigest()
    except OSError:
        h = None
    return {"resolved_path": str(p), "sha256": h}


def _launch_flag_identity() -> dict:
    """Decision 13 P3: the effective launch flags (the prefix's random path is excluded)."""
    fl = _runtime_flags()
    return {"isolated": bool(fl.get("isolated")),
            "dont_write_bytecode": bool(fl.get("dont_write_bytecode")),
            "pycache_prefix_set": bool(fl.get("pycache_prefix"))}


def _dont_write_bytecode() -> bool:
    return bool(sys.flags.dont_write_bytecode)


def _runtime_flags() -> dict:
    """DEC-E1 launch evidence (stubbable in tests): -I, -B and -X pycache_prefix."""
    return {"isolated": bool(sys.flags.isolated),
            "dont_write_bytecode": bool(sys.flags.dont_write_bytecode),
            "pycache_prefix": sys.pycache_prefix}


def _pyc_flags(pyc: Path) -> int:
    return int.from_bytes(Path(pyc).read_bytes()[4:8], "little")


def _module_name_for(rel: str) -> Optional[str]:
    """Import name under which a REQUIRED repository file is loaded (None: not importable
    by name, e.g. a script run as __main__, which is matched by path instead)."""
    if rel == HISTORICAL_DRIVER_RELPATH:
        return HISTORICAL_DRIVER_ALIAS
    if rel.startswith("src/") and rel.endswith(".py"):
        parts = rel[:-3].split("/")
        return ".".join(parts[:-1] if parts[-1] == "__init__" else parts)
    return None


def verify_loaded_modules(root: Path, code: dict, required) -> None:
    """Decision 13 B2: in-process identity of the repository modules ACTUALLY LOADED by this
    process (it supplements, never replaces, the launcher's pre-Python file checks).
      * a loaded module whose source resolves inside the repository (outside the virtual
        environment and any site-packages) must be a required file and match its frozen
        hash -- a repository-local dependency outside the inventory fails closed;
      * a module loaded under the import name of a required file (`src.*`, and the
        historical driver alias `p9_dc_driver_frozen_ro`), wherever it was loaded from,
        must match that file's frozen hash;
      * any loaded `src.phase9_corrected.*` module must be in the freeze.
    Standard-library and site-package modules are not repository sources. Required files
    that were not imported are not demanded here (their on-disk hashes are checked)."""
    def refuse(msg):
        raise CodeIdentityRefusal(msg)
    root_res = Path(root).resolve()
    required = set(required)
    by_name = {n: rel for rel in required if (n := _module_name_for(rel))}
    env_roots = []
    for pfx in {sys.prefix, sys.base_prefix, sys.exec_prefix, sys.base_exec_prefix}:
        try:
            env_roots.append(Path(pfx).resolve())
        except OSError:
            pass

    def inside(p: Path, d: Path) -> bool:
        return p == d or d in p.parents
    for name, mod in sorted(list(sys.modules.items()), key=lambda x: x[0]):
        f = getattr(mod, "__file__", None)
        if name == HISTORICAL_DRIVER_ALIAS and not f:
            refuse(f"loaded module {name} has no source file")
        if not f:
            continue
        try:
            p = Path(f).resolve()
        except OSError:
            refuse(f"loaded module {name} ({f}) cannot be resolved")
        local = (inside(p, root_res) and "site-packages" not in p.parts
                 and not any(inside(p, e) for e in env_roots))
        rel_name = by_name.get(name)
        if name.startswith(PACKAGE_DIR.replace("/", ".") + "."):
            rel_name = f"{PACKAGE_DIR}/{p.name}"
        if local:
            rel = p.relative_to(root_res).as_posix()
            if rel not in required:
                refuse(f"loaded repository-local module {name} ({rel}) is not a required "
                       "inventory file")
            if rel_name is not None and rel_name != rel:
                refuse(f"loaded module {name} resolves to {rel}, not {rel_name}")
        elif rel_name is not None:
            rel = rel_name
        else:
            continue
        if rel not in code or p.suffix != ".py":
            refuse(f"loaded module {name} ({f}) does not match the freeze")
        try:
            h = _sha(p.read_bytes())
        except OSError:
            refuse(f"loaded module {name} ({f}) cannot be read")
        if h != code[rel]:
            refuse(f"loaded module {name} ({f}) does not match the freeze")


def verify_code_identity(config: CorrectedConfig, manifest: Optional[dict] = None,
                         extra_files=()) -> dict:
    """Decision 9 condition 4, enforced on every path (and by the future runner, which
    passes its own files in `extra_files`). Fails closed with CodeIdentityRefusal.

    DEC-B2 verification order (owner_decisions_outside_schema.D13_driver_freeze):
      1. approved / date / reference;
      2. `freeze_record_sha256` (the anchor: SHA-256 of the complete freeze-record bytes)
         and `sha256` (SHA-256 of the frozen driver at `driver_path`) are both 64-char
         lowercase hex;
      3. sha256(freeze-record bytes) == freeze_record_sha256;
      4. the freeze record parses and has the required schema; it never contains a
         configuration file (no circular trust: the owner config anchors the freeze record,
         never the reverse);
      5. freeze.code[driver_path] == sha256;
      6. sha256(driver on disk) == sha256;
      7. DEC-E1 launch flags, then the package file set, every frozen source (including the
         launcher), loaded modules, environment, bytecode policy and the manifest.
    These are SOURCE-FILE and launch-flag checks; they are not an attestation of the code
    executing in memory (residual TOCTOU and interpreter-binary limitations, DEC-E1)."""
    root = Path(config.root)

    def refuse(msg):
        raise CodeIdentityRefusal(f"code identity not verified: {msg}")
    d13 = ((config.raw.get("owner_decisions_outside_schema") or {}).get("D13_driver_freeze") or {})
    if not (d13.get("approved") is True and _nonblank(d13.get("date")) and _nonblank(d13.get("reference"))):
        refuse("Decision 13 freeze is not approved (owner_decisions_outside_schema.D13_driver_freeze)")
    anchor, driver_sha = d13.get("freeze_record_sha256"), d13.get("sha256")
    driver_path = d13.get("driver_path")                         # required (no default)
    for label, h in (("freeze_record_sha256", anchor), ("sha256 (driver)", driver_sha)):
        if not (isinstance(h, str) and _HEX64.match(h)):
            refuse(f"approved Decision 13 {label} is missing or not 64 hex characters")
    if driver_path != DRIVER_RELPATH:
        refuse(f"D13_driver_freeze.driver_path must be {DRIVER_RELPATH}")
    fp = root / FREEZE_RECORD_RELPATH
    try:
        fraw = fp.read_bytes()
    except OSError:
        refuse(f"freeze record {FREEZE_RECORD_RELPATH} is missing or unreadable")
    if _sha(fraw) != anchor:
        refuse("freeze record does not match the approved Decision 13 hash")
    try:
        freeze = json.loads(fraw)
    except ValueError:
        refuse("freeze record is not valid JSON")
    if not isinstance(freeze, dict):
        refuse("freeze record schema invalid")
    code, envf = freeze.get("code"), freeze.get("environment")
    if freeze.get("schema") != FREEZE_SCHEMA or not isinstance(code, dict) or not isinstance(envf, dict):
        refuse("freeze record schema invalid")
    if any(str(k).startswith("config/") for k in code):
        refuse("freeze record must not contain configuration files (circular trust)")
    if code.get(driver_path) != driver_sha:
        refuse("approved driver sha256 differs from the freeze inventory's driver entry")
    try:
        if _sha((root / driver_path).read_bytes()) != driver_sha:
            refuse("driver file differs from the approved driver sha256")
    except OSError:
        refuse(f"{driver_path} is missing")
    fl = _runtime_flags()                               # DEC-E1 launch evidence
    if not (fl.get("isolated") is True and fl.get("dont_write_bytecode") is True):
        refuse("not launched with -I -B (use scripts/42_phase9_corrected_launch.sh)")
    pfx = fl.get("pycache_prefix")
    if not (isinstance(pfx, str) and pfx):
        refuse("not launched with -X pycache_prefix (fresh external bytecode directory)")
    pfx_path = Path(pfx).resolve()
    root_res = root.resolve()
    if pfx_path == root_res or root_res in pfx_path.parents:
        refuse("pycache_prefix must be outside the repository")
    if not pfx_path.is_dir() or any(pfx_path.iterdir()):
        refuse("pycache_prefix must be an existing, empty directory")
    pkg_disk = {f"{PACKAGE_DIR}/{p.name}" for p in (root / PACKAGE_DIR).glob("*.py")}
    pkg_frozen = {k for k in code if k.startswith(PACKAGE_DIR + "/")}
    if pkg_disk != pkg_frozen:
        refuse(f"package file set differs from the freeze: {sorted(pkg_disk ^ pkg_frozen)}")
    required = pkg_disk | set(REQUIRED_FILES) | {str(x) for x in extra_files}
    for rel in sorted(required):
        if rel not in code:
            refuse(f"{rel} is not in the freeze record")
        try:
            h = _sha((root / rel).read_bytes())
        except OSError:
            refuse(f"{rel} is missing")
        if h != code[rel]:
            refuse(f"{rel} differs from the frozen hash")
    try:
        verify_loaded_modules(root, code, required)
    except CodeIdentityRefusal as ex:
        refuse(str(ex))
    cur = current_environment()
    for k in ENV_KEYS:
        if k not in envf:
            refuse(f"freeze environment lacks {k}")
        if _jsonable(cur[k]) != envf[k]:
            refuse(f"environment {k} differs from the freeze")
    recs = freeze.get("distribution_records")          # optional strongest check (item 19)
    if recs is not None:
        for name, rec_sha in recs.items():
            try:
                dist = importlib.metadata.distribution(name)
                rf = next(f for f in dist.files if f.name == "RECORD")
                if _sha(Path(dist.locate_file(rf)).read_bytes()) != rec_sha:
                    refuse(f"distribution {name} RECORD differs from the freeze")
            except (StopIteration, importlib.metadata.PackageNotFoundError, OSError, TypeError):
                refuse(f"distribution {name} RECORD cannot be verified")
    policy = freeze.get("bytecode_policy")
    if policy not in BYTECODE_POLICIES:
        refuse("freeze bytecode_policy missing or unknown")
    if policy == "no_pyc" and not _dont_write_bytecode():
        refuse("bytecode policy no_pyc requires PYTHONDONTWRITEBYTECODE / -B")
    for rel in sorted(required):
        if not rel.endswith(".py") or policy == E1_BYTECODE_POLICY:
            continue                     # under fresh_compile_prefix in-tree caches are never read
        pyc = Path(importlib.util.cache_from_source(str(root / rel)))
        if pyc.exists():
            if policy == "no_pyc":
                refuse(f"stale bytecode present: {pyc}")
            if _pyc_flags(pyc) & 0b11 != 0b11:
                refuse(f"bytecode {pyc} is not checked-hash based")
    if manifest is not None:
        mc = manifest.get("code") if isinstance(manifest, dict) else None
        if not isinstance(mc, dict) or not mc:
            refuse("manifest has no code hashes")
        for rel, h in sorted(mc.items()):
            if code.get(rel) != h:
                refuse(f"manifest code hash for {rel} differs from the freeze")
    return {"freeze_record_sha256": anchor, "driver_sha256": driver_sha, "freeze": freeze,
            "bytecode_policy": policy}


# ----------------------------------------------------------------- Decision 7B sizing (item 7)
def _records_by_k(raw: bytes) -> dict:
    _, _, recs, problems = _parse_checkpoint(raw)
    if problems:
        raise ValueError(problems)
    return {r["k"]: r for r in recs}


def compute_sizing_decision(config: CorrectedConfig, B_cur: int, records_by_run: dict) -> dict:
    """Deterministic Decision 5A/5B/6A/6B/7A/7B sizing computation at B_cur, from records.
    The single source of truth for the future runner and for Path C recomputation."""
    from .reliability import section_6_5
    from .statistics import (OwnerReviewRequired, combined_next_size, endpoint_mcse,
                             ep1_evaluate, ladder_assessment, w5_decision, w5_evaluate)
    lad, w5c = config.get("ladder"), config.get("w5") or {}
    per_run, run_values = {}, {}
    for r in sorted(records_by_run):
        recs = [records_by_run[r][k] for k in range(B_cur)]
        ok_vals = [x["dC_corrected"] for x in recs if x.get("status") == "ok"]
        # DEC-D: the §6.5 MC-SE population is every replicate with status "ok" (each has a
        # dC_corrected; apparent-failure and inner-shortfall replicates have none and are
        # excluded by §6.5 criteria 1-2). Per endpoint: MIN over {0.5h, h, 2h} (Decision 5B,
        # NOT the EP-1 Option B maximum); §6.5 uses the MIN of the two endpoints; the
        # threshold is computed BEFORE the bounding-width exclusion (criterion 3).
        mc = endpoint_mcse(ok_vals)
        mcse = {e: min(mc[e]["SE"].values()) for e in ("lower", "upper")}      # Decision 5B
        s65 = section_6_5(recs, min(mcse.values()))
        h1 = s65["H1"] == "UNRELIABLE"
        usable = set(s65["usable_k"])
        vbk = {k: (recs[k]["dC_corrected"] if k in usable else None) for k in range(B_cur)}
        run_values[r] = np.array([vbk[k] for k in range(B_cur) if vbk[k] is not None])
        w5 = w5_evaluate(vbk, B_cur, exclusion_rule=w5c.get("exclusion_rule"), h1_triggered=h1)
        per_run[r] = {"endpoint_mcse_min": mcse, "n_excluded_6_5": s65["n_excluded"], "H1": s65["H1"],
                      "w5": w5, "w5_decision": w5_decision(w5, w5c.get("action"),
                                                           w5c.get("corrective_rule"), B_cur, lad)}
    ep1 = ep1_evaluate(run_values, config.get("ep1_option"))
    ep1_dec = ladder_assessment(B_cur, ep1, lad)
    w5_decs = [per_run[r]["w5_decision"] for r in sorted(per_run)]
    combined, action = None, None
    if any(d["action"] == "stop_for_owner_review" for d in w5_decs):
        action = "stop_for_owner_review"
    else:
        w5_comb = max(w5_decs, key=lambda d: d.get("B_next") or 0)
        try:
            combined = combined_next_size(ep1_dec, w5_comb)
        except OwnerReviewRequired:
            action = "stop_for_owner_review"
        if action is None:
            if combined:
                action = "top_up"
            elif ep1_dec["action"] != "stop_ep1_pass":
                action = ep1_dec["action"]
            elif any(d["action"] == "report_persisting_breach" for d in w5_decs):
                action = "report_persisting_breach"
            else:
                action = "stop_ep1_pass"
    return _jsonable({"B_cur": B_cur, "ep1": ep1, "ep1_decision": ep1_dec, "per_run": per_run,
                      "combined_B_next": combined, "action": action})


def build_sizing_decision(config: CorrectedConfig, B_cur: int, manifest_sha256: str,
                          freeze_info: dict) -> dict:
    """The artifact CONTENT (p9c-sizing-1). Writing it is the future runner's job; it must
    be written with durable_create to _sizing_path(output_dir, B_cur)."""
    out = Path(config.output_dir)
    runs, records = {}, {}
    for r in config.run_ids:
        raw = (out / "runs" / f"run{r}.jsonl").read_bytes()
        st = invocation_state(out, r)
        latest = max(st["invocations"]) if st["invocations"] else None
        el = evaluate_eligibility(config, r, B_cur, manifest_sha256, max_invocation=latest,
                                  checkpoint_raw=raw,
                                  freeze_record_sha256=freeze_info["freeze_record_sha256"])
        if el["class"] != "ELIGIBLE":
            raise DispatchRefused(f"run {r} is {el['class']} at B_cur={B_cur}: no sizing decision")
        records[r] = _records_by_k(raw)
        runs[str(r)] = {"checkpoint_sha256": _sha(raw), "checkpoint_size_bytes": len(raw),
                        "n_records": len(records[r]), "latest_invocation": latest,
                        "eligibility_result_sha256": el["eligibility_result_sha256"]}
    code = freeze_info["freeze"]["code"]
    return _jsonable({
        "kind": "sizing_decision", "schema": "p9c-sizing-1", "B_cur": B_cur,
        "run_ids": list(config.run_ids), "config_sha256": config.sha256(),
        "manifest_sha256": manifest_sha256, "runs": runs,
        "decision": compute_sizing_decision(config, B_cur, records),
        "rules": {"ep1": "Decision 5A", "mcse_6_4_6_5_H2": "Decision 5B", "w5": "Decisions 6A/6B",
                  "ladder": "Decisions 7A/7B", "ladder_rungs": (config.get("ladder") or {}).get("rungs"),
                  "beyond_ladder_rule": (config.get("ladder") or {}).get("beyond_ladder_rule")},
        "code_identity": {"freeze_record_sha256": freeze_info["freeze_record_sha256"],
                          "package": {k: code[k] for k in sorted(code)}}})


def verify_sizing_decision(config: CorrectedConfig, run_id: int, L_cur: int, target: int,
                           freeze_info: dict, manifest_sha256: str) -> str:
    """Path C layer 2: accept the sizing file ONLY if every value is re-derived from the
    approved freeze anchor and the actual checkpoints. Returns its file SHA-256."""
    out = Path(config.output_dir)

    def refuse(msg):
        raise DispatchRefused(f"top-up refused: sizing decision for B_cur={L_cur}: {msg}")
    p = _sizing_path(out, L_cur)
    try:
        sraw = p.read_bytes()
        art = json.loads(sraw)
    except (OSError, ValueError):
        refuse(f"{p} missing or unreadable")
    try:
        if art.get("kind") != "sizing_decision" or art.get("schema") != "p9c-sizing-1":
            refuse("schema invalid")
        if art.get("B_cur") != L_cur or art.get("run_ids") != list(config.run_ids):
            refuse("B_cur or run_ids inconsistent")
        if art.get("config_sha256") != config.sha256() or art.get("manifest_sha256") != manifest_sha256:
            refuse("configuration or manifest identity differs")
        code = freeze_info["freeze"]["code"]
        ci = art.get("code_identity") or {}
        if ci.get("freeze_record_sha256") != freeze_info["freeze_record_sha256"] or \
                ci.get("package") != {k: code[k] for k in sorted(code)}:
            refuse("code identity does not resolve to the verified freeze")
        sizing_sha = _sha(sraw)
        records = {}
        for r in config.run_ids:
            e = art["runs"][str(r)]
            craw = (out / "runs" / f"run{r}.jsonl").read_bytes()
            size = e["checkpoint_size_bytes"]
            if len(craw) < size or _sha(craw[:size]) != e["checkpoint_sha256"]:
                refuse(f"run {r} checkpoint no longer contains the recorded state")
            if r == run_id and len(craw) != size:
                refuse(f"run {r} checkpoint changed since the sizing decision")
            # D11-D2: the COMPLETE CURRENT history of EVERY run is assessed, independent of
            # checkpoint growth. A later zero-growth disqualifying invocation is caught here.
            a = assess_run(config, r, manifest_sha256)
            bad = a["reasons"].codes("disqualifier") + a["reasons"].codes("blocking")
            if bad:
                refuse(f"run {r} history is not clean: {sorted(set(bad))}")
            why = _later_invocations_problem(out, r, e["latest_invocation"], L_cur, target,
                                             sizing_sha)
            if why:
                refuse(f"run {r}: {why}")
            snap = craw[:size]
            el = evaluate_eligibility(config, r, L_cur, manifest_sha256,
                                      max_invocation=e["latest_invocation"], checkpoint_raw=snap,
                                      freeze_record_sha256=freeze_info["freeze_record_sha256"])
            if el["class"] != "ELIGIBLE" or el["eligibility_result_sha256"] != e["eligibility_result_sha256"]:
                refuse(f"run {r} eligibility at B_cur is not reproduced ({el['class']})")
            records[r] = _records_by_k(snap)
            if e.get("n_records") != L_cur or len(records[r]) != L_cur:
                refuse(f"run {r} record count at B_cur is not {L_cur}")
        dec = compute_sizing_decision(config, L_cur, records)
    except DispatchRefused:
        raise
    except Exception as ex:                          # any recomputation failure fails closed
        refuse(f"recomputation failed: {type(ex).__name__}: {ex}")
    if _canon(dec) != _canon(art.get("decision")):
        refuse("recorded decision differs from the independent recomputation")
    if dec["action"] != "top_up" or dec["combined_B_next"] != target:
        refuse(f"recomputed decision is {dec['action']} -> {dec['combined_B_next']}, not {target}")
    return _sha(sraw)


def _later_invocations_problem(out: Path, r: int, latest: int, L_cur: int, target: int,
                               sizing_sha: str) -> Optional[str]:
    """D11-D2 transition model for invocations of run r recorded after the sizing snapshot.
    Validity depends on path, target, reference and outcome -- never on records appended.
    Allowed, in order: (a) benign no-op resumes at the sized level L_cur that completed
    with zero new records; then (b) this run's own top-up to `target`: exactly one
    `topup` invocation referencing THIS sizing file and starting from L_cur, followed by any
    number of `resume` continuations at `target` carrying the same top-up reference (after a
    runtime halt + authorisation or an acknowledged interruption, Decision 10).
    Disqualifying or blocking outcomes are rejected separately by assess_run."""
    st = invocation_state(out, r)
    invs = st["invocations"]
    opened, topup_ref = False, None
    for m in sorted(x for x in invs if x > latest):
        s, end = invs[m].get("start") or {}, invs[m].get("end")
        t, path = s.get("target_B_out"), s.get("path")
        ref = s.get("topup_authorization")
        ref = ref.strip() if isinstance(ref, str) else None
        if (not opened and t == L_cur and path == "resume" and end is not None
                and end.get("outcome") == "completed" and end.get("done_now") == 0):
            continue                                     # benign no-op at the sized level
        if t != target:
            return f"invocation {m} targets {t}, not the decided target {target}"
        if not opened:
            if (path != "topup" or s.get("sizing_decision_sha256") != sizing_sha
                    or s.get("n_start") != L_cur or not ref):
                return f"invocation {m} is not the top-up authorised by this sizing decision"
            opened, topup_ref = True, ref
        elif path != "resume" or ref != topup_ref:
            return (f"invocation {m} is not a continuation of the authorised top-up "
                    f"(path {path!r}, reference mismatch: {ref != topup_ref})")
    return None


def _refusal_message(R: _Reasons) -> str:
    codes = R.codes("disqualifier") + R.codes("blocking")
    msg = f"dispatch refused ({', '.join(sorted(set(codes)))})"
    if "chain_break" in codes or "region_mismatch" in codes:
        msg += ": checkpoint changed outside any recorded invocation - its content has changed"
    if "nonestimable_prefix_violation" in codes or "history_nonestimable_stop" in codes:
        msg += ": the non-estimable 5 % rule was violated (permanent under Option A)"
    if "preservation_failure_pending" in codes:
        msg += ": historical-preservation failure pending owner adjudication"
    if "window_recomputation_failed" in codes:
        msg += ": a window recomputation FAILED (permanent disqualifier)"
    return msg


# ======================================================================= run
def run_corrected(config: CorrectedConfig, run_id: int, target_B_out: int,
                  topup_authorization: str | None = None,
                  resume_authorization: str | None = None, clock=None) -> dict:
    """`clock` = (name, callable) for tests; default = the approved platform clock."""
    config.validate_for_execution()                       # fail fast on any owner item
    root, out = config.root, config.output_dir
    if run_id not in config.run_ids:
        raise ConfigurationError(f"run_id {run_id} not in configured run_ids")
    if target_B_out not in allowed_targets(config):
        raise ConfigurationError(f"target B_out {target_B_out} is not the initial size or a ladder rung")
    if topup_authorization is not None and not _nonblank(topup_authorization):
        raise ConfigurationError("--topup-authorization reference is blank")
    if target_B_out > config.B_out_initial and not topup_authorization:
        raise ConfigurationError("a top-up beyond B_out_initial needs an explicit authorisation "
                                 "reference (EP-1/W5 decision output)")
    if resume_authorization is not None and not str(resume_authorization).strip():
        raise InvocationStateError("--resume-authorization reference is blank")
    clock_name, clk = select_clock() if clock is None else clock
    assert_output_dir_new_or_ours(out)
    freeze_info = verify_code_identity(config)            # Decision 9 cond. 4, every path
    h4_check(root, "pre")
    pres = historical_preservation_check(root)
    if not pres["all_ok"]:
        _record_preservation_failure(out, preservation_evidence(pres, "pre_dispatch"))
        raise SystemExit("historical artifact changed -- refusing to run")
    manifest = build_manifest(config, root)
    msha = write_manifest_once(out / "manifest" / "input_manifest.json", manifest)
    verify_code_identity(config, manifest=manifest)       # manifest code == freeze
    streams = PlanStreams(config.base, run_id, config.convention, config.T)
    n_workers = int(config.get("n_workers", 8))
    budget = config.get("runtime_budget_hours_per_run")

    with CheckpointStore(out, run_id, config.sha256(), msha, root) as store:
        state = invocation_state(out, run_id)
        if state["state"] == "AMBIGUOUS":
            raise InvocationStateError(f"run {run_id}: invocation state AMBIGUOUS -- owner "
                                       f"inspection required: {state['problems']}")
        assessment = assess_run(config, run_id, msha)
        R = assessment["reasons"]
        if _dispatch(R, state["state"], False) in ("blocked_permanently", "blocked_pending_owner"):
            log_event(out, {"event": "dispatch_refused", "run_id": run_id,
                            "reasons": [i["code"] for i in R.items
                                        if i["category"] in ("disqualifier", "blocking")]})
            raise DispatchRefused(f"run {run_id}: " + _refusal_message(R))
        # ---- three-path classification (Decision 11 item 6)
        levels, ms = assessment["levels"], assessment["invocations"]
        L_cur = levels.get("L_cur")
        ref = topup_authorization.strip() if topup_authorization else None
        sizing_sha = None
        if not ms:
            if len(store.done) != 0 or target_B_out != config.B_out_initial or ref:
                raise DispatchRefused("a fresh run must target B_out_initial with no records "
                                      "and no top-up reference")
            path = "fresh"
        elif target_B_out == L_cur:
            want = levels["refs"].get(str(L_cur)) if L_cur > config.B_out_initial else None
            if ref != want:
                raise DispatchRefused(f"resume at level {L_cur} requires its original top-up "
                                      "reference (and none at the initial level)")
            path = "resume"
        elif target_B_out > L_cur:
            if state["state"] != "CLEAN":
                raise DispatchRefused("a top-up is not permitted while an invocation is pending")
            if L_cur not in levels["completed"] or len(store.done) != L_cur:
                raise DispatchRefused(f"top-up refused: level {L_cur} is not complete")
            if not ref or ref in levels["refs"].values():
                raise DispatchRefused("top-up refused: a new, non-blank reference is required")
            sizing_sha = verify_sizing_decision(config, run_id, L_cur, target_B_out,
                                                freeze_info, msha)
            path = "topup"
        else:
            raise DispatchRefused(f"target {target_B_out} is below the current level {L_cur}")
        # Decision 13: the complete gate, executed freshly for THIS attempt, before any data
        # is loaded, any authorisation is consumed or any start record is written.
        gate_ev = run_fresh_gate(config, out, run_id, state["next_invocation"], msha)
        ctx = load_context(root)          # after a genuine PASS, BEFORE any authorisation is used
        auth = resolve_authorization(out, run_id, state, resume_authorization, store,
                                     config.sha256(), msha)
        monitor = RuntimeMonitor(budget, len(store.done), target_B_out, n_workers, clk)
        m = state["next_invocation"]
        prior = state["cumulative_elapsed_s_known"]
        start = {"kind": "start", "run_id": run_id, "invocation": m,
                 "utc": datetime.now(timezone.utc).isoformat(), "target_B_out": target_B_out,
                 "clock": clock_name, **monitor.describe(),
                 "config_sha256": config.sha256(), "manifest_sha256": msha,
                 "checkpoint_sha256": store.file_sha256(),
                 "checkpoint_size_bytes": store.path.stat().st_size,
                 "topup_authorization": topup_authorization, "path": path,
                 "sizing_decision_sha256": sizing_sha,
                 "gate_evidence": gate_ev,
                 "freeze_record_sha256": freeze_info.get("freeze_record_sha256"),
                 "resolves_invocation": auth["invocation"] if auth else None,
                 "cumulative_elapsed_s_known_before": prior,
                 "n_invocations_elapsed_unknown_before": state["n_invocations_elapsed_unknown"]}
        durable_create(_inv_path(out, run_id, m, "start"), start)   # before any dispatch
        log_event(out, {"event": "run_start", "run_id": run_id, "invocation": m,
                        "target_B_out": target_B_out, "config_sha256": config.sha256(),
                        "manifest_sha256": msha, "topup_authorization": topup_authorization,
                        "clock": clock_name, "runtime_check": monitor.runtime_check,
                        "W": monitor.W, "n_new": monitor.n_new,
                        "threshold_s": monitor.threshold_s, "path": path})
        monitor.start()

        def on_record(rec, st):
            h4_check(root, f"during run {run_id} k={rec['k']}")
            ne = nonestimable_status(st.done.values())
            halt = monitor.observe(len(st.done))
            reasons = ([r for r, hit in (("nonestimable_rate", ne["result"]),
                                         ("runtime_projection", bool(halt))) if hit])
            if halt:
                log_event(out, {"event": "halt", "run_id": run_id, "invocation": m,
                                "reason": "projected runtime > 2x declared budget (plan 9.5)",
                                "stop_reasons": reasons, **halt})
                raise RuntimeBudgetHalt("STOP: projected runtime exceeds twice the declared "
                                        "budget (plan 9.5)", halt, nonestimable=ne, reasons=reasons)
            if ne["result"]:
                log_event(out, {"event": "halt", "run_id": run_id,
                                "reason": "outer non-estimable rate > 5 % (plan 9.5)"})
                raise NonEstimableStop("STOP: outer non-estimable rate exceeds 5 % (plan 9.5)",
                                       nonestimable=ne, reasons=reasons)

        try:
            execute_replicates(store, ctx, streams, RealPipeline(), config.B_in, target_B_out,
                               n_workers, bool(config.as_run_variant), on_record)
            h4_check(root, "post")
            post = historical_preservation_check(root)
            if not post["all_ok"]:
                raise HistoricalPreservationFailure("historical artifact changed during the run",
                                                    evidence=preservation_evidence(post, "post"))
        except RuntimeBudgetHalt as ex:
            _write_end(out, run_id, m, "runtime_halt", monitor, store, prior,
                       {**ex.details, "clock": clock_name}, stop_reasons=ex.reasons,
                       nonestimable=ex.nonestimable)
            raise
        except NonEstimableStop as ex:
            _write_end(out, run_id, m, "stop_nonestimable", monitor, store, prior,
                       stop_reasons=ex.reasons, nonestimable=ex.nonestimable)
            raise
        except HardFailureH4 as ex:
            _write_end(out, run_id, m, "hard_failure_h4", monitor, store, prior, {"error": str(ex)})
            raise
        except HistoricalPreservationFailure as ex:
            _write_end(out, run_id, m, "historical_preservation_failure", monitor, store, prior,
                       ex.evidence or preservation_evidence(None, "post"))
            raise
        except ImplementationDefectH6 as ex:
            _write_end(out, run_id, m, "hard_failure_h6", monitor, store, prior, {"error": str(ex)})
            raise
        # any other exception (worker error, KeyboardInterrupt, kill) leaves no end record:
        # the invocation is UNTERMINATED and needs explicit acknowledgement.
        _write_end(out, run_id, m, "completed", monitor, store, prior)
        res = {"run_id": run_id, "invocation": m, "n_records": len(store.done),
               "checkpoint_sha256": store.file_sha256(),
               "minutes": (monitor.elapsed()) / 60}
        log_event(out, {"event": "run_complete", **res})
        return res
