"""Execution configuration with explicit, unresolved owner fields.

No [OWNER] field has a default. `validate_for_execution` fails fast, listing every
missing or unapproved item, so the pipeline cannot run on a silently chosen
interpretation.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

UINT32_MAX = 2 ** 32 - 1

# [PLAN] fixed values (plan line numbers)
PLAN_BASE = 2026110300            # line 354 `BASE = 20261103_00`
PLAN_RUN_IDS = (0, 1, 2, 3, 4)    # line 355
PLAN_B_OUT = 500                  # line 237
PLAN_B_IN = 100                   # line 238
PLAN_LADDER = (750, 1000, 1500, 2000, 3000)  # line 119
PLAN_EP1_FRACTION = 0.05          # line 84

EP1_OPTIONS = ("A", "B")                       # A: Silverman h alone; B: max over {0.5h, h, 2h}
W5_ACTIONS = ("W5-A", "W5-B")                  # A: W5 mandates top-up; B: investigate/report only
W5_CORRECTIVE_RULES = ("next_rung_above_current",)   # the only rule drafted (adjudication D5)
W5_EXCLUSION_RULES = ("usable_within_k_positions_v1",)   # owner Decision 6B, 2026-11-06
CONFIG_SCHEMA = "p9c-config-2"     # -2 adds the required w5.exclusion_rule (Decision 6B)
BEYOND_LADDER_RULES = ("top_up_to_max_rung_then_reassess", "stop_and_label_imprecise")
APPARENT_CONVENTIONS = ("independent", "as_run")

# Every owner decision that must be explicitly approved before execution.
OWNER_FIELDS = (
    "cv_seed_deviation_adjudicated",   # D1
    "apparent_cv",                     # D2 convention (+ T if independent)
    "B_in_max",                        # D2 collision domain
    "as_run_variant",                  # D2 optional paired variant
    "ep1_option",                      # D4
    "w5",                              # D5 action (+ corrective rule if W5-A)
    "ladder",                          # D4/D5 ladder + beyond-ladder rule
    "output_dir",                      # D8
    "verification_gate",               # D8
    "runtime_budget_hours_per_run",    # plan 9.5 budget (Decision 10A): expected wall-clock
                                       # hours PER RUN PER 500 NEW REPLICATES at the configured
                                       # worker count; the halt threshold is scaled per invocation
                                       # as 2 x budget x n_new / 500 (orchestration.RuntimeMonitor)
)

# Directories the corrected pipeline may never write into (or contain).
PROTECTED_DIRS = (
    "reports/phase9/dc_interval",
    "reports/phase9/dc_interval_audit",
    "reports/phase9/execution",
    "reports/phase9/execution_s7",
    "reports/phase9/h5_spotcheck",
    "reports/phase9/final_project",
    "docs",
    "data",
    "config",
    "src",
    "scripts",
    "tests",
    ".git",
    ".venv",
)


class ConfigurationError(ValueError):
    """Configuration is structurally invalid or not approved for execution."""


class UnresolvedDecisionError(RuntimeError):
    """A function was asked to apply an owner decision that has not been supplied."""


def canonical_json(obj: Any) -> bytes:
    """Canonical serialisation for every record, record hash, eligibility hash and sizing
    file (sorted keys, compact separators).

    DEC-A1 (owner decision 2026-11-06): `allow_nan=True` is retained deliberately. Non-finite
    floats are written as the Python JSON tokens `NaN`, `Infinity` and `-Infinity`, which
    are NOT RFC 8259-compliant JSON values; strict JSON parsers reject them, while this
    project's loaders (`json.loads`) accept them and round-trip them exactly, so hashes
    are stable. `-Infinity` can occur as a degenerate tuning score (every CV fold skipped:
    `inner_cv_c` on success, `tuning_score_at_failure` on failure). `NaN` is not expected in
    any record; if one appears, V5's field comparison reports a mismatch (NaN != NaN) and
    fails closed. Reports must render non-finite values explicitly (e.g. "-inf: degenerate
    tuning, all CV folds skipped") and must never use them in arithmetic.
    """
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=True).encode()


def _is_int(x) -> bool:
    return isinstance(x, int) and not isinstance(x, bool)


@dataclass
class CorrectedConfig:
    raw: dict
    root: Path
    problems: list = field(default_factory=list)

    # ------------------------------------------------------------------ access
    def get(self, key, default=None):
        return self.raw.get(key, default)

    @property
    def base(self) -> int:
        return self.raw["base"]

    @property
    def run_ids(self) -> tuple:
        return tuple(self.raw["run_ids"])

    @property
    def B_out_initial(self) -> int:
        return self.raw["B_out_initial"]

    @property
    def B_in(self) -> int:
        return self.raw["B_in"]

    @property
    def convention(self):
        return (self.raw.get("apparent_cv") or {}).get("convention")

    @property
    def T(self):
        return (self.raw.get("apparent_cv") or {}).get("T")

    @property
    def B_in_max(self):
        return self.raw.get("B_in_max")

    @property
    def as_run_variant(self):
        return self.raw.get("as_run_variant")

    @property
    def output_dir(self):
        v = self.raw.get("output_dir")
        if v is None:
            return None
        p = Path(v)
        return p if p.is_absolute() else (self.root / p)

    def sha256(self) -> str:
        """Hash of the scientific configuration. n_workers is excluded: results do not
        depend on the worker count (plan line 371; verified by the integration tests)."""
        return hashlib.sha256(canonical_json({k: v for k, v in self.raw.items()
                                              if k != "n_workers"})).hexdigest()

    # ------------------------------------------------------------- validation
    def validate_structure(self) -> list:
        """Type/range checks that hold regardless of approval. Returns problems."""
        p = []
        r = self.raw
        if r.get("schema") != CONFIG_SCHEMA:
            p.append(f"schema must be '{CONFIG_SCHEMA}' (older schema versions are refused)")
        if r.get("base") != PLAN_BASE:
            p.append(f"base must equal the plan value {PLAN_BASE} (plan line 354)")
        rid = r.get("run_ids")
        if not (isinstance(rid, list) and rid and all(_is_int(x) and x in PLAN_RUN_IDS for x in rid)
                and len(set(rid)) == len(rid)):
            p.append("run_ids must be a non-empty list of distinct values from {0..4} (plan line 355)")
        if r.get("B_out_initial") != PLAN_B_OUT:
            p.append(f"B_out_initial must equal the plan value {PLAN_B_OUT} (plan line 237); "
                     "a different size is an owner decision with no approval field")
        if r.get("B_in") != PLAN_B_IN:
            p.append(f"B_in must equal the plan value {PLAN_B_IN} (plan line 238); raising B_in "
                     "(W3) has no defined procedure")
        ac = r.get("apparent_cv")
        if ac is not None:
            conv = ac.get("convention")
            if conv not in APPARENT_CONVENTIONS + (None,):
                p.append(f"apparent_cv.convention must be one of {APPARENT_CONVENTIONS}")
            T = ac.get("T")
            if conv == "independent":
                if T is None:
                    p.append("apparent_cv.T is required for the independent convention")
                elif not (_is_int(T) and 1 <= T <= UINT32_MAX):
                    p.append("apparent_cv.T must be an integer in [1, 2**32-1] "
                             "(non-zero, one uint32 word)")
            if conv == "as_run" and T is not None:
                p.append("apparent_cv.T must be null under the as_run convention")
        bmax = r.get("B_in_max")
        if bmax is not None:
            if not (_is_int(bmax) and 1 <= bmax <= UINT32_MAX):
                p.append("B_in_max must be an integer in [1, 2**32-1]")
            elif _is_int(r.get("B_in")) and r["B_in"] > bmax:
                p.append("B_in exceeds the declared B_in_max")
        if r.get("as_run_variant") not in (None, True, False):
            p.append("as_run_variant must be true, false or null")
        if r.get("as_run_variant") is True and (ac or {}).get("convention") == "as_run":
            p.append("as_run_variant is redundant when the primary convention is as_run")
        if r.get("ep1_option") not in EP1_OPTIONS + (None,):
            p.append(f"ep1_option must be one of {EP1_OPTIONS}")
        w5 = r.get("w5")
        if w5 is not None:
            if w5.get("action") not in W5_ACTIONS:
                p.append(f"w5.action must be one of {W5_ACTIONS}")
            if w5.get("action") == "W5-A" and w5.get("corrective_rule") not in W5_CORRECTIVE_RULES:
                p.append(f"w5.corrective_rule must be one of {W5_CORRECTIVE_RULES} under W5-A")
            if w5.get("definition") != "adjudication_D5_as_proposed":
                p.append("w5.definition must be 'adjudication_D5_as_proposed' (no other "
                         "definition is implemented)")
            if w5.get("exclusion_rule") not in W5_EXCLUSION_RULES:
                p.append(f"w5.exclusion_rule must be one of {W5_EXCLUSION_RULES} "
                         "(required whenever w5 is supplied; Decision 6B)")
        lad = r.get("ladder")
        if lad is not None:
            rungs = lad.get("rungs")
            if not (isinstance(rungs, list) and rungs and all(_is_int(x) for x in rungs)
                    and rungs == sorted(set(rungs))
                    and _is_int(r.get("B_out_initial")) and rungs[0] > r["B_out_initial"]
                    and rungs[-1] <= UINT32_MAX):
                p.append("ladder.rungs must be strictly increasing integers above B_out_initial")
            if lad.get("beyond_ladder_rule") not in BEYOND_LADDER_RULES:
                p.append(f"ladder.beyond_ladder_rule must be one of {BEYOND_LADDER_RULES}")
        rb = r.get("runtime_budget_hours_per_run")
        if rb is not None and not (isinstance(rb, (int, float)) and not isinstance(rb, bool) and rb > 0):
            p.append("runtime_budget_hours_per_run must be a positive number")
        vg = r.get("verification_gate")
        if vg is not None and vg != "V1-V5_H4_as_implemented_p9c-1":
            p.append("verification_gate must be 'V1-V5_H4_as_implemented_p9c-1'")
        od = r.get("output_dir")
        if od is not None:
            try:
                assert_output_dir_allowed(self.output_dir, self.root)
            except ConfigurationError as ex:
                p.append(str(ex))
        if not (_is_int(r.get("n_workers", 1)) and r.get("n_workers", 1) >= 1):
            p.append("n_workers must be a positive integer")
        return p

    def missing_owner_items(self) -> list:
        """Every owner field that is absent or not explicitly approved."""
        r = self.raw
        out = []
        appr = r.get("approvals") or {}
        for f in OWNER_FIELDS:
            if f == "cv_seed_deviation_adjudicated":
                if r.get(f) is not True:
                    out.append(f"{f}: Decision 1 not recorded as adjudicated")
            elif f == "apparent_cv":
                if not r.get(f) or r[f].get("convention") is None:
                    out.append("apparent_cv.convention: not supplied (Decision 2)")
                elif r[f]["convention"] == "independent" and r[f].get("T") is None:
                    out.append("apparent_cv.T: not supplied (Decision 2)")
            elif r.get(f) is None:
                out.append(f"{f}: not supplied")
            a = appr.get(f) or {}
            if not (a.get("approved") is True and a.get("date") and a.get("reference")):
                out.append(f"{f}: no explicit owner approval (approved/date/reference)")
        ea = r.get("execution_authorization") or {}
        if not (ea.get("authorized") is True and ea.get("date") and ea.get("reference")):
            out.append("execution_authorization: rerun not authorised")
        return out

    def validate_for_execution(self) -> None:
        problems = self.validate_structure() + self.missing_owner_items()
        if problems:
            raise ConfigurationError("configuration not executable:\n  - " + "\n  - ".join(problems))


def assert_output_dir_allowed(path, root: Path) -> None:
    """The output directory must be explicit, inside the repository, and disjoint
    from every protected (historical / source) directory."""
    if path is None:
        raise ConfigurationError("output_dir: not supplied")
    p = Path(path).resolve()
    root = Path(root).resolve()
    try:
        p.relative_to(root)
    except ValueError:
        raise ConfigurationError(f"output_dir {p} is outside the repository root") from None
    if p == root:
        raise ConfigurationError("output_dir may not be the repository root")
    for d in PROTECTED_DIRS:
        q = (root / d).resolve()
        if p == q or q in p.parents or p in q.parents:
            raise ConfigurationError(f"output_dir {p} overlaps protected directory {d}")


def load_config(path, root: Path) -> CorrectedConfig:
    raw = json.loads(Path(path).read_text())
    return CorrectedConfig(raw=raw, root=Path(root))


def template_config() -> dict:
    """A template with every owner field UNRESOLVED (null). It is not executable."""
    return {
        "schema": CONFIG_SCHEMA,
        "base": PLAN_BASE,
        "run_ids": list(PLAN_RUN_IDS),
        "B_out_initial": PLAN_B_OUT,
        "B_in": PLAN_B_IN,
        "n_workers": 8,
        "cv_seed_deviation_adjudicated": None,
        "apparent_cv": {"convention": None, "T": None},
        "B_in_max": None,
        "as_run_variant": None,
        "ep1_option": None,
        "w5": None,
        "ladder": None,
        "output_dir": None,
        "verification_gate": None,
        "runtime_budget_hours_per_run": None,
        "approvals": {},
        "execution_authorization": {"authorized": False, "date": None, "reference": None},
    }
