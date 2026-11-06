"""Decision 13 P1/P4 -- post-run verification runner (src/phase9_corrected/postrun.py).

SYNTHETIC ONLY. Campaigns are produced by run_corrected with a fake replicate function and
stubbed data loading, gate, H4 and code identity; the runner is given a synthetic
VerificationEnv (no real data, no real pipeline, no historical results). Record-level
verification functions are replaced by doubles where a test exercises the runner's
control flow and verdict semantics; the real functions run on synthetic records elsewhere.
These tests demonstrate software behaviour only, never scientific validity."""
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import types
import uuid
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from src.phase9_corrected import orchestration as orch
from src.phase9_corrected import postrun
from src.phase9_corrected import verification as V
from src.phase9_corrected.config import CorrectedConfig, canonical_json
from src.phase9_corrected.data import LAM_GRID
from src.phase9_corrected.engine import compute_outer_replicate
from src.phase9_corrected.rng import PlanStreams
from tests.p9c_support import FakePipeline, approved_test_config, synthetic_context
from tests.test_phase9_corrected_sizing_expected import expected as sizing_oracle

ROOT = Path(__file__).resolve().parents[1]
ORIGINAL_V4 = V.v4_as_run_apparent
ORIGINAL_VERIFY_LOADED = orch.verify_loaded_modules
SYNTH_FREEZE = {"freeze_record_sha256": "f" * 64,
                "freeze": {"code": {"x": "0" * 64}, "environment_record_only": {"jq_version": "synthetic"}}}
W5A = {"action": "W5-A", "corrective_rule": "next_rung_above_current",
       "definition": "adjudication_D5_as_proposed", "exclusion_rule": "usable_within_k_positions_v1"}
LADDER_750 = {"rungs": [750], "beyond_ladder_rule": "top_up_to_max_rung_then_reassess"}
LADDER_FULL = {"rungs": [750, 1000, 1500, 2000, 3000], "beyond_ladder_rule": "top_up_to_max_rung_then_reassess"}

S = {"drift": 0.0, "raise_at": None, "perturb": set(), "nan": False}


def value(run_id, k):
    return float(np.random.default_rng([run_id, k]).normal(0.0, 0.01) + S["drift"] * k / 500)


def fake_record(ctx, streams, pipeline, run_id, k, B_in, lam_grid, as_run_variant=False):
    if S["raise_at"] is not None and k == S["raise_at"]:
        raise RuntimeError(f"simulated crash at k={k}")
    v = value(run_id, k) + (1e-9 if (run_id, k) in S["perturb"] else 0.0)
    rec = {"run_id": int(run_id), "k": int(k), "status": "ok", "B_in_attempted": 0, "inner": [],
           "apparent": {"M1": {"ok": True}, "M2": {"ok": True}},
           "dC_corrected": v, "dC_app": v}
    if as_run_variant:
        rec["dC_app_as_run"] = v
        rec["apparent_as_run"] = {m: {"ok": True, "C_fit": 0.6, "lam": 0.0} for m in ("M1", "M2")}
    if S["nan"]:
        rec["diag_list"] = [1.0, float("nan")]
        rec["diag_dict"] = {"a": float("nan")}
        rec["diag_inf"] = float("-inf")
        if S["nan"] == "top":
            rec["diag"] = float("nan")
    return rec


class StepClock:
    def __init__(self, step=0.001):
        self.step, self.n = step, 0

    def __call__(self):
        v = self.n * self.step
        self.n += 1
        return v


class Campaign:
    def __init__(self, out_rel):
        self.out_rel, self.out = out_rel, ROOT / out_rel
        self.over = {}

    def cfg(self, **extra):
        return CorrectedConfig(approved_test_config(ROOT, self.out_rel, n_workers=1,
                                                    runtime_budget_hours_per_run=4.0,
                                                    **{**self.over, **extra}), ROOT)

    def run(self, run_id=0, target=500, **kw):
        return orch.run_corrected(self.cfg(), run_id, target, clock=("TEST_CLOCK", StepClock()), **kw)

    def msha(self):
        return hashlib.sha256((self.out / "manifest/input_manifest.json").read_bytes()).hexdigest()

    def all_runs(self, target=500, **kw):
        for r in range(5):
            self.run(run_id=r, target=target, **kw)

    def sizing(self, B_cur):
        art = orch.build_sizing_decision(self.cfg(), B_cur, self.msha(), SYNTH_FREEZE)
        orch.durable_create(orch._sizing_path(self.out, B_cur), art)
        return art

    def snapshot(self):
        return {p.relative_to(self.out).as_posix(): p.read_bytes()
                for p in sorted(self.out.rglob("*")) if p.is_file()}


# Hand-calculated V2-domain state counts (NOT the runner's formula): per run and k there are
# 1 outer + 1 apparent_cv (independent) + 100 b x 3 (inner, inner_cv, inner_reserve) = 302
# states, plus the 5 Phase 9 primary seeds.
#   5 runs x 3000 x 302 + 5 = 4,530,005   5 runs x 750 x 302 + 5 = 1,132,505
#   2 runs x 3000 x 302 + 5 = 1,812,005   1 run  x 3000 x 302 + 5 =   906,005
GATE_STATES = {(5, 3000, 100, "independent"): 4_530_005, (5, 750, 100, "independent"): 1_132_505,
               (2, 3000, 100, "independent"): 1_812_005, (1, 3000, 100, "independent"): 906_005}


def gate_body(config, msha):
    """A SYNTHETIC gate cache with the structure gate.py writes (the gate is never run)."""
    k_max = max([config.B_out_initial] + list(config.get("ladder")["rungs"]))
    n = GATE_STATES[(len(config.run_ids), k_max, config.B_in_max, config.convention)]
    return {"key": {"config_sha256": config.sha256(), "manifest_sha256": msha, "numpy": np.__version__},
            "H4": "PASS", "historical_preservation": "PASS",
            "V2_numpy": {"status": "PASS", "detail": {"all_pass": True}},
            "V2_structural": {"status": "PASS", "detail": {"holds": True}},
            "V2_domain": {"status": "PASS", "pass": True, "collisions": [], "n_states": n, "n_distinct": n,
                          "domain": {"run_ids": list(config.run_ids), "k_max": k_max,
                                     "b_max": config.B_in_max, "convention": config.convention,
                                     "T": config.T}},
            "V1_outer_historical": {"status": "PASS", "n_compared": 500 * len(config.run_ids),
                                    "mismatches": []},
            "all_pass": True}


GATE_CALLS = []


def synth_gate(config, msha, progress=None, out_path=None):
    """Writes the (synthetic, hand-counted) evidence to the requested path, as gate.py does."""
    p = Path(out_path) if out_path else Path(config.output_dir) / "gate" / "pre_execution_gate.json"
    body = gate_body(config, msha)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(canonical_json(body))
    GATE_CALLS.append(p.name)
    return body


def synth_env(recompute=fake_record):
    """Synthetic VerificationEnv: no real data, pipeline, history or frozen driver."""
    return postrun.VerificationEnv(
        load_context=lambda root: SimpleNamespace(e=np.zeros(8)),
        pipeline=lambda: "synthetic-pipeline",
        recompute=recompute,
        load_historical=lambda root, r: {},
        load_frozen_driver=lambda root: SimpleNamespace(name="synthetic-frozen-driver"),
        streams=lambda c, r: SimpleNamespace(run_id=r))


class Checks:
    """Doubles for verification.py's record checks; statuses configurable per check."""
    def __init__(self, monkeypatch):
        self.status = {}
        self.extra = {}
        self.calls = []
        for fn, name in (("v1_self_consistency", "V1-self"), ("v1_outer_vs_historical", "V1-outer-historical"),
                         ("v1_inner_vs_frozen_path", "V1-inner"), ("v2_record_assignment", "V2-record"),
                         ("v3_m1_invariance", "V3"), ("v4_as_run_apparent", "V4")):
            monkeypatch.setattr(V, fn, self.make(name))

    def make(self, name):
        def f(records, *a, **kw):
            self.calls.append(name)
            st = self.status.get(name, "PASS")
            if isinstance(st, Exception):
                raise st
            n = len(records)
            if name == "V4" and name not in self.status and a[1] != "as_run" and a[2] is not True:
                return {"check": name, "status": "NOT_APPLICABLE",         # as the real V4
                        "reason": "independent primary without the as-run variant"}
            if name == "V4":                   # count comparable records as the real V4 does
                hist, field = a[0], ("dC_app" if a[1] == "as_run" else "dC_app_as_run")
                n = sum(1 for r in records if (hist.get(r["k"]) or {}).get("status") == "ok" and field in r)
            return {"check": name, "status": st, "n_compared": n, **self.extra.get(name, {})}
        return f


@pytest.fixture
def camp(monkeypatch):
    c = Campaign(f"scratch_p9c_post_{uuid.uuid4().hex[:10]}")
    S.update(drift=0.0, raise_at=None, perturb=set(), nan=False)
    monkeypatch.setattr(orch, "h4_check", lambda root, stage, raise_on_fail=True: {"all_ok": True})
    monkeypatch.setattr(orch, "historical_preservation_check", lambda root: {"all_ok": True, "status": "PASS"})
    monkeypatch.setattr(orch, "build_manifest", lambda config, root: {
        "schema": "test-manifest", "config": config.sha256(), "code": {"x": "0" * 64}})
    monkeypatch.setattr(orch, "run_pre_execution_gate", synth_gate)
    monkeypatch.setattr(orch, "load_context", lambda root: None)
    monkeypatch.setattr(orch, "RealPipeline", lambda: None)
    monkeypatch.setattr(orch, "compute_outer_replicate", fake_record)
    c.identity_calls = []

    def identity(config, manifest=None, extra_files=()):
        c.identity_calls.append(manifest)
        if manifest is not None and manifest.get("code") != SYNTH_FREEZE["freeze"]["code"]:
            raise orch.CodeIdentityRefusal("code identity not verified: manifest code hash differs")
        return SYNTH_FREEZE
    monkeypatch.setattr(orch, "verify_code_identity", identity)
    monkeypatch.setattr(orch, "verify_loaded_modules", lambda root, code, required: None)

    def no_real(*a, **k):
        raise AssertionError("a test attempted to use the REAL verification environment")
    monkeypatch.setattr(postrun.VerificationEnv, "real", classmethod(no_real))
    yield c
    shutil.rmtree(c.out, ignore_errors=True)
    assert not c.out.exists()


@pytest.fixture
def checks(monkeypatch):
    return Checks(monkeypatch)


def verify(c, env=None, **extra):
    return postrun.run_verification(c.cfg(**extra), env or synth_env())


def by_check(rep, name):
    return [r for r in rep["results"] if r["check"] == name]


# ================================================================ end-to-end PASS and outputs
def test_clean_campaign_pass_order_counts_and_create_only_outputs(camp, checks):
    camp.all_runs()
    rep = verify(camp)
    assert rep["verdict"] == "PASS" and rep["reporting_permitted"] and not rep["hard_failure_h5"]
    assert rep["level"] == 500 and rep["status_counts"]["FAIL"] == 0
    assert checks.calls[:6] == ["V1-self", "V1-outer-historical", "V1-inner", "V2-record", "V3", "V4"]
    order = [r["check"] for r in rep["results"] if r["run_id"] == 0 and r["check"] in
             postrun.RECORD_CHECKS + ("V5",)]
    assert order == list(postrun.RECORD_CHECKS) + ["V5"]                    # ordered sequence
    assert all(set(r) >= {"status", "counts", "reasons"} for r in rep["results"])
    sz = by_check(rep, "sizing")
    assert len(sz) == 1 and sz[0]["status"] == "PASS" and sz[0]["evidence"]["action"] == "stop_ep1_pass"
    j, t = postrun.report_paths(camp.out, 500)
    assert json.loads(j.read_bytes()) == json.loads(canonical_json(rep))
    assert "verdict: PASS" in t.read_text()
    # manifest-to-freeze comparison happened, with the stored manifest
    stored = json.loads((camp.out / "manifest/input_manifest.json").read_text())
    assert stored in camp.identity_calls


def test_verification_never_starts_scientific_execution(camp, checks, monkeypatch):
    camp.all_runs()
    before = camp.snapshot()

    def forbidden(name):
        def f(*a, **k):
            raise AssertionError(f"verification called {name}")
        return f
    for name in ("run_corrected", "execute_replicates", "build_sizing_decision", "build_manifest",
                 "write_manifest_once", "run_pre_execution_gate", "compute_outer_replicate",
                 "resolve_authorization", "CheckpointStore", "load_context", "RealPipeline"):
        monkeypatch.setattr(orch, name, forbidden(name))
    rep = verify(camp)
    assert rep["verdict"] == "PASS"
    after = camp.snapshot()
    changed = {k for k in before if after.get(k) != before[k]}
    assert not changed                                                     # nothing modified
    new = set(after) - set(before)
    assert new and all(k.startswith(("verification/", "windows/")) for k in new)
    assert not any(k.startswith(("runs/", "halts/", "manifest/", "sizing/")) for k in new)


def test_runner_source_has_no_dispatch_or_write_paths():
    src = Path(postrun.__file__).read_text()
    for name in ("run_corrected(", "execute_replicates(", "build_sizing_decision(", "CheckpointStore(",
                 "write_manifest_once(", "build_manifest(", "run_pre_execution_gate("):
        assert name not in src, name


def test_report_collision_is_refused_and_never_replaced(camp, checks):
    camp.all_runs()
    verify(camp)
    j, t = postrun.report_paths(camp.out, 500)
    saved = j.read_bytes(), t.read_bytes()
    with pytest.raises(postrun.VerificationRefused, match="never replaced"):
        verify(camp)
    assert (j.read_bytes(), t.read_bytes()) == saved
    t.unlink()                                       # an incomplete prior output: still refused
    with pytest.raises(postrun.VerificationRefused):
        verify(camp)


# ================================================================ P4 normalisation
def test_v4_not_applicable_only_under_the_approved_condition(camp, checks):
    camp.all_runs()
    checks.status["V4"] = "NOT_APPLICABLE"
    rep = verify(camp)                                # as_run_variant false, independent
    v4 = by_check(rep, "V4")
    assert {x["status"] for x in v4} == {"NOT_APPLICABLE"} and rep["verdict"] == "PASS"
    assert "as_run_variant = false" in v4[0]["evidence"]["approved_condition"]


def test_v4_not_applicable_is_fail_when_as_run_variant_true(camp, checks):
    camp.over = {"as_run_variant": True}             # the campaign's own approved configuration
    camp.all_runs()
    checks.status["V4"] = "NOT_APPLICABLE"           # claims inapplicability although comparable
    env = synth_env()
    env.load_historical = lambda root, r: hist_rows()
    rep = verify(camp, env=env)
    assert {x["status"] for x in by_check(rep, "V4")} == {"FAIL"}
    assert rep["verdict"] == "FAIL" and not rep["reporting_permitted"]


@pytest.mark.parametrize("check,raw,why", [
    ("V1-outer-historical", "NOT_APPLICABLE", "inapplicability is not established"),
    ("V3", "NOT_APPLICABLE", "inapplicability is not established"),
    ("V1-inner", "BLOCKED", "evidence missing or malformed"),
    ("V2-record", ValueError("broken record"), "evidence missing or malformed"),
    ("V1-self", "SOMETHING_ELSE", "evidence missing or malformed")])
def test_missing_malformed_or_unapproved_evidence_is_fail(camp, checks, check, raw, why):
    camp.all_runs()
    checks.status[check] = raw
    rep = verify(camp)
    res = by_check(rep, check)
    assert {x["status"] for x in res} == {"FAIL"} and why in res[0]["reasons"][0]
    assert rep["verdict"] == "FAIL"


def test_pass_with_no_evidence_is_fail(camp, checks):
    camp.all_runs()
    checks.extra["V1-inner"] = {"n_compared": 0}
    rep = verify(camp)
    assert {x["status"] for x in by_check(rep, "V1-inner")} == {"FAIL"}


def test_not_comparable_is_never_pass(camp, checks):
    camp.all_runs()
    checks.extra["V3"] = {"n_not_comparable": 3,
                          "not_comparable": [(1, "paired set differs"), (2, "x"), (3, "y")]}
    rep = verify(camp)
    nc = by_check(rep, "V3:not_comparable")
    assert len(nc) == 5 and all(x["status"] == "NOT_COMPARABLE" and x["counts"]["n_not_comparable"] == 3
                                for x in nc)
    assert rep["verdict"] == "OWNER_REVIEW_REQUIRED" and not rep["reporting_permitted"]


def test_normalise_unit_rules():
    cfg = CorrectedConfig(approved_test_config(ROOT, "scratch_unused"), ROOT)
    n = postrun.normalise
    assert n({"check": "V3", "status": "PASS", "n_compared": 4}, cfg)[0]["status"] == "PASS"
    assert n({"check": "V3", "status": "PASS", "n_compared": 0}, cfg)[0]["status"] == "FAIL"
    assert n({"check": "D12-lambda", "status": "PASS", "n_fits": 0}, cfg)[0]["status"] == "FAIL"
    assert n({"check": "V4", "status": "BLOCKED", "reason": "x"}, cfg)[0]["status"] == "FAIL"
    assert n({"check": "V5", "status": "NOT_APPLICABLE"}, cfg)[0]["status"] == "FAIL"
    assert n({"check": "V4", "status": "NOT_APPLICABLE"}, cfg)[0]["status"] == "NOT_APPLICABLE"
    assert n({"check": "V1-self", "status": None}, cfg)[0]["status"] == "FAIL"
    assert set(postrun.APPROVED_NOT_COMPARABLE) == {"V4"}                 # D2, owner-approved
    assert set(postrun.APPROVED_NOT_APPLICABLE) == {"V4"}


# ================================================================ V5 / H5
def test_v5_failure_is_hard_h5_and_blocks_reporting(camp, checks):
    camp.all_runs()

    def recompute(*a):
        rec = fake_record(*a)
        if rec["run_id"] == 0 and rec["k"] == 25:
            rec["dC_corrected"] += 1e-12
        return rec
    with pytest.raises(postrun.HardFailureH5):
        verify(camp, env=synth_env(recompute))
    rep = json.loads(postrun.report_paths(camp.out, 500)[0].read_bytes())
    assert rep["hard_failure_h5"] and rep["verdict"] == "FAIL" and not rep["reporting_permitted"]
    assert by_check(rep, "V5")[0]["status"] == "FAIL"


def test_v5_missing_selected_replicate_is_fail_not_na(camp, checks, monkeypatch):
    camp.all_runs()
    monkeypatch.setattr(V, "V5_SELECTED_K", (0, 25, 9999))
    real = V.v5_spot_check
    monkeypatch.setattr(V, "v5_spot_check", lambda stored, rc, selected_k=(0, 25, 9999): real(stored, rc, selected_k))
    rep = verify(camp)
    assert by_check(rep, "V5")[0]["status"] == "FAIL" and rep["verdict"] == "FAIL"


# ================================================================ Stage A / Stage B (DEC-C1)
def windowed(camp):
    """Run 2 crashes at k=200 (UNTERMINATED), is acknowledged and completed: one window."""
    for r in (0, 1, 3, 4):
        camp.run(run_id=r)
    S["raise_at"] = 200
    with pytest.raises(RuntimeError):
        camp.run(run_id=2)
    S["raise_at"] = None
    camp.run(run_id=2, resume_authorization="ACK-2")
    w = orch.assess_run(camp.cfg(), 2, camp.msha())["windows"]
    assert len(w) == 1
    return w[0]


def test_stage_a_recomputes_window_writes_bound_result_and_enables_eligibility(camp, checks):
    w = windowed(camp)
    assert orch.evaluate_eligibility(camp.cfg(), 2, 500, camp.msha())["class"] == "UNVERIFIABLE"
    rep = verify(camp)
    sa = by_check(rep, "StageA-window")
    assert len(sa) == 1 and sa[0]["status"] == "PASS" and sa[0]["counts"]["n_identical"] == len(w["k_list"])
    res = json.loads(orch.window_result_path(camp.out, 2, w["invocation"]).read_text())
    assert res["status"] == "PASS" and res["region_sha256"] == w["region_sha256"]
    assert orch.evaluate_eligibility(camp.cfg(), 2, 500, camp.msha())["class"] == "ELIGIBLE"
    assert rep["verdict"] == "PASS"


def test_stage_a_mismatch_fails_disqualifies_and_skips_record_checks(camp, checks):
    w = windowed(camp)
    S["perturb"] = {(2, w["k_list"][3])}
    rep = verify(camp)
    assert by_check(rep, "StageA-window")[0]["status"] == "FAIL"
    elig = [x for x in by_check(rep, "StageB-eligibility") if x["run_id"] == 2][0]
    assert elig["status"] == "FAIL" and elig["evidence"]["class"] == "INELIGIBLE"
    assert checks.calls == []                                             # no record checks ran
    assert {x["status"] for x in by_check(rep, "V1-self")} == {"NOT_COMPARABLE"}
    assert rep["verdict"] == "FAIL"


def test_stored_pass_window_file_is_not_trusted(camp, checks):
    w = windowed(camp)
    forged = {"schema": "p9c-window-1", "run_id": 2, "invocation": w["invocation"],
              "region_sha256": w["region_sha256"], "k_list": w["k_list"],
              "code_identity_verified": True, "status": "PASS",
              "per_k": {str(k): {"record_identical": True, "dC_hex_equal": True} for k in w["k_list"]}}
    p = orch.window_result_path(camp.out, 2, w["invocation"])
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(forged))
    saved = p.read_bytes()
    S["perturb"] = {(2, w["k_list"][0])}                                   # recomputation differs
    rep = verify(camp)
    sa = by_check(rep, "StageA-window")[0]
    assert sa["status"] == "FAIL" and "not replaced" in sa["reasons"][0]
    assert p.read_bytes() == saved and rep["verdict"] == "FAIL"


def test_identical_existing_window_result_is_accepted(camp, checks):
    windowed(camp)
    verify(camp)
    for q in postrun.report_paths(camp.out, 500):
        q.unlink()                                   # scratch only: allow a second verification
    rep = verify(camp)
    assert by_check(rep, "StageA-window")[0]["status"] == "PASS" and rep["verdict"] == "PASS"


def test_incomplete_run_is_not_eligible_and_fails(camp, checks):
    for r in range(4):
        camp.run(run_id=r)
    with pytest.raises(orch.RuntimeBudgetHalt):
        orch.run_corrected(camp.cfg(), 4, 500, clock=("TEST_CLOCK", StepClock(600)))
    rep = verify(camp)
    e4 = [x for x in by_check(rep, "StageB-eligibility") if x["run_id"] == 4][0]
    assert e4["status"] == "FAIL" and e4["evidence"]["class"] == "INCOMPLETE"
    assert rep["verdict"] == "FAIL"


def test_identity_binding_of_every_invocation(camp, checks):
    camp.all_runs()
    p = camp.out / "halts/run3.inv0001.start.json"
    o = json.loads(p.read_text())
    o["freeze_record_sha256"] = "e" * 64
    p.write_text(json.dumps(o))
    rep = verify(camp)
    ib = [x for x in by_check(rep, "identity-binding") if x["run_id"] == 3][0]
    assert ib["status"] == "FAIL" and rep["verdict"] == "FAIL"


# ================================================================ sizing (approved logic + independent oracle)
def test_sizing_requiring_topup_is_fail_and_matches_independent_oracle(camp, checks):
    camp.over = {"w5": W5A, "ladder": LADDER_FULL}
    S["drift"] = 0.05
    camp.all_runs()
    rep = verify(camp)
    sz = by_check(rep, "sizing")[0]
    by = {r: {x["k"]: x for x in orch.assess_run(camp.cfg(), r, camp.msha())["records"]} for r in range(5)}
    exp = sizing_oracle(by, 500)                                          # independent derivation
    assert sz["evidence"]["action"] == exp["action"] == "top_up"
    assert sz["evidence"]["combined_B_next"] == exp["combined"] == 750
    assert math.isclose(sz["evidence"]["rho_max"], exp["rho_max"], rel_tol=1e-10)
    assert sz["status"] == "FAIL" and "not complete" in sz["reasons"][0] and rep["verdict"] == "FAIL"


def test_topped_up_campaign_reproduces_every_level(camp, checks):
    camp.over = {"w5": W5A, "ladder": LADDER_750}
    S["drift"] = 0.05
    camp.all_runs()
    camp.sizing(500)
    for r in range(5):
        camp.run(run_id=r, target=750, topup_authorization="OWNER-750")
    rep = verify(camp)
    assert rep["level"] == 750
    sz = by_check(rep, "sizing")
    assert [x["counts"]["level"] for x in sz] == [500, 750] and all(x["status"] == "PASS" for x in sz)
    assert sz[1]["evidence"]["action"] == "report_persisting_breach"
    assert by_check(rep, "sizing-file")[0]["status"] == "PASS"
    assert rep["verdict"] == "PASS"


def test_runs_at_different_levels_fail_cross_run_compliance(camp, checks):
    camp.over = {"w5": W5A, "ladder": LADDER_750}
    S["drift"] = 0.05
    camp.all_runs()
    camp.sizing(500)
    camp.run(run_id=0, target=750, topup_authorization="OWNER-750")
    rep = verify(camp)
    assert by_check(rep, "campaign-level")[0]["status"] == "FAIL" and rep["level"] is None
    assert rep["verdict"] == "FAIL"


def test_doctored_sizing_file_is_detected(camp, checks, tmp_path):
    camp.over = {"w5": W5A, "ladder": LADDER_750}
    S["drift"] = 0.05
    camp.all_runs()
    art = camp.sizing(500)
    cfg = camp.cfg()
    cps = {r: (camp.out / f"runs/run{r}.jsonl").read_bytes() for r in range(5)}
    dec = orch.compute_sizing_decision(cfg, 500, {r: orch._records_by_k(cps[r]) for r in range(5)})
    ok = postrun._verify_sizing_file(cfg, 500, orch._sizing_path(camp.out, 500), dec, camp.msha(),
                                     SYNTH_FREEZE, cps)
    assert ok[0]["status"] == "PASS"
    bad = json.loads(json.dumps(art))
    bad["decision"]["combined_B_next"] = 1000
    p = tmp_path / "doctored.json"
    p.write_text(json.dumps(bad))
    res = postrun._verify_sizing_file(cfg, 500, p, dec, camp.msha(), SYNTH_FREEZE, cps)
    assert res[0]["status"] == "FAIL" and "differs from the recomputation" in " ".join(res[0]["reasons"])
    cps_bad = dict(cps)
    cps_bad[1] = b"x" + cps[1][1:]
    res2 = postrun._verify_sizing_file(cfg, 500, orch._sizing_path(camp.out, 500), dec, camp.msha(),
                                       SYNTH_FREEZE, cps_bad)
    assert "run 1 snapshot is not a prefix" in " ".join(res2[0]["reasons"])


# ================================================================ preconditions
def test_refusals_before_any_computation(camp, checks):
    called = []
    env = synth_env()
    env.load_context = lambda root: called.append("ctx")
    with pytest.raises(postrun.VerificationRefused, match="does not exist"):
        verify(camp, env=env)                                            # no output_dir
    camp.all_runs()
    with pytest.raises(postrun.VerificationRefused, match="not executable"):
        verify(camp, env=env, execution_authorization={"authorized": False, "date": None, "reference": None})
    (camp.out / "manifest/input_manifest.json").rename(camp.out / "manifest/moved.json")
    with pytest.raises(postrun.VerificationRefused, match="manifest"):
        verify(camp, env=env)
    (camp.out / "manifest/moved.json").rename(camp.out / "manifest/input_manifest.json")
    assert called == []
    assert not (camp.out / "verification").exists()


def test_manifest_not_matching_freeze_is_refused(camp, checks, monkeypatch):
    camp.all_runs()
    p = camp.out / "manifest/input_manifest.json"
    m = json.loads(p.read_text())
    m["code"] = {"x": "1" * 64}
    p.write_text(json.dumps(m))
    with pytest.raises(orch.CodeIdentityRefusal, match="manifest"):
        verify(camp)


# ================================================================ Decision 12 statuses
def fit_fail(model, status, lam, score, ftype, consumed=True):
    return {"ok": False, "lambda_status": status, "lambda_at_failure": lam,
            "tuning_score_at_failure": score, "failure_type": ftype, "cv_rng_consumed": consumed}


def d12_record(fd, model="M2", stage="apparent", mirror=True):
    rec = {"k": 0, "apparent": {"M1": {"ok": True}, model: fd} if model == "M2" else {model: fd, "M2": {"ok": True}},
           "inner": [], "failures": []}
    if mirror:
        rec["failures"].append({"stage": stage, "model": model, "b": None,
                                "lambda_status": fd.get("lambda_status"),
                                "lambda_at_failure": fd.get("lambda_at_failure"),
                                "tuning_score_at_failure": fd.get("tuning_score_at_failure")})
    return rec


G0 = float(LAM_GRID[0])


@pytest.mark.parametrize("fd,model,mirror,ok", [
    (fit_fail("M2", "captured", G0, float("-inf"), "ConvergenceFailure"), "M2", True, True),
    (fit_fail("M2", "captured", 0.123, 0.6, "ConvergenceFailure"), "M2", True, False),      # off grid
    (fit_fail("M2", "captured", G0, 0.6, "ConvergenceFailure", consumed=False), "M2", True, False),
    (fit_fail("M1", "captured", G0, 0.6, "ConvergenceFailure"), "M1", True, False),         # M1 never tunes
    (fit_fail("M2", "not_applicable_before_tuning", None, None, "NonEstimableReplicate"), "M2", True, True),
    (fit_fail("M2", "not_applicable_before_tuning", 0.0, None, "NonEstimableReplicate"), "M2", True, False),
    (fit_fail("M1", "fixed_zero_no_molecular_block", 0.0, None, "ConvergenceFailure"), "M1", True, True),
    (fit_fail("M1", "fixed_zero_no_molecular_block", 0.0, None, "NonEstimableReplicate"), "M1", True, False),
    (fit_fail("M2", None, None, None, "ConvergenceFailure"), "M2", True, False),            # capture missing
    (fit_fail("M2", "captured", G0, 0.6, "ConvergenceFailure"), "M2", False, False)])       # not mirrored
def test_d12_lambda_status_consistency(fd, model, mirror, ok):
    res = postrun.d12_lambda_check([d12_record(fd, model, mirror=mirror)])
    assert (res["status"] == "PASS") is ok, res["mismatches"]


def test_d12_success_with_status_and_orphan_failure_entry_fail():
    rec = {"k": 1, "apparent": {"M1": {"ok": True, "lambda_status": "captured"}, "M2": {"ok": True}},
           "inner": [], "failures": []}
    assert postrun.d12_lambda_check([rec])["status"] == "FAIL"
    rec2 = {"k": 1, "apparent": {"M1": {"ok": True}, "M2": {"ok": True}}, "inner": [],
            "failures": [{"stage": "inner", "model": "M2", "b": 3}]}
    assert postrun.d12_lambda_check([rec2])["status"] == "FAIL"


# ================================================================ real check functions on synthetic records
def test_real_record_checks_on_synthetic_records_pass_and_detect_tampering():
    ctx = synthetic_context(n=40)
    cfg = CorrectedConfig(approved_test_config(ROOT, "scratch_unused"), ROOT)
    streams = PlanStreams(cfg.base, 0, cfg.convention, cfg.T)
    recs = [json.loads(canonical_json(compute_outer_replicate(ctx, streams, FakePipeline(), 0, k, 2, LAM_GRID)))
            for k in range(3)]
    for raw in (V.v1_self_consistency(recs, ctx.e, streams), V.v2_record_assignment(recs, streams)):
        assert postrun.normalise(raw, cfg, run_id=0)[0]["status"] == "PASS"
    recs[1]["outer"]["n_distinct"] += 1
    assert postrun.normalise(V.v1_self_consistency(recs, ctx.e, streams), cfg)[0]["status"] == "FAIL"
    # the fake pipeline has no Decision 12 capture: failed fits carry no status -> FAIL
    fails = compute_outer_replicate(ctx, streams, FakePipeline(fail_when=lambda ids, panel: "x" if panel else None),
                                    0, 7, 2, LAM_GRID)
    assert postrun.d12_lambda_check([fails])["status"] == "FAIL"


# ================================================================ DEC-A1 rendering
def test_nonfinite_values_rendered_explicitly(camp, checks):
    camp.all_runs()
    checks.extra["V3"] = {"worst_gap": float("-inf"), "ratio": float("nan"), "top": float("inf")}
    rep = verify(camp)
    assert rep["nonfinite_values"]["count"] == 15                         # 3 fields x 5 runs
    j, t = postrun.report_paths(camp.out, 500)
    raw = j.read_bytes()
    assert b"-Infinity" in raw and b"NaN" in raw
    txt = t.read_text()
    assert "-Infinity" in txt and "+Infinity" in txt and "NaN" in txt
    assert "non-finite values in this report (DEC-A1): 15" in txt
    assert [postrun.render_value(x) for x in (float("nan"), float("inf"), float("-inf"), None, 1.5)] == \
        ["NaN", "+Infinity", "-Infinity", "null", "1.5"]


# ================================================================ entry point exit codes
@pytest.mark.parametrize("outcome,code", [("PASS", 0), ("FAIL", 1), ("OWNER_REVIEW_REQUIRED", 3),
                                          ("H5", 5), ("H4", 4), ("REFUSED", 2)])
def test_main_exit_codes(monkeypatch, tmp_path, outcome, code):
    from src.phase9_corrected.integrity import HardFailureH4
    cfgp = tmp_path / "c.json"
    cfgp.write_text(json.dumps(approved_test_config(ROOT, "scratch_unused")))

    def fake(cfg, env=None):
        if outcome == "H5":
            raise postrun.HardFailureH5("h5")
        if outcome == "H4":
            raise HardFailureH4("h4")
        if outcome == "REFUSED":
            raise postrun.VerificationRefused("no")
        return {"verdict": outcome, "level": 500, "status_counts": {}}
    monkeypatch.setattr(postrun, "run_verification", fake)
    assert postrun.main([str(cfgp)], ROOT) == code
    assert postrun.main([], ROOT) == 2


# ================================================================ launcher `verify` (temporary repository, stubs)
STUB = r'''import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
(ROOT / "{marker}").write_text(json.dumps({{"argv": sys.argv[1:], "isolated": sys.flags.isolated,
    "no_bytecode": sys.flags.dont_write_bytecode, "prefix": bool(sys.pycache_prefix)}}))
'''
FROZEN_NONPKG = list(orch.REQUIRED_FILES)


@pytest.fixture
def vrepo(tmp_path):
    root = tmp_path / "repo"
    for rel in [f"src/phase9_corrected/{p.name}" for p in (ROOT / "src/phase9_corrected").glob("*.py")] + FROZEN_NONPKG:
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, root / rel)
    (root / "scripts/41_phase9_corrected_driver.py").write_text(STUB.format(marker="DRIVER_STARTED"))
    (root / "scripts/43_phase9_corrected_verify.py").write_text(STUB.format(marker="VERIFY_STARTED"))
    os.symlink(ROOT / ".venv", root / ".venv")

    def freeze(drop=None):
        files = sorted([f"src/phase9_corrected/{p.name}" for p in (root / "src/phase9_corrected").glob("*.py")]
                       + FROZEN_NONPKG)
        fr = {"schema": orch.FREEZE_SCHEMA, "bytecode_policy": "fresh_compile_prefix", "environment": {},
              "code": {f: hashlib.sha256((root / f).read_bytes()).hexdigest() for f in files if f != drop}}
        fp = root / orch.FREEZE_RECORD_RELPATH
        fp.parent.mkdir(parents=True, exist_ok=True)
        fp.write_bytes(json.dumps(fr, sort_keys=True).encode())
        cfg = approved_test_config(root, "out")
        cfg["owner_decisions_outside_schema"] = {"D13_driver_freeze": {
            "approved": True, "date": "2026-11-06", "reference": "TEST FIXTURE",
            "driver_path": orch.DRIVER_RELPATH,
            "sha256": fr["code"][orch.DRIVER_RELPATH],
            "freeze_record_sha256": hashlib.sha256(fp.read_bytes()).hexdigest()}}
        (tmp_path / "owner.json").write_text(json.dumps(cfg))
    freeze()
    return tmp_path, root, freeze


def launch(base, root, *args):
    for m in ("DRIVER_STARTED", "VERIFY_STARTED"):
        (root / m).unlink(missing_ok=True)
    t = base / f"t_{uuid.uuid4().hex[:6]}"
    t.mkdir()
    return subprocess.run(["/bin/bash", str(root / "scripts/42_phase9_corrected_launch.sh"), *args],
                          capture_output=True, text=True, cwd=str(base),
                          env={"PATH": "/usr/bin:/bin", "HOME": str(base), "TMPDIR": str(t)})


def test_launcher_verify_starts_only_the_runner_with_the_same_flags(vrepo):
    base, root, _ = vrepo
    r = launch(base, root, "verify", str(base / "owner.json"))
    assert r.returncode == 0, r.stderr
    assert not (root / "DRIVER_STARTED").exists()
    o = json.loads((root / "VERIFY_STARTED").read_text())
    assert o["argv"] == [str((base / "owner.json").resolve())]
    assert o["isolated"] == 1 and o["no_bytecode"] == 1 and o["prefix"]
    r2 = launch(base, root, "validate", str(base / "owner.json"))
    assert r2.returncode == 0 and (root / "DRIVER_STARTED").exists() and not (root / "VERIFY_STARTED").exists()


def test_launcher_verify_refusals_before_python(vrepo):
    base, root, freeze = vrepo
    r = launch(base, root, "verify", str(base / "owner.json"), "--extra")
    assert r.returncode == 2 and not (root / "VERIFY_STARTED").exists()
    with open(root / "scripts/43_phase9_corrected_verify.py", "a") as fh:
        fh.write("# modified after the freeze\n")
    r = launch(base, root, "verify", str(base / "owner.json"))
    assert r.returncode == 16 and not (root / "VERIFY_STARTED").exists()
    freeze(drop="scripts/43_phase9_corrected_verify.py")                   # runner not frozen
    r = launch(base, root, "verify", str(base / "owner.json"))
    assert r.returncode == 15 and "required file is not in the freeze inventory" in r.stderr
    freeze(drop="scripts/40_phase9_dc_interval.py")
    r = launch(base, root, "validate", str(base / "owner.json"))
    assert r.returncode == 15 and not (root / "DRIVER_STARTED").exists()


def test_launcher_required_list_equals_required_files():
    txt = (ROOT / "scripts/42_phase9_corrected_launch.sh").read_text()
    block = txt[txt.index("REQUIRED_REL=("):].split(")", 1)[0]
    names = {"$DRIVER_REL": orch.DRIVER_RELPATH, "$SELF_REL": orch.LAUNCHER_RELPATH,
             "$VERIFY_REL": orch.VERIFY_RELPATH}
    items = [names.get(x.strip('"'), x.strip('"')) for x in block[len("REQUIRED_REL=("):].split()]
    assert sorted(items) == sorted(orch.REQUIRED_FILES)



# ================================================================ review remediation (second cycle)
# ---- D1: identity re-check before every output write
def test_d1_identity_guard_precedes_every_output_write(camp, checks, monkeypatch):
    windowed(camp)
    log = []
    real_guard = postrun.identity_guard
    monkeypatch.setattr(postrun, "identity_guard", lambda root, code, what: (log.append(("guard", what)),
                                                                              real_guard(root, code, what)))
    real_dc, real_cx = orch.durable_create, postrun._create_exclusive
    monkeypatch.setattr(orch, "durable_create", lambda path, obj: (log.append(("write", Path(path).name)),
                                                                   real_dc(path, obj)))
    monkeypatch.setattr(postrun, "_create_exclusive", lambda path, data: (log.append(("write", Path(path).name)),
                                                                          real_cx(path, data)))
    verify(camp)
    writes = [i for i, e in enumerate(log) if e[0] == "write"]
    assert len(writes) == 3                                   # one window result + two reports
    for i in writes:
        assert log[i - 1][0] == "guard", log                  # immediately preceded by a re-check
    assert [log[i][1] for i in writes][1:] == ["L00500.report.json", "L00500.report.txt"]


def test_d1_module_loaded_during_recomputation_blocks_the_window_write(camp, checks, monkeypatch):
    w = windowed(camp)
    marker = "zz_p9c_lazily_loaded_shadow"

    def recompute(*a):
        sys.modules[marker] = types.ModuleType(marker)        # e.g. a lazily imported shadow
        return fake_record(*a)

    def lm(root, code, required):
        if marker in sys.modules:
            raise orch.CodeIdentityRefusal(f"loaded repository-local module {marker} is not a required inventory file")
    monkeypatch.setattr(orch, "verify_loaded_modules", lm)
    try:
        with pytest.raises(postrun.VerificationRefused, match="before writing Stage A window result"):
            verify(camp, env=synth_env(recompute))
    finally:
        sys.modules.pop(marker, None)
    assert not orch.window_result_path(camp.out, 2, w["invocation"]).exists()  # never written as verified
    assert not (camp.out / "verification").exists()


def test_d1_identity_check_errors_fail_closed(camp, checks, monkeypatch):
    camp.all_runs()

    def broken(root, code, required):
        raise RuntimeError("unexpected internal error")
    monkeypatch.setattr(orch, "verify_loaded_modules", broken)
    with pytest.raises(postrun.VerificationRefused, match="RuntimeError"):
        verify(camp)
    assert not (camp.out / "verification").exists()


def test_d1_root_shadow_detection_and_namespace_limitation(tmp_path):
    root = tmp_path / "r"
    for d, init in (("src", True), ("tests", True), ("data", False), ("multiprocessing", True)):
        (root / d).mkdir(parents=True)
        if init:
            (root / d / "__init__.py").write_text("")
    (root / "lifelines.py").write_text("")
    (root / "fast.cpython-311-darwin.so").write_bytes(b"")
    (root / "README.md").write_text("")
    (root / "pytest.ini").write_text("")
    assert postrun.root_shadow_entries(root) == ["fast.cpython-311-darwin.so", "lifelines.py", "multiprocessing"]
    with pytest.raises(postrun.VerificationRefused, match="shadow imports"):
        postrun.identity_guard(root, {}, "test output")
    # documented limitation: a namespace directory (no __init__) is not reported -- it cannot
    # shadow a regular module that exists later on sys.path
    assert "data" not in postrun.root_shadow_entries(root)
    assert postrun.root_shadow_entries(ROOT) == []                       # the real repository is clean


# ---- D3: gate evidence (inspection only)
def gate_file(camp, r=0, m=1):
    st = json.loads((camp.out / f"halts/run{r}.inv{m:04d}.start.json").read_text())
    return camp.out / st["gate_evidence"]["path"], st


def damage_body(g, damage):
    if damage == "extra_check":
        g["V9"] = "PASS"
    elif damage == "drop_check":
        del g["V2_numpy"]
    elif damage == "numpy_key":
        g["key"]["numpy"] = "0.0.0"
    elif damage == "config_key":
        g["key"]["config_sha256"] = "0" * 64
    elif damage == "manifest_key":
        g["key"]["manifest_sha256"] = "0" * 64
    elif damage == "h4_fail":
        g["H4"] = "FAIL"
    elif damage == "domain_fail":
        g["V2_domain"]["status"] = "FAIL"
    elif damage == "domain_bmax":
        g["V2_domain"]["domain"]["b_max"] = 50
    elif damage == "domain_count":
        g["V2_domain"]["n_states"] -= 1
    elif damage == "collision":
        g["V2_domain"]["collisions"] = [["a", "b"]]
    elif damage == "v1_short":
        g["V1_outer_historical"]["n_compared"] = 499
    elif damage == "all_pass":
        g["all_pass"] = False
    elif damage == "structural":
        g["V2_structural"]["detail"]["holds"] = False
    elif damage == "status_missing":
        del g["V1_outer_historical"]["status"]
    return g


GATE_DAMAGE = [("extra_check", "exactly the approved checks"), ("drop_check", "exactly the approved checks"),
               ("numpy_key", "key does not match"), ("config_key", "key does not match"),
               ("manifest_key", "key does not match"), ("h4_fail", "H4 is not PASS"),
               ("domain_fail", "V2_domain is not PASS"), ("domain_bmax", "V2_domain domain"),
               ("domain_count", "V2_domain counts"), ("collision", "V2_domain counts"),
               ("v1_short", "V1_outer_historical does not cover"), ("all_pass", "all_pass is not true"),
               ("structural", "V2_structural detail"), ("status_missing", "V1_outer_historical is not PASS")]


@pytest.mark.parametrize("damage,match", GATE_DAMAGE)
def test_gate_evidence_content_rules(camp, checks, damage, match):
    """Content validation, independent of the SHA-256 binding."""
    camp.all_runs()
    p, _ = gate_file(camp)
    raw = canonical_json(damage_body(json.loads(p.read_bytes()), damage))
    probs = orch.gate_evidence_problems(raw, camp.cfg(), camp.msha())
    assert probs and any(match in x for x in probs), probs
    assert orch.gate_evidence_problems(p.read_bytes(), camp.cfg(), camp.msha()) == []


@pytest.mark.parametrize("raw,match", [(None, "missing or unreadable"), (b"{not json", "not valid JSON"),
                                       (b"[]", "not a JSON object"), (b"", "not valid JSON")])
def test_gate_evidence_malformed_is_structured_never_raises(camp, checks, raw, match):
    camp.all_runs()
    probs = orch.gate_evidence_problems(raw, camp.cfg(), camp.msha())
    assert any(match in x for x in probs)
    cfg_other = camp.cfg(verification_gate="something_else")
    assert "verification_gate" in orch.gate_evidence_problems(b"{}", cfg_other, camp.msha())[0]


def test_every_invocation_binds_its_own_verified_gate_evidence(camp, checks):
    camp.all_runs()
    rep = verify(camp)
    ge = by_check(rep, "gate-evidence")
    assert len(ge) == 5 and all(x["status"] == "PASS" for x in ge)
    for r in range(5):
        p, st = gate_file(camp, r)
        assert st["gate_evidence"]["path"] == f"gate/run{r}.inv0001.attempt01.json"
        assert st["gate_evidence"]["sha256"] == hashlib.sha256(p.read_bytes()).hexdigest()
        assert st["gate_evidence"]["all_pass"] is True
    assert "not cryptographically authenticated" in ge[0]["reasons"][0]


# ---- D2: V4 skipped records -- owner-approved NOT_COMPARABLE
def hist_rows(n=500, not_ok=(), missing=()):
    return {k: {"status": "not_ok" if k in not_ok else "ok"} for k in range(n) if k not in missing}


def test_d2_in_domain_skips_are_not_comparable_with_runs_k_and_reasons(camp, checks):
    camp.over = {"as_run_variant": True}
    camp.all_runs()
    env = synth_env()
    env.load_historical = lambda root, r: hist_rows(not_ok={3, 4}, missing={7})
    rep = verify(camp, env=env)
    v4 = [x for x in by_check(rep, "V4") if x["run_id"] == 0][0]
    assert v4["status"] == "PASS" and v4["counts"]["n_compared"] == 497     # the 497 real comparisons
    nc = [x for x in by_check(rep, "V4:not_comparable") if x["run_id"] == 0][0]
    assert nc["status"] == "NOT_COMPARABLE" and nc["counts"]["n_not_comparable"] == 3
    assert nc["evidence"]["skipped_k"]["historical_status_not_ok"] == [3, 4]
    assert nc["evidence"]["skipped_k"]["historical_row_missing"] == [7]
    assert any("historical_status_not_ok: k=[3, 4]" == r for r in nc["reasons"])
    assert {x["run_id"] for x in by_check(rep, "V4:not_comparable")} == set(range(5))
    assert rep["verdict"] == "OWNER_REVIEW_REQUIRED" and not rep["reporting_permitted"]
    assert rep["status_counts"]["FAIL"] == 0                                 # never a scientific FAIL
    assert "NOT_COMPARABLE" in rep["policies"]["D2_V4_skipped_records"]


def test_d2_every_in_domain_record_skipped_keeps_v4_not_comparable(camp, checks, monkeypatch):
    camp.over = {"as_run_variant": True}
    camp.all_runs()
    monkeypatch.setattr(V, "v4_as_run_apparent", ORIGINAL_V4)               # the real V4
    env = synth_env()
    env.load_historical = lambda root, r: hist_rows(not_ok=set(range(500)))
    rep = verify(camp, env=env)
    v4 = [x for x in by_check(rep, "V4") if x["run_id"] == 2][0]
    assert v4["status"] == "NOT_COMPARABLE" and v4["raw_status"] == "NOT_APPLICABLE"
    assert rep["verdict"] == "OWNER_REVIEW_REQUIRED" and rep["status_counts"]["FAIL"] == 0


def test_d2_comparison_mismatch_remains_fail(camp, checks, monkeypatch):
    camp.over = {"as_run_variant": True}
    camp.all_runs()
    monkeypatch.setattr(V, "v4_as_run_apparent", ORIGINAL_V4)
    env = synth_env()
    env.load_historical = lambda root, r: {k: {"status": "ok", "dC_app": 9.0,
                                               "C_app": {"M1": 0.6, "M2": 0.6},
                                               "lambda_apparent": {"M1": 0.0, "M2": 0.0}}
                                           for k in range(500)}
    rep = verify(camp, env=env)
    assert {x["status"] for x in by_check(rep, "V4")} <= {"FAIL"} and rep["verdict"] == "FAIL"


def test_d2_beyond_domain_disclosed_separately_and_pass_without_skips(camp, checks):
    camp.over = {"w5": W5A, "ladder": LADDER_750, "as_run_variant": True}
    S["drift"] = 0.05
    camp.all_runs()
    camp.sizing(500)
    for r in range(5):
        camp.run(run_id=r, target=750, topup_authorization="OWNER-750")
    env = synth_env()
    env.load_historical = lambda root, r: hist_rows()
    rep = verify(camp, env=env)
    v4 = [x for x in by_check(rep, "V4") if x["run_id"] == 1][0]
    assert v4["status"] == "PASS" and v4["counts"]["n_beyond_historical_domain"] == 250
    assert v4["evidence"]["skip_accounting"]["n_skipped"]["beyond_historical_range"] == 250
    assert by_check(rep, "V4:not_comparable") == [] and rep["verdict"] == "PASS"


def test_d2_accounting_unit_conditions():
    cfg_var = CorrectedConfig(approved_test_config(ROOT, "scratch_unused", as_run_variant=True), ROOT)
    recs = [{"k": 0, "dC_app_as_run": 0.1}, {"k": 1}, {"k": 2, "dC_app_as_run": 0.1},
            {"k": 3, "dC_app_as_run": 0.1}, {"k": 600, "dC_app_as_run": 0.1},
            {"k": 4, "apparent_as_run": {"M1": {"ok": True}, "M2": {"ok": False}}}]
    hist = {0: {"status": "ok"}, 1: {"status": "ok"}, 2: {"status": "unreliable"}, 4: {"status": "ok"}}
    a = postrun.v4_skip_accounting(recs, hist, cfg_var)
    assert a["n_skipped"] == {"beyond_historical_range": 1, "historical_row_missing": 1,
                              "historical_status_not_ok": 1, "corrected_fit_failed": 1,
                              "corrected_field_missing": 1}
    assert a["skipped_k"]["corrected_fit_failed"] == [4] and a["skipped_k"]["historical_row_missing"] == [3]
    assert a["n_compared_expected"] == 1 and a["n_policy_relevant_skips"] == 4
    cfg_na = CorrectedConfig(approved_test_config(ROOT, "scratch_unused"), ROOT)   # variant off
    assert postrun.v4_skip_accounting(recs, {}, cfg_na) is None


def test_d2_accounting_exception_is_structured_not_pass(camp, checks, monkeypatch):
    camp.over = {"as_run_variant": True}
    camp.all_runs()

    def broken(*a, **k):
        raise KeyError("k")
    monkeypatch.setattr(postrun, "v4_skip_accounting", broken)
    env = synth_env()
    env.load_historical = lambda root, r: hist_rows()                    # genuine raw V4 PASS
    rep = verify(camp, env=env)
    v4 = by_check(rep, "V4")
    assert {x["raw_status"] for x in v4} == {"PASS"}
    assert {x["status"] for x in v4} == {"NOT_COMPARABLE"}                # raw PASS, but never PASS
    sa = by_check(rep, "V4:skip-accounting")
    assert sa[0]["evidence"]["stage"] == "V4 skip accounting" and sa[0]["evidence"]["exception_type"] == "KeyError"
    assert rep["verdict"] == "OWNER_REVIEW_REQUIRED" and not rep["reporting_permitted"]


def test_d2_accounting_exception_never_hides_a_real_v4_fail(camp, checks, monkeypatch):
    camp.over = {"as_run_variant": True}
    camp.all_runs()
    checks.status["V4"] = "FAIL"

    def broken(*a, **k):
        raise AttributeError("non-dict apparent_as_run")
    monkeypatch.setattr(postrun, "v4_skip_accounting", broken)
    rep = verify(camp)
    assert {x["status"] for x in by_check(rep, "V4")} == {"FAIL"}
    assert {x["status"] for x in by_check(rep, "V4:skip-accounting")} == {"NOT_COMPARABLE"}
    assert rep["verdict"] == "FAIL"


def test_d2_count_disagreement_is_fail(camp, checks):
    camp.over = {"as_run_variant": True}
    camp.all_runs()
    env = synth_env()
    env.load_historical = lambda root, r: hist_rows(missing={1})
    checks.extra["V4"] = {"n_compared": 500}               # claims more than were comparable
    rep = verify(camp, env=env)
    assert {x["status"] for x in by_check(rep, "V4")} == {"FAIL"}


# ---- D5: run selections with and without run 0
def test_d5_run_selection_without_run_0_is_owner_review_not_a_crash(camp, checks):
    camp.over = {"run_ids": [1, 2]}
    camp.run(run_id=1)
    camp.run(run_id=2)
    rep = verify(camp)
    v5 = by_check(rep, "V5")
    assert len(v5) == 1 and v5[0]["status"] == "NOT_COMPARABLE" and "run 0" in v5[0]["reasons"][0]
    assert rep["verdict"] == "OWNER_REVIEW_REQUIRED" and not rep["hard_failure_h5"]
    assert {x["run_id"] for x in by_check(rep, "V1-self")} == {1, 2}


def test_d5_run_selection_with_run_0(camp, checks):
    camp.over = {"run_ids": [0, 3]}
    camp.run(run_id=0)
    camp.run(run_id=3)
    rep = verify(camp)
    assert by_check(rep, "V5")[0]["status"] == "PASS" and rep["verdict"] == "PASS"


def test_d5_selected_run_missing_or_unselected_run_present_refuses(camp, checks):
    camp.over = {"run_ids": [0, 1]}
    camp.run(run_id=0)
    with pytest.raises(postrun.VerificationRefused, match=r"no campaign records: \[1\]"):
        verify(camp)
    camp.run(run_id=1)
    with pytest.raises(postrun.VerificationRefused, match=r"outside the selected run_ids: \[1\]"):
        verify(camp, run_ids=[0])
    assert not (camp.out / "verification").exists()


# ---- D6 / D4: V5 exceptions and NaN under the approved comparison
def test_d6_v5_recomputation_exception_is_fail_not_h5(camp, checks):
    camp.all_runs()

    def recompute(*a):
        if a[3] == 0 and a[4] == 50:
            raise RuntimeError("synthetic recomputation crash")
        return fake_record(*a)
    rep = verify(camp, env=synth_env(recompute))              # no HardFailureH5 raised
    v5 = by_check(rep, "V5")[0]
    assert v5["status"] == "FAIL" and not rep["hard_failure_h5"] and rep["verdict"] == "FAIL"
    assert v5["evidence"]["exception_type"] == "RuntimeError" and v5["evidence"]["run_id"] == 0
    assert v5["evidence"]["stage"] == "V5 recomputation" and "synthetic recomputation crash" in v5["evidence"]["message"]


def test_d4_top_level_nan_keeps_approved_fail_and_is_classified(camp, checks):
    S["nan"] = "top"                                           # every record: top-level and nested NaN, -inf
    camp.all_runs()
    with pytest.raises(postrun.HardFailureH5):                 # approved `!=`: top-level NaN != NaN -> FAIL
        verify(camp)
    rep = json.loads(postrun.report_paths(camp.out, 500)[0].read_bytes())
    v5 = by_check(rep, "V5")[0]
    assert v5["status"] == "FAIL" and rep["hard_failure_h5"]
    nf = v5["evidence"]["nonfinite"]
    assert nf["fields_by_kind"] == {"nested_nan_in_dict": ["diag_dict"], "nested_nan_in_list": ["diag_list"],
                                    "other_nonfinite": ["diag_inf"], "top_level_nan": ["diag"]}
    assert {tuple(x["fields"]) for x in nf["nan_only_differences"]} == {("diag",)}   # only the top-level one
    row = [x for x in v5["evidence"]["rows"] if x["k"] == 0][0]
    assert row["differing_fields"] == ["diag"]                 # nested NaN and -inf compare equal
    assert "top-level NaN" in " ".join(v5["reasons"])


def test_d4_nested_nan_and_infinity_compare_equal_and_are_disclosed(camp, checks):
    S["nan"] = "nested"
    camp.all_runs()
    rep = verify(camp)                                         # no H5
    v5 = by_check(rep, "V5")[0]
    assert v5["status"] == "PASS" and not rep["hard_failure_h5"]
    assert set(v5["evidence"]["nonfinite"]["fields_by_kind"]) == {"nested_nan_in_dict", "nested_nan_in_list",
                                                                  "other_nonfinite"}
    assert v5["evidence"]["nonfinite"]["nan_only_differences"] == []



# ================================================================ final remediation cycle
# ---- D3: independent expected state counts
def test_d3_expected_state_count_matches_the_real_enumerator_on_small_domains():
    from src.phase9_corrected.rng import enumerate_stream_states
    for conv, T in (("independent", 7), ("as_run", None)):
        for run_ids, k_max, b_max in (((0,), 3, 2), ((1, 4), 2, 3), ((0, 1, 2, 3, 4), 2, 1)):
            dom = enumerate_stream_states(PlanStreams(2026110300, run_ids[0], conv, T), run_ids, k_max, b_max)
            assert dom["n_states"] == postrun.expected_v2_domain_states(run_ids, k_max, b_max, conv)
            assert dom["n_distinct"] == dom["n_states"] and dom["collisions"] == []


@pytest.mark.parametrize("run_ids,k_max,b_max,conv,literal", [
    ((0, 1, 2, 3, 4), 3000, 100, "independent", 4_530_005),
    ((0, 1, 2, 3, 4), 3000, 100, "as_run", 4_515_005),          # 1 + 300 per k
    ((0, 1, 2, 3, 4), 750, 100, "independent", 1_132_505),
    ((1, 2), 1000, 120, "independent", 724_005),                # 2 x 1000 x 362 + 5
    ((0,), 500, 100, "independent", 151_005)])
def test_d3_expected_state_count_hand_literals(run_ids, k_max, b_max, conv, literal):
    assert postrun.expected_v2_domain_states(run_ids, k_max, b_max, conv) == literal


@pytest.mark.parametrize("damage", ["k_max_low", "b_max_high", "run_missing", "convention",
                                    "states_missing", "duplicate_states", "count_high"])
def test_d3_gate_domain_and_count_damage(camp, checks, damage):
    camp.all_runs()
    p, _ = gate_file(camp)
    g = json.loads(p.read_bytes())
    d = g["V2_domain"]
    if damage == "k_max_low":
        d["domain"]["k_max"] = 2000
    elif damage == "b_max_high":
        d["domain"]["b_max"] = 101
    elif damage == "run_missing":
        d["domain"]["run_ids"] = [0, 1, 2, 3]
    elif damage == "convention":
        d["domain"]["convention"] = "as_run"
    elif damage == "states_missing":
        d["n_states"] = d["n_distinct"] = 4_530_004
    elif damage == "duplicate_states":
        d["n_distinct"] = 4_530_004                              # one state seen twice
    elif damage == "count_high":
        d["n_states"] = d["n_distinct"] = 4_530_006
    probs = orch.gate_evidence_problems(canonical_json(g), camp.cfg(), camp.msha())
    assert any("V2_domain" in x for x in probs), probs


# ---- D5: run artifacts
@pytest.mark.parametrize("artifact,match", [
    ("runs/run1.lock", "incomplete run artifacts"),               # stale or active lock, selected run
    ("runs/run2.jsonl.tmp", "incomplete run artifacts"),
    ("halts/run3.inv0002.end.json.tmp", "incomplete run artifacts"),
    ("runs/run4.jsonl.bak", "unrecognised run artifacts"),
    ("halts/run0.invXX.start.json", "unrecognised run artifacts")])
def test_d5_incomplete_or_unrecognised_artifacts_refuse(camp, checks, artifact, match):
    camp.all_runs()
    (camp.out / artifact).write_text("x")
    with pytest.raises(postrun.VerificationRefused, match=match):
        verify(camp)
    assert not (camp.out / "verification").exists()


@pytest.mark.parametrize("artifact", ["runs/run3.lock", "runs/run4.jsonl.tmp", "halts/run2.inv0001.start.json.tmp"])
def test_d5_unselected_run_seen_only_through_incomplete_artifacts_refuses(camp, checks, artifact):
    camp.over = {"run_ids": [0, 1]}
    camp.run(run_id=0)
    camp.run(run_id=1)
    (camp.out / artifact).write_text("x")
    with pytest.raises(postrun.VerificationRefused, match="incomplete run artifacts|outside the selected"):
        verify(camp)


def test_d5_unrelated_files_are_not_runs_and_are_disclosed(camp, checks):
    camp.all_runs()
    (camp.out / "runs/.DS_Store").write_text("x")
    (camp.out / "halts/README").write_text("x")
    rep = verify(camp)
    art = by_check(rep, "campaign-artifacts")[0]
    assert art["status"] == "PASS" and art["evidence"]["unrelated_files"] == ["runs/.DS_Store", "halts/README"]
    assert rep["verdict"] == "PASS"


def test_d5_lock_left_by_a_live_writer_is_detected(camp, checks):
    from src.phase9_corrected.checkpoint import CheckpointStore
    camp.all_runs()
    cfg = camp.cfg()
    with CheckpointStore(camp.out, 0, cfg.sha256(), camp.msha(), ROOT):    # a writer holding run 0
        with pytest.raises(postrun.VerificationRefused, match="run0.lock"):
            verify(camp)
    rep = verify(camp)                                                     # released: verifies
    assert rep["verdict"] == "PASS"


# ---- D1: real module loading, packages, extension packages
@pytest.fixture
def shadow_repo(tmp_path, monkeypatch):
    root = tmp_path / "repo"
    pkg = [f"src/phase9_corrected/{q.name}" for q in (ROOT / "src/phase9_corrected").glob("*.py")]
    for rel in pkg + list(orch.REQUIRED_FILES):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, root / rel)
    code = {rel: hashlib.sha256((root / rel).read_bytes()).hexdigest() for rel in pkg + list(orch.REQUIRED_FILES)}
    monkeypatch.syspath_prepend(str(root))
    return root, code


def _import_fresh(monkeypatch, name):
    import importlib
    monkeypatch.delitem(sys.modules, name, raising=False)
    mod = importlib.import_module(name)
    monkeypatch.setitem(sys.modules, name, mod)
    return mod


def test_d1_real_loaded_module_shadow_is_refused(shadow_repo, monkeypatch):
    root, code = shadow_repo
    postrun.identity_guard(root, code, "baseline")                        # clean: passes (real checks)
    (root / "zz_p9c_shadow_mod.py").write_text("LOADED = True\n")
    with pytest.raises(postrun.VerificationRefused, match="shadow imports"):
        postrun.identity_guard(root, code, "output")                      # caught before import
    mod = _import_fresh(monkeypatch, "zz_p9c_shadow_mod")                  # actually loaded
    assert mod.LOADED and Path(mod.__file__).parent == root
    with pytest.raises(orch.CodeIdentityRefusal, match="not a required inventory file"):
        orch.verify_loaded_modules(root, code, postrun.required_set(root))  # real check, loaded module


def test_d1_real_loaded_regular_package_shadow_is_refused(shadow_repo, monkeypatch):
    root, code = shadow_repo
    (root / "zz_p9c_shadow_pkg").mkdir()
    (root / "zz_p9c_shadow_pkg/__init__.py").write_text("X = 1\n")
    _import_fresh(monkeypatch, "zz_p9c_shadow_pkg")
    assert postrun.root_shadow_entries(root) == ["zz_p9c_shadow_pkg"]
    with pytest.raises(orch.CodeIdentityRefusal, match="zz_p9c_shadow_pkg"):
        orch.verify_loaded_modules(root, code, postrun.required_set(root))


def test_d1_versioned_extension_package_detected_by_the_root_scan(tmp_path):
    root = tmp_path / "r"
    for name, init in (("zz_ext", "__init__.cpython-313-darwin.so"), ("zz_abi", "__init__.abi3.so"),
                       ("zz_ns", None)):
        (root / name).mkdir(parents=True)
        if init:
            (root / name / init).write_bytes(b"")
    assert postrun.root_shadow_entries(root) == ["zz_abi", "zz_ext"]      # namespace dir: limitation


def test_d1_real_shadow_loaded_during_stage_a_prevents_all_output(camp, checks, monkeypatch, shadow_repo):
    root, code = shadow_repo
    w = windowed(camp)
    real_lm = ORIGINAL_VERIFY_LOADED

    def lm(_root, _code, _required):                                      # the REAL check, temp repo
        real_lm(root, code, postrun.required_set(root))
    monkeypatch.setattr(orch, "verify_loaded_modules", lm)
    (root / "zz_p9c_lazy.py").write_text("V = 1\n")
    monkeypatch.setattr(postrun, "root_shadow_entries", lambda r: [])      # isolate the post-import check

    def recompute(*a):
        _import_fresh(monkeypatch, "zz_p9c_lazy")                         # imported lazily, mid-run
        return fake_record(*a)
    with pytest.raises(postrun.VerificationRefused, match="zz_p9c_lazy"):
        verify(camp, env=synth_env(recompute))
    assert not orch.window_result_path(camp.out, 2, w["invocation"]).exists()
    assert not (camp.out / "verification").exists()



# ================================================================ gate remediation: adversarial tests
def test_gate_is_rerun_for_every_invocation_including_after_a_pass(camp, checks):
    GATE_CALLS.clear()
    with pytest.raises(orch.RuntimeBudgetHalt):
        orch.run_corrected(camp.cfg(), 0, 500, clock=("TEST_CLOCK", StepClock(600)))   # inv1 PASS
    camp.run(run_id=0, resume_authorization="RT-1")                                  # inv2
    assert GATE_CALLS == ["run0.inv0001.attempt01.json", "run0.inv0002.attempt01.json"]
    p1, s1 = gate_file(camp, 0, 1)
    p2, s2 = gate_file(camp, 0, 2)
    assert p1 != p2 and s1["gate_evidence"]["sha256"] == hashlib.sha256(p1.read_bytes()).hexdigest()
    assert not (camp.out / "gate/pre_execution_gate.json").exists()                   # no shared cache


def test_failed_gate_dispatches_nothing_and_keeps_its_evidence(camp, checks, monkeypatch):
    from tests.p9c_support import make_synthetic_gate
    monkeypatch.setattr(orch, "run_pre_execution_gate", make_synthetic_gate(fail_check="V2_domain"))
    with pytest.raises(SystemExit, match="gate FAILED.*nothing dispatched"):
        camp.run(run_id=0)
    ev = camp.out / "gate/run0.inv0001.attempt01.json"
    assert json.loads(ev.read_text())["V2_domain"]["status"] == "FAIL"               # preserved
    assert not (camp.out / "halts").exists() or not list((camp.out / "halts").glob("run0.*start*"))
    ck = camp.out / "runs/run0.jsonl"
    assert not ck.exists() or orch._records_by_k(ck.read_bytes()) == {}               # zero replicates


def test_failed_gate_edited_to_pass_cannot_authorise_dispatch(camp, checks, monkeypatch):
    """The demonstrated attack: FAIL evidence rewritten to a fully consistent PASS."""
    from tests.p9c_support import make_synthetic_gate, synthetic_gate_body
    monkeypatch.setattr(orch, "run_pre_execution_gate", make_synthetic_gate(fail_check="V2_numpy"))
    with pytest.raises(SystemExit):
        camp.run(run_id=0)
    ev1 = camp.out / "gate/run0.inv0001.attempt01.json"
    forged = canonical_json(synthetic_gate_body(camp.cfg(), camp.msha()))             # perfect forgery
    ev1.write_bytes(forged)
    with pytest.raises(SystemExit, match="gate FAILED"):                              # gate re-executed
        camp.run(run_id=0)
    assert ev1.read_bytes() == forged                                                 # never read as a cache
    assert json.loads((camp.out / "gate/run0.inv0001.attempt02.json").read_text())["all_pass"] is False
    assert not list((camp.out / "halts").glob("run0.*")) if (camp.out / "halts").exists() else True


@pytest.mark.parametrize("lie", ["return_value_only", "file_all_pass_flag"])
def test_forged_all_pass_true_cannot_authorise_dispatch(camp, checks, monkeypatch, lie):
    from tests.p9c_support import synthetic_gate_body

    def lying_gate(config, msha, progress=None, out_path=None):
        body = synthetic_gate_body(config, msha, fail_check="V1_outer_historical")
        if lie == "file_all_pass_flag":
            body["all_pass"] = True                                 # flag says PASS, a check says FAIL
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        Path(out_path).write_bytes(canonical_json(body))
        return {**body, "all_pass": True}                           # the returned value lies too
    monkeypatch.setattr(orch, "run_pre_execution_gate", lying_gate)
    with pytest.raises(SystemExit, match="V1_outer_historical is not PASS"):
        camp.run(run_id=0)
    assert not (camp.out / "halts").exists() or not list((camp.out / "halts").glob("run0.*start*"))


def test_malformed_gate_output_fails_closed_with_a_structured_reason(camp, checks, monkeypatch):
    def broken_gate(config, msha, progress=None, out_path=None):
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        Path(out_path).write_bytes(b"{truncated")
        return {"all_pass": True}
    monkeypatch.setattr(orch, "run_pre_execution_gate", broken_gate)
    with pytest.raises(SystemExit, match="not valid JSON"):
        camp.run(run_id=0)
    monkeypatch.setattr(orch, "run_pre_execution_gate", lambda c, m, progress=None, out_path=None: {"all_pass": True})
    with pytest.raises(SystemExit, match="missing or unreadable"):                    # no file written
        camp.run(run_id=1)


@pytest.mark.parametrize("tamper", ["modified", "missing", "malformed_rebound", "invalid_rebound"])
def test_tampered_evidence_after_dispatch_is_never_eligible(camp, checks, tamper):
    from tests.p9c_support import synthetic_gate_body
    camp.all_runs()
    p, st = gate_file(camp, 2)
    sp = camp.out / "halts/run2.inv0001.start.json"
    if tamper == "modified":
        p.write_bytes(p.read_bytes().replace(b'"PASS"', b'"PASS" ', 1))
    elif tamper == "missing":
        p.unlink()
    else:                                    # coordinated rewrite of evidence AND start record hash
        new = b"{bad json" if tamper == "malformed_rebound" else \
            canonical_json(damage_body(synthetic_gate_body(camp.cfg(), camp.msha()), "domain_fail"))
        p.write_bytes(new)
        st["gate_evidence"]["sha256"] = hashlib.sha256(new).hexdigest()
        sp.write_text(json.dumps(st))
    el = orch.evaluate_eligibility(camp.cfg(), 2, 500, camp.msha())
    assert el["class"] != "ELIGIBLE"
    codes = {x["code"] for x in el["reasons"]}
    assert codes & {"gate_evidence_invalid", "gate_evidence_unverifiable"}
    if tamper == "missing":
        assert "gate_evidence_unverifiable" in codes
    rep = verify(camp)
    g2 = [x for x in by_check(rep, "gate-evidence") if x["run_id"] == 2][0]
    assert g2["status"] == "FAIL" and rep["verdict"] == "FAIL" and not rep["reporting_permitted"]
    assert [x for x in by_check(rep, "StageB-eligibility") if x["run_id"] == 2][0]["status"] == "FAIL"


def test_tampered_evidence_blocks_further_dispatch(camp, checks):
    with pytest.raises(orch.RuntimeBudgetHalt):
        orch.run_corrected(camp.cfg(), 0, 500, clock=("TEST_CLOCK", StepClock(600)))
    p, _ = gate_file(camp, 0, 1)
    p.write_bytes(p.read_bytes() + b" ")
    with pytest.raises(orch.DispatchRefused, match="gate_evidence_invalid"):
        camp.run(run_id=0, resume_authorization="RT-1")


def test_gate_failure_on_resume_consumes_no_authorisation_then_genuine_pass_continues(camp, checks, monkeypatch):
    from tests.p9c_support import make_synthetic_gate
    with pytest.raises(orch.RuntimeBudgetHalt):
        orch.run_corrected(camp.cfg(), 0, 500, clock=("TEST_CLOCK", StepClock(600)))   # inv1 halts
    good = orch.run_pre_execution_gate
    monkeypatch.setattr(orch, "run_pre_execution_gate", make_synthetic_gate(fail_check="H4"))
    with pytest.raises(SystemExit, match="gate FAILED"):
        camp.run(run_id=0, resume_authorization="RT-1")
    assert not (camp.out / "halts/run0.inv0002.start.json").exists()
    assert not (camp.out / "halts/run0.inv0001.auth.json").exists()                   # not consumed
    assert orch.invocation_state(camp.out, 0)["state"] == "HALTED"                     # still pending
    monkeypatch.setattr(orch, "run_pre_execution_gate", good)
    assert camp.run(run_id=0, resume_authorization="RT-1")["n_records"] == 500          # same reference
    st = json.loads((camp.out / "halts/run0.inv0002.start.json").read_text())
    assert st["gate_evidence"]["path"] == "gate/run0.inv0002.attempt02.json"           # attempt 1 kept
    assert json.loads((camp.out / "gate/run0.inv0002.attempt01.json").read_text())["H4"] == "FAIL"
    assert orch.evaluate_eligibility(camp.cfg(), 0, 500, camp.msha())["class"] == "ELIGIBLE"


# ================================================================ gate remediation review fixes (D1, D2, D4)
def test_d1_data_load_failure_after_gate_consumes_no_authorisation(camp, checks, monkeypatch):
    with pytest.raises(orch.RuntimeBudgetHalt):
        orch.run_corrected(camp.cfg(), 0, 500, clock=("TEST_CLOCK", StepClock(600)))   # inv1 halts

    def failing_load(root):
        raise OSError("synthetic data-load failure")
    monkeypatch.setattr(orch, "load_context", failing_load)
    with pytest.raises(OSError, match="synthetic data-load failure"):
        camp.run(run_id=0, resume_authorization="RT-1")                 # gate PASSes, then load fails
    assert (camp.out / "gate/run0.inv0002.attempt01.json").exists()     # the gate did run
    assert not (camp.out / "halts/run0.inv0001.auth.json").exists()     # authorisation NOT consumed
    assert not (camp.out / "halts/run0.inv0002.start.json").exists()
    assert orch.invocation_state(camp.out, 0)["state"] == "HALTED"
    monkeypatch.setattr(orch, "load_context", lambda root: None)
    assert camp.run(run_id=0, resume_authorization="RT-1")["n_records"] == 500   # same reference works


@pytest.mark.parametrize("malformed,match", [
    ("detail_list", "V2_numpy detail does not pass"), ("structural_detail_str", "V2_structural detail"),
    ("deep_nesting", "malformed"), ("domain_not_dict", "malformed|V2_domain")])
def test_d4_malformed_evidence_types_never_raise(camp, checks, malformed, match):
    import re as _re
    from tests.p9c_support import synthetic_gate_body
    camp.all_runs()
    g = synthetic_gate_body(camp.cfg(), camp.msha())
    if malformed == "detail_list":
        g["V2_numpy"]["detail"] = [1]
        raw = canonical_json(g)
    elif malformed == "structural_detail_str":
        g["V2_structural"]["detail"] = "holds"
        raw = canonical_json(g)
    elif malformed == "deep_nesting":
        raw = b"[" * 100000 + b"]" * 100000
    else:
        g["V2_domain"]["domain"] = ["not", "a", "dict"]
        g["V2_domain"]["n_states"] = {"x": 1}
        raw = canonical_json(g)
    probs = orch.gate_evidence_problems(raw, camp.cfg(), camp.msha())        # must not raise
    assert probs and any(_re.search(match, x) for x in probs), probs


@pytest.mark.parametrize("ref", [{"path": "gate/run0.inv0001.attempt01.json\x00x", "sha256": "a" * 64,
                                  "all_pass": True, "attempt": 1},
                                 {"path": ["gate"], "sha256": "a" * 64}, "not-a-dict"])
def test_d4_malformed_start_record_references_are_structured(camp, checks, ref):
    camp.all_runs()
    sp = camp.out / "halts/run1.inv0001.start.json"
    st = json.loads(sp.read_text())
    st["gate_evidence"] = ref
    sp.write_text(json.dumps(st))
    cat, probs = orch.start_gate_problems(camp.out, st, camp.cfg(), camp.msha())   # must not raise
    assert cat in ("blocking", "disqualifier") and probs
    el = orch.evaluate_eligibility(camp.cfg(), 1, 500, camp.msha())
    assert el["class"] != "ELIGIBLE"
    rep = verify(camp)
    assert [x for x in by_check(rep, "gate-evidence") if x["run_id"] == 1][0]["status"] == "FAIL"


@pytest.mark.parametrize("target", ["other_run", "other_invocation", "attempt_mismatch"])
def test_d2_evidence_bound_to_its_own_invocation(camp, checks, target):
    with pytest.raises(orch.RuntimeBudgetHalt):
        orch.run_corrected(camp.cfg(), 0, 500, clock=("TEST_CLOCK", StepClock(600)))
    camp.run(run_id=0, resume_authorization="RT-1")                     # run 0: inv1, inv2
    camp.run(run_id=1)
    sp = camp.out / "halts/run0.inv0002.start.json"
    st = json.loads(sp.read_text())
    if target == "other_run":
        rel = "gate/run1.inv0001.attempt01.json"                       # a genuine PASS, wrong run
    elif target == "other_invocation":
        rel = "gate/run0.inv0001.attempt01.json"                       # a genuine PASS, wrong invocation
    else:
        rel = st["gate_evidence"]["path"]
        st["gate_evidence"]["attempt"] = 7
    st["gate_evidence"]["path"] = rel
    st["gate_evidence"]["sha256"] = hashlib.sha256((camp.out / rel).read_bytes()).hexdigest()
    sp.write_text(json.dumps(st))
    cat, probs = orch.start_gate_problems(camp.out, st, camp.cfg(), camp.msha())
    assert cat == "disqualifier" and "does not belong" in probs[0]
    assert orch.evaluate_eligibility(camp.cfg(), 0, 500, camp.msha())["class"] == "INELIGIBLE"
