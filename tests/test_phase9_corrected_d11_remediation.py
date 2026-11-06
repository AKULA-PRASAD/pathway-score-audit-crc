"""Decision 11 remediation (owner decisions 2026-11-06): D11-D1 (window FAIL blocks dispatch;
DEC-C1 pending-window rules), D11-D2 (Path C assesses the complete invocation history of
every run), D11-T2 (code-identity gaps), DEC-B2 (freeze / driver hash order), DEC-E1
(runtime launch flags and the external launcher), DEC-A1 (non-finite serialisation).

Synthetic only. The scientific fit, data loading, gate, H4 and (unless a test says
otherwise) code identity are stubbed; temporary repositories live under pytest's tmp_path
or in scratch directories that are removed. The launcher tests run a HARMLESS STUB driver
in a temporary repository copy -- never the real scientific driver."""
import hashlib
import importlib.util
import json
import marshal
import math
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import numpy as np
import pytest

from src.phase9_corrected import orchestration as orch
from src.phase9_corrected.checkpoint import CheckpointStore
from src.phase9_corrected.config import CorrectedConfig, canonical_json
from src.phase9_corrected.integrity import HardFailureH4
from src.phase9_corrected.verification import v5_spot_check
from tests.p9c_support import approved_test_config, synthetic_gate

ROOT = Path(__file__).resolve().parents[1]
REAL_VERIFY = orch.verify_code_identity
LAUNCHER = ROOT / "scripts/42_phase9_corrected_launch.sh"
SYNTH_FREEZE = {"freeze_record_sha256": "f" * 64, "freeze": {"code": {"x": "0" * 64}}}
W5A = {"action": "W5-A", "corrective_rule": "next_rung_above_current",
       "definition": "adjudication_D5_as_proposed", "exclusion_rule": "usable_within_k_positions_v1"}
LADDER = {"rungs": [750, 1000, 1500, 2000, 3000], "beyond_ladder_rule": "top_up_to_max_rung_then_reassess"}

# ---------------------------------------------------------------- synthetic run environment
B = {"bad": {}, "raise_at": None, "calls": 0, "drift": 0.0}


def fake_record(ctx, streams, pipeline, run_id, k, B_in, lam_grid, as_run_variant=False):
    B["calls"] += 1
    if B["raise_at"] is not None and k == B["raise_at"]:
        raise RuntimeError(f"simulated crash at k={k}")
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
        self.h4_fail_stage, self.pres_fail_next = None, False
        self.cfg_over = {}

    def cfg(self, n_workers=1):
        return CorrectedConfig(approved_test_config(ROOT, self.out_rel, n_workers=n_workers,
                                                    runtime_budget_hours_per_run=4.0,
                                                    **self.cfg_over), ROOT)

    def run(self, run_id=0, target=500, clock=None, **kw):
        return orch.run_corrected(self.cfg(), run_id, target,
                                  clock=("TEST_CLOCK", clock or StepClock()), **kw)

    def msha(self):
        return hashlib.sha256((self.out / "manifest/input_manifest.json").read_bytes()).hexdigest()

    def elig(self, r=0, target=500):
        return orch.evaluate_eligibility(self.cfg(), r, target, self.msha())

    def halts(self, r, m, kind):
        return self.out / "halts" / f"run{r}.inv{m:04d}.{kind}.json"

    def rewrite(self, r, m, kind, fn):
        p = self.halts(r, m, kind)
        o = json.loads(p.read_text())
        fn(o)
        p.write_text(json.dumps(o))

    def ckpt(self, r=0):
        return self.out / "runs" / f"run{r}.jsonl"


@pytest.fixture
def env(monkeypatch):
    e = Env(f"scratch_p9c_rem_{uuid.uuid4().hex[:10]}")
    B.update(bad={}, raise_at=None, calls=0, drift=0.0)

    def h4(root, stage, raise_on_fail=True):
        if stage == e.h4_fail_stage:
            raise HardFailureH4(f"HARD FAILURE H4 ({stage}) -- simulated")
        return {"status": "PASS", "all_ok": True}

    def preservation(root):
        if e.pres_fail_next:
            e.pres_fail_next = False
            return {"all_ok": False, "status": "FAIL", "files": {
                "reports/phase9/dc_interval/run0_outer.jsonl":
                    {"expected_prefix": "885bbd1bd891b570", "sha256": "ab" * 32, "ok": False}}}
        return {"all_ok": True, "status": "PASS"}

    monkeypatch.setattr(orch, "h4_check", h4)
    monkeypatch.setattr(orch, "historical_preservation_check", preservation)
    monkeypatch.setattr(orch, "build_manifest", lambda config, root: {
        "schema": "test-manifest", "config": config.sha256(), "code": {"x": "0" * 64}})
    monkeypatch.setattr(orch, "run_pre_execution_gate", synthetic_gate)
    monkeypatch.setattr(orch, "load_context", lambda root: None)
    monkeypatch.setattr(orch, "RealPipeline", lambda: None)
    monkeypatch.setattr(orch, "compute_outer_replicate", fake_record)
    monkeypatch.setattr(orch, "verify_code_identity",
                        lambda config, manifest=None, extra_files=(): SYNTH_FREEZE)
    yield e
    shutil.rmtree(e.out, ignore_errors=True)
    assert not e.out.exists()


def halt(env, **kw):
    with pytest.raises(orch.RuntimeBudgetHalt):
        env.run(clock=StepClock(600), **kw)


def crash(env, at, **kw):
    B["raise_at"] = at
    with pytest.raises(RuntimeError):
        env.run(**kw)
    B["raise_at"] = None


def window_result(env, r, m, status="PASS", bind=True):
    w = [x for x in orch.assess_run(env.cfg(), r, env.msha())["windows"] if x["invocation"] == m][0]
    res = {"schema": "p9c-window-1", "run_id": r, "invocation": m,
           "region_sha256": w["region_sha256"] if bind else "0" * 64, "k_list": w["k_list"],
           "code_identity_verified": True, "status": status,
           "per_k": {str(k): {"record_identical": status == "PASS", "dC_hex_equal": status == "PASS"}
                     for k in w["k_list"]}}
    p = orch.window_result_path(env.out, r, m)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(res))


def codes(e):
    return {x["code"] for x in e["reasons"]}


# ================================================================ D11-D1 / DEC-C1
def incomplete_with_window(env):
    """inv1 crashes at k=120 (UNTERMINATED); inv2 is the acknowledged continuation and
    halts on runtime: the run is incomplete (HALTED) with one pending window (inv1)."""
    crash(env, 120)
    halt(env, resume_authorization="ACK-1")
    st = orch.invocation_state(env.out, 0)
    assert st["state"] == "HALTED"


def test_d1_incomplete_run_with_pending_window_may_continue(env):
    incomplete_with_window(env)
    e = env.elig()
    assert e["class"] == "INCOMPLETE" and e["dispatch"] == "permitted_with_authorization"
    assert "window_not_recomputed" in codes(e)
    assert env.run(resume_authorization="R-1")["n_records"] == 500       # continuation permitted
    e2 = env.elig()
    assert e2["class"] == "UNVERIFIABLE" and e2["dispatch"] == "not_applicable"  # never ELIGIBLE


def test_d1_window_fail_refuses_resume_and_is_ineligible(env):
    incomplete_with_window(env)
    window_result(env, 0, 1, status="FAIL")
    e = env.elig()
    assert e["class"] == "INELIGIBLE" and e["dispatch"] == "blocked_permanently"
    assert "window_recomputation_failed" in codes(e)
    with pytest.raises(orch.DispatchRefused, match="FAILED"):            # actual dispatch agrees
        env.run(resume_authorization="R-1")
    assert not env.halts(0, 2, "auth").exists()                           # nothing consumed
    a = orch.assess_run(env.cfg(), 0, env.msha())
    assert "window_recomputation_failed" in a["reasons"].codes("disqualifier")


def test_d1_unbound_result_is_not_pass_but_permits_incomplete_continuation(env):
    incomplete_with_window(env)
    window_result(env, 0, 1, status="PASS", bind=False)
    e = env.elig()
    assert "window_result_not_bound" in codes(e) and e["class"] == "INCOMPLETE"
    assert e["dispatch"] == "permitted_with_authorization"
    env.run(resume_authorization="R-1")
    assert env.elig()["class"] == "UNVERIFIABLE"                          # not ELIGIBLE


def test_d1_one_pass_and_one_pending_window_is_still_not_eligible(env):
    crash(env, 100)
    crash(env, 300, resume_authorization="ACK-1")
    env.run(resume_authorization="ACK-2")
    window_result(env, 0, 1, status="PASS")
    e = env.elig()
    assert e["class"] == "UNVERIFIABLE"
    assert "window_recomputed_bit_identical" in {d["code"] for d in e["disclosures"]}
    assert "window_not_recomputed" in codes(e)
    window_result(env, 0, 2, status="PASS")
    assert env.elig()["class"] == "ELIGIBLE"                              # all required PASS


# ================================================================ D11-D2: Path C history model
def complete_five(env):
    for r in range(5):
        env.run(run_id=r)


@pytest.fixture
def topup(env):
    env.cfg_over = {"w5": W5A, "ladder": LADDER}
    B["drift"] = 0.05                                                     # W5-A breach -> 750
    complete_five(env)
    art = orch.build_sizing_decision(env.cfg(), 500, env.msha(), SYNTH_FREEZE)
    orch.durable_create(orch._sizing_path(env.out, 500), art)
    assert art["decision"]["action"] == "top_up" and art["decision"]["combined_B_next"] == 750
    return env


def topup_run(env, r, **kw):
    return env.run(run_id=r, target=750, topup_authorization="OWNER-750", **kw)


def test_d2_later_zero_growth_h4_in_another_run_refuses(topup):
    env = topup
    env.h4_fail_stage = "post"
    with pytest.raises(HardFailureH4):
        env.run(run_id=2)                                                  # no-op, zero growth
    env.h4_fail_stage = None
    assert json.loads(env.halts(2, 2, "end").read_text())["done_now"] == 0
    with pytest.raises(orch.DispatchRefused, match="run 2 history is not clean"):
        topup_run(env, 0)


@pytest.mark.parametrize("outcome", ["stop_nonestimable", "hard_failure_h6"])
def test_d2_later_zero_growth_disqualifying_outcome_refuses(topup, outcome):
    env = topup
    env.run(run_id=2)                                                      # benign no-op ...
    env.rewrite(2, 2, "end", lambda o: o.update(outcome=outcome))          # ... recorded as
    with pytest.raises(orch.DispatchRefused, match="run 2"):
        topup_run(env, 0)


def test_d2_benign_noop_resumes_do_not_block(topup):
    env = topup
    env.run(run_id=1)                                                      # another run
    env.run(run_id=0)                                                      # the run itself
    assert topup_run(env, 0)["n_records"] == 750


def test_d2_runtime_halted_topup_in_another_run_with_multiple_continuations(topup):
    env = topup
    with pytest.raises(orch.RuntimeBudgetHalt):
        topup_run(env, 1, clock=StepClock(600))
    with pytest.raises(orch.RuntimeBudgetHalt):
        topup_run(env, 1, clock=StepClock(600), resume_authorization="RT-1")
    topup_run(env, 1, resume_authorization="RT-2")                         # target reached
    assert topup_run(env, 0)["n_records"] == 750                           # Path C allowed
    for r in (2, 3, 4):
        topup_run(env, r)
    for r in range(5):
        assert env.elig(r, 750)["class"] == "ELIGIBLE"


def test_d2_pending_halt_in_another_run_is_not_disqualifying_but_not_eligible(topup):
    env = topup
    with pytest.raises(orch.RuntimeBudgetHalt):
        topup_run(env, 1, clock=StepClock(600))
    e1 = env.elig(1, 750)
    assert e1["class"] == "INCOMPLETE" and not e1["permanent"]
    assert topup_run(env, 0)["n_records"] == 750                           # snapshot at B_cur ELIGIBLE


def test_d2_acknowledged_interruption_then_same_target_continuation(topup):
    env = topup
    B["raise_at"] = 600
    with pytest.raises(RuntimeError):
        topup_run(env, 1)
    B["raise_at"] = None
    topup_run(env, 1, resume_authorization="ACK-1")
    assert topup_run(env, 0)["n_records"] == 750


def test_d2_wrong_target_path_or_reference_refuses(topup):
    env = topup
    with pytest.raises(orch.RuntimeBudgetHalt):
        topup_run(env, 1, clock=StepClock(600))
    topup_run(env, 1, resume_authorization="RT-1")                         # inv3 = resume
    for field, value in (("topup_authorization", "OTHER-REF"), ("path", "topup"),
                         ("target_B_out", 1000)):
        saved = json.loads(env.halts(1, 3, "start").read_text())
        env.rewrite(1, 3, "start", lambda o: o.update({field: value}))
        with pytest.raises(orch.DispatchRefused, match="run 1"):
            topup_run(env, 0)
        env.halts(1, 3, "start").write_text(json.dumps(saved))
    assert topup_run(env, 0)["n_records"] == 750                           # restored: allowed


def test_d2_preservation_failure_refuses_globally(topup):
    env = topup
    env.pres_fail_next = True
    with pytest.raises(SystemExit):
        env.run(run_id=4)
    with pytest.raises(orch.DispatchRefused, match="preservation"):
        topup_run(env, 0)


def test_d2_complete_run_with_pending_window_blocks_new_topup(env):
    env.cfg_over = {"w5": W5A, "ladder": LADDER}
    B["drift"] = 0.05
    for r in (0, 1, 2, 4):
        env.run(run_id=r)
    crash(env, 200, run_id=3)
    env.run(run_id=3, resume_authorization="ACK-3")                        # run 3: one window
    window_result(env, 3, 1, status="PASS")
    art = orch.build_sizing_decision(env.cfg(), 500, env.msha(), SYNTH_FREEZE)
    orch.durable_create(orch._sizing_path(env.out, 500), art)
    orch.window_result_path(env.out, 3, 1).unlink()                        # verification removed
    assert env.elig(3)["class"] == "UNVERIFIABLE"
    with pytest.raises(orch.DispatchRefused, match="run 3 eligibility at B_cur"):
        topup_run(env, 0)


# ================================================================ D11-T2 / DEC-B2 / DEC-E1 (Python side)
FROZEN = [f"src/phase9_corrected/{p.name}" for p in (ROOT / "src/phase9_corrected").glob("*.py")] + \
    ["scripts/41_phase9_corrected_driver.py", "scripts/42_phase9_corrected_launch.sh",
     "src/phase9_execution.py", "src/pathway_scoring.py", "src/score_standardisation.py",
     # Decision 13 (B1): the expanded mandatory inventory
     "scripts/43_phase9_corrected_verify.py", "scripts/44_phase9_freeze_prepare.sh",
     "scripts/44_phase9_freeze_capture.py", "scripts/40_phase9_dc_interval.py",
     "src/__init__.py", "src/paths.py", "src/reproducibility.py"]


def make_repo(tmp_path, policy="fresh_compile_prefix", edit=None, freeze_edit=None):
    root = tmp_path / "repo"
    for rel in FROZEN:
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, root / rel)
    if edit:
        edit(root)
    freeze = {"schema": orch.FREEZE_SCHEMA,
              "code": {rel: hashlib.sha256((root / rel).read_bytes()).hexdigest()
                       for rel in FROZEN if (root / rel).exists()},
              "environment": orch._jsonable(orch.current_environment()), "bytecode_policy": policy}
    if freeze_edit:
        freeze_edit(freeze)
    fp = root / orch.FREEZE_RECORD_RELPATH
    fp.parent.mkdir(parents=True, exist_ok=True)
    fp.write_bytes(json.dumps(freeze, sort_keys=True).encode())
    return root, hashlib.sha256(fp.read_bytes()).hexdigest(), freeze


def owner_cfg(root, anchor, out_rel="out", **d13):
    raw = approved_test_config(root, out_rel)
    drv = hashlib.sha256((root / "scripts/41_phase9_corrected_driver.py").read_bytes()).hexdigest()
    raw["owner_decisions_outside_schema"] = {"D13_driver_freeze": {
        "approved": True, "date": "2026-11-06", "reference": "TEST FIXTURE",
        "driver_path": "scripts/41_phase9_corrected_driver.py",
        "freeze_record_sha256": anchor, "sha256": drv, **d13}}
    return CorrectedConfig(raw, root)


@pytest.fixture
def e1_flags(monkeypatch, tmp_path_factory):
    prefix = tmp_path_factory.mktemp("pyc_prefix")
    state = {"isolated": True, "dont_write_bytecode": True, "pycache_prefix": str(prefix)}
    monkeypatch.setattr(orch, "_runtime_flags", lambda: dict(state))
    return state


def test_t2_freeze_chain_passes_with_e1_flags(tmp_path, e1_flags):
    root, anchor, freeze = make_repo(tmp_path)
    info = REAL_VERIFY(owner_cfg(root, anchor))
    assert info["freeze_record_sha256"] == anchor and info["bytecode_policy"] == "fresh_compile_prefix"
    assert info["driver_sha256"] == freeze["code"]["scripts/41_phase9_corrected_driver.py"]


def test_t2_missing_package_file_refused(tmp_path, e1_flags):
    root, anchor, _ = make_repo(tmp_path)
    (root / "src/phase9_corrected/statistics.py").unlink()
    with pytest.raises(orch.CodeIdentityRefusal, match="package file set"):
        REAL_VERIFY(owner_cfg(root, anchor))


def test_t2_loaded_module_mismatch_refused(tmp_path, e1_flags):
    def edit(root):                         # disk copy (and freeze) differ from the LOADED module
        with open(root / "src/phase9_corrected/statistics.py", "a") as fh:
            fh.write("\n# differs from the module loaded in this process\n")
    root, anchor, _ = make_repo(tmp_path, edit=edit)
    with pytest.raises(orch.CodeIdentityRefusal, match="loaded module src.phase9_corrected.statistics"):
        REAL_VERIFY(owner_cfg(root, anchor))


@pytest.mark.parametrize("case,match", [
    ("unapproved", "not approved"), ("frh_bad", "freeze_record_sha256 is missing or not 64"),
    ("drv_bad", "sha256 \\(driver\\) is missing or not 64"), ("driver_path", "driver_path must be"),
    ("bytes_before_parse", "does not match the approved"), ("inventory_driver", "freeze inventory's driver entry"),
    ("driver_disk", "driver file differs"), ("config_in_inventory", "configuration files")])
def test_b2_verification_order(tmp_path, e1_flags, case, match):
    fe = {"inventory_driver": lambda f: f["code"].update({"scripts/41_phase9_corrected_driver.py": "1" * 64}),
          "config_in_inventory": lambda f: f["code"].update({"config/analysis/owner.json": "2" * 64})}.get(case)
    root, anchor, _ = make_repo(tmp_path, freeze_edit=fe)
    kw = {}
    if case == "unapproved":
        kw = {"approved": False}
    elif case == "frh_bad":
        kw = {"freeze_record_sha256": "A" * 64}                        # uppercase: not lowercase hex
    elif case == "drv_bad":
        kw = {"sha256": "xyz"}
    elif case == "driver_path":
        kw = {"driver_path": "scripts/other.py"}
    elif case == "bytes_before_parse":
        (root / orch.FREEZE_RECORD_RELPATH).write_text("not json at all")   # step 3 precedes 4
    cfg = owner_cfg(root, anchor, **kw)
    if case == "driver_disk":                                          # owner + inventory agree,
        with open(root / "scripts/41_phase9_corrected_driver.py", "a") as fh:   # disk differs
            fh.write("\n# changed\n")
    with pytest.raises(orch.CodeIdentityRefusal, match=match):
        REAL_VERIFY(cfg)


@pytest.mark.parametrize("change,match", [
    ({"isolated": False}, "-I -B"), ({"dont_write_bytecode": False}, "-I -B"),
    ({"pycache_prefix": None}, "pycache_prefix"), ("inside", "outside the repository"),
    ("nonempty", "existing, empty directory"), ("missing", "existing, empty directory")])
def test_e1_runtime_flag_checks(tmp_path, e1_flags, change, match):
    root, anchor, _ = make_repo(tmp_path)
    if isinstance(change, dict):
        e1_flags.update(change)
    elif change == "inside":
        d = root / "pyc_inside"
        d.mkdir()
        e1_flags["pycache_prefix"] = str(d)
    elif change == "nonempty":
        (Path(e1_flags["pycache_prefix"]) / "stray.pyc").write_bytes(b"x")
    else:
        e1_flags["pycache_prefix"] = str(tmp_path / "does_not_exist")
    with pytest.raises(orch.CodeIdentityRefusal, match=match):
        REAL_VERIFY(owner_cfg(root, anchor))


def test_t2_successful_run_corrected_with_real_code_identity(tmp_path, e1_flags, monkeypatch):
    root, anchor, _ = make_repo(tmp_path)
    cfg = owner_cfg(root, anchor, out_rel="corrected_out")
    B.update(bad={}, raise_at=None, calls=0, drift=0.0)
    monkeypatch.setattr(orch, "h4_check", lambda r, s, raise_on_fail=True: {"all_ok": True})
    monkeypatch.setattr(orch, "historical_preservation_check", lambda r: {"all_ok": True})
    monkeypatch.setattr(orch, "run_pre_execution_gate", synthetic_gate)
    monkeypatch.setattr(orch, "load_context", lambda r: None)
    monkeypatch.setattr(orch, "RealPipeline", lambda: None)
    monkeypatch.setattr(orch, "compute_outer_replicate", fake_record)
    res = orch.run_corrected(cfg, 0, 500, clock=("TEST_CLOCK", StepClock()))   # REAL verify + manifest
    assert res["n_records"] == 500
    start = json.loads((root / "corrected_out/halts/run0.inv0001.start.json").read_text())
    assert start["freeze_record_sha256"] == anchor and start["path"] == "fresh"
    manifest = json.loads((root / "corrected_out/manifest/input_manifest.json").read_text())
    assert manifest["code"]["scripts/41_phase9_corrected_driver.py"] == \
        hashlib.sha256((root / "scripts/41_phase9_corrected_driver.py").read_bytes()).hexdigest()


# ================================================================ DEC-A1: non-finite serialisation
def test_a1_nonfinite_tokens_round_trip_and_hash_stably():
    rec = {"a": float("nan"), "b": float("inf"), "c": float("-inf"), "d": 1.5}
    b = canonical_json(rec)
    assert b == b'{"a":NaN,"b":Infinity,"c":-Infinity,"d":1.5}'           # not RFC 8259
    back = json.loads(b)
    assert math.isnan(back["a"]) and back["b"] == math.inf and back["c"] == -math.inf
    assert canonical_json(back) == b and hashlib.sha256(canonical_json(back)).digest() == \
        hashlib.sha256(b).digest()
    with pytest.raises(ValueError):                                       # strict parse rejects
        json.loads(b, parse_constant=lambda c: (_ for _ in ()).throw(ValueError(c)))


def test_a1_v5_fails_closed_on_nan_and_passes_on_infinity():
    nan_rec = {"k": 0, "x": float("nan"), "dC_corrected": 0.1}
    res = v5_spot_check({0: nan_rec}, lambda k: dict(nan_rec), selected_k=(0,))
    assert res["status"] == "FAIL" and res["rows"][0]["differing_fields"] == ["x"]
    inf_rec = {"k": 0, "x": float("-inf"), "dC_corrected": 0.1}
    assert v5_spot_check({0: json.loads(canonical_json(inf_rec))}, lambda k: dict(inf_rec),
                         selected_k=(0,))["status"] == "PASS"


# ================================================================ DEC-E1: launcher (temporary repo, stub driver)
STUB_DRIVER = r'''"""HARMLESS STUB standing in for the scientific driver in launcher tests."""
import json, os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
(ROOT / "STUB_STARTED").write_text("1")
mode = os.environ.get("P9C_STUB_MODE", "ok")
if mode == "import_fail":
    raise ImportError("synthetic import failure")
import numpy, pandas, scipy, lifelines, pyarrow
from src.phase9_corrected import zz_canary
from src.phase9_corrected.config import CorrectedConfig
from src.phase9_corrected import orchestration as orch
cfg = CorrectedConfig(json.loads(Path(sys.argv[2]).read_text()), ROOT)
info = orch.verify_code_identity(cfg)
print(json.dumps({"cmd": sys.argv[1], "canary": zz_canary.VALUE, "policy": info["bytecode_policy"],
                  "isolated": sys.flags.isolated, "no_bytecode": sys.flags.dont_write_bytecode,
                  "prefix": sys.pycache_prefix, "numpy": numpy.__version__,
                  "lifelines": lifelines.__version__, "pyarrow": pyarrow.__version__}))
if mode == "write_pyc":
    Path(sys.pycache_prefix, "intruder.pyc").write_bytes(b"x")
'''


def plant_stale_pyc(src: Path, stale_source: str):
    """A TIMESTAMP-valid in-tree .pyc whose code differs from the source on disk."""
    st = src.stat()
    code = compile(stale_source, str(src), "exec")
    data = (importlib.util.MAGIC_NUMBER + (0).to_bytes(4, "little")
            + (int(st.st_mtime) & 0xFFFFFFFF).to_bytes(4, "little")
            + (st.st_size & 0xFFFFFFFF).to_bytes(4, "little") + marshal.dumps(code))
    pyc = Path(importlib.util.cache_from_source(str(src)))
    pyc.parent.mkdir(parents=True, exist_ok=True)
    pyc.write_bytes(data)
    return pyc


def py_env(tmpdir, mode="ok"):
    e = {"PATH": "/usr/bin:/bin", "HOME": os.environ.get("HOME", "/tmp"), "TMPDIR": str(tmpdir),
         "P9C_STUB_MODE": mode}
    return e


@pytest.fixture(scope="module")
def launch_repo(tmp_path_factory):
    """Temporary repository: package + frozen modules + launcher + STUB driver + canary,
    .venv -> the real virtual environment, a synthetic freeze generated under the exact
    launch flags, and a synthetic approved owner configuration."""
    base = tmp_path_factory.mktemp("launch")
    root = base / "repo"
    for rel in FROZEN:
        if rel == "scripts/41_phase9_corrected_driver.py":
            continue
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, root / rel)
    (root / "scripts/41_phase9_corrected_driver.py").write_text(STUB_DRIVER)
    (root / "src/phase9_corrected/zz_canary.py").write_text('VALUE = "source"\n')
    os.symlink(ROOT / ".venv", root / ".venv")
    py = str(root / ".venv/bin/python")
    gen_pfx = base / "gen_pfx"
    gen_pfx.mkdir()
    out = subprocess.run([py, "-I", "-B", "-X", f"pycache_prefix={gen_pfx}", "-c",
                          "import json,sys; sys.path.insert(0, sys.argv[1]);"
                          "from src.phase9_corrected import orchestration as o;"
                          "print(json.dumps(o._jsonable(o.current_environment())))", str(root)],
                         capture_output=True, text=True, env=py_env(base), check=True)
    files = [p.relative_to(root).as_posix() for p in (root / "src/phase9_corrected").glob("*.py")] + \
        [rel for rel in FROZEN if not rel.startswith("src/phase9_corrected/")]
    freeze = {"schema": orch.FREEZE_SCHEMA, "bytecode_policy": "fresh_compile_prefix",
              "environment": json.loads(out.stdout),
              "code": {f: hashlib.sha256((root / f).read_bytes()).hexdigest() for f in sorted(files)}}
    fp = root / "config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json"
    fp.parent.mkdir(parents=True)
    fp.write_bytes(json.dumps(freeze, sort_keys=True).encode())
    cfg = approved_test_config(root, "out")
    cfg["owner_decisions_outside_schema"] = {"D13_driver_freeze": {
        "approved": True, "date": "2026-11-06", "reference": "TEST FIXTURE",
        "driver_path": "scripts/41_phase9_corrected_driver.py",
        "sha256": freeze["code"]["scripts/41_phase9_corrected_driver.py"],
        "freeze_record_sha256": hashlib.sha256(fp.read_bytes()).hexdigest()}}
    (base / "owner.json").write_text(json.dumps(cfg))
    plant_stale_pyc(root / "src/phase9_corrected/zz_canary.py", 'VALUE = "stale"\n')
    return base, root


def launch(base, root, *args, cwd=None, tmpdir=None, mode="ok", config=None):
    (root / "STUB_STARTED").unlink(missing_ok=True)
    t = tmpdir or (base / f"tmp_{uuid.uuid4().hex[:6]}")
    if tmpdir is None:
        t.mkdir()
    r = subprocess.run(["/bin/bash", str(root / "scripts/42_phase9_corrected_launch.sh"),
                        *args[:1], str(config or base / "owner.json"), *args[1:]],
                       capture_output=True, text=True, cwd=str(cwd or base),
                       env=py_env(t, mode))
    return r, t


def pycache_files(root):
    return sorted(str(p) for p in root.rglob("*.pyc") if ".venv" not in p.parts)


def test_e1_control_stale_pyc_is_read_without_prefix(launch_repo):
    """Control: without pycache_prefix, Python DOES read the planted in-tree stale .pyc
    (so the next test is meaningful)."""
    base, root = launch_repo
    r = subprocess.run([str(root / ".venv/bin/python"), "-I", "-B", "-c",
                        "import sys; sys.path.insert(0, sys.argv[1]);"
                        "import src.phase9_corrected.zz_canary as c; print(c.VALUE)", str(root)],
                       capture_output=True, text=True, env=py_env(base))
    assert r.stdout.strip() == "stale", r.stderr


def test_e1_launcher_positive_run_flags_imports_and_fresh_compile(launch_repo):
    base, root = launch_repo
    before = pycache_files(root)
    r, t = launch(base, root, "validate")
    assert r.returncode == 0, r.stderr
    o = json.loads(r.stdout.strip().splitlines()[-1])
    assert o["canary"] == "source"                    # in-tree stale .pyc NOT used
    assert o["policy"] == "fresh_compile_prefix" and o["isolated"] == 1 and o["no_bytecode"] == 1
    assert o["prefix"] and not o["prefix"].startswith(str(root))
    assert o["numpy"] and o["lifelines"] and o["pyarrow"]             # imports work under -I
    assert pycache_files(root) == before              # no repository bytecode created
    assert list(t.iterdir()) == []                    # temporary directories cleaned up


def test_e1_launcher_independent_of_caller_cwd(launch_repo, tmp_path):
    base, root = launch_repo
    rel_cfg = os.path.relpath(base / "owner.json", tmp_path)
    r, _ = launch(base, root, "validate", cwd=tmp_path, config=rel_cfg)
    assert r.returncode == 0, r.stderr


def test_e1_launcher_preserves_python_failure_status(launch_repo):
    base, root = launch_repo
    r, t = launch(base, root, "validate", mode="import_fail")
    assert r.returncode == 1 and "synthetic import failure" in r.stderr
    assert list(t.iterdir()) == []


def test_e1_launcher_refuses_bytecode_written_despite_policy(launch_repo):
    base, root = launch_repo
    r, _ = launch(base, root, "validate", mode="write_pyc")
    assert r.returncode == 22 and "bytecode was written" in r.stderr


def mutate(base, root, how):
    """Apply a pre-Python failure condition; return an undo callable."""
    fp = root / "config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json"
    cfg_p = base / "owner.json"
    saved_fp, saved_cfg = fp.read_bytes(), cfg_p.read_text()
    extra = root / "src/phase9_corrected/zz_unexpected.py"
    launcher = root / "scripts/42_phase9_corrected_launch.sh"
    saved_launcher = launcher.read_bytes()
    pipe = root / "src/pathway_scoring.py"
    saved_pipe = pipe.read_bytes()

    def set_d13(**kw):
        c = json.loads(saved_cfg)
        c["owner_decisions_outside_schema"]["D13_driver_freeze"].update(kw)
        cfg_p.write_text(json.dumps(c))
    if how == "freeze_absent":
        fp.unlink()
    elif how == "owner_hash_missing":
        set_d13(freeze_record_sha256=None)
    elif how == "malformed_hash":
        set_d13(freeze_record_sha256="ABC")
    elif how == "unapproved":
        set_d13(approved=False)
    elif how == "freeze_mismatch":
        fp.write_bytes(saved_fp + b" ")
    elif how == "launcher_modified":
        launcher.write_bytes(saved_launcher + b"\n# modified after freeze\n")
    elif how == "driver_mismatch":
        set_d13(sha256="0" * 64)
    elif how == "unexpected_file":
        extra.write_text("X = 1\n")
    elif how == "missing_file":
        pipe.unlink()

    def undo():
        fp.write_bytes(saved_fp)
        cfg_p.write_text(saved_cfg)
        launcher.write_bytes(saved_launcher)
        pipe.write_bytes(saved_pipe)
        extra.unlink(missing_ok=True)
    return undo


@pytest.mark.parametrize("how,code", [
    ("freeze_absent", 11), ("owner_hash_missing", 12), ("malformed_hash", 13), ("unapproved", 12),
    ("freeze_mismatch", 14), ("launcher_modified", 16), ("driver_mismatch", 17),
    ("unexpected_file", 18), ("missing_file", 16)])
def test_e1_launcher_fails_closed_before_python(launch_repo, how, code):
    base, root = launch_repo
    undo = mutate(base, root, how)
    try:
        r, t = launch(base, root, "validate")
    finally:
        undo()
    assert r.returncode == code, (how, r.returncode, r.stderr)
    assert "REFUSED" in r.stderr
    assert not (root / "STUB_STARTED").exists()        # Python never started
    assert list(t.iterdir()) == []


def test_e1_launcher_mktemp_failure_and_missing_python(launch_repo, tmp_path):
    base, root = launch_repo
    r, _ = launch(base, root, "validate", tmpdir=tmp_path / "does_not_exist")
    assert r.returncode == 19 and not (root / "STUB_STARTED").exists()
    venv = root / ".venv"
    venv.unlink()
    try:
        r2, _ = launch(base, root, "validate")
    finally:
        os.symlink(ROOT / ".venv", venv)
    assert r2.returncode == 21 and not (root / "STUB_STARTED").exists()


def test_e1_launcher_rejects_unknown_command(launch_repo):
    base, root = launch_repo
    r, _ = launch(base, root, "template")
    assert r.returncode == 2 and not (root / "STUB_STARTED").exists()


# ================================================================ review follow-ups (Stage 4)
@pytest.mark.parametrize("per_k_kind,expected", [
    ("missing_keys_all_true", "pending_verification"),      # incomplete evidence, not FAIL
    ("one_false", "disqualifier"),                           # PASS contradicted by per-k
    ("not_a_dict", "pending_verification"),                  # malformed -> incomplete
    ("value_not_dict", "pending_verification"),
    ("empty", "pending_verification"),
    ("complete", "disclosure")])
def test_assess_windows_per_k_classification(per_k_kind, expected):
    ks = list(range(120, 130))
    full = {str(k): {"record_identical": True, "dC_hex_equal": True} for k in ks}
    per_k = {"missing_keys_all_true": dict(list(full.items())[:-1]),
             "one_false": {**full, "125": {"record_identical": False, "dC_hex_equal": True}},
             "not_a_dict": [1, 2], "value_not_dict": {**full, "121": [True]},
             "empty": {}, "complete": full}[per_k_kind]
    w = {"invocation": 1, "region_sha256": "a" * 64, "k_list": ks}
    res = {"run_id": 0, "invocation": 1, "region_sha256": "a" * 64, "k_list": ks,
           "status": "PASS", "code_identity_verified": True, "per_k": per_k}
    R = orch._Reasons()
    orch._assess_windows(R, 0, [w], {1: res})
    cats = {c for c in ("disqualifier", "pending_verification", "disclosure") if R.codes(c)}
    assert cats == {expected}


def _inv(path, t, n_start=750 - 250, sizing="S" * 64, ref="OWNER-750", end=None):
    return {"start": {"path": path, "target_B_out": t, "n_start": n_start,
                      "sizing_decision_sha256": sizing, "topup_authorization": ref}, "end": end}


NOOP = {"outcome": "completed", "done_now": 0}


@pytest.mark.parametrize("later,ok", [
    ({}, True),
    ({2: _inv("resume", 500, end=NOOP)}, True),                                       # benign no-op
    ({2: _inv("topup", 750, n_start=500)}, True),
    ({2: _inv("topup", 750, n_start=500), 3: _inv("resume", 750, n_start=600)}, True),
    ({2: _inv("topup", 750, n_start=500, sizing="T" * 64)}, False),                    # other sizing
    ({2: _inv("topup", 750, n_start=499)}, False),                                     # n_start != L
    ({2: _inv("topup", 750, n_start=500, ref="  ")}, False),                           # blank ref
    ({2: _inv("topup", 750, n_start=500), 3: _inv("topup", 750, n_start=600)}, False), # 2nd top-up
    ({2: _inv("topup", 750, n_start=500), 3: _inv("resume", 500, end=NOOP)}, False),   # L_cur after
    ({2: _inv("topup", 750, n_start=500), 3: _inv("resume", 750, ref="OTHER")}, False),
    ({2: _inv("resume", 500, end={"outcome": "completed", "done_now": 3})}, False),    # growth
    ({2: _inv("topup", 1000, n_start=500)}, False)])                                   # target
def test_later_invocations_transition_model_direct(monkeypatch, later, ok):
    invs = {1: _inv("fresh", 500, n_start=0, end={"outcome": "completed", "done_now": 500}), **later}
    monkeypatch.setattr(orch, "invocation_state", lambda out, r: {"invocations": invs})
    why = orch._later_invocations_problem(Path("/nonexistent"), 0, 1, 500, 750, "S" * 64)
    assert (why is None) is ok, why


def test_b2_missing_driver_path_and_non_object_freeze_refused(tmp_path, e1_flags):
    root, anchor, _ = make_repo(tmp_path)
    cfg = owner_cfg(root, anchor)
    del cfg.raw["owner_decisions_outside_schema"]["D13_driver_freeze"]["driver_path"]
    with pytest.raises(orch.CodeIdentityRefusal, match="driver_path must be"):
        REAL_VERIFY(cfg)
    fp = root / orch.FREEZE_RECORD_RELPATH
    fp.write_bytes(b"[]")
    with pytest.raises(orch.CodeIdentityRefusal, match="schema invalid"):
        REAL_VERIFY(owner_cfg(root, hashlib.sha256(b"[]").hexdigest()))
