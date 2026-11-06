"""Decision 13 P2/P3/B2 -- freeze-candidate preparation, environment identity and the
in-process loaded-module identity check.

SYNTHETIC ONLY. Every subprocess runs in a TEMPORARY repository copy whose scientific
driver is a harmless stub; candidates are written to temporary directories outside it.
The real freeze candidate is never generated, the repository's freeze record is never
created, and no scientific operation is invoked."""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import types
import uuid
from pathlib import Path

import pytest

from src.phase9_corrected import orchestration as orch
from src.phase9_corrected import verification as V
from tests.p9c_support import approved_test_config

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_REQUIRED = {
    "scripts/41_phase9_corrected_driver.py", "scripts/42_phase9_corrected_launch.sh",
    "scripts/43_phase9_corrected_verify.py", "scripts/44_phase9_freeze_prepare.sh",
    "scripts/44_phase9_freeze_capture.py", "scripts/40_phase9_dc_interval.py",
    "src/__init__.py", "src/paths.py", "src/reproducibility.py", "src/phase9_execution.py",
    "src/pathway_scoring.py", "src/score_standardisation.py"}
PKG = [f"src/phase9_corrected/{p.name}" for p in (ROOT / "src/phase9_corrected").glob("*.py")]


# ================================================================ inventory and environment (P3)
def test_required_files_are_explicit_and_complete():
    assert set(orch.REQUIRED_FILES) == EXPECTED_REQUIRED
    assert len(orch.REQUIRED_FILES) == len(EXPECTED_REQUIRED)
    assert V.FROZEN_DRIVER == orch.HISTORICAL_DRIVER_RELPATH
    assert orch.FREEZE_SCHEMA == "p9c-freeze-2"
    for rel in EXPECTED_REQUIRED:
        assert (ROOT / rel).is_file(), rel
    assert not any(r.startswith("config/") for r in orch.REQUIRED_FILES)


def test_environment_blocking_fields_and_record_only_separation(monkeypatch):
    monkeypatch.setattr(orch, "_runtime_flags", lambda: {"isolated": True, "dont_write_bytecode": True,
                                                         "pycache_prefix": "/tmp/x"})
    env = orch.current_environment()
    assert set(orch.ENV_KEYS) == set(env)
    exe = Path(sys.executable).resolve()
    assert env["executable"] == {"resolved_path": str(exe),
                                 "sha256": hashlib.sha256(exe.read_bytes()).hexdigest()}
    assert env["launch_flags"] == {"isolated": True, "dont_write_bytecode": True, "pycache_prefix_set": True}
    assert "macos_build" not in env and "jq_version" not in env           # record-only: never compared


# ================================================================ B2: loaded-module identity
@pytest.fixture
def b2(tmp_path):
    root = tmp_path / "repo"
    for rel in PKG + list(orch.REQUIRED_FILES):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, root / rel)
    code = {rel: hashlib.sha256((root / rel).read_bytes()).hexdigest() for rel in PKG + list(orch.REQUIRED_FILES)}
    return root, code, set(code)


def fake_module(monkeypatch, name, path):
    m = types.ModuleType(name)
    if path is not None:
        m.__file__ = str(path)
    monkeypatch.setitem(sys.modules, name, m)
    return m


def test_b2_baseline_passes_and_unimported_required_files_are_not_demanded(b2):
    root, code, req = b2
    orch.verify_loaded_modules(root, code, req)


def test_b2_historical_driver_alias_checked_against_script_40(b2, monkeypatch):
    root, code, req = b2
    fake_module(monkeypatch, orch.HISTORICAL_DRIVER_ALIAS, root / orch.HISTORICAL_DRIVER_RELPATH)
    orch.verify_loaded_modules(root, code, req)                            # identical: passes
    with open(root / orch.HISTORICAL_DRIVER_RELPATH, "a") as fh:
        fh.write("\n# modified\n")
    with pytest.raises(orch.CodeIdentityRefusal, match="loaded module p9_dc_driver_frozen_ro"):
        orch.verify_loaded_modules(root, code, req)


def test_b2_alias_loaded_from_elsewhere_or_without_source(b2, monkeypatch, tmp_path):
    root, code, req = b2
    other = tmp_path / "elsewhere.py"
    other.write_text("X = 1\n")
    fake_module(monkeypatch, orch.HISTORICAL_DRIVER_ALIAS, other)           # outside the repo
    with pytest.raises(orch.CodeIdentityRefusal, match="does not match the freeze"):
        orch.verify_loaded_modules(root, code, req)
    fake_module(monkeypatch, orch.HISTORICAL_DRIVER_ALIAS, root / "scripts/41_phase9_corrected_driver.py")
    with pytest.raises(orch.CodeIdentityRefusal, match="resolves to scripts/41"):
        orch.verify_loaded_modules(root, code, req)
    fake_module(monkeypatch, orch.HISTORICAL_DRIVER_ALIAS, None)
    with pytest.raises(orch.CodeIdentityRefusal, match="has no source file"):
        orch.verify_loaded_modules(root, code, req)


def test_b2_repository_local_module_outside_inventory_fails_closed(b2, monkeypatch):
    root, code, req = b2
    extra = root / "src/unlisted_helper.py"
    extra.write_text("Y = 2\n")
    fake_module(monkeypatch, "src.unlisted_helper", extra)
    with pytest.raises(orch.CodeIdentityRefusal, match="not a required inventory file"):
        orch.verify_loaded_modules(root, code, req)


def test_b2_repository_local_module_hash_mismatch_and_missing_inventory_entry(b2, monkeypatch):
    root, code, req = b2
    fake_module(monkeypatch, "src.paths", root / "src/paths.py")
    orch.verify_loaded_modules(root, code, req)
    original = (root / "src/paths.py").read_bytes()
    with open(root / "src/paths.py", "a") as fh:
        fh.write("\n# changed\n")
    with pytest.raises(orch.CodeIdentityRefusal, match="loaded module src.paths"):
        orch.verify_loaded_modules(root, code, req)
    (root / "src/paths.py").write_bytes(original)
    code2 = {k: v for k, v in code.items() if k != "src/reproducibility.py"}
    fake_module(monkeypatch, "src.reproducibility", root / "src/reproducibility.py")
    with pytest.raises(orch.CodeIdentityRefusal, match="loaded module src.reproducibility"):
        orch.verify_loaded_modules(root, code2, req)


def test_b2_script_run_as_main_is_matched_by_path(b2, monkeypatch):
    root, code, req = b2
    fake_module(monkeypatch, "__p9c_test_main__", root / "scripts/43_phase9_corrected_verify.py")
    orch.verify_loaded_modules(root, code, req)
    fake_module(monkeypatch, "__p9c_test_main__", root / "scripts/99_unlisted.py")
    (root / "scripts/99_unlisted.py").write_text("pass\n")
    with pytest.raises(orch.CodeIdentityRefusal, match="not a required inventory file"):
        orch.verify_loaded_modules(root, code, req)


def test_b2_site_packages_and_stdlib_are_not_repository_sources(b2, monkeypatch):
    root, code, req = b2
    sp = root / ".venv/lib/python3.x/site-packages/thirdparty.py"
    sp.parent.mkdir(parents=True)
    sp.write_text("Z = 3\n")
    fake_module(monkeypatch, "thirdparty", sp)
    orch.verify_loaded_modules(root, code, req)                            # not repository code
    import json as stdlib_json
    assert stdlib_json.__file__                                           # loaded stdlib ignored


def test_b2_package_module_loaded_from_elsewhere_is_compared(b2, monkeypatch):
    import src.phase9_corrected.statistics  # noqa: F401  (imported lazily by the package)
    root, code, req = b2
    with open(root / "src/phase9_corrected/statistics.py", "a") as fh:
        fh.write("\n# differs from the loaded module\n")
    code = {**code, "src/phase9_corrected/statistics.py":
            hashlib.sha256((root / "src/phase9_corrected/statistics.py").read_bytes()).hexdigest()}
    with pytest.raises(orch.CodeIdentityRefusal, match="loaded module src.phase9_corrected.statistics"):
        orch.verify_loaded_modules(root, code, req)


def test_b2_is_part_of_verify_code_identity(monkeypatch):
    src = Path(orch.__file__).read_text()
    body = src[src.index("def verify_code_identity"):src.index("def _records_by_k")]
    assert "verify_loaded_modules(root, code, required)" in body
    assert "required = pkg_disk | set(REQUIRED_FILES)" in body


# ================================================================ P2: capture entry point (in-process)
def load_capture(root):
    spec = __import__("importlib.util").util.spec_from_file_location(
        f"p9c_capture_{uuid.uuid4().hex[:6]}", root / "scripts/44_phase9_freeze_capture.py")
    mod = __import__("importlib.util").util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_capture_is_strictly_non_scientific_and_writes_outside(b2, monkeypatch, tmp_path):
    root, code, _ = b2
    cap = load_capture(root)
    monkeypatch.setattr(cap, "ROOT", root)

    def forbidden(name):
        def f(*a, **k):
            raise AssertionError(f"capture called {name}")
        return f
    for name in ("run_corrected", "execute_replicates", "run_pre_execution_gate", "load_context",
                 "RealPipeline", "compute_outer_replicate", "build_manifest", "h4_check",
                 "historical_preservation_check", "verify_code_identity", "durable_create"):
        monkeypatch.setattr(orch, name, forbidden(name))
    monkeypatch.setattr(orch, "_runtime_flags", lambda: {"isolated": True, "dont_write_bytecode": True,
                                                         "pycache_prefix": str(tmp_path)})
    monkeypatch.setattr(orch, "verify_loaded_modules", lambda r, c, req: None)  # pytest process
    out = tmp_path / "cand.json"
    assert cap.main([str(out)]) == 0
    cand = json.loads(out.read_text())
    assert cand["code"] == code and cand["schema"] == orch.FREEZE_SCHEMA
    assert cand["bytecode_policy"] == "fresh_compile_prefix" and "CANDIDATE" in cand["status"]
    assert set(cand["environment"]) == set(orch.ENV_KEYS)
    assert set(cand["environment_record_only"]) == {"macos_build", "jq_version", "platform_string"}
    assert not any(k.startswith("config/") for k in cand["code"])
    assert not (root / orch.FREEZE_RECORD_RELPATH).exists()
    with pytest.raises(SystemExit) as ex:                                   # never overwritten
        cap.main([str(out)])
    assert ex.value.code == 5
    src = (ROOT / "scripts/44_phase9_freeze_capture.py").read_text()
    for name in ("run_corrected(", "load_context(", "RealPipeline(", "run_pre_execution_gate(",
                 "compute_outer_replicate(", "fit_and_score(", "verify_code_identity("):
        assert name not in src


# ================================================================ P2: preparation wrapper (subprocess)
STUB_DRIVER = r'''"""HARMLESS STUB standing in for the scientific driver."""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
(ROOT / "STUB_STARTED").write_text("1")
from src.phase9_corrected.config import CorrectedConfig
from src.phase9_corrected import orchestration as orch
cfg = CorrectedConfig(json.loads(Path(sys.argv[2]).read_text()), ROOT)
info = orch.verify_code_identity(cfg)
print(json.dumps({"verified": info["freeze_record_sha256"], "policy": info["bytecode_policy"]}))
'''


def penv(tmpdir):
    return {"PATH": "/usr/bin:/bin", "HOME": str(tmpdir), "TMPDIR": str(tmpdir)}


@pytest.fixture(scope="module")
def prep_repo(tmp_path_factory):
    base = tmp_path_factory.mktemp("prep")
    root = base / "repo"
    for rel in PKG + list(orch.REQUIRED_FILES):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, root / rel)
    (root / "scripts/41_phase9_corrected_driver.py").write_text(STUB_DRIVER)
    os.symlink(ROOT / ".venv", root / ".venv")
    return base, root


def prepare(base, root, out, script="scripts/44_phase9_freeze_prepare.sh", cwd=None):
    t = base / f"t_{uuid.uuid4().hex[:6]}"
    t.mkdir()
    r = subprocess.run(["/bin/bash", str(root / script), str(out)], capture_output=True, text=True,
                       cwd=str(cwd or base), env=penv(t))
    return r, t


def repo_pyc(root):
    return sorted(str(p) for p in root.rglob("*.pyc") if ".venv" not in p.parts)


def test_prepare_bootstrap_end_to_end_with_stub_driver(prep_repo):
    base, root = prep_repo
    out = base / f"cand_{uuid.uuid4().hex[:6]}.json"
    before = repo_pyc(root)
    r, t = prepare(base, root, out)
    assert r.returncode == 0, r.stderr
    printed = [l for l in r.stdout.splitlines() if l.strip().startswith("sha256:")][0].split()[-1]
    assert printed == hashlib.sha256(out.read_bytes()).hexdigest()
    assert "NOT approved" in r.stdout and "all re-checked with shasum" in r.stdout
    cand = json.loads(out.read_text())
    assert set(cand["code"]) == set(PKG) | set(orch.REQUIRED_FILES)
    assert cand["environment"]["launch_flags"] == {"isolated": True, "dont_write_bytecode": True,
                                                   "pycache_prefix_set": True}
    ro = cand["environment_record_only"]
    assert ro["jq_version"]["stdout"].startswith("jq") and "ProductVersion" in ro["macos_build"]["stdout"]
    assert repo_pyc(root) == before and list(t.iterdir()) == []
    assert not (root / orch.FREEZE_RECORD_RELPATH).exists()                 # never installed
    # owner steps simulated in the TEMPORARY repository: install, anchor, launch the stub
    fp = root / orch.FREEZE_RECORD_RELPATH
    fp.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(out, fp)
    cfg = approved_test_config(root, "out")
    cfg["owner_decisions_outside_schema"] = {"D13_driver_freeze": {
        "approved": True, "date": "2026-11-06", "reference": "TEST FIXTURE",
        "driver_path": orch.DRIVER_RELPATH, "sha256": cand["code"][orch.DRIVER_RELPATH],
        "freeze_record_sha256": printed}}
    (base / "owner.json").write_text(json.dumps(cfg))
    t2 = base / f"t_{uuid.uuid4().hex[:6]}"
    t2.mkdir()
    try:
        lr = subprocess.run(["/bin/bash", str(root / "scripts/42_phase9_corrected_launch.sh"), "validate",
                             str(base / "owner.json")], capture_output=True, text=True, cwd=str(base),
                            env=penv(t2))
        assert lr.returncode == 0, lr.stderr                                  # environment matches
        assert json.loads(lr.stdout.strip().splitlines()[-1])["verified"] == printed
    finally:
        fp.unlink()


@pytest.mark.parametrize("case,code", [("inside_repo", 19), ("exists", 11), ("no_dir", 19),
                                       ("missing_required", 4)])
def test_prepare_refusals(prep_repo, case, code):
    base, root = prep_repo
    out = base / f"cand_{uuid.uuid4().hex[:6]}.json"
    restore = None
    if case == "inside_repo":
        out = root / "cand.json"
    elif case == "exists":
        out.write_text("{}")
    elif case == "no_dir":
        out = base / "missing_dir" / "cand.json"
    elif case == "missing_required":
        p = root / "src/paths.py"
        restore = (p, p.read_bytes())
        p.unlink()
    try:
        r, t = prepare(base, root, out)
    finally:
        if restore:
            restore[0].write_bytes(restore[1])
    assert r.returncode == code, (case, r.returncode, r.stderr)
    assert "REFUSED" in r.stderr
    if case != "exists":
        assert not out.exists()
    assert not (root / "cand.json").exists() and list(t.iterdir()) == []


def test_capture_refuses_without_launch_flags_and_relative_output(prep_repo):
    base, root = prep_repo
    py = str(root / ".venv/bin/python")
    cap = str(root / "scripts/44_phase9_freeze_capture.py")
    r = subprocess.run([py, cap, str(base / "x.json")], capture_output=True, text=True, env=penv(base))
    assert r.returncode == 3 and not (base / "x.json").exists()
    pfx = base / f"pfx_{uuid.uuid4().hex[:6]}"
    pfx.mkdir()
    r = subprocess.run([py, "-I", "-B", "-X", f"pycache_prefix={pfx}", cap, "relative.json"],
                       capture_output=True, text=True, cwd=str(base), env=penv(base))
    assert r.returncode == 5 and not (base / "relative.json").exists()


def test_prepare_never_touches_the_real_repository():
    assert not (ROOT / orch.FREEZE_RECORD_RELPATH).exists()



# ================================================================ review remediation: record-only and thread environment
def test_platform_blocking_field_has_no_os_version_and_full_string_is_record_only(monkeypatch):
    import platform as pf
    env = orch.current_environment()
    assert env["platform"] == {"system": pf.system(), "machine": pf.machine(),
                               "architecture": pf.architecture()[0]}
    ver = pf.mac_ver()[0]
    if ver:
        assert ver not in json.dumps(env["platform"])
    src = (ROOT / "scripts/44_phase9_freeze_capture.py").read_text()
    assert '"platform_string": platform.platform()' in src


def _identity_repo(tmp_path, monkeypatch):
    root = tmp_path / "repo"
    for rel in PKG + list(orch.REQUIRED_FILES):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, root / rel)
    pfx = tmp_path / "pfx"
    pfx.mkdir()
    monkeypatch.setattr(orch, "_runtime_flags", lambda: {"isolated": True, "dont_write_bytecode": True,
                                                         "pycache_prefix": str(pfx)})

    def write_freeze(record_only=None):
        fr = {"schema": orch.FREEZE_SCHEMA, "bytecode_policy": "fresh_compile_prefix",
              "code": {r: hashlib.sha256((root / r).read_bytes()).hexdigest() for r in PKG + list(orch.REQUIRED_FILES)},
              "environment": orch._jsonable(orch.current_environment()),
              "environment_record_only": record_only or {"macos_build": "synthetic", "jq_version": "synthetic"}}
        fp = root / orch.FREEZE_RECORD_RELPATH
        fp.parent.mkdir(parents=True, exist_ok=True)
        fp.write_bytes(json.dumps(fr, sort_keys=True).encode())
        raw = approved_test_config(root, "out")
        raw["owner_decisions_outside_schema"] = {"D13_driver_freeze": {
            "approved": True, "date": "2026-11-06", "reference": "TEST FIXTURE", "driver_path": orch.DRIVER_RELPATH,
            "sha256": fr["code"][orch.DRIVER_RELPATH], "freeze_record_sha256": hashlib.sha256(fp.read_bytes()).hexdigest()}}
        from src.phase9_corrected.config import CorrectedConfig
        return CorrectedConfig(raw, root)
    return root, write_freeze


def test_thread_environment_must_match_between_capture_and_launch(tmp_path, monkeypatch):
    real_verify = orch.verify_code_identity
    for v in orch.THREAD_ENV_VARS:
        monkeypatch.delenv(v, raising=False)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    root, write_freeze = _identity_repo(tmp_path, monkeypatch)
    cfg = write_freeze()                                       # "captured" with OMP_NUM_THREADS=1
    real_verify(cfg)                                           # same shell settings: accepted
    monkeypatch.setenv("OMP_NUM_THREADS", "2")
    with pytest.raises(orch.CodeIdentityRefusal, match="environment thread_env differs"):
        real_verify(cfg)
    monkeypatch.delenv("OMP_NUM_THREADS")
    with pytest.raises(orch.CodeIdentityRefusal, match="environment thread_env differs"):
        real_verify(cfg)


def test_record_only_fields_never_block(tmp_path, monkeypatch):
    real_verify = orch.verify_code_identity
    root, write_freeze = _identity_repo(tmp_path, monkeypatch)
    cfg = write_freeze(record_only={"macos_build": {"stdout": "ProductVersion: 1.0 (a different build)"},
                                    "jq_version": {"stdout": "jq-0.0"}, "platform_string": "macOS-1.0-x"})
    assert real_verify(cfg)["bytecode_policy"] == "fresh_compile_prefix"



# ================================================================ manifest platform identity (integrity.py)
def _manifest_cfg():
    from src.phase9_corrected.config import CorrectedConfig
    return CorrectedConfig(approved_test_config(ROOT, "scratch_unused"), ROOT)


def _fake_manifest(monkeypatch, system="Darwin", machine="arm64", arch="64bit", full="macOS-15.3-arm64-arm-64bit"):
    from src.phase9_corrected import integrity as I
    monkeypatch.setattr(I.platform, "system", lambda: system)
    monkeypatch.setattr(I.platform, "machine", lambda: machine)
    monkeypatch.setattr(I.platform, "architecture", lambda *a, **k: (arch, ""))
    monkeypatch.setattr(I.platform, "platform", lambda *a, **k: full)
    monkeypatch.setattr(I, "check_hashes", lambda files, root: {"files": {"in.csv": {"sha256": "1" * 64}}})
    return I.build_manifest(_manifest_cfg(), ROOT)


def test_manifest_blocking_platform_is_os_family_and_architecture(monkeypatch):
    m = _fake_manifest(monkeypatch)
    assert m["schema"] == "p9c-manifest-2"
    assert m["environment"]["platform"] == {"system": "Darwin", "machine": "arm64", "architecture": "64bit"}
    assert m["record_only"] == {"platform_string": "macOS-15.3-arm64-arm-64bit"}
    assert "15.3" not in json.dumps(m["environment"])


def test_macos_point_update_alone_does_not_change_blocking_identity(monkeypatch, tmp_path):
    from src.phase9_corrected import integrity as I
    p = tmp_path / "manifest/input_manifest.json"
    h1 = I.write_manifest_once(p, _fake_manifest(monkeypatch))
    stored = p.read_bytes()
    h2 = I.write_manifest_once(p, _fake_manifest(monkeypatch, full="macOS-15.3.1-arm64-arm-64bit"))
    assert h1 == h2 == hashlib.sha256(stored).hexdigest() and p.read_bytes() == stored   # unchanged on disk


@pytest.mark.parametrize("change", [{"system": "Linux"}, {"machine": "x86_64"}, {"arch": "32bit"}])
def test_os_family_or_architecture_change_is_h4(monkeypatch, tmp_path, change):
    from src.phase9_corrected import integrity as I
    p = tmp_path / "m.json"
    I.write_manifest_once(p, _fake_manifest(monkeypatch))
    with pytest.raises(I.HardFailureH4):
        I.write_manifest_once(p, _fake_manifest(monkeypatch, **change))


def test_other_manifest_fields_and_stored_file_remain_strict(monkeypatch, tmp_path):
    from src.phase9_corrected import integrity as I
    p = tmp_path / "m.json"
    m = _fake_manifest(monkeypatch)
    I.write_manifest_once(p, m)
    for mutate in (lambda x: x["environment"].update(numpy="0.0"), lambda x: x["code"].update(extra="0" * 64),
                   lambda x: x.update(config_sha256="0" * 64), lambda x: x["inputs"].update({"in.csv": "2" * 64})):
        m2 = json.loads(json.dumps(m))
        mutate(m2)
        with pytest.raises(I.HardFailureH4):
            I.write_manifest_once(p, m2)
    stored = p.read_bytes()
    p.write_bytes(stored + b" ")                                           # stored file tampered
    with pytest.raises(I.HardFailureH4):
        I.write_manifest_once(p, m)
    p.write_bytes(b"[]")
    with pytest.raises(I.HardFailureH4):
        I.write_manifest_once(p, m)
