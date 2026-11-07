"""Refusal paths of the Path C sizing-decision writer (scripts/45_phase9_sizing_decision.py).

The decision content and its verification are the frozen build_sizing_decision /
verify_sizing_decision, already covered end to end with synthetic campaigns in
test_phase9_corrected_nonestimable_stop.py (write_sizing, test_valid_topup_and_progression)
and test_phase9_corrected_postrun.py. These tests cover only the writer's own guards and
never read study data or write inside the repository.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "45_phase9_sizing_decision.py"
PY = sys.executable


def run(args, flags=(), cwd=ROOT):
    return subprocess.run([PY, *flags, str(SCRIPT), *args], capture_output=True, text=True,
                          cwd=cwd, timeout=120)


def test_usage_is_refused_with_code_2(tmp_path):
    pfx = tmp_path / "pyc"
    pfx.mkdir()
    flags = ("-I", "-B", "-X", f"pycache_prefix={pfx}")
    for args in ([], ["only-one"], ["cfg.json", "five-hundred"], ["a", "500", "extra"]):
        r = run(args, flags)
        assert r.returncode == 2, (args, r.stderr)
        assert "usage" in r.stderr


def test_missing_launch_flags_are_refused_before_anything_is_read(tmp_path):
    r = run(["does-not-exist.json", "500"], flags=("-B",))
    assert r.returncode == 4
    assert "pycache_prefix" in r.stderr


def test_unapproved_or_missing_configuration_is_refused(tmp_path):
    pfx = tmp_path / "pyc"
    pfx.mkdir()
    r = run([str(tmp_path / "missing.json"), "500"], ("-I", "-B", "-X", f"pycache_prefix={pfx}"))
    assert r.returncode != 0
    assert not list(pfx.iterdir()), "bytecode was written despite -B"


def test_entry_point_exclusion_restores_main_module(monkeypatch):
    import importlib.util
    spec = importlib.util.spec_from_file_location("p9c_sizing_writer_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    fake_main = type(sys)("__main__")
    fake_main.__file__ = str(SCRIPT)
    monkeypatch.setitem(sys.modules, "__main__", fake_main)
    with mod._entry_point_excluded():
        assert sys.modules["__main__"] is not fake_main
        assert getattr(sys.modules["__main__"], "__file__", None) is None
    assert sys.modules["__main__"] is fake_main


def test_entry_point_exclusion_refuses_a_foreign_main_module(monkeypatch):
    import importlib.util
    import pytest
    spec = importlib.util.spec_from_file_location("p9c_sizing_writer_under_test2", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    other = type(sys)("__main__")
    other.__file__ = str(ROOT / "scripts" / "43_phase9_corrected_verify.py")
    monkeypatch.setitem(sys.modules, "__main__", other)
    with pytest.raises(SystemExit) as ei:
        with mod._entry_point_excluded():
            pass
    assert ei.value.code == 4
    assert sys.modules["__main__"] is other
