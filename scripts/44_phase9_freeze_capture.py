"""Decision 13 P2 -- CAPTURE-ONLY freeze-candidate generator (NON-SCIENTIFIC).

Run ONLY through scripts/44_phase9_freeze_prepare.sh, which starts it as
    .venv/bin/python -I -B -X pycache_prefix=<fresh, empty, external directory> \
        scripts/44_phase9_freeze_capture.py <candidate.json>
i.e. with exactly the interpreter and launch flags of the scientific launcher, so that the
environment it records is the environment the launcher will later be compared against.

It reads repository SOURCE FILES and interpreter/package metadata, and nothing else. It
never calls the driver, the pre-execution gate, run_corrected, data loading, fitting,
tuning, resampling or any other scientific operation; it never reads data/ or reports/,
never writes inside the repository and never reads or writes the owner configuration.
Importing the package (for `current_environment`) only defines functions.

Output: a CANDIDATE freeze record (schema orchestration.FREEZE_SCHEMA), created
exclusively (O_EXCL) at a path OUTSIDE the repository, and its SHA-256 on stdout. The
candidate is NOT a trust anchor: it becomes usable only after the owner independently
re-hashes and reviews it, copies it to config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json
and records both D13 hashes with an approval (docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md).

Exit status: 0 written; 2 usage; 3 not launched with the approved flags; 4 inventory
problem; 5 output path refused; 6 loaded-module identity refused.
"""
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))                    # -I removes the script directory

from src.phase9_corrected import orchestration as orch  # noqa: E402

CANDIDATE_STATUS = ("CANDIDATE - NOT APPROVED. Not a trust anchor until the owner has "
                    "independently re-hashed and reviewed it (Decision 13).")


def fail(code: int, msg: str):
    print(f"p9c-freeze-capture: REFUSED (code {code}): {msg}", file=sys.stderr)
    raise SystemExit(code)


def _tool_output(argv) -> dict:
    """Record-only metadata from a sealed system tool; absence is recorded, not fatal."""
    try:
        r = subprocess.run(argv, capture_output=True, text=True, timeout=30,
                           env={"PATH": "/usr/bin:/bin", "LC_ALL": "C"})
        return {"argv": list(argv), "returncode": r.returncode, "stdout": r.stdout.strip()}
    except (OSError, subprocess.SubprocessError) as ex:
        return {"argv": list(argv), "returncode": None, "error": type(ex).__name__}


def record_only_environment() -> dict:
    """Decision 13 P3 record-only fields: visible in the freeze and in review output, never
    compared at launch."""
    import platform
    return {"macos_build": _tool_output(["/usr/bin/sw_vers"]),
            "jq_version": _tool_output(["/usr/bin/jq", "--version"]),
            "platform_string": platform.platform()}


def inventory(root: Path) -> dict:
    files = sorted({f"{orch.PACKAGE_DIR}/{p.name}" for p in (root / orch.PACKAGE_DIR).glob("*.py")}
                   | set(orch.REQUIRED_FILES))
    code = {}
    for rel in files:
        if rel.startswith("config/") or rel.startswith("/") or ".." in rel.split("/"):
            fail(4, f"inventory path not permitted: {rel}")
        p = root / rel
        if not p.is_file() or p.is_symlink():
            fail(4, f"required file missing or not a regular file: {rel}")
        code[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
    return code


def build_candidate(root: Path, code: dict, utc: str) -> dict:
    return {"schema": orch.FREEZE_SCHEMA,
            "status": CANDIDATE_STATUS,
            "code": code,
            "environment": orch._jsonable(orch.current_environment()),
            "environment_record_only": record_only_environment(),
            "bytecode_policy": orch.E1_BYTECODE_POLICY,
            "launch": {"launcher": orch.LAUNCHER_RELPATH,
                       "interpreter": ".venv/bin/python",
                       "flags": ["-I", "-B", "-X", "pycache_prefix=<fresh, empty, external directory>"]},
            "generated": {"utc": utc, "by": orch.FREEZE_CAPTURE_RELPATH},
            "excludes": "the owner configuration and its hash (the owner config anchors this "
                        "record, never the reverse)"}


def main(argv) -> int:
    if len(argv) != 1:
        fail(2, "usage: 44_phase9_freeze_capture.py <candidate-output.json>")
    fl = orch._runtime_flags()
    if not (fl.get("isolated") is True and fl.get("dont_write_bytecode") is True
            and isinstance(fl.get("pycache_prefix"), str) and fl.get("pycache_prefix")):
        fail(3, "not launched with -I -B -X pycache_prefix (use scripts/44_phase9_freeze_prepare.sh)")
    out = Path(argv[0])
    if not out.is_absolute():
        fail(5, "the output path must be absolute")
    parent = out.parent.resolve()
    target = parent / out.name
    root = ROOT.resolve()
    if target == root or root in target.parents:
        fail(5, "the candidate must be written OUTSIDE the repository")
    if not parent.is_dir():
        fail(5, f"output directory does not exist: {parent}")
    code = inventory(root)
    try:
        orch.verify_loaded_modules(root, code, set(code))
    except orch.CodeIdentityRefusal as ex:
        fail(6, str(ex))
    cand = build_candidate(root, code, datetime.now(timezone.utc).isoformat())
    data = (json.dumps(cand, sort_keys=True, indent=1) + "\n").encode()
    try:
        fd = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except OSError as ex:
        fail(5, f"cannot create {target} (it must not exist): {ex}")
    try:
        os.write(fd, data)
        os.fsync(fd)
    finally:
        os.close(fd)
    print(json.dumps({"candidate": str(target), "sha256": hashlib.sha256(data).hexdigest(),
                      "n_files": len(code), "status": CANDIDATE_STATUS}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
