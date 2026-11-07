"""Phase 9 corrected pipeline -- Path C sizing-decision writer
(docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md, section 7).

Usage (from the repository root, with the BLAS/OpenMP thread variables unset):
    P=$(mktemp -d /tmp/p9c_pyc.XXXXXX)
    .venv/bin/python -I -B -X pycache_prefix=$P scripts/45_phase9_sizing_decision.py \
        <owner-config.json> <B_cur>

What it does
  1. Requires the launch flags of the scientific launcher (-I -B -X pycache_prefix=<fresh,
     empty, external directory>) and an approved configuration.
  2. Runs the frozen code-identity check (orchestration.verify_code_identity): the freeze
     anchor, every frozen source file, the loaded package modules and the environment.
     This script is NOT part of the Decision 13 freeze (adding it would change the freeze
     record and therefore the configuration hash that every checkpoint is bound to), so the
     check is run with this entry-point module excluded; its own SHA-256 is printed instead.
  3. Confirms that the stored input manifest still matches the current inputs, code and
     environment, and that the historical artifacts are intact.
  4. Computes the sizing decision at B_cur with the frozen
     orchestration.build_sizing_decision (the single source of truth for Decisions 5A, 5B,
     6A, 6B, 7A and 7B), which refuses unless all five runs are ELIGIBLE at B_cur.
  5. Writes <output_dir>/sizing/level<B_cur>.json ONCE (orchestration.durable_create; an
     existing file is never overwritten) and ONLY when the decision is `top_up`. Any other
     decision is printed and nothing is written.
  6. Re-reads the file and checks it for every run with the frozen
     orchestration.verify_sizing_decision, exactly as a top-up dispatch would.

Why the file needs no further trust: the sizing file is an input, not an authority. Before
any top-up is dispatched, run_corrected re-derives every value in it from the checkpoints
with the frozen code (verify_sizing_decision), and the post-run runner does so again. A
file that differs in any respect is refused. A top-up additionally requires the owner's
explicit --topup-authorization reference for each run.

It never dispatches, resumes or tops up replicates, and never writes a checkpoint, manifest,
configuration or invocation record.

Exit status: 0 top-up decision written and verified; 2 usage; 3 decision is not a top-up
(nothing written); 4 refused (identity, manifest, eligibility or existing file).
"""
import contextlib
import hashlib
import json
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))                    # -I removes the script directory

from src.phase9_corrected import orchestration as orch  # noqa: E402
from src.phase9_corrected.config import ConfigurationError, load_config  # noqa: E402
from src.phase9_corrected.integrity import (HardFailureH4, build_manifest,  # noqa: E402
                                            write_manifest_once)

SELF_RELPATH = "scripts/45_phase9_sizing_decision.py"


def fail(code: int, msg: str):
    print(f"p9c-sizing: REFUSED (code {code}): {msg}", file=sys.stderr)
    raise SystemExit(code)


@contextlib.contextmanager
def _entry_point_excluded():
    """Hide ONLY this entry-point module from the loaded-module identity check. Every
    package module, every other loaded repository module and every frozen file is still
    checked. The original module is restored afterwards."""
    main = sys.modules.get("__main__")
    if main is None or Path(getattr(main, "__file__", "")).resolve() != Path(__file__).resolve():
        fail(4, "unexpected entry-point module")
    sys.modules["__main__"] = types.ModuleType("__main__")
    try:
        yield
    finally:
        sys.modules["__main__"] = main


def main(argv) -> int:
    if len(argv) != 2 or not argv[1].isdigit():
        fail(2, "usage: 45_phase9_sizing_decision.py <owner-config.json> <B_cur>")
    B_cur = int(argv[1])
    fl = orch._runtime_flags()
    if not (fl.get("isolated") is True and fl.get("dont_write_bytecode") is True
            and isinstance(fl.get("pycache_prefix"), str) and fl.get("pycache_prefix")):
        fail(4, "not launched with -I -B -X pycache_prefix=<fresh, empty, external directory>")
    self_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    try:
        config = load_config(argv[0], ROOT)
        config.validate_for_execution()
    except ConfigurationError as ex:
        fail(4, f"configuration: {ex}")
    out = Path(config.output_dir)
    try:
        with _entry_point_excluded():
            freeze_info = orch.verify_code_identity(config)
        mpath = out / "manifest" / "input_manifest.json"
        mraw = mpath.read_bytes()
        manifest = json.loads(mraw)
        with _entry_point_excluded():
            orch.verify_code_identity(config, manifest=manifest)
        msha = write_manifest_once(mpath, build_manifest(config, ROOT))   # compare-only: it exists
        if msha != orch._sha(mraw):
            fail(4, "stored manifest bytes changed during the check")
        orch.h4_check(ROOT, "sizing")
        if not orch.historical_preservation_check(ROOT).get("all_ok"):
            fail(4, "historical preservation check failed")
    except (orch.CodeIdentityRefusal, HardFailureH4, OSError, ValueError) as ex:
        fail(4, f"identity or manifest: {ex}")
    path = orch._sizing_path(out, B_cur)
    if path.exists():
        fail(4, f"{path} already exists; it is never overwritten")
    try:
        art = orch.build_sizing_decision(config, B_cur, msha, freeze_info)
    except (orch.DispatchRefused, ValueError, KeyError) as ex:
        fail(4, f"sizing decision refused: {ex}")
    dec = art["decision"]
    summary = {"B_cur": B_cur, "action": dec["action"], "combined_B_next": dec["combined_B_next"],
               "rho_max": (dec.get("ep1") or {}).get("rho_max"),
               "freeze_record_sha256": freeze_info["freeze_record_sha256"],
               "manifest_sha256": msha, "config_sha256": config.sha256(),
               "writer": {"path": SELF_RELPATH, "sha256": self_sha}}
    if dec["action"] != "top_up":
        print(json.dumps({**summary, "written": False,
                          "note": "only a top_up decision is written (procedure section 7, item 4)"},
                         indent=1))
        return 3
    path.parent.mkdir(parents=True, exist_ok=True)
    orch.durable_create(path, art)
    target = dec["combined_B_next"]
    checks = {}
    try:
        for r in config.run_ids:
            checks[str(r)] = orch.verify_sizing_decision(config, r, B_cur, target, freeze_info, msha)
    except orch.DispatchRefused as ex:
        fail(4, f"written file did not pass the frozen verification: {ex}")
    print(json.dumps({**summary, "written": True, "path": str(path.relative_to(ROOT)),
                      "sizing_decision_sha256": orch._sha(path.read_bytes()),
                      "verified_for_runs": sorted(checks)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
