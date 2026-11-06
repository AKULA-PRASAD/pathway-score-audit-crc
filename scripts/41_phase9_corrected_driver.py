"""Phase 9 corrected nested-bootstrap driver (new file; the frozen driver
scripts/40_phase9_dc_interval.py is untouched).

Usage:
  python scripts/41_phase9_corrected_driver.py validate <config.json>
  python scripts/41_phase9_corrected_driver.py template
  python scripts/41_phase9_corrected_driver.py run <config.json> <run_id> <target_B_out> [--topup-authorization REF] [--resume-authorization REF]

`run` refuses to start unless every owner decision is supplied and approved in the
configuration, execution is authorised, H4 passes and all historical artifacts are intact.

--resume-authorization REF (owner Decision 10): required to continue a run whose latest
invocation ended in a runtime-budget halt, or has no terminal record (interrupted, cause
unknown). It is refused when nothing is pending, single-use, and recorded durably.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.phase9_corrected.config import (ConfigurationError, load_config,  # noqa: E402
                                         template_config)


def main(argv) -> int:
    if not argv or argv[0] not in ("validate", "template", "run"):
        print(__doc__)
        return 2
    if argv[0] == "template":
        print(json.dumps(template_config(), indent=1))
        return 0
    resume_ref = None
    if argv[0] == "run" and "--resume-authorization" in argv:
        i = argv.index("--resume-authorization")
        if i + 1 >= len(argv) or argv[i + 1].startswith("--") or not argv[i + 1].strip():
            print("--resume-authorization requires a non-blank REF")
            print(__doc__)
            return 2
        resume_ref = argv[i + 1]
    cfg = load_config(argv[1], ROOT)
    if argv[0] == "validate":
        probs = cfg.validate_structure() + cfg.missing_owner_items()
        print(json.dumps({"executable": not probs, "problems": probs}, indent=1))
        return 0 if not probs else 3
    from src.phase9_corrected.orchestration import run_corrected
    ref = None
    if "--topup-authorization" in argv:
        i = argv.index("--topup-authorization")
        if i + 1 >= len(argv) or argv[i + 1].startswith("--") or not argv[i + 1].strip():
            print("--topup-authorization requires a non-blank REF")
            print(__doc__)
            return 2
        ref = argv[i + 1]
    try:
        res = run_corrected(cfg, int(argv[2]), int(argv[3]), topup_authorization=ref,
                            resume_authorization=resume_ref)
    except ConfigurationError as ex:
        print(f"REFUSED: {ex}")
        return 3
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
