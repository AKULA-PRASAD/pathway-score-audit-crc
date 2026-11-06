"""Phase 9 corrected pipeline -- post-run verification runner entry point (Decision 13 P1).

Start ONLY through the launcher, which verifies the frozen package before Python starts:
    scripts/42_phase9_corrected_launch.sh verify <owner-config.json>
The runner never dispatches, resumes or tops up bootstrap replicates and never writes a
checkpoint, manifest or sizing file (src/phase9_corrected/postrun.py). Its Stage A, V1-inner
and V5 steps perform real-data recomputation: executing it requires Decision 14.

Exit status: 0 PASS; 1 FAIL; 2 refused (nothing computed); 3 OWNER_REVIEW_REQUIRED;
4 hard failure H4; 5 hard failure H5 (V5).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.phase9_corrected.postrun import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:], ROOT))
