#!/usr/bin/env python3
"""Phase 9 theta-hat alignment check (docs/PHASE9_THETA_ALIGNMENT_AUDIT.md). Read-only.

Usage: .venv/bin/python -I -B -X pycache_prefix=<empty dir> scripts/47_phase9_theta_alignment_check.py

Check: the corrected engine's fit_and_score on D, under the historical
apparent convention (fresh default_rng(20261102) per model, scripts/38:90-92), must
reproduce the apparent C values stored in phase9_execution_results.json."""
import sys, json, hashlib
from pathlib import Path
ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from src.phase9_corrected.data import load_context, LAM_GRID
from src.phase9_corrected.fitting import RealPipeline, fit_and_score
ex = json.loads((ROOT/"reports/phase9/execution/phase9_execution_results.json").read_text())
ctx = load_context(ROOT); p = RealPipeline()
print("n", len(ctx.ids), "events", int(ctx.e.sum()))
print("grid equal", [float(x) for x in LAM_GRID] == ex["lambda_grid"])
out = {}
for name, hist in (("M1", "M1"), ("M2", "M2_primary")):
    fo = fit_and_score(p, ctx, name, "apparent", ctx.ids, None, ctx.t, ctx.e, None, None,
                       np.random.default_rng(20261102), LAM_GRID)
    out[name] = fo.C_fit
    print(name, repr(fo.C_fit), repr(ex["apparent"][hist]["C"]), "lambda", fo.lam, ex["apparent"][hist]["lambda"],
          "identical", fo.C_fit == ex["apparent"][hist]["C"])
d = out["M2"] - out["M1"]
print("dC_app(D) corrected engine", repr(d))
mo = {m: np.mean([r["mean_optimism_C"] for r in ex["monte_carlo"][h]["per_seed"]]) for m, h in (("M1","M1"),("M2","M2_primary"))}
print("theta reconstructed = dC_app - dOpt(2500) =", repr(d - (mo["M2"] - mo["M1"])))
print("theta stored          =", repr(ex["monte_carlo"]["M2_primary"]["corrected_C_mean"] - ex["monte_carlo"]["M1"]["corrected_C_mean"]))
print("dOpt(D, B=2500) =", repr(mo["M2"] - mo["M1"]))
print("input_hashes:")
for k, v in ex["input_hashes"].items():
    pth = ROOT / k if (ROOT / k).exists() else None
    if pth is None:
        cands = list(ROOT.rglob(Path(k).name))[:1]; pth = cands[0] if cands else None
    cur = hashlib.sha256(pth.read_bytes()).hexdigest() if pth else None
    print("  ", k, "match" if cur == v else f"MISMATCH/unknown ({pth})")
