"""Phase 9 -- nested bootstrap CI for the optimism-corrected Delta-C.

Implements docs/PHASE9_DELTAC_INTERVAL_PLAN.md v2.1 exactly.
Estimand: DeltaC_corrected = C_corrected(M2 four-set) - C_corrected(M1).

Reuses src/phase9_execution.py UNCHANGED. Writes only to reports/phase9/dc_interval/.

Usage:  python scripts/40_phase9_dc_interval.py <run_id 0..4> [B_out] [B_in]
"""
from __future__ import annotations

import hashlib, json, os, sys, time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import paths                                                  # noqa: E402
from src.phase9_execution import (AGE_COL, ConvergenceFailure, c_index,  # noqa: E402
                                  linear_predictor, pipeline_fit)
from src.score_standardisation import NonEstimableReplicate             # noqa: E402
from src.reproducibility import environment_report                      # noqa: E402

ANALYSIS = paths.PROCESSED_DIR / "analysis"
PANELS = paths.ROOT / "config" / "pathways"
OUT = paths.ROOT / "reports" / "phase9" / "dc_interval"
EV, TM = "rfs_event_full", "rfs_time_full"
M1_VARS = [AGE_COL, "sex", "tnm.stage", "tumor.location", "chemotherapy.adjuvant"]
LAM_GRID = np.geomspace(1e-3, 10.0, 12)

BASE = 2026110300            # v2.1 section 8.1 -- disjoint from the Phase 9 seeds
N_WORKERS = int(os.environ.get("P9_WORKERS", "8"))

FROZEN = {  # safeguard 1 / hard failure H4
 "docs/PHASE9_MODELING_PROTOCOL_DRAFT.md": "ab56708919f3746f",
 "config/pathways/m2_panel_4set.json": "43ce15b2f0e06f75",
 "config/pathways/M2_FREEZE_RECORD.json": "64a86fa8ed0d6a8a",
 "config/analysis/A19_FREEZE_RECORD.json": "a77de7b37fe09eeb",
 "src/pathway_scoring.py": "38cac035c0b2d450",
 "src/score_standardisation.py": "41b5c091125e1b11",
 "src/phase9_execution.py": "571d0c233f367be6",
 "data/processed/analysis/expression_median_S0.parquet": "a1ca41807c534584",
 "data/processed/analysis/covariates.csv": "47b2d9829a4fccec",
 "data/processed/analysis/endpoint.csv": "802a68575a34a3b6",
 "reports/phase9/execution/phase9_execution_results.json": "05d13b8b339fa7d4",
 "reports/phase9/execution_s7/phase9_s7_results.json": "b9202305bb6416e1",
}


def check_frozen():
    bad = [k for k, v in FROZEN.items()
           if not hashlib.sha256((ROOT / k).read_bytes()).hexdigest().startswith(v)]
    if bad:
        raise SystemExit(f"HARD FAILURE H4 -- frozen artifact changed: {bad}")


def strat_boot(ev_mask, rng):
    """Event-stratified patient bootstrap. Returns POSITIONAL indices with repeats."""
    idx = np.arange(len(ev_mask))
    out = [rng.choice(idx[ev_mask == lab], size=int((ev_mask == lab).sum()), replace=True)
           for lab in (1.0, 0.0)]
    return np.sort(np.concatenate(out))


_G = {}


def init_worker():
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    cov = pd.read_csv(ANALYSIS / "covariates.csv", keep_default_na=False,
                      dtype=str).set_index("sample_id")
    ids = [s for s in s3 if all(cov.loc[s, v] != "" for v in M1_VARS)]
    ep = pd.read_csv(ANALYSIS / "endpoint.csv").set_index("sample_id").reindex(ids)
    _G["ids"] = ids
    _G["cov"] = cov
    _G["t"] = ep[TM].to_numpy(float)
    _G["e"] = ep[EV].to_numpy(float)
    _G["X"] = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[ids].astype(np.float64)
    _G["panel"] = json.loads((PANELS / "m2_panel_4set.json").read_text())


def outer_replicate(args):
    """One outer replicate: apparent on D*_k + B_in paired inner optimism replicates."""
    run_id, k, B_in = args
    if not _G:
        init_worker()
    ids, cov, X = _G["ids"], _G["cov"], _G["X"]
    t, e, panel = _G["t"], _G["e"], _G["panel"]
    MODELS = {"M1": None, "M2": panel}

    ss_root = np.random.SeedSequence(entropy=[BASE, run_id])
    ss_outer = ss_root.spawn(k + 1)[k]                 # deterministic in k alone
    rng_o = np.random.default_rng(ss_outer)

    I = strat_boot(e, rng_o)                            # OUTER indices, with repeats
    tO, eO = t[I], e[I]
    oids = [ids[i] for i in I]
    mult = np.bincount(I, minlength=len(ids)).astype(int)

    rec = {"run_id": run_id, "k": int(k), "B_in_attempted": int(B_in),
           "entropy": [BASE, run_id], "multiplicity_sha": hashlib.sha256(
               mult.tobytes()).hexdigest()[:16],
           "n_distinct_outer": int((mult > 0).sum()),
           "events_outer": int(eO.sum()), "censored_outer": int((eO == 0).sum()),
           "failures": [], "status": "ok"}

    # ---- apparent on D*_k, both models, SAME outer sample ----
    app, knots, lams = {}, {}, {}
    for name, pan in MODELS.items():
        try:
            m, d_tr, _, dg = pipeline_fit(X, cov, pan, oids, oids, tO, eO,
                                          LAM_GRID, np.random.default_rng(ss_outer))
            app[name] = c_index(tO, eO, linear_predictor(m, d_tr))
            knots[name] = dg["age_knots"]; lams[name] = dg["lambda"]
        except (ConvergenceFailure, NonEstimableReplicate) as ex:
            rec["failures"].append({"stage": "apparent", "model": name,
                                    "reason": getattr(ex, "reason", str(ex))})
    if len(app) < 2:
        rec["status"] = "unreliable_apparent_failed"
        return rec
    rec["C_app"] = app; rec["age_knots"] = knots; rec["lambda_apparent"] = lams
    rec["dC_app"] = app["M2"] - app["M1"]

    # ---- inner optimism, evaluated on D*_k ----
    ss_inner = ss_outer.spawn(B_in)
    inner = {"M1": [], "M2": []}
    paired_ok = 0
    for b in range(B_in):
        rng_i = np.random.default_rng(ss_inner[b])
        J = strat_boot(eO, rng_i)                       # INNER indices into D*_k
        jids = [oids[i] for i in J]
        tJ, eJ = tO[J], eO[J]
        got = {}
        for name, pan in MODELS.items():
            try:
                m, d_tr, d_ap, _ = pipeline_fit(X, cov, pan, jids, oids, tJ, eJ,
                                                LAM_GRID, np.random.default_rng(ss_inner[b]))
                cb = c_index(tJ, eJ, linear_predictor(m, d_tr))
                co = c_index(tO, eO, linear_predictor(m, d_ap))
                got[name] = (cb, co)
            except (ConvergenceFailure, NonEstimableReplicate) as ex:
                rec["failures"].append({
                    "stage": "inner", "b": int(b), "model": name,
                    "reason": getattr(ex, "reason", str(ex)),
                    "n_distinct_inner": int(len(set(jids))),
                    "events_inner": int(eJ.sum())})
        # PAIRING-PRESERVING EXCLUSION (v2.1 section 6.2)
        if len(got) == 2:
            paired_ok += 1
            for name in MODELS:
                inner[name].append([got[name][0], got[name][1]])
    rec["B_in_paired_complete"] = int(paired_ok)
    if paired_ok < 0.95 * B_in:
        rec["status"] = "unreliable_inner_shortfall"
        return rec

    o1 = np.array(inner["M1"]); o2 = np.array(inner["M2"])
    opt1 = o1[:, 0] - o1[:, 1]; opt2 = o2[:, 0] - o2[:, 1]
    rec["opt_M1_mean"] = float(opt1.mean()); rec["opt_M2_mean"] = float(opt2.mean())
    rec["dOpt"] = float(opt2.mean() - opt1.mean())
    rec["dC_corrected"] = float(rec["dC_app"] - rec["dOpt"])
    rec["inner_C"] = {"M1": o1.round(8).tolist(), "M2": o2.round(8).tolist()}
    rec["dOpt_per_replicate_min"] = float((opt2 - opt1).min())
    rec["dOpt_per_replicate_max"] = float((opt2 - opt1).max())
    return rec


def main() -> int:
    check_frozen()
    run_id = int(sys.argv[1])
    B_out = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    B_in = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    OUT.mkdir(parents=True, exist_ok=True)
    ckpt = OUT / f"run{run_id}_outer.jsonl"
    done = set()
    if ckpt.exists():
        for line in ckpt.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["k"])
    todo = [(run_id, k, B_in) for k in range(B_out) if k not in done]
    print(f"[run {run_id}] B_out={B_out} B_in={B_in} workers={N_WORKERS} "
          f"resume: {len(done)} done, {len(todo)} to do", flush=True)

    t0 = time.time()
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    init_worker()                      # load once, before fork -> shared copy-on-write
    with ctx.Pool(N_WORKERS) as pool, open(ckpt, "a") as fh:
        for i, rec in enumerate(pool.imap_unordered(outer_replicate, todo, chunksize=1), 1):
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            if i % 50 == 0:
                el = time.time() - t0
                print(f"[run {run_id}] {i}/{len(todo)}  {el/60:.1f} min  "
                      f"eta {el/i*(len(todo)-i)/60:.1f} min", flush=True)
    check_frozen()
    print(f"[run {run_id}] COMPLETE in {(time.time()-t0)/60:.1f} min -> {ckpt.name}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
