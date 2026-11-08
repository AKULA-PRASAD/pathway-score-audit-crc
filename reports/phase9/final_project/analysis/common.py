"""Shared constants and READ-ONLY loaders for the Phase 9 final-project analysis.

Nothing in this module writes to disk. All paths outside reports/phase9/final_project/
are opened read-only.
"""
from __future__ import annotations

import hashlib, json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
CKPT_DIR = ROOT / "reports" / "phase9" / "dc_interval"
FP = ROOT / "reports" / "phase9" / "final_project"

BASE = 2026110300          # plan v2.1 section 8.1
B_OUT, B_IN, N_RUNS = 500, 100, 5

# 17-file verification set (driver FROZEN set + cohort membership + driver + plan + audit records)
EXPECTED_FILES = {
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
    "data/processed/analysis/cohort_membership.csv": "bd3671f8df64c594",
    "scripts/40_phase9_dc_interval.py": "4b1738a4d0ee0e39",
    "docs/PHASE9_DELTAC_INTERVAL_PLAN.md": "d27c52a29930a338",
    "reports/phase9/dc_interval/HALT_RECORD.json": "c1b404309bd63744",
    "reports/phase9/dc_interval/RESUME_AUTHORIZATION.json": "0571c10e02e9059d",
}
EXPECTED_CKPT = {0: "885bbd1bd891b570", 1: "1bb27c56612a5f9d", 2: "0f12dce12bfce95e",
                 3: "4115bd4cd014d345", 4: "2202d559b0418392"}


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_run(r: int) -> list[dict]:
    """Records of run r, sorted by k (the canonical replicate index)."""
    p = CKPT_DIR / f"run{r}_outer.jsonl"
    rows = [json.loads(l) for l in p.read_text().splitlines() if l.strip()]
    return sorted(rows, key=lambda x: x["k"])


def run_values(r: int) -> np.ndarray:
    """dC_corrected(k), k = 0..499 in k order, for run r."""
    rows = load_run(r)
    assert [x["k"] for x in rows] == list(range(B_OUT)), f"run {r}: k not 0..499"
    return np.array([x["dC_corrected"] for x in rows], dtype=np.float64)


def theta_hat() -> dict:
    """theta_hat = DeltaC_corrected on the original data, as fixed in plan v2.1 header (+0.0043).

    Source: reports/phase9/execution/phase9_execution_results.json (frozen),
    monte_carlo.M2_primary.corrected_C_mean - monte_carlo.M1.corrected_C_mean,
    which is the definition used by docs/PHASE9_EXECUTION_REPORT.md section 'Optimism-corrected
    Delta-C' (+0.0043). It is itself a Monte-Carlo estimate (5 seeds x B=500).
    """
    d = json.loads((ROOT / "reports/phase9/execution/phase9_execution_results.json").read_text())
    m2 = d["monte_carlo"]["M2_primary"]; m1 = d["monte_carlo"]["M1"]
    per_seed = [p2["corrected_C"] - p1["corrected_C"] for p1, p2 in zip(m1["per_seed"], m2["per_seed"])]
    opt_seed = [p2["mean_optimism_C"] - p1["mean_optimism_C"] for p1, p2 in zip(m1["per_seed"], m2["per_seed"])]
    sd = float(np.std(per_seed, ddof=1))
    dC_app_D = d["apparent"]["M2_primary"]["C"] - d["apparent"]["M1"]["C"]
    return {"value": m2["corrected_C_mean"] - m1["corrected_C_mean"],
            "C_corrected_M2_primary": m2["corrected_C_mean"], "C_corrected_M1": m1["corrected_C_mean"],
            "MC_SE_M2_primary": m2["monte_carlo_SE"], "MC_SE_M1": m1["monte_carlo_SE"],
            "per_seed_dC_corrected": per_seed, "between_seed_SD": sd,
            "MC_SE_of_mean": sd / len(per_seed) ** 0.5,
            "dC_app_D": dC_app_D, "dOpt_D": float(np.mean(opt_seed)), "per_seed_dOpt": opt_seed,
            "source": "reports/phase9/execution/phase9_execution_results.json "
                      "monte_carlo.M2_primary.corrected_C_mean - monte_carlo.M1.corrected_C_mean"}
