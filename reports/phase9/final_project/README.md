# Phase 9 final project: nested-bootstrap interval for the optimism-corrected ΔC

All ΔC_corr-derived results here are **PROVISIONAL — calculated using the as-run implementation; CV-seed protocol deviation not yet adjudicated.**

**Start here:** `PROJECT_COMPLETION_REPORT.md`, then `documentation/PHASE9_FINAL_REPORT.md`.

| Directory | Contents |
|---|---|
| `analysis/` | stage scripts: `s0` (verify inputs), `s1` (§9.1 statistics), `s2` (figures), `s3` (render documents), `s9` (final QC); `stats_core.py` (pure statistical functions); `common.py` (constants and read-only loaders) |
| `results/` | `s0_starting_state.json`, `s1_interval_results.json` (every computed quantity) |
| `tables/` | `T1–T6*.csv` (machine-readable), `*.md` (rendered into the report) |
| `figures/` | `F1–F6` as PNG and SVG, plus `FIGURE_CAPTIONS.md` |
| `validation/` | `independent_check.py` and its results, `cv_seed_state_evidence.py` and its evidence, and `tests/` (unit, integration, document tests) |
| `documentation/` | final report (rendered and template), CV-seed deviation record **v2 (DRAFT — NOT YET ADJUDICATED)**, statistical conventions |
| `presentation/` | self-contained HTML technical presentation (open locally; arrow keys to navigate) |
| `reproducibility/` | `README_REPRODUCE.md`, `run_all.sh`, `environment.json`, `MANIFEST_SHA256.txt` |
| `audit/` | `EXECUTION_LEDGER.md`, stage run logs, pytest logs, `final_qc.json`, superseded outputs kept for audit |

Original artifacts are referenced by hash, never copied or modified:
- `reports/phase9/dc_interval/` (five checkpoints, halt and resume records, driver log);
- `reports/phase9/h5_spotcheck/` (H5 evidence, deviation draft v1);
- `scripts/40_phase9_dc_interval.py` (frozen driver);
- `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` (plan v2.1).

Nothing here has been committed, pushed or published.
