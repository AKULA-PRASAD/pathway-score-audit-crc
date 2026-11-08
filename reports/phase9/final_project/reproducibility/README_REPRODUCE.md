# Reproducing the Phase 9 final-project outputs

## What is reproduced

Every number, table, figure and document in `reports/phase9/final_project/`, computed from the **five frozen run checkpoints**. The nested bootstrap itself is **not** re-run. Its reproducibility was established separately by H5 (`reports/phase9/h5_spotcheck/`).

## Inputs (read-only; all verified by hash before use)

| Input | sha256/16 |
|---|---|
| `reports/phase9/dc_interval/run0_outer.jsonl` … `run4_outer.jsonl` | `885bbd1bd891b570`, `1bb27c56612a5f9d`, `0f12dce12bfce95e`, `4115bd4cd014d345`, `2202d559b0418392` |
| `reports/phase9/execution/phase9_execution_results.json` (source of θ̂) | `05d13b8b339fa7d4` |
| 17-file frozen/audit set (driver, plan, data, panels, freeze records, audit records) | listed in `analysis/common.py::EXPECTED_FILES` |

Checkpoint record structure: one JSON object per outer replicate, with fields `run_id, k, B_in_attempted, entropy, multiplicity_sha, n_distinct_outer, events_outer, censored_outer, failures, status, C_app{M1,M2}, age_knots, lambda_apparent{M1,M2}, dC_app, B_in_paired_complete, opt_M1_mean, opt_M2_mean, dOpt, dC_corrected, inner_C{M1,M2: [[C_boot, C_outer] × 100]}` (rounded to 8 dp), `dOpt_per_replicate_min/max`.

## Environment

The project interpreter `.venv/bin/python` is used (Python 3.13.13; numpy 2.5.4, scipy 1.18.1, pandas 2.3.3, matplotlib 3.11.2, pytest 9.2.0, lifelines 0.30.3). The full capture is in `environment.json`. **No package was installed or upgraded.**

## Command

From the repository root:

```
zsh reports/phase9/final_project/reproducibility/run_all.sh
```

It refuses to start if a Phase 9 driver process is running. It then runs these stages in order, each writing only inside `reports/phase9/final_project/`:

| Stage | Script | Output |
|---|---|---|
| S0 | `analysis/s0_verify_starting_state.py` | `results/s0_starting_state.json` |
| S1 | `analysis/s1_interval_analysis.py` (pure functions in `analysis/stats_core.py`) | `results/s1_interval_results.json`, `tables/T1–T6*.csv` |
| V1 | `validation/independent_check.py` (separate implementation; doesn't import `stats_core`) | `validation/independent_check_results.json` |
| V2 | `validation/cv_seed_state_evidence.py` | `validation/cv_seed_state_evidence.json` |
| S2 | `analysis/s2_figures.py` | `figures/F1–F6.{png,svg}`, `figures/FIGURE_CAPTIONS.md` |
| S3 | `analysis/s3_build_docs.py` (renders `*.template.*`; no hand-typed numbers) | report, presentation, completion report, `tables/*.md`, `environment.json` |
| T | `pytest validation/tests` | unit (`test_stats_core.py`), integration (`test_integration.py`), documents (`test_documents.py`) |
| QC | `analysis/s9_final_qc.py` | `audit/final_qc.json`, `reproducibility/MANIFEST_SHA256.txt` |

The project's own regression suite (`tests/`, 526 tests) can be run with `.venv/bin/python -B -m pytest -p no:cacheprovider -q`.

## Expected results

- S0 `all_ok: true`
- independent check: all values agree (166 / 166 at the time of writing)
- pytest: 0 failures
- final QC: `all_ok: true`

All outputs are deterministic given the inputs. The only fields that change between reruns are timestamps (`generated_utc`, `captured_utc`) and, through them, the hashes of the rendered documents.

## Validation procedure for a reviewer

1. Check the input hashes above (S0 does this, and so does `MANIFEST_SHA256.txt`).
2. Compare `results/s1_interval_results.json` with `validation/independent_check_results.json`. The independent check recomputes every interval, z₀, kernel MC-SE, between-run SD, EP-1 ratio, W3/W5/W6 value and table entry by a different code path.
3. Read `documentation/STATISTICAL_CONVENTIONS.md`: every choice the plan leaves open is listed there.
4. Read `documentation/CV_SEED_DEVIATION_DRAFT_v2.md` before relying on any result. All ΔC_corr-derived results are **PROVISIONAL**.

## Random-seed documentation

- No random numbers are drawn by any final-project stage. Every quantity is a deterministic function of the stored checkpoints.
- The test suite's synthetic data use fixed `default_rng(seed)` seeds 1–8.
- The bootstrap's seed hierarchy, and its deviation from plan §8.2, are documented in report §8 and §14 and in the deviation record.
