# PROJECT COMPLETION REPORT: Phase 9 nested-bootstrap interval

Rendered by `analysis/s3_build_docs.py` from `results/`, `validation/` and `audit/final_qc.json`; no number is hand-typed.

> **Bottom line.**
> - **Computation:** complete and verified.
> - **Statistics:** computed, but **PROVISIONAL**, because every ΔC_corr-derived result depends on M2 λ tuning carried out under an unadjudicated CV-seed protocol deviation.
> - **Methodological and scientific completion:** **not reached**. Four decisions remain (§15).

## 1. Initial project status (verified, not assumed)

- **Five checkpoints:** each has 500 unique records with k = 0–499, all with status `ok`; 100 / 100 paired-complete; 0 inner failures; 136 / 316 events / censored; hashes as recorded.
- **Frozen artifacts:** the driver's `check_frozen()` returned {{ S0['driver_check_frozen'] }}, and {{ sum(v['ok'] for v in S0['files_17'].values()) }} / 17 frozen/audit hashes matched.
- **H5 evidence:** hashes intact; independently rechecked at {{ S0['h5_evidence']['independent_recheck_bit_identical'] }} / 20 bit-identical, so H5 was **not** rerun.
- **CV-seed deviation:** draft v1, unadjudicated.

{{ TABLE:CKPT }}

## 2. Completed computational work

- **Starting-state verification (S0)**, including a recheck of record arithmetic: ΔC_app = C_app(M2) − C_app(M1) exactly; ΔOpt from `inner_C` within {{ S0['record_arithmetic_max_abs_error']['dOpt_from_inner'] | .1e }}, from the 8-dp rounding; ΔC_corr = ΔC_app − ΔOpt exactly.
- **§9.1 analysis (S1):** for every run, the percentile, basic and z₀ / BC intervals; kernel MC-SE at three bandwidths; the §6.2–6.5 accounting; the §6.4 bounding analysis; W3, W5, W6 and W7. Across runs, the between-run SD, EP-1 and W4.
- **Independent re-computation:** {{ IC['n_ok'] }} / {{ IC['n_checks'] }} values agree.
- **Seed-state evidence** for the deviation.
- **Six figures, tables, the report, the presentation, the reproducibility package and final QC.**
- **Not computed:** no new bootstrap replicate was computed, and no original artifact was modified.

## 3. H5 validation and evidence

**PASS**, {{ H5['n_bit_identical'] }} / 20 bit-identical, {{ H5['n_full_record_identical'] }} / 20 full records (`reports/phase9/h5_spotcheck/h5_comparison.json`). This establishes the reproducibility of the **as-run** implementation. **It is not evidence of conformance with plan §8.2.**

## 4. CV-seed deviation and adjudication status

**UNADJUDICATED.** The reviewed record is `documentation/CV_SEED_DEVIATION_DRAFT_v2.md` (DRAFT — NOT YET ADJUDICATED); v1 is retained unchanged.
- **Mismatch:** the plan's `cvseed = inner[b].spawn(2)[0]` (line 364) against the driver's fresh `default_rng(ss_inner[b])` and `default_rng(ss_outer)` (lines 142, 117).
- **Shared starting state:** confirmed in {{ sum(c['executed_cv_state_equals_resample_state'] for c in CV['cases']) }} / {{ len(CV['cases']) }} cases.
- **Affected:** M2 λ.
- **Not affected:** M1 (stored apparent λ = 0 in 2 500 / 2 500; the inner level is inferred from the code path).
- **Not performed:** no alternative-seed recomputation, which was not authorised.

## 5. §9.1 statistical results (PROVISIONAL)

θ̂ = {{ TH['value'] | +.5f }} (MC-SE {{ TH['MC_SE_of_mean'] | .5f }}), taken from Phase 9 primary.

{{ TABLE:T1 }}

{{ TABLE:T3 }}

- **H1:** 0 / 2 500.
- **H2:** vacuous pass (0 failures).
- **H6:** 0 mismatches.

## 6. EP-1 result

**UNRESOLVED (decision required).** The threshold of 0.05 was not changed.

{{ TABLE:T4 }}

- **The one failure:** {{ EP['failing_cases'][0]['run'] if EP['failing_cases'] else '—' }} / {{ EP['failing_cases'][0]['endpoint'] if EP['failing_cases'] else '—' }} / {{ EP['failing_cases'][0]['estimator'] if EP['failing_cases'] else '—' }}, ratio {{ EP['max_ratio'] | .4f }}.
- **If every estimator must pass:** the §2.4 ladder requires B_out = {{ EP['plan_2_4_ladder_if_strict']['ladder_rung'] }}, which was **not executed**.
- **H3:** not triggered.
- **Scope:** EP-1 concerns Monte-Carlo precision only. It does not validate the interval.

## 7. Robustness and sensitivity results

{{ TABLE:T5 }}

- **Warnings breached:** W6 / A7 and W7 / A8 in every run, and W5 in run 1 (upper).
- **No breach:** W3 and W4.
- **Not applicable:** W1 and W2 (0 failures), §6.3 (no failures to compare), and the S-7 nested interval (not in plan v2.1).
- **Blocked:** the §2.4 top-up.
- **Not authorised:** the CV-seed sensitivity recomputation.

## 8. Scientific interpretation (summary; full text in report §24)

- **Observed:**
  - every one of the 15 intervals (5 runs × 3 methods) includes 0;
  - the outer distribution is right-skewed and centred above θ̂, with mean {{ BD['mean_dC_corrected_k'] | +.4f }} against θ̂ {{ TH['value'] | +.4f }};
  - z₀ = {{ rng(None, "x['z0_diagnostic']['z0']", "+.3f") }}.
- **Supported:**
  - there is no evidence that M2 improves optimism-corrected discrimination over M1;
  - the bootstrap distribution is not centred on θ̂;
  - endpoint Monte-Carlo error is small next to the interval width.
- **Plausible but untested:** duplicate rows inflate M2's apparent ΔC on D*ₖ, and inner optimism evaluated within D*ₖ is smaller than on D. Plan §7.2 predicted the opposite direction.
- **Unresolved:** see §15.
- **Not established:** no clinical utility, mechanism, external validity or transportability.

## 9. Tests performed

| Suite | Result |
|---|---|
| Unit, integration and document tests (`validation/tests`) | {{ QC['checks'].get('pytest_final_project_log', 'see audit/pytest_final_project.log') }} |
| Project regression suite (`tests/`) | {{ QC['checks'].get('project_regression_log', 'see audit/pytest_project_regression.log') }} |
| Independent re-computation | {{ IC['n_ok'] }} / {{ IC['n_checks'] }} agree |
| End-to-end | `reproducibility/run_all.sh` regenerates every output |

## 10. Errors detected and corrected

1. **H3 flag:** set to True whenever EP-1 was mixed, but H3 requires the ladder to be exhausted. Corrected; the v1 output is preserved in `audit/`.
2. **H6:** was a literal value; it is now computed from the records.
3. **Figures:** label and legend collisions in five figures and a stray label in a sixth. Fixed and re-inspected.
4. **Report prose:** six hand-typed numbers converted to computed placeholders, and one internal inconsistency fixed ("two" against four open decisions).
5. **Deviation record:** one line-range citation corrected in v2 (plan 367–372).

None of these affected a computed statistic. See `audit/EXECUTION_LEDGER.md`.

## 11. Independent verification

`validation/independent_check.py` shares no code with `stats_core.py`. It recomputes θ̂, all intervals, z₀ and its ties, the BC levels, the Silverman h, the 30 kernel SEs, the W3 / W5 / W6 values, the between-run SDs, the 40 EP-1 widths and the pass count, plus the CSV table entries: {{ IC['n_ok'] }} / {{ IC['n_checks'] }} agree.

## 12. Original artifact integrity (final QC)

**Overall QC: {{ 'ALL CHECKS PASS' if QC['all_ok'] else ('NOT YET RUN' if QC['all_ok'] is None else 'FAILED: ' + ', '.join(QC['failed'])) }}** (`audit/final_qc.json`)

| Check | Result |
|---|---|
| Driver `check_frozen()` | {{ QC['checks'].get('driver_check_frozen') }} |
| 17 frozen/audit hashes unchanged | {{ QC['checks'].get('files_17_unchanged') }} |
| 5 checkpoint hashes unchanged | {{ QC['checks'].get('checkpoints_5_unchanged') }} |
| H5 evidence unchanged | {{ QC['checks'].get('h5_evidence_unchanged') }} |
| No new files in `dc_interval/` | {{ QC['checks'].get('dc_interval_no_new_files') }} |
| No file outside `final_project/` modified since {{ QC['checks'].get('work_order_start_utc', '?') }} | {{ QC['checks'].get('no_outside_modification') }} |
| Results computed on the verified checkpoints | {{ QC['checks'].get('results_on_verified_checkpoints') }} |
| Documents free of placeholders and labelled PROVISIONAL | {{ QC['checks'].get('docs_no_placeholders') }} / {{ QC['checks'].get('docs_labelled_provisional') }} |

## 13. Final deliverable locations

All under `reports/phase9/final_project/` (see `README.md`):
- `documentation/PHASE9_FINAL_REPORT.md`
- `documentation/CV_SEED_DEVIATION_DRAFT_v2.md`
- `documentation/STATISTICAL_CONVENTIONS.md`
- `presentation/PHASE9_TECHNICAL_PRESENTATION.html`
- `figures/`, `tables/`, `results/`, `validation/`, `reproducibility/`
- `audit/EXECUTION_LEDGER.md`

## 14. Outstanding limitations

- **The interval itself:** the CV-seed deviation; θ̂ from a different seeding convention; EP-1's governing estimator unspecified; the W5 breach (run 1); the A7 / A8 breaches in every run (by plan §7, a warning that neither interval should be trusted).
- **Conventions adopted where the plan is silent:** listed in `STATISTICAL_CONVENTIONS.md`.
- **The plan's standing limitations (§13–14):** conditional on 136 events; asymptotic bootstrap; dilution; df = 4; internal validation only.

## 15. Remaining user decisions (the only genuinely unresolved items)

1. **Adjudicate the CV-seed deviation.** Is plan §8.2 or the hash-locked driver binding, and is any sensitivity investigation wanted? Until this is decided, every result stays PROVISIONAL.
2. **EP-1 governing estimator.** If every prescribed MC-SE estimator must pass, EP-1 fails narrowly and §2.4 requires topping up each run to B_out = {{ EP['plan_2_4_ladder_if_strict']['ladder_rung'] }} (≈ 250 × 5 new outer replicates, appended to the checkpoints). The run-1 W5 breach also points to "raise B_out". Authorising a top-up is your decision.
3. **Reporting an interval that the diagnostics flag.** A7 and A8 fire in every run. Which interval, if any, should be presented as *the* result, and with what qualification? The plan also doesn't designate a single headline across the five runs.
4. **θ̂ centring.** Should θ̂ remain the Phase 9 primary value, or be re-derived by applying the nested driver's functional to D? Re-deriving it would be a new computation.

## Task-completion matrix

| # | Work-order stage | Status |
|---|---|---|
| 0 | Verify starting state | **PASS** |
| A | Deviation documentation (v2) | **PROVISIONAL** (draft; adjudication is the user's) |
| B 5.1 | Percentile intervals | **PROVISIONAL** (computed, verified) |
| B 5.2 | Basic intervals | **PROVISIONAL** (computed, verified) |
| B 5.3 | z₀ / BC intervals | **PROVISIONAL** (computed, verified) |
| B 5.4 | Kernel endpoint MC-SE | **PROVISIONAL** (computed, verified) |
| B 5.5 | Bounding analysis | **PASS** (vacuous) |
| B 5.6 | Across-run variability | **PROVISIONAL** (computed, verified) |
| B 5.7 | Deviation labelling | **PASS** |
| C | EP-1 | **BLOCKED**: 39 / 40, governing estimator undecided |
| D | Robustness inventory | **PASS** / NOT APPLICABLE / BLOCKED per item (report §23) |
| E | Scientific interpretation | **PROVISIONAL** |
| F | Final report | **PASS** (content PROVISIONAL) |
| G | Figures and tables | **PASS** |
| H | Reproducibility package | **PASS** |
| I | Technical presentation | **PASS** (content PROVISIONAL) |
| 14 | Final QC | {{ 'PASS' if QC['all_ok'] else 'see §12' }} |
| — | Computational completion | **COMPLETE** |
| — | Methodological completion | **NOT COMPLETE**: decisions 1–2 |
| — | Scientific completion | **NOT COMPLETE**: decisions 1–4 |
