# PROJECT COMPLETION REPORT: Phase 9 nested-bootstrap interval

Rendered by `analysis/s3_build_docs.py` from `results/`, `validation/` and `audit/final_qc.json`; no number is hand-typed.

> **Bottom line.**
> - **Computation:** complete and verified.
> - **Statistics:** computed, but **PROVISIONAL**, because every ΔC_corr-derived result depends on M2 λ tuning carried out under an unadjudicated CV-seed protocol deviation.
> - **Methodological and scientific completion:** **not reached**. Four decisions remain (§15).

## 1. Initial project status (verified, not assumed)

- **Five checkpoints:** each has 500 unique records with k = 0–499, all with status `ok`; 100 / 100 paired-complete; 0 inner failures; 136 / 316 events / censored; hashes as recorded.
- **Frozen artifacts:** the driver's `check_frozen()` returned PASS, and 17 / 17 frozen/audit hashes matched.
- **H5 evidence:** hashes intact; independently rechecked at 20 / 20 bit-identical, so H5 was **not** rerun.
- **CV-seed deviation:** draft v1, unadjudicated.

| run | records | unique k | k = 0–499 | status ok | paired-complete min | inner failures | events / censored | sha256/16 |
|---|---|---|---|---|---|---|---|---|
| 0 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `885bbd1bd891b570` |
| 1 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `1bb27c56612a5f9d` |
| 2 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `0f12dce12bfce95e` |
| 3 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `4115bd4cd014d345` |
| 4 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `2202d559b0418392` |

## 2. Completed computational work

- **Starting-state verification (S0)**, including a recheck of record arithmetic: ΔC_app = C_app(M2) − C_app(M1) exactly; ΔOpt from `inner_C` within 2.7e-09, from the 8-dp rounding; ΔC_corr = ΔC_app − ΔOpt exactly.
- **§9.1 analysis (S1):** for every run, the percentile, basic and z₀ / BC intervals; kernel MC-SE at three bandwidths; the §6.2–6.5 accounting; the §6.4 bounding analysis; W3, W5, W6 and W7. Across runs, the between-run SD, EP-1 and W4.
- **Independent re-computation:** 166 / 166 values agree.
- **Seed-state evidence** for the deviation.
- **Six figures, tables, the report, the presentation, the reproducibility package and final QC.**
- **Not computed:** no new bootstrap replicate was computed, and no original artifact was modified.

## 3. H5 validation and evidence

**PASS**, 20 / 20 bit-identical, 20 / 20 full records (`reports/phase9/h5_spotcheck/h5_comparison.json`). This establishes the reproducibility of the **as-run** implementation. **It is not evidence of conformance with plan §8.2.**

## 4. CV-seed deviation and adjudication status

**UNADJUDICATED.** The reviewed record is `documentation/CV_SEED_DEVIATION_DRAFT_v2.md` (DRAFT — NOT YET ADJUDICATED); v1 is retained unchanged.
- **Mismatch:** the plan's `cvseed = inner[b].spawn(2)[0]` (line 364) against the driver's fresh `default_rng(ss_inner[b])` and `default_rng(ss_outer)` (lines 142, 117).
- **Shared starting state:** confirmed in 9 / 9 cases.
- **Affected:** M2 λ.
- **Not affected:** M1 (stored apparent λ = 0 in 2 500 / 2 500; the inner level is inferred from the code path).
- **Not performed:** no alternative-seed recomputation, which was not authorised.

## 5. §9.1 statistical results (PROVISIONAL)

θ̂ = +0.00432 (MC-SE 0.00019), taken from Phase 9 primary.

| run | median | mean | SD | % < 0 | percentile 95 % (primary) | basic 95 % (diag.) | BC z₀-only 95 % (diag.) | z₀ | #{< θ̂}/500 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | +0.01086 | +0.01400 | 0.01634 | 21.2 | [-0.00886, +0.05163] | [-0.04300, +0.01749] | [-0.01115, +0.03028] | -0.468 | 160 |
| 1 | +0.01137 | +0.01335 | 0.01572 | 20.0 | [-0.00906, +0.04687] | [-0.03823, +0.01769] | [-0.01146, +0.02728] | -0.519 | 151 |
| 2 | +0.01292 | +0.01421 | 0.01561 | 19.2 | [-0.00927, +0.04838] | [-0.03975, +0.01791] | [-0.01143, +0.02587] | -0.571 | 142 |
| 3 | +0.01180 | +0.01363 | 0.01641 | 23.2 | [-0.00975, +0.05133] | [-0.04269, +0.01839] | [-0.01252, +0.03270] | -0.429 | 167 |
| 4 | +0.01208 | +0.01490 | 0.01615 | 17.8 | [-0.00921, +0.05151] | [-0.04288, +0.01784] | [-0.01197, +0.02788] | -0.577 | 141 |

| quantity | run 0 | run 1 | run 2 | run 3 | run 4 | mean | SD (ddof = 1, df = 4) | range |
|---|---|---|---|---|---|---|---|---|
| percentile lower | -0.00886 | -0.00906 | -0.00927 | -0.00975 | -0.00921 | -0.00923 | 0.000334 | 0.000898 |
| percentile upper | +0.05163 | +0.04687 | +0.04838 | +0.05133 | +0.05151 | +0.04994 | 0.002186 | 0.004764 |
| basic lower | -0.04300 | -0.03823 | -0.03975 | -0.04269 | -0.04288 | -0.04131 | 0.002186 | 0.004764 |
| basic upper | +0.01749 | +0.01769 | +0.01791 | +0.01839 | +0.01784 | +0.01786 | 0.000334 | 0.000898 |
| BC z₀ lower | -0.01115 | -0.01146 | -0.01143 | -0.01252 | -0.01197 | -0.01171 | 0.000543 | 0.001374 |
| BC z₀ upper | +0.03028 | +0.02728 | +0.02587 | +0.03270 | +0.02788 | +0.02880 | 0.002699 | 0.006831 |
| z₀ | -0.46770 | -0.51866 | -0.57100 | -0.42889 | -0.57691 | -0.51263 | 0.064432 | 0.148016 |
| percentile width | +0.06049 | +0.05592 | +0.05766 | +0.06108 | +0.06072 | +0.05917 | 0.002269 | 0.005157 |
| median ΔC_corr(k) | +0.01086 | +0.01137 | +0.01292 | +0.01180 | +0.01208 | +0.01181 | 0.000776 | 0.002060 |
| mean ΔC_corr(k) | +0.01400 | +0.01335 | +0.01421 | +0.01363 | +0.01490 | +0.01402 | 0.000595 | 0.001553 |

- **H1:** 0 / 2 500.
- **H2:** vacuous pass (0 failures).
- **H6:** 0 mismatches.

## 6. EP-1 result

**UNRESOLVED (decision required).** The threshold of 0.05 was not changed.

| MC-SE estimator | evaluations | pass (≤ 0.05) | max MC-SE / width |
|---|---|---|---|
| kernel_bw_x0.5 | 10 | 9 | 0.0535 |
| kernel_bw_x1.0 | 10 | 10 | 0.0476 |
| kernel_bw_x2.0 | 10 | 10 | 0.0375 |
| between_run_SD | 10 | 10 | 0.0391 |

- **The one failure:** 1 / upper / kernel_bw_x0.5, ratio 0.0535.
- **If every estimator must pass:** the §2.4 ladder requires B_out = 750, which was **not executed**.
- **H3:** not triggered.
- **Scope:** EP-1 concerns Monte-Carlo precision only. It does not validate the interval.

## 7. Robustness and sensitivity results

| run | inner failures | unreliable outer (6.5) | H1 | bounding (6.4) | H2 | W3 SD inflation | W5 spread lo / hi (bound) | W6 A7 lo / hi | W7 |z₀| |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | no | VACUOUS | no | 0.35 % | 0.00053 / 0.00160 (0.00302) | 56.4 % / 56.4 % **breach** | 0.468 **breach** |
| 1 | 0 | 0 | no | VACUOUS | no | 0.36 % | 0.00026 / 0.00299 (0.00280) **breach** | 52.2 % / 52.2 % **breach** | 0.519 **breach** |
| 2 | 0 | 0 | no | VACUOUS | no | 0.36 % | 0.00089 / 0.00107 (0.00288) | 52.9 % / 52.9 % **breach** | 0.571 **breach** |
| 3 | 0 | 0 | no | VACUOUS | no | 0.34 % | 0.00001 / 0.00133 (0.00305) | 53.9 % / 53.9 % **breach** | 0.429 **breach** |
| 4 | 0 | 0 | no | VACUOUS | no | 0.34 % | 0.00040 / 0.00291 (0.00304) | 55.5 % / 55.5 % **breach** | 0.577 **breach** |

- **Warnings breached:** W6 / A7 and W7 / A8 in every run, and W5 in run 1 (upper).
- **No breach:** W3 and W4.
- **Not applicable:** W1 and W2 (0 failures), §6.3 (no failures to compare), and the S-7 nested interval (not in plan v2.1).
- **Blocked:** the §2.4 top-up.
- **Not authorised:** the CV-seed sensitivity recomputation.

## 8. Scientific interpretation (summary; full text in report §24)

- **Observed:**
  - every one of the 15 intervals (5 runs × 3 methods) includes 0;
  - the outer distribution is right-skewed and centred above θ̂, with mean +0.0140 against θ̂ +0.0043;
  - z₀ = -0.577 to -0.429.
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
| Unit, integration and document tests (`validation/tests`) | 40 passed in 0.67s |
| Project regression suite (`tests/`) | ======================= 526 passed, 1 warning in 41.34s ======================== |
| Independent re-computation | 166 / 166 agree |
| End-to-end | `reproducibility/run_all.sh` regenerates every output |

## 10. Errors detected and corrected

1. **H3 flag:** set to True whenever EP-1 was mixed, but H3 requires the ladder to be exhausted. Corrected; the v1 output is preserved in `audit/`.
2. **H6:** was a literal value; it is now computed from the records.
3. **Figures:** label and legend collisions in five figures and a stray label in a sixth. Fixed and re-inspected.
4. **Report prose:** six hand-typed numbers converted to computed placeholders, and one internal inconsistency fixed ("two" against four open decisions).
5. **Deviation record:** one line-range citation corrected in v2 (plan 367–372).

None of these affected a computed statistic. See `audit/EXECUTION_LEDGER.md`.

## 11. Independent verification

`validation/independent_check.py` shares no code with `stats_core.py`. It recomputes θ̂, all intervals, z₀ and its ties, the BC levels, the Silverman h, the 30 kernel SEs, the W3 / W5 / W6 values, the between-run SDs, the 40 EP-1 widths and the pass count, plus the CSV table entries: 166 / 166 agree.

## 12. Original artifact integrity (final QC)

**Overall QC: ALL CHECKS PASS** (`audit/final_qc.json`)

| Check | Result |
|---|---|
| Driver `check_frozen()` | True |
| 17 frozen/audit hashes unchanged | True |
| 5 checkpoint hashes unchanged | True |
| H5 evidence unchanged | True |
| No new files in `dc_interval/` | True |
| No file outside `final_project/` modified since 2026-11-04T16:24:03.655430+00:00 | True |
| Results computed on the verified checkpoints | True |
| Documents free of placeholders and labelled PROVISIONAL | True / True |

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
2. **EP-1 governing estimator.** If every prescribed MC-SE estimator must pass, EP-1 fails narrowly and §2.4 requires topping up each run to B_out = 750 (≈ 250 × 5 new outer replicates, appended to the checkpoints). The run-1 W5 breach also points to "raise B_out". Authorising a top-up is your decision.
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
| 14 | Final QC | PASS |
| — | Computational completion | **COMPLETE** |
| — | Methodological completion | **NOT COMPLETE**: decisions 1–2 |
| — | Scientific completion | **NOT COMPLETE**: decisions 1–4 |
