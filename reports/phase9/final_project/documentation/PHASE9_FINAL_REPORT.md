# Phase 9 — Nested-Bootstrap Confidence Interval for the Optimism-Corrected ΔC

**Final research report** · generated 2026-11-04 from `results/s1_interval_results.json` ·
plan `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1

> **Status of every ΔC_corr-derived result in this report: PROVISIONAL — calculated using the as-run
> implementation; CV-seed protocol deviation not yet adjudicated.**
>
> *Conditional on the observed 136 recurrence events and 316 censored observations in the 452-patient A-19
> population* (mandatory description, plan §4).

---

## 1. Title

Nested-bootstrap 95 % confidence interval for the optimism-corrected difference in Harrell's C between the M2
four-pathway model and the clinical M1 model (Phase 9, median representation, A-19 cohort).

## 2. Executive summary

| Area | Status | Evidence |
|---|---|---|
| Five nested-bootstrap runs (5 × 500 outer × 100 inner; 2 500 outer / 250 000 inner replicates) | **COMPLETE** | §10; checkpoints re-verified (§12) |
| H1: unreliable outer replicates > 5 % | **PASS**: 0 / 2 500 | §11 |
| H2: bounding analysis | **PASS, vacuous**: 0 inner failures, so the bounds equal the primary values | §20 |
| H4: frozen-artifact hashes | **PASS**: 17 / 17 files and 5 / 5 checkpoints unchanged, before and after | §12 |
| H5: reproducibility spot-check | **PASS**: 20 / 20 bit-identical (reproducibility of the **as-run** implementation only) | §13 |
| H6: paired-complete count M1 = M2 | **PASS**: 0 / 2 500 records unequal | §11 |
| CV-seed protocol deviation (plan §8.2 vs driver) | **UNADJUDICATED** | §14; deviation record v2 |
| §9.1 intervals, MC-SE, bounding, between-run SD | **COMPUTED: PROVISIONAL** | §16–§21 |
| EP-1 (MC-SE ≤ 0.05 × width) | **UNRESOLVED**: 39 / 40 evaluations pass; it passes under the named Silverman bandwidth and under the between-run SD, and fails once at bandwidth × 0.5 | §22 |
| H3: EP-1 unmet after exhausting the B_out ladder | **NOT TRIGGERED**: ladder not run | §22 |
| W6 / A7: percentile vs basic divergence > 25 % of width | **BREACH, all 5 runs** | §23 |
| W7 / A8: \|z₀\| > 0.25 | **BREACH, all 5 runs** (z₀ -0.577 to -0.429) | §18, §23 |
| W5: MC-1 final-40 % trace spread > EP-1 bound | **BREACH, run 1, upper endpoint** | §23 |
| W1 / W2 / W3 / W4 | no breach | §23 |

**Main numerical result (provisional).**
- The point estimate is θ̂ = **+0.0043**, carried forward from Phase 9 primary.
- The five per-run 95 % **percentile** intervals (primary) have lower endpoints -0.0098 to -0.0089 and upper endpoints +0.0469 to +0.0516. The between-run means are [-0.0092, +0.0499].
- **Every** run's percentile, basic and bias-corrected (z₀) interval **includes zero.**

**What the diagnostics say.**
- The outer distribution of ΔC_corr(k) is **right-skewed** and **centred well above θ̂**. Its pooled mean is +0.0140, a shift of +0.0097 that is larger than θ̂ itself.
- That shift is why A7 and A8 are breached.
- The plan (§7) says a material percentile-vs-basic divergence "is evidence of skewness or bootstrap bias and a warning that **neither** interval should be trusted".

**What can be concluded.**
- The analysis provides **no evidence** that the M2 four-pathway model improves optimism-corrected discrimination over M1.
- Nor does it establish the *size* of any difference with interval methods whose assumptions are met.

**Completion status.** The project is **computationally complete** but **not methodologically or scientifically complete**. Four decisions remain open (§29):
1. adjudication of the CV-seed deviation;
2. which MC-SE estimator governs EP-1, and with it the B_out top-up;
3. how to report an interval that the plan's A7 / A8 diagnostics flag;
4. whether θ̂ should be re-derived under the nested driver's functional.

## 3. Scientific background

Phase 9 compared a clinical Cox model (M1: age as a restricted cubic spline, sex, TNM stage, tumour location and adjuvant chemotherapy) with the same model plus four frozen pathway scores (M2 four-set; ridge-penalised molecular block, unpenalised clinical block, P6B-A1). The endpoint was relapse-free survival (`rfs_time_full`, `rfs_event_full`) in the A-19 population (n = 452, 136 events), using the median expression representation.

Phase 9 primary estimated:
- ΔC_app(D) = +0.0140;
- ΔOpt(D) = +0.0096 (Harrell bootstrap optimism, 5 seeds × B = 500);
- hence ΔC_corrected = **+0.0043**.

**No interval was available for ΔC_corrected** (Phase 9 execution report). The completed runs did not retain the replicate-level arrays, so a nested bootstrap was the only valid way to obtain one (plan §1; `docs/PHASE9_DELTAC_INTERVAL_AUDIT.md`).

## 4. Research question

How precisely is the optimism-corrected ΔC = C_corr(M2 four-set) − C_corr(M1) determined *in a cohort with this event composition*? And does its 95 % interval support a non-zero difference?

## 5. Objectives

1. Compute the plan v2.1 percentile interval (primary) with basic and z₀ diagnostics, per run.
2. Quantify the endpoints' Monte-Carlo error (kernel quantile SE at three bandwidths; between-run SD) separately from sampling uncertainty.
3. Evaluate EP-1, H1–H6 and W1–W7 exactly as written.
4. Verify reproducibility (H5) and document every deviation from the written protocol.

## 6. Computational design

- **Runs:** five independent full-size runs, `run_id` 0–4, B_out = 500, B_in = 100, 8 fork-based workers (plan §5).
- **Fits:** 500 000 in total.
- **Outer level:** each outer replicate is an event-stratified patient bootstrap D*ₖ of D.
- **Apparent fits:** M1 and M2 are fitted on D*ₖ with the full pipeline refitted, giving ΔC_app(k).
- **Inner level:** 100 event-stratified bootstraps of D*ₖ, with optimism evaluated on D*ₖ.
- **Inner failures:** pairing-preserving exclusion.
- **Per replicate:** ΔOpt(k) = mean_b opt(M2) − mean_b opt(M1), and ΔC_corr(k) = ΔC_app(k) − ΔOpt(k) (plan §3.1).

## 7. Data and inputs

| Input | Role | sha256/16 |
|---|---|---|
| `data/processed/analysis/expression_median_S0.parquet` | expression (median representation) | `a1ca41807c534584` |
| `data/processed/analysis/covariates.csv` | clinical covariates | `47b2d9829a4fccec` |
| `data/processed/analysis/endpoint.csv` | RFS endpoint | `802a68575a34a3b6` |
| `data/processed/analysis/cohort_membership.csv` | S3 / A-19 membership | `bd3671f8df64c594` |
| `config/pathways/m2_panel_4set.json` | frozen M2 panel | `43ce15b2f0e06f75` |
| `reports/phase9/execution/phase9_execution_results.json` | source of θ̂ | `05d13b8b339fa7d4` |
| `reports/phase9/dc_interval/run{0..4}_outer.jsonl` | the five checkpoints (analysis input) | §12 |

The cohort was originally published as GEO GSE39582 (Marisa et al., 2013).

## 8. Random-seed architecture

- **Base seed:** `BASE = 2026110300`, disjoint from all earlier seeds (plan §8.1).
- **Outer and inner streams:** `SeedSequence([BASE, run_id]).spawn(k+1)[k]` for outer k (equivalent to the plan's `spawn(B_out)[k]`, verified), and `.spawn(B_in)[b]` for inner b.
- **What that guarantees:** runs are independent; within a run, every replicate is a pure function of (BASE, run_id, k, b); resume and scheduling invariance hold.
- **Deviation:** the plan specifies `cvseed = inner[b].spawn(2)[0]` for the λ-tuning CV. The driver instead gives the CV a fresh generator on the *same* seed as the resample (§14).

## 9. Nested-bootstrap methodology

As in §6. Every data-estimated quantity is refitted in every fit, and nothing is cached across layers (plan §3.4):
- age-spline knots;
- gene-level μ/σ (ddof = 1);
- pathway-score standardisation (M2-A2);
- ridge λ by 5-fold CV;
- Cox fit and convergence checks.

The C-index is Harrell's C on `−lp` (lifelines).

## 10. Original five-run execution

| Run | Status | Notes |
|---|---|---|
| 0 | complete | halted at 144 / 500 under the runtime stop condition. Resumed under a **user-authorised exception to the runtime ceiling only** (`RESUME_AUTHORIZATION.json`). Paused in place while the Mac slept (≈ 15:47–18:15, 3 Nov), with no process signalled. Completed 18:31 |
| 1–4 | complete | launched automatically by the original orchestrator; finished 22:55 (run 1), 02:30 (run 2), ≈ 05:50 (run 3), 09:38 (run 4) |

The orchestrator log ends `=== ALL FIVE RUNS COMPLETE ===`. Every run printed `COMPLETE` only after its closing frozen-hash check (H4) passed. Runtimes in the log include sleep time.

## 11. H1 reliability results (and H6, §6.2–6.5 accounting)

| run | records | unique k | k = 0–499 | status ok | paired-complete min | inner failures | events / censored | sha256/16 |
|---|---|---|---|---|---|---|---|---|
| 0 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `885bbd1bd891b570` |
| 1 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `1bb27c56612a5f9d` |
| 2 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `0f12dce12bfce95e` |
| 3 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `4115bd4cd014d345` |
| 4 | 500 | 500 | ✓ | 500 | 100 | 0 | 136 / 316 | `2202d559b0418392` |

- **Failures:** there are 0 inner failures in 250 000 inner replicates, so W1 and W2 do not apply.
- **§6.3:** the representativeness test is **vacuous**, because there are no failed replicates to compare.
- **§6.5:** no outer replicate meets any unreliability criterion, so **H1 is 0 / 2 500**.
- **H6:** paired-complete counts are equal for M1 and M2 in all 2 500 records.

## 12. H4 frozen-artifact verification

The driver's own `check_frozen()` returned **PASS**. The 17-file set gave:

| file | sha256/16 | status |
|---|---|---|
| `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` | `ab56708919f3746f` | OK |
| `config/pathways/m2_panel_4set.json` | `43ce15b2f0e06f75` | OK |
| `config/pathways/M2_FREEZE_RECORD.json` | `64a86fa8ed0d6a8a` | OK |
| `config/analysis/A19_FREEZE_RECORD.json` | `a77de7b37fe09eeb` | OK |
| `src/pathway_scoring.py` | `38cac035c0b2d450` | OK |
| `src/score_standardisation.py` | `41b5c091125e1b11` | OK |
| `src/phase9_execution.py` | `571d0c233f367be6` | OK |
| `data/processed/analysis/expression_median_S0.parquet` | `a1ca41807c534584` | OK |
| `data/processed/analysis/covariates.csv` | `47b2d9829a4fccec` | OK |
| `data/processed/analysis/endpoint.csv` | `802a68575a34a3b6` | OK |
| `reports/phase9/execution/phase9_execution_results.json` | `05d13b8b339fa7d4` | OK |
| `reports/phase9/execution_s7/phase9_s7_results.json` | `b9202305bb6416e1` | OK |
| `data/processed/analysis/cohort_membership.csv` | `bd3671f8df64c594` | OK |
| `scripts/40_phase9_dc_interval.py` | `4b1738a4d0ee0e39` | OK |
| `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` | `d27c52a29930a338` | OK |
| `reports/phase9/dc_interval/HALT_RECORD.json` | `c1b404309bd63744` | OK |
| `reports/phase9/dc_interval/RESUME_AUTHORIZATION.json` | `0571c10e02e9059d` | OK |

The verification was run at the start of this work (`results/s0_starting_state.json`) and again in final QC (`audit/final_qc.json`).

## 13. H5 reproducibility validation

The H5 spot-check recomputed run 0 at k = 0, 25, …, 475 with the **frozen driver**. It loaded `outer_replicate` through importlib without calling `main()`, under the same configuration as the original run (`reports/phase9/h5_spotcheck/`).

**Result: 20 / 20 bit-identical ΔC_corr (`float.hex`), 20 / 20 identical full records.** These 20 indices span all three execution segments: before the halt, after the first resume, and after the sleep pause.

| k | original ΔC_corr (float.hex) | recomputed | bit-identical | multiplicity | full record |
|---|---|---|---|---|---|
| 0 | `0x1.2ce4ed45c7998p-5` | `0x1.2ce4ed45c7998p-5` | ✓ | ✓ | ✓ |
| 25 | `0x1.4639dcc6610acp-5` | `0x1.4639dcc6610acp-5` | ✓ | ✓ | ✓ |
| 50 | `0x1.579716e3a5209p-5` | `0x1.579716e3a5209p-5` | ✓ | ✓ | ✓ |
| 75 | `0x1.52771508a92bbp-6` | `0x1.52771508a92bbp-6` | ✓ | ✓ | ✓ |
| 100 | `0x1.2b003f81a5184p-5` | `0x1.2b003f81a5184p-5` | ✓ | ✓ | ✓ |
| 125 | `0x1.41bafcbdc4f14p-8` | `0x1.41bafcbdc4f14p-8` | ✓ | ✓ | ✓ |
| 150 | `0x1.54c2ad7bb2e66p-7` | `0x1.54c2ad7bb2e66p-7` | ✓ | ✓ | ✓ |
| 175 | `0x1.2383e8297f7fcp-8` | `0x1.2383e8297f7fcp-8` | ✓ | ✓ | ✓ |
| 200 | `0x1.33d316327c1b0p-7` | `0x1.33d316327c1b0p-7` | ✓ | ✓ | ✓ |
| 225 | `0x1.0b4663235fc70p-4` | `0x1.0b4663235fc70p-4` | ✓ | ✓ | ✓ |
| 250 | `0x1.4eb9398647a34p-8` | `0x1.4eb9398647a34p-8` | ✓ | ✓ | ✓ |
| 275 | `0x1.0c7260bd2a407p-5` | `0x1.0c7260bd2a407p-5` | ✓ | ✓ | ✓ |
| 300 | `0x1.10c21e14c935cp-8` | `0x1.10c21e14c935cp-8` | ✓ | ✓ | ✓ |
| 325 | `0x1.dc41e317010b4p-7` | `0x1.dc41e317010b4p-7` | ✓ | ✓ | ✓ |
| 350 | `0x1.b72378b65435ep-6` | `0x1.b72378b65435ep-6` | ✓ | ✓ | ✓ |
| 375 | `0x1.7c86b55ca9a16p-8` | `0x1.7c86b55ca9a16p-8` | ✓ | ✓ | ✓ |
| 400 | `0x1.3bccfb486a11cp-7` | `0x1.3bccfb486a11cp-7` | ✓ | ✓ | ✓ |
| 425 | `-0x1.b6fdfe4ec3300p-10` | `-0x1.b6fdfe4ec3300p-10` | ✓ | ✓ | ✓ |
| 450 | `0x1.053a483415426p-6` | `0x1.053a483415426p-6` | ✓ | ✓ | ✓ |
| 475 | `0x1.f29c4a07e911fp-6` | `0x1.f29c4a07e911fp-6` | ✓ | ✓ | ✓ |

**What H5 establishes:** the as-run implementation is deterministic and resumable.
**What H5 does not establish:** conformance with plan §8.2, or correctness.

## 14. CV-seed protocol deviation

Full record: `documentation/CV_SEED_DEVIATION_DRAFT_v2.md` (**DRAFT — NOT YET ADJUDICATED**). Version 1 is retained at `reports/phase9/h5_spotcheck/`.

- **Plan.** §8.2 line 364 gives the λ-tuning CV its own child stream, `inner[b].spawn(2)[0]`.
- **Driver.** Lines 117 and 142 pass `default_rng(ss_outer)` (apparent fit) and `default_rng(ss_inner[b])` (inner fit). Each starts in **the same PCG64 state** as the generator that drew that fit's bootstrap resample: 9 / 9 cases in `validation/cv_seed_state_evidence.json`, with identical raw outputs. A plan-conformant seed never matches.
- **Affected:** M2 λ selection, and through it every ΔC_corr-derived quantity.
- **Not affected (stored evidence):** M1 apparent λ = 0 in all 2 500 records.
- **Not affected (code path, not measured):** M1 inner λ = 0.
- **θ̂ uses a third convention.** It comes from Phase 9 primary, where the CV draws from one advancing generator after the resample.
- **Not attributed:** the size and direction of any effect are unknown. This report does not attribute any observed feature, including the bias shift in §24, to the deviation.

## 15. Statistical methods

Conventions not fixed by the plan are listed in `documentation/STATISTICAL_CONVENTIONS.md`.

| Quantity | Definition used | Plan reference |
|---|---|---|
| Percentile interval (primary) | [q₀.₀₂₅, q₀.₉₇₅] of the 500 ΔC_corr(k), type-7 quantiles | §3.1, §7 |
| Basic interval (diagnostic) | [2θ̂ − q₀.₉₇₅, 2θ̂ − q₀.₀₂₅] | §7 |
| z₀ (diagnostic) | Φ⁻¹(#{ΔC_corr(k) < θ̂} / 500), strict "<", ties counted separately (there were 0) | §7 |
| BC z₀-only interval (diagnostic) | percentiles at Φ(2z₀ ± 1.96), acceleration 0; BCa excluded | §7; audit line 98 |
| SE(z₀) | √(p(1−p)/B) / φ(z₀) | §7.1 |
| Kernel endpoint MC-SE | √(p(1−p)/B_out) / f̂(q_p), f̂ a Gaussian KDE at the endpoint | §2.2 |
| Bandwidth | Silverman's rule of thumb, 0.9·min(s, IQR/1.34)·n^(−1/5), × {0.5, 1, 2} | §2.2 |
| Between-run endpoint variability | sample SD (ddof = 1, df = 4) of the five endpoints | §2.1 #3, §2.3 |
| EP-1 | MC-SE(endpoint) ≤ 0.05 × (upper − lower), both endpoints; evaluated for **every** prescribed MC-SE estimator | §2 |
| Bounding | missing paired optimism differences imputed at the within-replicate min / max | §6.4 |
| Inner-noise inflation (W3) | mean_k var(d_b)/n_b relative to Var(observed) − that | §1.2 |
| MC-1 trace (W5) | endpoints over the first B ∈ {300, 400, 500} replicates, k order | §11.2; Phase 6B MC-1 |

**Confidence level.** The plan never states the level explicitly. 95 % is used because:
1. Phase 9 primary reports "bootstrap 95 % pct" intervals;
2. the plan's predecessor audit uses the 2.5–97.5 percentile spread;
3. the plan's own arithmetic in §2.2 (expected MC-SE 0.00073 ≈ 3.0 % of width, with a proxy outer SD of 0.0061) implies a width of about 4 SD, i.e. ≈ 2 × 1.96 SD.

## 16. Percentile interval results (primary, PROVISIONAL)

| run | median | mean | SD | % < 0 | percentile 95 % (primary) | basic 95 % (diag.) | BC z₀-only 95 % (diag.) | z₀ | #{< θ̂}/500 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | +0.01086 | +0.01400 | 0.01634 | 21.2 | [-0.00886, +0.05163] | [-0.04300, +0.01749] | [-0.01115, +0.03028] | -0.468 | 160 |
| 1 | +0.01137 | +0.01335 | 0.01572 | 20.0 | [-0.00906, +0.04687] | [-0.03823, +0.01769] | [-0.01146, +0.02728] | -0.519 | 151 |
| 2 | +0.01292 | +0.01421 | 0.01561 | 19.2 | [-0.00927, +0.04838] | [-0.03975, +0.01791] | [-0.01143, +0.02587] | -0.571 | 142 |
| 3 | +0.01180 | +0.01363 | 0.01641 | 23.2 | [-0.00975, +0.05133] | [-0.04269, +0.01839] | [-0.01252, +0.03270] | -0.429 | 167 |
| 4 | +0.01208 | +0.01490 | 0.01615 | 17.8 | [-0.00921, +0.05151] | [-0.04288, +0.01784] | [-0.01197, +0.02788] | -0.577 | 141 |

## 17. Basic interval results (diagnostic, PROVISIONAL)

The basic interval is the reflection of the percentile interval about θ̂ = +0.00432. Because the outer distribution sits mostly above θ̂, the reflection moves it far to the left:
- lower endpoints -0.0430 to -0.0382;
- upper endpoints +0.0175 to +0.0184.

At each endpoint it differs from the percentile interval by 52.2 to 56.4 % of the percentile width (the W6 / A7 breach). The two endpoint differences are equal by construction: both equal |q₀.₀₂₅ + q₀.₉₇₅ − 2θ̂|.

## 18. z₀ interval results (diagnostic, PROVISIONAL)

- **z₀ per run:** -0.577 to -0.429. Only 28.2 to 33.4 % of outer replicates fall below θ̂.
- **Precision:** SE(z₀) ≈ 0.058 per run (plan §7.1). The pooled z₀ over 2 500 replicates is -0.512, with SE 0.026.
- **Bias-corrected (z₀-only) intervals:**
  - percentile levels α₁ ≈ 0.0019, α₂ ≈ 0.847 (run 0);
  - lower endpoints -0.0125 to -0.0112;
  - upper endpoints +0.0259 to +0.0327.
- **Boundary handling:** none needed, since no run had p = 0 or 1.
- **θ̂'s own Monte-Carlo error cannot explain z₀.** θ̂ has MC-SE 0.00019 (SD across the 5 Phase 9 seeds / √5). Moving θ̂ by ± 2 MC-SE changes z₀ by at most 0.034:

| run | z₀ at θ̂ − 2 MC-SE | z₀ at θ̂ | z₀ at θ̂ + 2 MC-SE |
|---|---|---|---|
| 0 | -0.479 | -0.468 | -0.445 |
| 1 | -0.530 | -0.519 | -0.485 |
| 2 | -0.577 | -0.571 | -0.553 |
| 3 | -0.473 | -0.429 | -0.418 |
| 4 | -0.607 | -0.577 | -0.553 |

## 19. Kernel endpoint MC-SE results (PROVISIONAL)

| run | endpoint | q | Silverman h | SE h×0.5 | SE h×1 | SE h×2 | f̂ h×0.5 | f̂ h×1 | f̂ h×2 | SE (scipy 1.06σ rule, cross-check) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | lower | -0.00886 | 0.00414 | 0.000458 | 0.000559 | 0.000574 | 15.24 | 12.48 | 12.17 | 0.000572 |
| 0 | upper | +0.05163 | 0.00414 | 0.002195 | 0.002365 | 0.002270 | 3.18 | 2.95 | 3.08 | 0.002376 |
| 1 | lower | -0.00906 | 0.00392 | 0.000435 | 0.000541 | 0.000589 | 16.07 | 12.92 | 11.86 | 0.000566 |
| 1 | upper | +0.04687 | 0.00392 | 0.002993 | 0.002663 | 0.001969 | 2.33 | 2.62 | 3.55 | 0.002433 |
| 2 | lower | -0.00927 | 0.00400 | 0.000511 | 0.000582 | 0.000630 | 13.67 | 12.01 | 11.09 | 0.000599 |
| 2 | upper | +0.04838 | 0.00400 | 0.001685 | 0.002214 | 0.002023 | 4.14 | 3.15 | 3.45 | 0.002238 |
| 3 | lower | -0.00975 | 0.00426 | 0.000459 | 0.000522 | 0.000576 | 15.22 | 13.39 | 12.11 | 0.000537 |
| 3 | upper | +0.05133 | 0.00426 | 0.002165 | 0.002505 | 0.002265 | 3.23 | 2.79 | 3.08 | 0.002523 |
| 4 | lower | -0.00921 | 0.00419 | 0.000487 | 0.000603 | 0.000653 | 14.33 | 11.59 | 10.69 | 0.000628 |
| 4 | upper | +0.05151 | 0.00419 | 0.001527 | 0.001853 | 0.001899 | 4.57 | 3.77 | 3.68 | 0.001919 |

- **Relative to width:**
  - the lower endpoint has a high density and a small MC-SE (0.00052 to 0.00060, 0.9 to 1.0 % of width);
  - the upper endpoint lies in the sparse right tail (0.00185 to 0.00266, 3.1 to 4.8 % of width).
- **Bandwidth sensitivity:** across bandwidth × {0.5, 1, 2}, the kernel MC-SE of a single endpoint varies by up to a factor of 1.52. This is the plan's stated reason for reporting three bandwidths.
- **Cross-check:** the scipy "Silverman" variant (1.06 σ n^(−1/5)) is given for comparison only.
- **Agreement with the between-run SD:** kernel and between-run estimates agree within 2× for every run and endpoint, so W4 is not breached. Agreement is *consistency*, not validation (plan §2.3).

## 20. Bounding analysis (§6.4)

- **Status:** **vacuous.** No outer replicate has a missing paired inner replicate (0 failures), so ΔC_corr^min(k) = ΔC_corr^max(k) = ΔC_corr(k) for every k, and both bounded intervals equal the primary interval.
- **Endpoint movement:** 0. It is therefore below every endpoint MC-SE, so **H2 is not breached**.
- **§6.5 criterion 3** (bounding width > MC-SE) is met by no replicate.
- **Testing:** the imputation code is unit-tested on synthetic replicates with missing pairs (`validation/tests/test_stats_core.py`).

## 21. Across-run variability

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

- **Precision of the SD estimates:** each SD rests on 5 values (df = 4), so its relative SE is about 35 % (plan §2.3).
- **Spread:** the upper percentile endpoint varies about 6.5× more between runs than the lower endpoint.
- **Scope:** these SDs bound Monte-Carlo variability only. They do not bound any systematic error that all five runs share.

## 22. EP-1 evaluation

**Criterion (plan §2, approved decision 1):** MC-SE(endpoint) ≤ 0.05 × (upper − lower), for **both** endpoints, on the percentile interval of each run. The plan prescribes two MC-SE estimators: kernel at Silverman × {0.5, 1, 2} (§2.2) and the between-run SD (§2.1 #3). **It does not say which governs**, so all four are evaluated.

| MC-SE estimator | evaluations | pass (≤ 0.05) | max MC-SE / width |
|---|---|---|---|
| kernel_bw_x0.5 | 10 | 9 | 0.0535 |
| kernel_bw_x1.0 | 10 | 10 | 0.0476 |
| kernel_bw_x2.0 | 10 | 10 | 0.0375 |
| between_run_SD | 10 | 10 | 0.0391 |

The only failing evaluation is **run 1, upper endpoint, kernel at bandwidth × 0.5**. Its MC-SE / width is 0.0535, against a bound of 0.05. At the Silverman bandwidth itself the same endpoint gives 0.0476.

**EP-1 status: UNRESOLVED, decision required.**
- **What it shows:** under the plan's named bandwidth and under the between-run SD, every endpoint of every run meets EP-1.
- **If every prescribed estimator must pass:** EP-1 fails narrowly. Plan §2.4 would then require a top-up to B_out′ = 500 × (0.0535 / 0.05)² = 573, i.e. ladder rung **750**.
- **Why the top-up was not run:** it is a resource and scientific decision, and it would add records to the run checkpoints. **The threshold was not relaxed.**
- **Overlap with W5:** run 1's upper-endpoint W5 breach (§23) points the same way. The plan's consequence for W5 is "raise B_out".
- **What EP-1 does not establish:** EP-1 concerns Monte-Carlo precision only. Even a clean pass would say nothing about the validity of the interval, which is the subject of the A7 / A8 diagnostics.

The full 40-row grid:

| run | endpoint | estimator | MC-SE | width | MC-SE / width | result |
|---|---|---|---|---|---|---|
| 0 | lower | kernel_bw_x0.5 | 0.000458 | 0.06049 | 0.0076 | PASS |
| 0 | lower | kernel_bw_x1.0 | 0.000559 | 0.06049 | 0.0092 | PASS |
| 0 | lower | kernel_bw_x2.0 | 0.000574 | 0.06049 | 0.0095 | PASS |
| 0 | lower | between_run_SD | 0.000334 | 0.06049 | 0.0055 | PASS |
| 0 | upper | kernel_bw_x0.5 | 0.002195 | 0.06049 | 0.0363 | PASS |
| 0 | upper | kernel_bw_x1.0 | 0.002365 | 0.06049 | 0.0391 | PASS |
| 0 | upper | kernel_bw_x2.0 | 0.002270 | 0.06049 | 0.0375 | PASS |
| 0 | upper | between_run_SD | 0.002186 | 0.06049 | 0.0361 | PASS |
| 1 | lower | kernel_bw_x0.5 | 0.000435 | 0.05592 | 0.0078 | PASS |
| 1 | lower | kernel_bw_x1.0 | 0.000541 | 0.05592 | 0.0097 | PASS |
| 1 | lower | kernel_bw_x2.0 | 0.000589 | 0.05592 | 0.0105 | PASS |
| 1 | lower | between_run_SD | 0.000334 | 0.05592 | 0.0060 | PASS |
| 1 | upper | kernel_bw_x0.5 | 0.002993 | 0.05592 | 0.0535 | **FAIL** |
| 1 | upper | kernel_bw_x1.0 | 0.002663 | 0.05592 | 0.0476 | PASS |
| 1 | upper | kernel_bw_x2.0 | 0.001969 | 0.05592 | 0.0352 | PASS |
| 1 | upper | between_run_SD | 0.002186 | 0.05592 | 0.0391 | PASS |
| 2 | lower | kernel_bw_x0.5 | 0.000511 | 0.05766 | 0.0089 | PASS |
| 2 | lower | kernel_bw_x1.0 | 0.000582 | 0.05766 | 0.0101 | PASS |
| 2 | lower | kernel_bw_x2.0 | 0.000630 | 0.05766 | 0.0109 | PASS |
| 2 | lower | between_run_SD | 0.000334 | 0.05766 | 0.0058 | PASS |
| 2 | upper | kernel_bw_x0.5 | 0.001685 | 0.05766 | 0.0292 | PASS |
| 2 | upper | kernel_bw_x1.0 | 0.002214 | 0.05766 | 0.0384 | PASS |
| 2 | upper | kernel_bw_x2.0 | 0.002023 | 0.05766 | 0.0351 | PASS |
| 2 | upper | between_run_SD | 0.002186 | 0.05766 | 0.0379 | PASS |
| 3 | lower | kernel_bw_x0.5 | 0.000459 | 0.06108 | 0.0075 | PASS |
| 3 | lower | kernel_bw_x1.0 | 0.000522 | 0.06108 | 0.0085 | PASS |
| 3 | lower | kernel_bw_x2.0 | 0.000576 | 0.06108 | 0.0094 | PASS |
| 3 | lower | between_run_SD | 0.000334 | 0.06108 | 0.0055 | PASS |
| 3 | upper | kernel_bw_x0.5 | 0.002165 | 0.06108 | 0.0354 | PASS |
| 3 | upper | kernel_bw_x1.0 | 0.002505 | 0.06108 | 0.0410 | PASS |
| 3 | upper | kernel_bw_x2.0 | 0.002265 | 0.06108 | 0.0371 | PASS |
| 3 | upper | between_run_SD | 0.002186 | 0.06108 | 0.0358 | PASS |
| 4 | lower | kernel_bw_x0.5 | 0.000487 | 0.06072 | 0.0080 | PASS |
| 4 | lower | kernel_bw_x1.0 | 0.000603 | 0.06072 | 0.0099 | PASS |
| 4 | lower | kernel_bw_x2.0 | 0.000653 | 0.06072 | 0.0108 | PASS |
| 4 | lower | between_run_SD | 0.000334 | 0.06072 | 0.0055 | PASS |
| 4 | upper | kernel_bw_x0.5 | 0.001527 | 0.06072 | 0.0251 | PASS |
| 4 | upper | kernel_bw_x1.0 | 0.001853 | 0.06072 | 0.0305 | PASS |
| 4 | upper | kernel_bw_x2.0 | 0.001899 | 0.06072 | 0.0313 | PASS |
| 4 | upper | between_run_SD | 0.002186 | 0.06072 | 0.0360 | PASS |

## 23. Robustness findings: hard-failure and warning inventory

| run | inner failures | unreliable outer (6.5) | H1 | bounding (6.4) | H2 | W3 SD inflation | W5 spread lo / hi (bound) | W6 A7 lo / hi | W7 |z₀| |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | no | VACUOUS | no | 0.35 % | 0.00053 / 0.00160 (0.00302) | 56.4 % / 56.4 % **breach** | 0.468 **breach** |
| 1 | 0 | 0 | no | VACUOUS | no | 0.36 % | 0.00026 / 0.00299 (0.00280) **breach** | 52.2 % / 52.2 % **breach** | 0.519 **breach** |
| 2 | 0 | 0 | no | VACUOUS | no | 0.36 % | 0.00089 / 0.00107 (0.00288) | 52.9 % / 52.9 % **breach** | 0.571 **breach** |
| 3 | 0 | 0 | no | VACUOUS | no | 0.34 % | 0.00001 / 0.00133 (0.00305) | 53.9 % / 53.9 % **breach** | 0.429 **breach** |
| 4 | 0 | 0 | no | VACUOUS | no | 0.34 % | 0.00040 / 0.00291 (0.00304) | 55.5 % / 55.5 % **breach** | 0.577 **breach** |

| Plan item | Required analysis | Status |
|---|---|---|
| §2.2 bandwidth × {0.5, 1, 2} | kernel MC-SE sensitivity | done (§19) |
| §2.1 #3, §2.3 | between-run endpoint SD | done (§21) |
| §6.3 | failed-vs-successful representativeness | **NOT APPLICABLE**: 0 failures |
| §6.4 / H2 | bounding analysis | done; vacuous (§20) |
| §6.5 / H1 | unreliable outer replicates | done; 0 (§11) |
| §1.2 / W3 | inner-noise SD inflation | done; 0.34 to 0.36 % (plan's sizing anticipated ≈ 1.2 %) |
| W4 | kernel vs between-run, > 2× | done; no breach |
| W5 / MC-1 | running-trace spread, final 40 % | done; **breach in run 1 (upper)** |
| W6 / A7 | percentile vs basic | done; **breach in all runs** |
| W7 / A8 | \|z₀\| > 0.25 | done; **breach in all runs** |
| §2.4 | B_out ladder top-up if EP-1 fails | **BLOCKED**: EP-1 governing estimator undecided; top-up not authorised |
| S-7 representation | nested interval for S-7 | **NOT APPLICABLE**: plan v2.1 specifies the median representation only; S-7 is sensitivity-only for the point estimate |
| CV-seed sensitivity recomputation | alternative seeding | **NOT AUTHORISED**: listed for adjudication only |

## 24. Scientific interpretation

### 24.1 Observed results

- **The estimand.** ΔC_corr is a difference in Harrell's C between two optimism-corrected models. A positive value means M2 ranks relapse risk better than M1 after optimism correction.
- **The outer distribution:**
  - per-run medians +0.0109 to +0.0129;
  - means +0.0134 to +0.0149;
  - SD 0.0156 to 0.0164;
  - right skew 0.64 to 1.02;
  - 17.8 to 23.2 % of replicates below zero.
- **All 15 intervals include zero:** 5 runs × 3 methods.
- **Decomposition:**

  | quantity | original data D (Phase 9 primary) | mean over 2 500 outer replicates | shift |
|---|---|---|---|
| ΔC_app | +0.01396 | +0.02224 | +0.00828 |
| ΔOpt | +0.00965 | +0.00822 | -0.00142 |
| ΔC_corr = ΔC_app − ΔOpt | +0.00432 | +0.01402 | +0.00970 |

  The outer mean of ΔC_corr(k) exceeds θ̂ by +0.0097. Of that:
  - +0.0083 comes from a larger apparent ΔC on D*ₖ;
  - +0.0014 comes from a smaller ΔOpt within D*ₖ.

| run | ΔC_app(k) mean | ΔOpt(k) mean | opt M1 mean | opt M2 mean | ΔC_corr skew | ΔC_corr excess kurtosis | distinct patients (mean) | M2 λ median | M2 λ max |
|---|---|---|---|---|---|---|---|---|---|
| 0 | +0.02233 | +0.00833 | +0.02117 | +0.02950 | 0.89 | 0.75 | 286.0 | 0.152 | 10 |
| 1 | +0.02167 | +0.00831 | +0.02067 | +0.02898 | 1.02 | 2.00 | 286.1 | 0.152 | 10 |
| 2 | +0.02224 | +0.00802 | +0.02147 | +0.02950 | 0.75 | 0.88 | 286.3 | 0.152 | 10 |
| 3 | +0.02183 | +0.00820 | +0.02097 | +0.02917 | 0.69 | 0.24 | 286.4 | 0.152 | 10 |
| 4 | +0.02314 | +0.00824 | +0.02103 | +0.02927 | 0.64 | 0.01 | 285.9 | 0.152 | 10 |

### 24.2 Supported conclusions

1. At the reported Monte-Carlo precision, **the 95 % interval for ΔC_corr is compatible with zero under every interval method and every run.** The data do not support a claim that the four pathway scores improve optimism-corrected discrimination over the clinical model. That conclusion does not depend on which interval method is chosen, or on the A7 / A8 outcome.
2. **The bootstrap distribution is not centred on θ̂.** |z₀| = 0.43 to 0.58 is 7.4 to 9.7 SE(z₀) from 0, and θ̂'s own Monte-Carlo error can't account for it. The plan's diagnostic for this (A8) fires in every run.
3. **Method choice changes the *location and upper limit* of the interval substantially:**
   - percentile upper ≈ +0.050;
   - basic upper ≈ +0.018;
   - BC upper ≈ +0.029.

   Under plan §7, a divergence this large is "a warning that neither interval should be trusted". **No method is preferred here because of its result.**
4. **Endpoint Monte-Carlo error is small next to sampling uncertainty.** at the Silverman bandwidth, MC-SE is 0.9–4.8 % of the interval width, and inner-noise inflation of the outer SD is 0.34 to 0.36 %. Numerical precision is not what limits the conclusions.

### 24.3 Plausible explanations (inferential, not tested here)

- **A larger apparent ΔC on bootstrap samples.** On D*ₖ, 36.7 % of rows duplicate other rows, which plausibly favours the more flexible M2 (four penalised scores) in apparent C.
- **A smaller inner ΔOpt.** The inner optimism is evaluated on D*ₖ, which shares most of its patients with each inner sample. That could make the estimated optimism *smaller* than on D.
- **Plan §7.2 predicted the opposite.** It expected the outer distribution to sit *below* θ̂, because inner optimism would be larger. The observed direction is the reverse. Note also that §7.2 says a distribution below θ̂ would give a *negative* z₀, which is inconsistent with §7's own definition (a distribution below θ̂ gives z₀ > 0). The observed z₀ is negative because the distribution sits *above* θ̂.
- **θ̂ and the outer replicates come from different implementations.** θ̂ is the Phase 9 primary estimate (different seeding, §14), not the frozen nested driver applied to D. How comparable they are is unverified.
- **The CV-seed deviation is not offered as an explanation.** No evidence links it to the shift.

### 24.4 Unresolved questions

1. Adjudication of the CV-seed deviation (record v2, §10).
2. Which MC-SE estimator governs EP-1, and so whether the §2.4 top-up (to B_out = 750) is required. Also how to act on the run-1 W5 breach.
3. How to report an interval when A7 and A8 both fire in every run: which interval, if any, is presented as *the* result, and with what qualification. The plan doesn't designate a single headline interval across the five runs.
4. Whether θ̂ from Phase 9 primary is the right centring value for the basic interval and z₀, or whether the frozen driver's functional should be evaluated on D. That evaluation would be a new computation and is not authorised.

## 25. Limitations

- **This analysis:**
  - the CV-seed protocol deviation (unadjudicated);
  - θ̂ produced under a different seeding convention;
  - EP-1's governing estimator unspecified;
  - the run-1 W5 breach;
  - the A7 / A8 breaches;
  - conventions adopted where the plan is silent (95 % level, Silverman 1986 rule, type-7 quantiles).
- **Standing limitations (plan §13–14), unchanged:**
  - conditional on 136 events;
  - bootstrap validity asymptotic at n = 452;
  - distinct-patient dilution (≈ 286 of 452 distinct per outer sample);
  - between-run SD on df = 4;
  - A7 has no statistical justification and A8 is a loose detector;
  - no deconvolution of inner noise;
  - Harrell's C depends on the censoring distribution;
  - asymmetric penalisation (P6B-A1) makes the comparison conservative by design;
  - internal validation only; U1 unresolved.
- **What the interval can't establish:** it quantifies uncertainty and nothing more. It establishes **no clinical utility, biological mechanism, external validity or transportability.**

## 26. Reproducibility instructions

See `reproducibility/README_REPRODUCE.md`. In summary: `reproducibility/run_all.sh` runs stages S0 → S3, the independent check, the tests and final QC, writing only inside `reports/phase9/final_project/`. Outputs are deterministic given the five checkpoints. No stage re-runs the bootstrap.

## 27. Software and environment

| | |
|---|---|
| Python | 3.13.13 |
| Platform | macOS-27.0.1-arm64-arm-64bit-Mach-O |
| numpy / scipy / pandas / lifelines / matplotlib | 2.5.4 / 1.18.1 / 2.3.3 / 0.30.3 / ? |

The full capture is in `reproducibility/environment.json`. No dependency was installed or upgraded.

## 28. Audit trail

The project-level record is `audit/EXECUTION_LEDGER.md` (every task, command, test, error and disposition), with run logs in `audit/`. The earlier records are retained unchanged: `HALT_RECORD.json`, `RESUME_AUTHORIZATION.json`, `resume_driver.log`, `reports/phase9/h5_spotcheck/*`.

Two corrections were made during this work and are preserved in `audit/`:
1. the H3 flag was wrongly set when EP-1 was mixed;
2. figure label collisions.

## 29. Conclusions

1. **Computation is complete and verified.** 5 × 500 × 100 replicates; H1, H2, H4, H5 and H6 pass; 0 inner failures; 166 / 166 independent re-computation checks agree.
2. **Provisional statistical result.** Every 95 % interval for the optimism-corrected ΔC includes zero. The primary percentile interval averages [-0.0092, +0.0499] across runs. The analysis does **not** support an improvement of M2 over M1.
3. **Interval validity is in question by the plan's own diagnostics.** A7 and A8 breach in every run: the bootstrap distribution sits well above θ̂ and is right-skewed. No interval should be quoted without that qualification.
4. **Open decisions:**
   - (a) CV-seed deviation;
   - (b) EP-1 governing estimator, and with it the §2.4 / W5 B_out top-up;
   - (c) how to report an interval that the plan's diagnostics flag;
   - (d) whether θ̂ should be re-derived under the nested driver's functional.

   Until these are settled, the project is **computationally complete but methodologically and scientifically provisional.**

## 30. References

1. Efron B, Tibshirani RJ. *An Introduction to the Bootstrap.* Chapman & Hall, 1993.
2. Davison AC, Hinkley DV. *Bootstrap Methods and their Application.* Cambridge University Press, 1997.
3. Efron B. Better bootstrap confidence intervals. *J Am Stat Assoc* 1987;82:171–185.
4. Harrell FE, Lee KL, Mark DB. Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. *Stat Med* 1996;15:361–387.
5. Steyerberg EW. *Clinical Prediction Models*, 2nd ed. Springer, 2019.
6. Silverman BW. *Density Estimation for Statistics and Data Analysis.* Chapman & Hall, 1986.
7. Hyndman RJ, Fan Y. Sample quantiles in statistical packages. *Am Stat* 1996;50:361–365.
8. van der Vaart AW. *Asymptotic Statistics.* Cambridge University Press, 1998 (asymptotic variance of sample quantiles).
9. Davidson-Pilon C. lifelines: survival analysis in Python. *J Open Source Softw* 2019;4(40):1317.
10. Marisa L, et al. Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. *PLoS Med* 2013;10(5):e1001453.
11. Project documents: `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1; `docs/PHASE9_DELTAC_INTERVAL_AUDIT.md`; `docs/PHASE9_EXECUTION_REPORT.md`; `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` (MC-1…MC-4).

---

### Figures

![F1](../figures/F1_outer_distributions.png)
![F2](../figures/F2_interval_comparison.png)
![F3](../figures/F3_endpoint_variability.png)
![F4](../figures/F4_MC1_running_trace.png)
![F5](../figures/F5_EP1_grid.png)
![F6](../figures/F6_shift_decomposition.png)

Captions: `figures/FIGURE_CAPTIONS.md`. Every figure is generated by `analysis/s2_figures.py` from `results/s1_interval_results.json`.
