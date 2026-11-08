# Phase 9 — Execution Report (v1.0)

**Protocol:** v3.1 (`ab56708919f3746f`) · **Executed:** 2026-11-03T03:50:56Z · **Runtime:** 4494.5 s
**Script:** `scripts/38_phase9_execute.py` · **Engine:** `src/phase9_execution.py`
**Results:** `reports/phase9/execution/phase9_execution_results.json`
**Tests:** `tests/test_phase9_execution.py` · **Smoke archive:** `reports/phase9/execution/phase9_smoketest_B5_ARCHIVE.json`

> **Three authorised analyses only.** M3 was **not** fitted and no transcriptome-wide selection
> was performed. No fixed-time, absolute-risk, CIF, time-dependent AUC, Brier or fixed-horizon
> cloglog quantity was computed (U1).

---

## 1. Execution status

PID 42843, `B=500`, **exit code 0**, completed 2026-11-03T03:50:56Z after 4494.5 s. No traceback,
no truncation; JSON valid with all sections present. Recorded configuration matches the run:
**B = 500**, 5 seeds, `B_effective = 500` in every seed for every model (**7 500 fits total**),
models `{M1, M2_primary, M2_sensitivity}`, `M3_fitted: false`.

**The B=5 smoke test was overwritten** at the same path by the B=500 run. Its values are archived
separately and labelled *"SMOKE TEST ONLY — NOT AN AUTHORISED RESULT"*; at B=5 MC-3 **failed** for
all three models (ratios 1.08–1.33), as expected from five replicates.

## 2. Analysis population

| | n | events | censored |
|---|---:|---:|---:|
| S3 cohort | 453 | **136** | 317 |
| **A-19 complete-case (model-fitting)** | **452** | **136** | 316 |

Sole exclusion **GSM972262** — missing `chemotherapy.adjuvant`, reconciled against the A-19
record. It is **censored**, so no event is lost. S3's 136 events match Register A-01 exactly, and
Phase 9's pre-execution bound of [135, 136] is honoured. Endpoint columns `rfs_event_full` /
`rfs_time_full` (A-07, full follow-up). Time scale and censoring unaltered; U1/U2 unresolved.

**Parameter-count note.** The protocol's planning table counted age as one covariate (M2 p = 9,
EPV 15.1). §4.1 specifies a **2-df age spline**, so the fitted models carry **6 / 10 / 11 df** and
EPV is **22.7 / 13.6 / 12.4**. The df-based figures are used throughout.

## 3. M1 — clinical baseline

6 df, **λ = 0 (no penalised term)**, knots [49.0, 68.5, 82.9], no linear fallback. logLik −777.640.

| covariate | β | SE | HR | 95 % CI |
|---|---|---|---|---|
| age_rcs1 | +0.0048 | 0.0140 | 1.005 | 0.978–1.033 |
| age_rcs2 | +0.0100 | 0.0148 | 1.010 | 0.981–1.040 |
| sex = Male | +0.3212 | 0.1800 | 1.379 | 0.969–1.962 |
| **stage III** | **+0.6761** | 0.2125 | **1.966** | **1.296–2.982** |
| location = proximal | −0.2180 | 0.1854 | 0.804 | 0.559–1.156 |
| adjuvant chemo = Y | +0.1552 | 0.2185 | 1.168 | 0.761–1.792 |

Converged; separation sentinel not tripped (max |β| 0.676, max SE 0.219). **Stage is the only
covariate whose CI excludes 1.** Chemotherapy is a prognostic covariate only — **never** a
treatment effect (A-12).

## 4. M2 primary — frozen four-set panel

10 df, **λ = 0.3511** (molecular block only), same clinical block and knots as M1. logLik −775.271.

Clinical terms essentially unchanged from M1 (stage III HR 1.939, 1.277–2.944). Molecular terms,
**per 1 training-fold SD** (M2-A2):

| pathway | β | SE | HR | 95 % CI |
|---|---|---|---|---|
| EMT | +0.0672 | 0.0597 | 1.070 | 0.951–1.202 |
| E2F_TARGETS | −0.0991 | 0.0588 | 0.906 | 0.807–1.016 |
| WNT_BETA_CATENIN | +0.0154 | 0.0580 | 1.016 | 0.906–1.138 |
| P53_PATHWAY | −0.0047 | 0.0618 | 0.995 | 0.882–1.123 |

**Every molecular CI includes 1.** Coverage rules met for all sets (198/195/41/194 genes), 0 score
drops, 0 zero-variance exclusions.

## 5. M2 sensitivity — frozen five-set panel

11 df, **λ = 0.1520**, logLik −774.232. EMT +0.0899, E2F −0.1485 (HR 0.862, 0.746–**0.996**), WNT
+0.0243, P53 −0.0220, **PI3K +0.0417** (HR 1.043, 0.891–1.221).

**Collinearity, as D8D-2 requires:** P53 ↔ PI3K **r = +0.6360**, the only |r| > 0.5 pair; condition
number **2.909 → 7.084**; minimum eigenvalue 0.4848 → 0.2762; VIF rises for P53 (1.23 → 1.96) and
PI3K is 2.19. Every molecular SE inflates (×1.20–1.34, worst for P53). The foreseeable instability
**materialised**: E2F's coefficient moves −0.099 → −0.149 and its CI brushes 1.

## 6. Discrimination and calibration

| model | C apparent | optimism | **C corrected** | bootstrap 95 % pct |
|---|---|---|---|---|
| M1 | 0.6242 | +0.0226 | **0.6016** | 0.5873–0.6243 |
| M2 primary | 0.6382 | +0.0323 | **0.6059** | 0.5995–0.6393 |
| M2 sensitivity | 0.6413 | +0.0347 | **0.6066** | 0.6004–0.6412 |

Calibration slope (log-hazard scale, the only calibration permitted): M1 1.000 → **0.833**,
M2 primary 1.083 → **0.862**, M2 sensitivity 1.075 → **0.836**. **Slopes well below 1 indicate
real overfitting**; predictions are too extreme and would need shrinking by ~15 %.

Optimism is **larger for the molecular models** (+0.032/+0.035 vs +0.023) — the penalised
molecular block still buys more apparent fit than it retains.

## 7. ΔC — and an important distinction

| quantity | value |
|---|---|
| Apparent ΔC (M2 primary − M1) | +0.0140 |
| **Paired out-of-original ΔC** | **+0.0121**, bootstrap 95 % pct **[+0.0002, +0.0240]**, MC-SE 0.00010 |
| **Optimism-corrected ΔC** | **+0.0043** |

> **These are different quantities and the interval belongs to the middle one.** The paired ΔC is
> computed on the original sample from bootstrap-fitted models; its interval barely excludes zero
> (lower bound **+0.0002**). The corrected ΔC — the difference of the two optimism-corrected
> C-indices — is **+0.0043**, less than half as large, because M2 carries more optimism than M1.
> **No interval was computed for the corrected ΔC**, so the "excludes zero" result must **not** be
> transferred to it. This is a limitation of the executed analysis, not a finding.

## 8. Monte-Carlo stability (MC-1 … MC-4)

**MC-1 — running trace.** Corrected C at B = 50/100/200/300/400/500, seed 20261102:
M1 0.6054→0.6013, M2p 0.6115→0.6061, M2s 0.6132→0.6070. Spread over the final 40 % of replicates:
**0.00033 / 0.00022 / 0.00027** — visibly stable.

**MC-2 — five seeds.** Range of corrected C: M1 **0.00157**, M2 primary **0.00111**,
M2 sensitivity **0.00089**.

**MC-3 — `MC-SE ≤ 0.10 × sampling SE`**, the two quantities reported separately:

| model | MC-SE | sampling SE | threshold (0.10 × SE) | ratio | **MC-3** |
|---|---|---|---|---|---|
| M1 | 0.00058 | 0.00944 | 0.00094 | 0.062 | **PASS** |
| M2 primary | 0.00047 | 0.01015 | 0.00101 | 0.046 | **PASS** |
| M2 sensitivity | 0.00041 | 0.01042 | 0.00104 | 0.039 | **PASS** |

**All three pass at B = 500. No B-increase is required**, so B stays at 500 — it was never
lowered, and no seed, threshold or criterion was changed.

**MC-4 — reporting.** MC-SE accompanies every corrected estimate above and is never summed with
the sampling SE. *Caveat:* the sampling SE is approximated as (bootstrap 97.5th − 2.5th)/3.92 of
the out-of-original C distribution — a normal-approximation proxy, not an exact sampling SE.

## 9. Non-estimable replicates

**0 of 2 500 per model (0.0 %)** — M1, M2 primary and M2 sensitivity alike. Attempted 2 500,
successful 2 500, non-estimable 0. No failure reasons to report; **nothing was redrawn**
(`redrawn: false`). The 5 % threshold is not approached, so no estimate is flagged UNRELIABLE.
Uncertainty was computed on the full replicate set, not a surviving subset.

## 10. Batch-stratified sensitivity (D8D-4-R1)

Nine `combat_batch` strata, **individually, no pooling, no batch coefficients in any primary model**.

| stratum | n | events | risk-set comparisons | events with a comparator |
|---|---:|---:|---:|---:|
| batch01 | 33 | 8 | 195 | 7 |
| batch02 | 14 | 6 | 62 | 6 |
| batch03 | 102 | 24 | 1 923 | 24 |
| batch04 | 103 | 52 | 3 983 | 52 |
| batch05 | 48 | 6 | 217 | 6 |
| batch07 | 81 | 13 | 796 | 13 |
| batch08 | 41 | 16 | 498 | 16 |
| batch09 | 15 | 4 | 50 | 4 |
| batch10 | 15 | 7 | 78 | 7 |
| **total** | **452** | **136** | **7 802** | **135** |

**No zero-event stratum**; all nine contribute ≥ 1 event with a comparator.

**Risk-set comparison retention** = N_strat / N_unstrat = 7 802 / 47 368 = **0.1647**.
Formula: N = Σ over events of (|risk set| − 1), within-stratum for N_strat. **Descriptive only.**
It counts **comparisons, not Fisher information**, and is **not** a measure of precision; a
comparison between two subjects with identical covariates is counted yet contributes nothing.
**No threshold is applied.**

**Stratified fit:** converged; covariance **positive definite** (min eigenvalue 2.80 × 10⁻⁵);
separation sentinel not tripped (max |β| 0.520, max SE 0.242); **0 non-estimable fits**.
Molecular coefficients barely move (Δ ≤ 0.027); **clinical** terms move more — stage III
0.662 → 0.520, chemo 0.137 → 0.029.

> A positive-definite information matrix proves local identification and invertibility. It does
> **not** prove adequate precision, and a numerical failure is not the same thing as a wide
> interval.

**Leave-one-batch-out** (nine refits, 0 failures): apparent C spans **0.6143–0.6592** (range
0.0449). Molecular coefficients are **markedly unstable** — EMT +0.009 to +0.204, E2F −0.249 to
−0.014, and **WNT and P53 each change sign**. Omitting batch03 or batch04 removes ~23 % of the
cohort, so the spread conflates batch identity with sample size (pre-registered caveat).

## 11. Representation sensitivity (D8D-3 / estimand S-7)

**OUTSTANDING.** Only `expression_median_S0.parquet` was used; the representative-IQR matrix was
not referenced. The S-7 re-fit was **not part of this authorisation and was not executed and not
launched**. **Median-versus-representative-IQR robustness is therefore not established.** WNT —
the weakest-agreeing score (r = 0.858) — is also one of the two whose sign flips under
leave-one-batch-out, so this gap matters.

## 12. Leakage and specification audit

Clinical block **identical in all three models** (same six terms, same order); M1 **λ = 0**; in M2
the penalty vector is `[0,0,0,0,0,0,λ,…]` — clinical unpenalised, molecular penalised (test E07).
Knots identical across models at the full-sample fit and **recomputed per training fold** inside
resampling; no linear fallback triggered. Coverage rules met, 0 score drops. Tests prove the
pipeline is training-fold-only: poisoning held-out rows leaves the fitted model **bit-identical**;
held-out scores are **not** forced to unit variance.

## 13. Endpoint restrictions honoured

No five-year or named-horizon absolute risk, no CIF at named years, no time-dependent AUC, no
fixed-horizon Brier, no fixed-horizon cloglog calibration. Only relative-risk and
internal-validation measures permitted by v3.1 were computed. Cause-specific hazards are **not**
interpreted as absolute recurrence probabilities.

## 14. Integrity

**526 tests pass.** Bundle VERIFICATION PASSED; `27_freeze_m2_specification` and
`29_freeze_a19_missing_data` exit 0; Phase 9 final audit **36/36**; raw checksums OK. **14/14
frozen artifacts unchanged.** Amendment chain `M2-A1 → M2-A2` (spec rev 3) resolves at every
revision; P6B-A1 intact with the original Phase 6B wording preserved; D8D-4-R1 intact with the
squared-SE ratio still rejected; D8D-2 and D8D-3 approved. Protocol v3.1 hash `ab56708919f3746f`.
All 11 recorded input hashes still match disk. Python 3.13.13, numpy 2.5.4, pandas 2.3.3,
lifelines 0.30.3. **Zero commits.**

## 15. Assessment

**Does M2 primary improve on M1?** Marginally at best. Corrected C rises **0.6016 → 0.6059**, a
corrected ΔC of **+0.0043**. The paired out-of-original ΔC of +0.0121 has an interval whose lower
bound is +0.0002 — but that interval belongs to a different quantity, and no interval was computed
for the corrected ΔC. **Every individual molecular CI includes 1.**

**Is it stable?** *Computationally*, yes — MC-3 passes with large margin, the running traces are
flat, and 0/7 500 fits failed. *Statistically*, no. Leave-one-batch-out moves molecular
coefficients by an order of magnitude and flips the sign of WNT and P53. Corrected calibration
slopes of 0.83–0.86 show real overfitting. **Computational stability is not statistical stability,
and only the former is established.**

**Does the five-set change the conclusion?** No. Corrected C 0.6066 vs 0.6059 — a difference of
0.0007, **well inside the MC-SE** and declared **"no evidence of a difference"** under D8D-2. The
four-set panel **remains primary**, as pre-specified, regardless. PI3K's addition did what the
collinearity predicted: SEs inflate up to ×1.34 and E2F's estimate shifts by 50 %.

**What is insufficient?** A corrected ΔC of +0.0043 on an internally validated C of ~0.61, with
all molecular CIs crossing 1, unstable coefficients under batch perturbation, and the
representation sensitivity unexecuted. The molecular block is **penalised while the clinical block
is not**, so this comparison is **conservative by design and is not an unbiased estimate of the
molecular contribution** — but a conservative design does not convert a null-consistent result
into a positive one. **No claim of clinical utility or biological causality is made or supported.**
No external validation exists, and none is claimed.

**Ready for final interpretation?** **No.** The three authorised analyses are complete and clean,
but the S-7 representation sensitivity is outstanding, the corrected ΔC lacks an interval, and the
coefficient instability under batch perturbation is unexplained. Each requires separate
authorisation.
