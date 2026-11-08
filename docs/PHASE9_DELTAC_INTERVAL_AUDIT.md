# Phase 9 — Methodological Audit: a Confidence Interval for the Optimism-Corrected ΔC

**Version:** 1.0 · **Date:** 2026-11-03 · **Status:** ⛔ **AUDIT AND RECOMMENDATION ONLY — NO NEW COMPUTATION RUN**

> Nothing was executed beyond a 10-replicate timing probe. No frozen decision, primary output,
> S-7 output or report was modified. No commit.

---

## 1. The estimand and the estimator

**Estimand.** The population value of

> **ΔC_corrected = C_corrected(M2 four-set) − C_corrected(M1)**

where `C_corrected = C_apparent − E[optimism]`, optimism being Harrell's bootstrap bias term.

**Estimator, exactly as implemented** (`scripts/38_phase9_execute.py`):

```
ΔC_corrected = ΔC_apparent − ( mean_b optimism_M2(b) − mean_b optimism_M1(b) )
             = +0.0140     − ( +0.0323            − +0.0226            )
             = +0.0140     − +0.0096
             = +0.0043
```

It is a **functional of the original 452-patient sample**. Its sampling distribution therefore
requires resampling **that sample** — a layer the executed design never ran.

**Resampling units and pairing as implemented.** Patient-level, event-stratified bootstrap
(`boot_indices` fixes 136 events in every replicate). Within replicate *b*, M1 and M2 are fitted on
**the same** bootstrap sample and both evaluated on the original data, so `dC_orig(b)` is **paired
at replicate level**. λ is retuned by inner 5-fold CV in every fit; knots, gene-level μ/σ and
M2-A2 score standardisation are all recomputed per fit.

**Stored outputs.** Only per-seed *aggregates* survive: `mean_optimism_C`, `corrected_C`,
`orig_C_p2.5/p97.5`, `dC_orig_mean/p2.5/p97.5`, and six cumulative means for MC-1. **The
replicate-level arrays are local variables and were discarded.** Nothing further can be extracted
post hoc without re-running.

## 2. Why `[+0.0002, +0.0240]` does not quantify uncertainty in ΔC_corrected

That interval is the 2.5–97.5 percentile spread of `dC_orig(b)`. Four independent defects:

1. **Wrong quantity.** `dC_orig(b)` is the paired ΔC of the *b*-th **bootstrap-fitted** model evaluated on the original data. Its spread measures **variability of the fitted model across resamples**, not sampling variability of the estimator ΔC_corrected.
2. **Wrong centre.** It is centred on **+0.0121**; the estimand is **+0.0043** — a factor of **2.80**. An interval centred 2.8× away from its own point estimate is not a confidence interval for it.
3. **The evaluation-sample variability is entirely absent.** Every `dC_orig(b)` is evaluated on the *same* 452 patients. The uncertainty arising from having *these* 452 patients rather than another 452 appears nowhere in the computation.
4. **It ignores that optimism is itself estimated.** ΔC_corrected subtracts a quantity with its own uncertainty; the percentile spread of `dC_orig` does not propagate it.

The existing interval is therefore closer to a **model-stability diagnostic** than to a confidence
interval, and it is neither a subset nor a superset of the correct variance.

> ### ⛔ A tempting shortcut that is **not** legitimate
> Transplanting the half-width (±0.0119) onto the corrected estimate would give
> `[−0.0076, +0.0162]`. **Do not do this.** It assumes
> Var(ΔC_corrected) = Var(dC_orig), which is unjustified: the two differ in both the fitting and
> the evaluation layer. The number is reported here only so it is recognised and refused.

## 3. Can the existing B=500 replicates support a valid interval? **No.**

They estimate the **bias term only**. Formally, the executed bootstrap draws
`D* ~ D̂` and evaluates on fixed `D`; the inner expectation estimates `E[optimism]`. A confidence
interval needs the distribution of the *whole* estimator under resampling of `D`, which requires an
**outer** layer that was never run. Re-running the identical design with more replicates would
shrink Monte-Carlo error and change nothing about this.

## 4. Candidate methods

| Method | Targets the estimand? | Assumptions | Cost | Verdict |
|---|---|---|---|---|
| **Nested (double) bootstrap** — outer patient resample, full optimism correction inside each | **Yes, exactly** | bootstrap validity for a rank statistic; inner MC noise controlled | 12 h single-core at B_out 500 / B_in 100 (~1.5 h on 8 cores) | **RECOMMENDED** |
| Jackknife / influence function | approximately | **smoothness** — the C-index is a non-smooth rank statistic, for which the jackknife is known to be inconsistent | 5–11 h | **Rejected** — assumption fails |
| Treat optimism as constant; bootstrap apparent ΔC and shift | only if optimism is dataset-stable | **untested.** The observed optimism stability (SD 0.0005) is *Monte-Carlo* stability on **one** dataset — it says nothing about variation **across** datasets | ~1 h | **Rejected as a primary method**; testable as a by-product of the nested run |
| .632 / .632+ bootstrap | no — a *different* estimator | — | moderate | **Rejected** — changes the estimand |
| Analytic / DeLong-type SE for a difference of C-indices | no | exists for *apparent* correlated C-indices only; no standard variance for an optimism-corrected difference with internal tuning | trivial | **Rejected** — not defensible here |

## 5. Recommended method

**Paired nested (double) bootstrap**, targeting the estimator exactly:

```
for k = 1..B_out:                       # OUTER: sampling uncertainty
    D*_k = patient-level resample of the 452          (event-stratified)
    for b = 1..B_in:                    # INNER: optimism, as in the frozen protocol
        fit M1 and M2 on a bootstrap of D*_k, on the SAME resample  (pairing preserved)
        optimism_M(b) = C(boot) − C(D*_k)
    ΔC_corrected(k) = ΔC_apparent(D*_k) − (mean_b opt_M2 − mean_b opt_M1)
interval = percentile of { ΔC_corrected(k) }
```

**Design points, each tied to the frozen protocol**

- **Pairing:** M1 and M2 share every inner bootstrap sample within every outer resample, so the difference is paired at both layers.
- **Event-stratified outer resampling** matches the frozen protocol. **Consequence, which must be stated:** fixing 136 events in every outer resample makes the interval **conditional on the observed event count**, and therefore slightly narrower than an unconditional one. Changing to unconditional outer resampling would be a **new scientific decision requiring approval** — it is not proposed here.
- **Re-fitting and tuning:** knots, gene-level μ/σ, M2-A2 standardisation and λ are recomputed inside every inner fit, exactly as now. No shortcut.
- **Inner Monte-Carlo noise inflates the outer variance:** Var(observed) ≈ Var(true) + E[Var_inner]. With B_in = 100 the inner MC-SE of optimism is ≈ 0.0011 per outer replicate. This is **estimable from the run itself** and a variance correction can be reported alongside the uncorrected interval.
- **Non-estimable replicates:** recorded with seed, outer index, inner index and reason; excluded from that outer replicate's optimism mean and **counted**; the 5 % threshold applied **per outer replicate**. **No silent redraws, no selective exclusion.** An outer replicate whose optimism cannot be computed at all is reported, not dropped quietly.
- **Interval type:** **percentile** as primary, with **basic (reverse-percentile)** reported alongside, plus a **bias-corrected (z₀-only)** variant — z₀ is computable cheaply from the outer distribution. **Full BCa is not proposed:** the acceleration constant needs a jackknife over 452 patients (5–11 h extra) and rests on the same smoothness assumption that disqualified the jackknife above.

**Sampling vs Monte-Carlo uncertainty stay separate.** The outer spread is the **sampling**
uncertainty; the inner-layer and finite-B_out noise are **Monte-Carlo** uncertainty. They are
reported separately and never summed, exactly as MC-4 requires.

## 6. Cost (measured, not guessed)

Timing probe: **0.852 s per inner replicate** for M1 + M2 four-set (20/20 fits succeeded).

| B_out | B_in | inner replicates | single-core | 8 cores |
|---:|---:|---:|---:|---:|
| 200 | 100 | 20 000 | 4.7 h | 0.6 h |
| **500** | **100** | **50 000** | **11.8 h** | **1.5 h** |
| 500 | 200 | 100 000 | 23.7 h | 3.0 h |
| 1000 | 100 | 100 000 | 23.7 h | 3.0 h |

**Proposed: B_out = 500, B_in = 100** — ~1.5 h on 8 cores.

**Monte-Carlo checks.** MC-1 (running trace of both interval endpoints against B_out), MC-2
(endpoints across the declared seeds), MC-4 (MC-SE reported beside each endpoint).
**MC-3 needs a decision:** its frozen criterion (`MC-SE ≤ 0.10 × sampling SE`) is written for a
**point estimate**. Applying it to interval *endpoints* is an **extension of the frozen rule** and
requires your approval; I have not assumed it.

**Seed policy also needs a decision.** The frozen policy declares five seeds for MC-2 on a point
estimate. For a nested design the options are (a) one declared seed for the main run plus a reduced
B_out across the other four for MC-2, or (b) all five at full B_out (5× cost). **Not chosen here.**

## 7. Code changes required

A new `scripts/40_phase9_dc_interval.py`, reusing `src/phase9_execution.py` **unchanged** and
importing `boot_indices` from the primary driver so the resampling is the same object. It must
**store replicate-level arrays** (the gap that makes the current outputs unextendable), write to a
**new** directory `reports/phase9/dc_interval/`, and parallelise across outer replicates.
**No protocol amendment, no change to existing B or seeds, no existing output touched.**

## 8. Limitations that persist whatever interval is produced

1. **Conditional on 136 events** under event-stratified outer resampling.
2. **Bootstrap validity for a rank statistic** with internal tuning is asymptotic; n = 452 with 136 events is not large.
3. **Harrell's C is censoring-distribution dependent** — a relative measure on *this* cohort, not transportable.
4. **The molecular block is penalised and the clinical block is not** (P6B-A1), so the comparison is **conservative by design and is not an unbiased estimate of the molecular contribution**. An interval does not repair that.
5. **Internal validation only.** No external cohort exists; no transportability claim is possible.
6. **U1 unresolved** — nothing here licenses any fixed-time or absolute-risk statement.
7. An interval around **+0.0043** quantifies uncertainty; it establishes **no** clinical utility and **no** biological mechanism.

## 9. Recommendation

**A defensible interval for ΔC_corrected is obtainable, but not from the current design.** The
executed B=500 run estimates only the bias term, and its replicate-level arrays were not retained.

**Recommended:** a paired, event-stratified **nested bootstrap at B_out = 500, B_in = 100**
(~1.5 h on 8 cores), reporting percentile, basic and bias-corrected intervals with the inner-noise
variance correction, Monte-Carlo diagnostics separated from sampling uncertainty, and full
non-estimable accounting.

**Three decisions are yours before any execution:**
1. Whether to apply **MC-3** to interval endpoints (an extension of the frozen rule).
2. The **seed policy** for the outer layer.
3. Confirmation that **event-stratified** outer resampling — and hence a **conditional** interval — is accepted.

**Awaiting approval. Nothing further will run.**
