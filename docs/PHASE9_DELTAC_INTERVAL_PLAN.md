# Phase 9 — Revised Execution Plan: Confidence Interval for the Optimism-Corrected ΔC

**Version:** 2.1 · **Date:** 2026-11-03 · supersedes v2.0 (retained in §12)
**Status:** ⛔ **PLAN ONLY — NOTHING EXECUTED, NOTHING AMENDED, NO COMMIT**

**Estimand (unchanged):** `ΔC_corrected = C_corrected(M2 four-set) − C_corrected(M1)` = **+0.0043**.
Median representation **primary**; S-7 **sensitivity-only**.


### Revision history

| Version | Status | Summary |
|---|---|---|
| 1.0 | superseded | Methodological audit; nested bootstrap identified as the only valid route |
| **2.0** | **superseded** | Execution plan: staged design C, EP-1 proposed, A7/A8 presented as acceptance criteria |
| **2.1** | **current** | Decisions 1–4 adopted. **Design changed to five full-size runs** (not staged C). **Inner-failure treatment rebuilt** — pairing-preserving exclusion + bounding analysis replaces "ignorable below 5 %". **A7/A8 demoted to warnings**, with A7 declared to have *no* distributional justification. RNG hierarchy, checkpointing and storage specified. |

---

## 1. Inner Monte-Carlo error — derivation, consequence, and why no deconvolution is used

### 1.1 How the figure is obtained

The v1.0 audit quoted "≈0.0011". **That was an unpaired hand-calculation and it was wrong in a way
that matters.** The correct derivation uses the *paired* optimism difference, because M1 and M2 are
fitted on the same inner bootstrap sample.

From the stored per-seed aggregates (R = 5 seeds, B = 500 each):

| quantity | value |
|---|---|
| SD across seeds, mean optimism M1 | 0.000585 |
| SD across seeds, mean optimism M2 | 0.000465 |
| **SD across seeds, paired difference Δopt** | **0.000422** |
| naive √(σ₁² + σ₂²) ignoring pairing | 0.000747 |
| **reduction due to pairing** | **44 %** → implied ρ ≈ **0.68** |

Per-replicate SD: **σ_Δ = 0.000422 × √500 = 0.00943**, hence

| B_in | 50 | 100 | 200 | 400 |
|---|---|---|---|---|
| inner MC-SE of Δopt | 0.00133 | **0.00094** | 0.00067 | 0.00047 |

**Precision of this estimate itself:** it rests on df = 4, so its relative SE is ≈ **35 %**
(σ_Δ = 0.00943 ± 0.0033). It is used **only for sizing B_in**, never reported as a result.

### 1.2 Contribution to the outer variance

For outer replicate *k*: `ΔC_corr(k) = ΔC_corr,true(D*ₖ) + εₖ`, with εₖ pure simulation noise given
D*ₖ, so to good approximation

> **Var(observed outer) = Var(sampling) + σ_Δ²/B_in**

The outer distribution is therefore **over-dispersed**, making the percentile interval
**conservative (too wide)** — never anti-conservative. Sized against a proxy outer SD of ≈ 0.0061:

| B_in | inner MC-SE | inflation of **variance** | inflation of **SD** |
|---|---|---|---|
| 50 | 0.00133 | 4.8 % | 2.4 % |
| **100** | **0.00094** | **2.4 %** | **1.2 %** |
| 200 | 0.00067 | 1.2 % | 0.6 % |

### 1.3 Resolution: **control it, do not subtract it**

At **B_in = 100 the inflation of the reported SD is ≈ 1.2 %** — an order of magnitude below the
endpoint Monte-Carlo noise itself. **No variance deconvolution will be performed.** That is
deliberate, and the reasons are stated rather than assumed:

- A percentile interval is a function of the **whole distribution**, not of its variance. There is no valid operation that "subtracts a variance" from a set of quantiles.
- Variance subtraction can return a negative estimate.
- It assumes exact additivity and independence of εₖ, which is an approximation.
- A normal-theory width correction would discard precisely the skewness-respecting property that motivates the percentile method.

**Preventing Monte-Carlo noise from being read as sampling uncertainty.** Three mechanisms:

1. **Sizing** — B_in chosen so the inflation (1.2 % of SD) is negligible and *reported as a known conservative bias of stated magnitude*.
2. **Separate reporting** — the interval carries the **sampling** uncertainty; the endpoint MC-SE (§2) is reported *beside* it and **never summed with it**, per MC-4.
3. **Direction stated** — the contamination makes the interval **wider**, so it cannot manufacture a spuriously narrow or significant result.

---

## 2. Endpoint stability — EP-1 (approved) and how it is estimated

**EP-1 (approved, decision 1):** `MC-SE(endpoint) ≤ 0.05 × (upper − lower)` for **both** endpoints,
with the **absolute** endpoint MC-SE reported **separately from sampling uncertainty**.

It is **not** the point-estimate MC-3 rule. What matters for an interval is that its endpoints are
pinned down finely relative to **the width being reported**.

### 2.1 Five quantities, reported separately and never combined

| # | Quantity | What it is | How obtained |
|---|---|---|---|
| **1** | **Sampling uncertainty** | the percentile interval itself — the spread of `ΔC_corr(k)` over outer resamples | the reported interval |
| **2** | **Within-run endpoint MC-SE** | noise in each quantile from finite B_out | kernel quantile SE (§2.2) |
| **3** | **Between-run endpoint variability** | SD of each endpoint across the five independent runs | direct, at full size |
| **4** | Kernel density detail | the `f̂(q_p)` entering (2) | §2.2 |
| **5** | Empirical endpoint variation | the five endpoint values themselves, listed | reported in full |

Quantity 1 is **sampling**; 2, 3 and 5 are **Monte-Carlo**. They are tabulated side by side and
**never summed**, per MC-4.

### 2.2 Kernel-based quantile MC-SE — method stated

`SE(q_p) = √(p(1−p)/B_out) / f̂(q_p)`, with `f̂` a **Gaussian kernel density** evaluated at the
endpoint. Bandwidth by **Silverman's rule**, with the estimate **recomputed at bandwidth × {0.5,
1.0, 2.0}** and all three reported — tail density estimates are bandwidth-sensitive and a single
number would overstate precision. Expected ≈ **0.00073**, i.e. **≈ 3.0 % of width**, at B_out = 500.

### 2.3 Limitations of five runs — stated, not glossed

- An SD from **5 observations has df = 4**, so the between-run endpoint SD carries a relative SE of ≈ **35 %**. It is a coarse estimate and is reported as such.
- **Agreement between the kernel and empirical estimators does not validate either.** Both treat the outer draws within a run as i.i.d.; a violation common to both would not show up as disagreement. Agreement is *consistency*, which is weaker than *correctness*, and will be described in those words.
- The five runs differ **only** in seed, so they bound Monte-Carlo variability — **not** any systematic error shared by all five (bootstrap validity at n = 452, the conditional design, the dilution effect of §7.2).

### 2.4 If EP-1 fails — the exact rule

Since `SE ∝ 1/√B_out`, the required size is `B_out′ = 500 × (observed ratio / 0.05)²`, rounded **up**
to the next rung of the pre-set ladder **{750, 1000, 1500, 2000, 3000}**:

| observed SE/width | required | ladder rung |
|---|---|---|
| 6 % | 720 | **750** |
| 8 % | 1 280 | **1 500** |
| 10 % | 2 000 | **2 000** |

Each run is **topped up** — additional outer replicates are drawn from the *same* `SeedSequence`
spawn list at indices `B_out … B_out′−1`, so existing replicates are reused unchanged, the draws
remain i.i.d., and the result is bit-reproducible. **The interval is never narrowed, re-centred or
re-quantiled to satisfy EP-1. B_out may only increase.** If the top ladder rung is exhausted, the
interval is reported with its measured endpoint MC-SE and labelled **imprecise** — not quietly
accepted.

## 3. Resampling implementation audit

### 3.1 Procedure

```
theta_hat = DeltaC_corrected(D)                       # +0.0043 on the A-19 452 / 136 events

for k = 1 .. B_out:                                   # OUTER -- sampling uncertainty
    I*_k = event-stratified patient bootstrap of D    # ONE index array, WITH multiplicities
    D*_k = D[I*_k]                                    # 136 events, 316 censored, by construction

    for M in {M1, M2}:                                # SAME I*_k for both models
        fit_full(M, D*_k)                             # complete pipeline, refitted
        C_app(M,k) = C(D*_k, lp)
    dC_app(k) = C_app(M2,k) - C_app(M1,k)

    for b = 1 .. B_in:                                # INNER -- optimism, WITHIN D*_k
        J*_{k,b} = event-stratified bootstrap of D*_k # ONE index array, WITH multiplicities
        for M in {M1, M2}:                            # SAME J*_{k,b} for both models
            fit(M, D*_k[J*_{k,b}])                    # complete pipeline, refitted
            opt(M,k,b) = C(D*_k[J], lp) - C(D*_k, lp) # evaluated on D*_k, NOT on D
    # paired exclusion: b enters only if BOTH models succeeded (section 6.2)
    dOpt(k) = mean_b opt(M2,k,b) - mean_b opt(M1,k,b)

    dC_corrected(k) = dC_app(k) - dOpt(k)

interval = percentile{ dC_corrected(k) }
```

### 3.2 Confirmation that this is the intended estimator

The estimator applied to the original data is `θ̂(D) = ΔC_app(D) − Δopt̂(D)`. The nested bootstrap
applies **the identical functional** to each `D*ₖ`, so the distribution of `θ̂(D*) − θ̂(D)`
approximates that of `θ̂(D) − θ`. Two consequences are load-bearing and are confirmed explicitly:

- **Optimism inside an outer replicate is evaluated on `D*ₖ`, not on `D`.** Evaluating on `D` would compute a different functional and would leak the original sample into the outer layer.
- **`ΔC_app(k)` is recomputed on `D*ₖ`** from its own full pipeline fit — it is not the original apparent value.

**Confirmed: this is the correct implementation of the optimism-corrected estimator.**

### 3.3 Pairing — enforced by construction at both levels

`I*ₖ` is drawn **once** per outer replicate and `J*_{k,b}` **once** per inner replicate; both
models consume the identical index array, **including multiplicities**. M1 and M2 never see
different patients, different duplicate counts, or different event/censoring composition.

### 3.4 Refit requirements — nothing cached across layers

| operation | refitted at | cached? |
|---|---|---|
| Age spline knots (10/50/90 pct of the fitting sample, `numpy.percentile` type 7) | **every fit** | never |
| Gene-level μ and σ, **`sigma_ddof = 1` (M2-A1)** | **every fit** | never |
| **M2-A2** pathway-score standardisation, with zero-variance drop (≤ 1e-12) and non-finite → non-estimable | **every fit** | never |
| Ridge λ, **frozen nested 5-fold CV** (P6B-A1: clinical unpenalised, molecular penalised) | **every fit** | never |
| Model fit, convergence and separation checks (ΔlogL < 1e-9, ≤ 100 iterations, \|β\| ≤ 10, SE ≤ 10) | **every fit** | never |
| Frozen panel membership, gene→probe map, cohort, endpoint | outside all layers | fixed metadata, by design |

**No fitted quantity is cached across models, across inner replicates, or across outer replicates.**
The frozen protocol permits caching of nothing that is estimated from data.

### 3.5 C-index and risk sets

Harrell's C via `lifelines.utils.concordance_index` on the **frozen** orientation (`-lp`, since
higher hazard ⇒ shorter survival), with the deposited `rfs_time_full` / `rfs_event_full` and Efron
tie handling — the same objects and conventions as the completed primary and S-7 runs. Duplicate
patients in a resample contribute **duplicate rows**, which increases ties; this is the standard
bootstrap behaviour and is unchanged from the frozen implementation.

## 4. Event stratification — **APPROVED (decision 3)**, interval is conditional

**Approved and retained** at **both** outer and inner levels, matching the frozen protocol.

> **Mandatory description, to appear wherever the interval is reported:**
> *"Conditional on the observed 136 recurrence events and 316 censored observations in the
> 452-patient A-19 population."*

**The resulting interval is CONDITIONAL on the observed event count of 136.** The limitation,
quantified:

| | event count per replicate |
|---|---|
| unconditional patient bootstrap | mean 136, **SD 9.75**, 95 % range **117–155** |
| **event-stratified (frozen)** | **fixed at 136** |

Since C-index precision scales roughly as 1/√(events), a ±19-event swing is ±14 % in events and up
to ≈ 8 % in SE at the low end. **The conditional interval is therefore expected to be a few per
cent narrower than an unconditional one.** This is an order-of-magnitude argument, not a measured
correction, and it is stated as such.

**Implication to carry into reporting:** the interval answers *"given a cohort with 136 recurrence
events, how precisely is ΔC_corrected determined?"* — not *"how would ΔC_corrected vary across
repeat studies of 452 patients?"* **No change to the scheme is proposed.** Switching to an
unconditional outer bootstrap would be a new scientific decision and would need your approval.

---

## 5. Monte-Carlo design — **APPROVED: five independent full-size runs**

The staged design C proposed in v2.0 is **withdrawn**. Per decision 2:

| | |
|---|---|
| Runs | **5 independent, full size** |
| B_out | **500** per run |
| B_in | **100** per outer replicate |
| Total | 2 500 outer · 250 000 inner replicates · **500 000 model fits** |
| Cost | ≈ **59.2 h single-core / 7.4 h on 8 cores** |

Model definitions, resampling procedure and frozen preprocessing are identical across all five
runs; **only the seed differs**. This removes the √0.5 extrapolation that design B/C required:
between-run endpoint variability is measured **directly at the size actually reported**.

**No pilot and no timing probe will be run** (instruction §6). The 0.852 s/inner-replicate figure
used for budgeting was measured *before* this task and is carried forward unchanged.

## 6. Non-estimable replicates — rebuilt

**The v2.0 rule ("exclude failures below 5 %, continue") is withdrawn. It assumed failures were
ignorable, and that assumption is not safe.**

### 6.1 Why exclusion is not neutral

Optimism is `C_boot − C_outer`. Fits fail on *strained* resamples — separation, a near-empty
covariate cell, an ill-conditioned information matrix. Those are plausibly the same resamples where
**overfitting, and therefore optimism, is largest**. Dropping them biases the optimism mean
**downward**, hence `C_corrected` **upward** — the optimistic direction. Silence about this would
be the worst option.

A second hazard is specific to a *paired* estimand: M2 has more parameters and a tuning step, so it
is the likelier to fail. Averaging `opt_M2` over a different replicate subset than `opt_M1` would
**break the pairing** that the whole design exists to preserve.

### 6.2 Rules

| | |
|---|---|
| **Recording** | every failure logged as `(run_seed, outer_k, inner_b, model, reason, n_distinct_patients, min_binary_cell, n_events_inner, lambda_at_failure)`. Reasons are the frozen set: non-convergence, iteration cap, singular/non-PD information, separation sentinel, non-finite scaling (M2-A2) |
| **Pairing-preserving exclusion** | if **either** model fails at inner replicate *b*, **both are excluded for that b**. `opt_M1` and `opt_M2` are *always* averaged over the identical replicate set. **New in v2.1** |
| **No substitution** | no redraws, no replacement replicates, no imputation of a failed optimism value |
| **Per-outer accounting** | each outer replicate reports `B_in_attempted`, `B_in_paired_complete`, failures by model and reason |

### 6.3 Are the survivors representative? — a test, not an assumption

For each run, failed and successful inner replicates are compared on the recorded covariates
(distinct-patient count, minimum binary cell, selected λ, inner event count) using a
distribution comparison, and **failures are cross-tabulated by model and by outer replicate** to
detect concentration in one model or in identifiable outer samples. Concentration is reported
whether or not a threshold is crossed.

### 6.4 Bounding analysis — replaces the ignorability assumption

Where any outer replicate has failures, `ΔC_corr(k)` is recomputed with the missing paired
optimism values imputed at the **minimum** and at the **maximum** observed within that outer
replicate, giving `[ΔC_corr^min(k), ΔC_corr^max(k)]`. The whole percentile interval is then
recomputed under both bounds.

> **Decision rule:** if both endpoints move by **less than the endpoint MC-SE**, the failures are
> *demonstrably* non-influential and the primary interval stands, with the bounds reported. If
> either endpoint moves by more, the failures are influential and the interval is labelled
> **UNRELIABLE** — regardless of how far below 5 % the failure rate is.

**The 5 % figure is a warning boundary only. It never licenses the assumption that failed fits are
ignorable.** Precedent suggests zero failures (0 / 7 500 inner replicates in both completed runs),
in which case every rule here records "0" and the bounding analysis is vacuous — which is itself
the reportable result.

### 6.5 When an outer replicate is unreliable

An outer replicate is declared unreliable, recorded and excluded from the quantile computation if
**any** of: its apparent fit fails · its paired-complete inner count falls below 95 % of B_in · its
bounding interval width exceeds the endpoint MC-SE. **> 5 % of outer replicates excluded ⇒ the
interval is labelled UNRELIABLE overall.**

## 7. Interval construction

**Percentile — primary** (decision 4). Transformation-respecting, needs no variance estimate, and
respects skewness in a bias-corrected functional of a rank statistic.

**Basic (reverse-percentile) — diagnostic.** Its pivotality assumption is weaker for a doubly
bootstrapped, bias-corrected estimator. Reported for comparison only: a material divergence is
evidence of skewness or bootstrap bias and a warning that **neither** interval should be trusted.

**z₀ — bias diagnostic.** `z₀ = Φ⁻¹(#{ΔC_corr(k) < θ̂}/B_out)`, free from the outer draws.

**BCa — excluded** (decision 4). Acceleration requires a jackknife that assumes the smoothness the
C-index lacks. No independently defensible acceleration estimate exists here.

### 7.1 A7 and A8 are **warning triggers, not acceptance criteria**

Per instruction §2B, both are demoted. Their standing is stated honestly:

| | threshold | statistical justification | standing |
|---|---|---|---|
| **A7** percentile vs basic | 25 % of width | **NONE.** There is no distributional result behind 25 %; it is a round number chosen for visibility. The *direction* of the comparison is informative, the cut-point is not | **WARNING — investigate and report** |
| **A8** \|z₀\| | 0.25 | **PARTIAL.** `SE(z₀) = √(p(1−p)/B)/φ(z₀)` = **0.056** at B_out = 500 and **0.025** pooled over five runs, so 0.25 is ≈ **4.5 SE** (single run) or **10 SE** (pooled). It is therefore a *loose detection* boundary for substantial median bias — **not** a validity criterion, and **\|z₀\| < 0.25 does not establish absence of bias** | **WARNING — investigate and report** |

**A breach of either triggers investigation and explicit reporting. Neither a breach nor a pass
establishes that the interval is valid or invalid.**

### 7.2 An anticipated, mechanistic source of z₀ ≠ 0

A nested bootstrap dilutes distinct patients: an outer sample holds ≈ **286 distinct** of 452
(63.3 %), and an inner sample of it ≈ **181** (40.0 %). Inner fits therefore see **fewer distinct
patients than the primary analysis did**, which can make optimism estimated *within* `D*ₖ` larger
than on `D`, shifting the outer distribution of `ΔC_corr(k)` **below** θ̂ and producing **negative
z₀**. This is a known property of the method, not a defect of the data, and it is pre-registered
here so it is interpreted rather than discovered.

## 8. Random-number and reproducibility policy

### 8.1 Seeds — new, and deliberately disjoint from the completed runs

The five Phase 9 seeds `[20261102 … 20261106]` drove the **single-layer** bootstrap of the
completed primary and S-7 runs. **Reusing them for the outer layer would make the new outer
resamples the same index sequences as the old inner resamples**, creating an unintended dependence
between this analysis and the one it is quantifying. They are therefore **not reused**.

New run seeds, declared here and never changed:

```
BASE = 20261103_00                  # distinct from every seed used previously
run_id = 0, 1, 2, 3, 4              # the five independent full-size runs
```

### 8.2 Deterministic stream hierarchy

```python
root_k = numpy.random.SeedSequence(entropy=[BASE, run_id])
outer  = root_k.spawn(B_out)                 # outer[k]  -> outer replicate k
inner  = outer[k].spawn(B_in)                # inner[b]  -> inner replicate b of outer k
cvseed = inner[b].spawn(2)                   # [0] = inner 5-fold CV, [1] = reserve
```

**Each replicate's stream is a pure function of `(BASE, run_id, k, b)`.** Consequences:

- **Independence between runs** — different `run_id` ⇒ disjoint spawn trees.
- **Independence between layers** — an inner stream descends from its own outer stream; no sharing.
- **Scheduling-invariance** — a replicate's draws do not depend on how many replicates ran before it, so results are identical under any worker count, any chunk size, and any completion order. This is the property that makes the parallel run reproducible.
- **Resumability** — a replicate can be recomputed in isolation from its indices alone, which is what makes checkpointing (§9.3) exact rather than approximate.

### 8.3 Multiplicities

A resample is carried as an **index array with repeats**, never a set. A patient drawn 3× supplies
3 rows to the design matrix, the Cox fit and the C-index. For each outer replicate the full
multiplicity vector is stored (counts per original patient id, 452 integers), so any replicate can
be reconstructed exactly without re-running the RNG.

### 8.4 Reproducibility verification

After the run, **20 outer replicates are selected by a fixed rule** (indices 0, 25, 50, …, 475 of
run 0) and recomputed from their stored seeds. Their `ΔC_corr(k)` must match **bit-identically**.
A mismatch is a **hard failure** (§11, H5), not a tolerance to be widened.

---

## 9. Computational plan

### 9.1 Execution sequence

1. Verify frozen-artifact hashes; abort on any mismatch.
2. Load the median matrix, covariates, endpoint and panels **once**, before forking workers.
3. For `run_id` 0…4, in order: dispatch outer replicates 0…499 across 8 fork-based workers; each worker computes one outer replicate end to end (apparent fit + 100 inner replicates, both models) and returns a record.
4. Append each completed outer record to that run's checkpoint file immediately.
5. After each run: compute percentile / basic / z₀, kernel endpoint MC-SE at three bandwidths, non-estimable accounting, bounding analysis.
6. After all five: between-run endpoint SD, EP-1 evaluation, reproducibility spot-check, final report.

### 9.2 Runtime, memory, storage

| | |
|---|---|
| Fits | 5 × 500 × (2 apparent + 100 × 2 inner) = **~505 000** |
| Runtime | **≈ 59.2 h single-core / 7.4 h on 8 cores** |
| Memory | expression matrix 21 815 × 452 float64 ≈ **79 MB**, loaded **before fork** so workers share it copy-on-write; observed single-process RSS was 372–424 MB, so **≈ 3.4 GB peak at 8 workers**. Headroom required, not assumed |
| Storage | outer-level 2 500 rows (negligible) · **inner-level 1 000 000 float64 ≈ 8 MB** · multiplicity vectors 2 500 × 452 int16 ≈ 2.3 MB · **total < 20 MB compressed** |

Inner-level retention is cheap, so **it is retained in full** — this is precisely the gap that made
the completed runs unextendable, and it is not repeated.

### 9.3 Checkpointing and recovery

Each run writes `run{run_id}_outer.jsonl`, one line per completed outer replicate, flushed on
write. On restart the completed `k` values are read and **skipped**; because streams are
index-derived (§8.2), a resumed run is **bit-identical** to an uninterrupted one. No partial
replicate is ever written — a record appears only when its 100 inner replicates have all been
attempted and accounted for.

### 9.4 Retained for independent audit

Per outer replicate: `run_id, k, seed entropy, multiplicity vector, ΔC_app, per-model apparent C,
per-model mean optimism, Δopt, ΔC_corr, B_in attempted, B_in paired-complete, failures by model and
reason, λ selected per model, age knots`. Per inner replicate: `C_boot` and `C_outer` for both
models. Plus the full failure log of §6.2.

### 9.5 Stopping conditions — stop rather than continue

- Any frozen-artifact hash changes mid-run.
- Outer non-estimable rate exceeds 5 % at any checkpoint.
- The bounding analysis (§6.4) shows failures are influential.
- Peak memory approaches the machine limit, risking a partial write.
- Projected runtime exceeds twice the budget, indicating a mis-specification.

In each case the run **halts and reports**; it does not continue to a result that would have to be
labelled unreliable afterwards.

---

## 10. Output location

**`reports/phase9/dc_interval/`** — a new directory, separate from `reports/phase9/execution/`
(median primary) and `reports/phase9/execution_s7/` (S-7). Neither is read-write by this analysis;
both are hash-verified before and after.

---

## 11. Acceptance criteria — hard failures separated from warnings

### 11.1 Hard failure conditions — the interval is not reported as valid

| # | Condition | Threshold | Consequence |
|---|---|---|---|
| **H1** | Outer replicates excluded as unreliable | > 5 % | interval labelled **UNRELIABLE** |
| **H2** | Bounding analysis on inner failures (§6.4) moves either endpoint | > endpoint MC-SE | interval labelled **UNRELIABLE** |
| **H3** | EP-1 unmet after exhausting the B_out ladder {750 … 3000} | — | interval reported and labelled **IMPRECISE** |
| **H4** | Any frozen artifact hash changed during execution | any | **abort**, report |
| **H5** | Reproducibility spot-check (§8.4) not bit-identical | any | **abort**, report |
| **H6** | Paired-complete inner count differs between M1 and M2 for any outer replicate | any | implementation defect — **abort** |

### 11.2 Warning conditions — investigate and report; never auto-validate or auto-invalidate

| # | Condition | Threshold | Justification of the threshold |
|---|---|---|---|
| **W1** | Inner fit failures, any model | > 0 | any failure triggers §6.3 representativeness checks and §6.4 bounding |
| **W2** | Inner failures concentrated in one model or a few outer replicates | qualitative | pairing asymmetry is the specific hazard |
| **W3** | Inner-noise SD inflation | > 5 % | **raise B_in**; never deconvolve the variance |
| **W4** | Kernel vs between-run endpoint MC-SE | disagree > 2× | both reported; agreement is *consistency*, not validation (§2.3) |
| **W5** | MC-1 endpoint running trace over the final 40 % of B_out | spread > EP-1 bound | raise B_out |
| **W6** | **A7** \|percentile − basic\| at either endpoint | 25 % of width | **NO distributional justification** — a visibility heuristic only |
| **W7** | **A8** \|z₀\| | 0.25 | **partial** — ≈ 4.5 SE at B_out = 500, ≈ 10 SE pooled; a *loose detection* boundary, and passing it does **not** establish absence of bias |

**No threshold is relaxed because a result is inconvenient. B_out and B_in may only increase.**

---

## 12. Superseded: the v2.0 staged design

v2.0 recommended staged design C (pilot 200 → 500 → 4 × 250, ≈ 5.0 h) to avoid spending the full
budget on a possibly mis-sized design, accepting a √0.5 extrapolation for between-run variability.
**Decision 2 selected five full-size runs instead** (≈ 7.4 h), which removes that extrapolation and
measures between-run endpoint variability directly at the reported size. The v2.0 reasoning is
retained here for audit; it is **not** the plan.

## 13. Standing limitations the interval cannot repair

Conditional on 136 events (§4) · bootstrap validity is asymptotic at n = 452 / 136 events ·
Harrell's C is censoring-distribution dependent and not transportable · the molecular block is
penalised while the clinical block is not (P6B-A1), so the comparison stays **conservative by
design and is not an unbiased estimate of the molecular contribution** · internal validation only,
no external cohort exists · U1 unresolved, so no fixed-time or absolute-risk statement is licensed.

**An interval around +0.0043 quantifies uncertainty. It establishes no clinical utility, no
biological mechanism, no external validity and no transportability.**

---

## 14. Unresolved methodological limitations

These are **not** resolved by anything in this plan and will persist in the result:

1. **Conditional design.** The interval conditions on 136 events / 316 censored. It answers "how precisely is ΔC_corrected determined *in a cohort with this event composition*", not "how would it vary across repeat studies of 452 patients". An unconditional variant would be a new scientific decision.
2. **Bootstrap validity is asymptotic.** n = 452 with 136 events, a non-smooth rank statistic, and internal tuning inside every fit. No finite-sample guarantee is available.
3. **Distinct-patient dilution (§7.2).** Inner fits see ≈ 181 of 452 distinct patients (40 %). The estimated optimism *within* outer samples may exceed that on `D`, biasing the outer distribution downward. The direction is predicted; the magnitude is not estimable from this design.
4. **Between-run variability rests on 5 observations** (df = 4, ≈ 35 % relative SE). It bounds Monte-Carlo noise only, never systematic error shared by all five runs.
5. **W6 (A7) has no statistical justification.** The 25 % cut-point is a visibility heuristic. The percentile-vs-basic comparison is directionally informative; its threshold is not principled.
6. **W7 (A8) is a loose detector.** Passing \|z₀\| < 0.25 does **not** establish absence of median bias.
7. **No deconvolution of inner noise.** The reported interval carries a known ≈ 1.2 % conservative SD inflation at B_in = 100. It is documented, not removed.
8. **Harrell's C is censoring-distribution dependent** — a relative measure on this cohort, not transportable.
9. **The molecular block is penalised and the clinical block is not** (P6B-A1). The comparison is **conservative by design and is not an unbiased estimate of the molecular contribution**. An interval does not repair that.
10. **Internal validation only.** No external cohort exists (Phase 4 G4 failed). **U1 unresolved** — no fixed-time or absolute-risk statement is licensed.

## 15. Recommendation

**The plan is ready for execution**, subject to nothing further from me.

All four decisions are incorporated: EP-1 approved as the endpoint criterion; five independent
full-size runs (B_out = 500, B_in = 100, ≈ 7.4 h on 8 cores); event-stratified outer resampling
with the conditional description made mandatory; percentile primary, basic and z₀ as diagnostics,
BCa excluded.

Two substantive methodological changes were made in response to the review, not merely
acknowledged:

- **Inner-failure treatment was rebuilt.** Pairing-preserving exclusion (if either model fails, both are excluded for that replicate) plus a **bounding analysis** that *tests* influence rather than assuming ignorability. The 5 % figure is now explicitly a warning boundary, and an outer replicate can be declared unreliable **below** it.
- **A7 and A8 were demoted to warnings**, with A7 stated to have **no** distributional justification and A8 only a partial one (≈ 4.5 SE at B_out = 500).

**The estimand is unchanged:** `ΔC_corrected = C_corrected(M2 four-set) − C_corrected(M1)`, point
estimate **+0.0043** under both the median primary and the S-7 sensitivity analysis. Median remains
primary; representative-IQR remains sensitivity-only. **An interval around +0.0043 quantifies
uncertainty and nothing more — it establishes no clinical utility, no biological mechanism, no
external validity and no transportability.**

**Awaiting explicit authorisation. Nothing will run until it is given.**
