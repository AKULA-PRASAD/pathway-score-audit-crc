# Methodological decision review: the four open Phase 9 decisions

| | |
|---|---|
| Status | **REVIEW AND RECOMMENDATIONS ONLY. NOTHING HERE IS ADJUDICATED OR IMPLEMENTED.** |
| Prepared | 2026-11-04, for the project owner |
| Scope | Decisions 1–4 in `PROJECT_COMPLETION_REPORT.md` §15 |
| Sources read in full or at the cited lines | `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1 (sha256/16 `d27c52a29930a338`, all 534 lines); `scripts/40_phase9_dc_interval.py` (frozen driver, `4b1738a4d0ee0e39`, all 206 lines); `src/phase9_execution.py` lines 119–206; `scripts/38_phase9_execute.py` lines 1–200; `reports/phase9/dc_interval/HALT_RECORD.json`, `RESUME_AUTHORIZATION.json`, `resume_driver.log`; `PROJECT_COMPLETION_REPORT.md`; `documentation/CV_SEED_DEVIATION_DRAFT_v2.md`; `documentation/STATISTICAL_CONVENTIONS.md`; `audit/EXECUTION_LEDGER.md`; `audit/final_qc.json`; `results/s1_interval_results.json`; `reports/phase9/execution/phase9_execution_results.json` |
| Computation | **None.** No bootstrap was run and no result was recomputed. The only arithmetic in this document combines values already stored in `results/s1_interval_results.json`, and each such case is marked *derived*. |

Status labels follow the project convention: **PASS**, **FAIL**, **PROVISIONAL**, **BLOCKED**, **NOT YET ADJUDICATED**. Each decision section separates **Facts** (plan text, code, stored values), **Interpretation**, and **Recommendation**.

---

## 1. Executive assessment

**Computation: PASS.** The five runs are complete and reproducible:
- H5 gave 20 / 20 bit-identical replicates.
- Final QC returned `all_ok: true`.
- The independent re-computation agreed on 166 / 166 values.

None of this bears on the four decisions below. They concern what the as-run numbers *mean*.

| # | Decision | What the evidence establishes | Recommended option | Status |
|---|---|---|---|---|
| 1 | CV-seed deviation | The driver does not implement plan line 364. Its folds use the same random starting state as the resample. All four properties that §8.2 states (lines 367–372) still hold. The effect on λ and ΔC_corr has not been measured. | Treat the plan as binding and record a protocol deviation. Run a protocol-conforming rerun with fresh seeds-as-planned. Declare **in advance** that its results will be primary, whatever they show. Keep the as-run results as the deviation-affected historical analysis. | **NOT YET ADJUDICATED**; every ΔC_corr result **PROVISIONAL** |
| 2 | EP-1 governing estimator | The plan names one estimator for "MC-SE": the kernel quantile SE with Silverman's bandwidth. It requires the ×0.5 and ×2 bandwidths to be **recomputed and reported**. It does not require them to pass. The between-run SD has no stated role in EP-1. | Govern EP-1 with the kernel SE at Silverman ×1. On the as-run data that gives 10 / 10 pass, so §2.4 mandates no top-up. Report ×0.5, ×2 and the between-run SD as supplementary. Disclose the ×0.5 exceedance and the W5 breach. | **BLOCKED** until adjudicated, then PROVISIONAL through Decision 1 |
| 3 | Interval interpretation | The plan designates the **percentile** interval as primary (line 309). It gives no rule for combining five runs. Here A7 and A8 measure the same thing: how far the outer distribution sits from θ̂. By the plan's own rules, a breach is a warning and does not change the interval's label. | Report all five per-run percentile intervals with the between-run variability, and no single combined headline. Show the W6 / W7 warnings prominently. Allow a qualitative "no evidence of improvement" only, as PROVISIONAL. Make no quantitative claim until the §7.1 investigation is complete. | **PROVISIONAL** |
| 4 | θ̂ reference value | θ̂ (+0.0043) is the plan's pre-specified point estimate (lines 6, 139). It comes from the same mathematical functional but under a third CV-seeding convention. Its reported MC-SE leaves out the randomness of the apparent-fit CV on D. The primary percentile interval does not use θ̂. | **Keep θ̂ unchanged.** Add, as a labelled sensitivity analysis, a measurement of how much θ̂ varies with the apparent-fit CV seed. | **NOT YET ADJUDICATED** |

**Scientific position (PROVISIONAL).**
- **Established, for the as-run implementation:** all 15 intervals (5 runs × percentile, basic and BC) include 0. That gives **no evidence that M2 improves optimism-corrected discrimination over M1**.
- **Not established:**
  - any improvement;
  - the absence of an improvement (the percentile upper limits reach +0.047 to +0.052);
  - that the interval is trustworthy in its location (W6 / W7 fire in every run);
  - protocol conformance (Decision 1).

---

## 2. Verified facts and source evidence

Evidence type: **P** = plan text, **C** = code, **S** = stored value, **D** = derived arithmetic on stored values.

| # | Fact | Source | Type |
|---|---|---|---|
| F1 | Estimand: "`ΔC_corrected = C_corrected(M2 four-set) − C_corrected(M1)` = **+0.0043**" | plan line 6 | P |
| F2 | `theta_hat = DeltaC_corrected(D)  # +0.0043`; "the nested bootstrap applies **the identical functional** to each `D*ₖ`" | plan lines 139, 165–166 | P |
| F3 | EP-1: "`MC-SE(endpoint) ≤ 0.05 × (upper − lower)` for **both** endpoints" | plan lines 84–85 | P |
| F4 | §2.1 row 2 is "Within-run endpoint MC-SE — kernel quantile SE (§2.2)". Row 3 is "Between-run endpoint variability — SD of each endpoint across the five independent runs". | plan lines 95–96 | P |
| F5 | "Bandwidth by **Silverman's rule**, with the estimate **recomputed at bandwidth × {0.5, 1.0, 2.0} and all three reported** — tail density estimates are bandwidth-sensitive and a single number would overstate precision. Expected ≈ 0.00073, i.e. ≈ 3.0 % of width, at B_out = 500." | plan lines 105–108 | P |
| F6 | §2.4: "If EP-1 fails — the exact rule". Ladder {750, 1000, 1500, 2000, 3000}. "Each run is **topped up**" from the same spawn list. "B_out may only increase." | plan lines 116–132 | P |
| F7 | "**Percentile — primary** (decision 4)". Basic: "a material divergence is … a warning that **neither** interval should be trusted". z₀ = Φ⁻¹(#{ΔC_corr(k) < θ̂}/B_out). | plan lines 309, 312–316 | P |
| F8 | A7 / A8: "**A breach of either triggers investigation and explicit reporting. Neither a breach nor a pass establishes that the interval is valid or invalid.**" | plan lines 330–331 | P |
| F9 | §11.2 heading: "investigate and report; **never auto-validate or auto-invalidate**". W5: "spread > EP-1 bound — raise B_out". | plan lines 461, 469 | P |
| F10 | Only H1 and H2 produce **UNRELIABLE**. Only H3 ("EP-1 unmet *after exhausting the B_out ladder*") produces **IMPRECISE**. | plan lines 454–456 | P |
| F11 | §8.2: `cvseed = inner[b].spawn(2)  # [0] = inner 5-fold CV, [1] = reserve`. Stated consequences: pure function of `(BASE, run_id, k, b)`; independence between runs; independence between layers; scheduling invariance; resumability. The plan gives no apparent-level CV seed. | plan lines 361–372 | P |
| F12 | Driver: outer resample from `default_rng(ss_outer)` (lines 97–100). Apparent CV from a **fresh** `default_rng(ss_outer)` (116–117). Inner resample from `default_rng(ss_inner[b])` (134–135). Inner CV from a **fresh** `default_rng(ss_inner[b])` (141–142). No `spawn(2)` call anywhere. The docstring says "Implements … v2.1 exactly" (line 3). | driver | C |
| F13 | The CV folds come from `rng.permutation` of the fitting sample's **row positions**, by event stratum (`_event_stratified_folds`, `src/phase9_execution.py:131–138`). `tune_lambda` is called only when `p > n_clin` (196–201), so M1 never consumes the CV generator. | shared code | C |
| F14 | Executed CV generator state equals the resampling generator state in 9 / 9 tested cases. Planned `spawn(2)[0]` state equals it in 0 / 9. | `validation/cv_seed_state_evidence.json` | S |
| F15 | Phase 9 primary (the source of θ̂):<br>• apparent fit on D uses a fresh `default_rng(SEEDS[0]=20261102)` per model (`38:90–92`);<br>• the bootstrap uses **one advancing** generator per seed, for the resample and then for each model's CV in turn (`38:117–129`);<br>• the corrected C is `apparent C − mean optimism`, with the **same apparent C for all five seeds** (`38:145`);<br>• MC-SE is the SD of the five seed-level corrected values (`38:170–176`). | `scripts/38_phase9_execute.py` | C |
| F16 | θ̂ = 0.0043163 = ΔC_app(D) 0.0139635 − ΔOpt(D) 0.0096472 (2 500 single-layer replicates). Its MC-SE of 0.000189 is the between-seed SD of 0.000422 divided by √5. | `s1_interval_results.json#theta_hat` | S |
| F17 | Apparent M2 λ on D = 0.351. Apparent M2 λ on D*ₖ: median 0.152 in every run, range 0.001–10 (the full grid); 2.5 % and 97.5 % quantiles at the grid limits. | `phase9_execution_results.json#apparent`; `s1#runs[].distribution.lambda_apparent_M2` | S |
| F18 | Pooled outer mean ΔC_corr(k) = 0.01402 against θ̂ = 0.00432, a shift of +0.00970. Of that, ΔC_app contributes +0.00828 (mean 0.02224 vs 0.01396) and ΔOpt contributes −0.00142 (mean 0.00822 vs 0.00965). | `s1#bias_decomposition_observed` | S |
| F19 | z₀ ranges from −0.577 to −0.429 by run, and is −0.512 pooled. The plan's SE(z₀) is 0.059 for one run and 0.026 pooled. Moving θ̂ by ±2 MC-SE changes z₀ by at most 0.044. | `s1#between_run.z0`, `#pooled_descriptive.z0`, `#z0_sensitivity_to_theta_MC_error` | S |
| F20 | A7 is 52.2–56.4 % of the width, **identical at both endpoints in every run**. | `s1#runs[].warnings.W6_A7_fraction_of_width` | S |
| F21 | Run 1, upper endpoint, width 0.055924, kernel SE by bandwidth:<br>• ×0.5: 0.002993 (ratio 0.0535, **FAIL**);<br>• ×1: 0.002663 (0.0476, PASS);<br>• ×2: 0.001969 (0.0352, PASS);<br>• scipy-variant Silverman cross-check: 0.002433 (0.0435).<br>The between-run SD of the upper endpoint is 0.002186 (0.0391, PASS). | `s1#runs[1].kernel_MCSE_2_2`, `#EP1` | S / D |
| F22 | W5, run 1 upper: the spread over B ∈ {300, 400, 500} is 0.002986, against a bound of 0.002796 (**breach**). The trace values are 0.04525, 0.04388 and 0.04687. W5 holds in every other run and at every other endpoint. The grid is a convention (Convention 12). | `s1#runs[1].MC1_running_trace`; `STATISTICAL_CONVENTIONS.md` #12 | S |
| F23 | The observed throughput was ≈ 2.13 outer replicates/min on 8 cores. Run 0 halted on the §9.5 runtime stop condition and resumed under a user-authorised exception to that ceiling only. Runs 1–4 took 199–264 min each. | `HALT_RECORD.json`; `RESUME_AUTHORIZATION.json`; `resume_driver.log` | S |
| F24 | The Phase 9 modelling protocol declares the five seeds (line 694). A keyword search found no CV-stream hierarchy there. A full conformance audit of Phase 9 primary has **not** been done. | `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` | P (search only) |

---

## 3. Decision 1: CV-seed deviation

### 3.1 Evidence

1. **What the plan specifies (F11).**
   - The inner 5-fold CV draws from `inner[b].spawn(2)`, element `[0]`: a child stream separate from the stream that generates the inner resample.
   - The plan says **nothing** about the CV seed for the apparent fit on D*ₖ.
   - Its stated rationale (lines 367–372) concerns reproducibility and independence *between runs and between layers*. It does not explain why the CV needs its own child stream *within* a fit.

2. **What each implementation does.**

| | outer / inner resample | inner CV | apparent CV on D*ₖ | apparent CV on D |
|---|---|---|---|---|
| Plan §8.2 | `outer[k]`, `outer[k].spawn(B_in)[b]` | separate child `inner[b].spawn(2)[0]` | not specified | not in scope |
| Frozen driver | same (equivalence verified, F14) | fresh generator on `ss_inner[b]`: **same starting state as the inner resample** | fresh generator on `ss_outer`: same starting state as the outer resample | not computed by the driver |
| Phase 9 primary | one advancing generator per seed | the bits **after** the resample (and after earlier models' CV); not a pure function of `(seed, b)` | n/a (single layer) | fresh `default_rng(20261102)` |

3. **Which calculations are affected.**
   - **Inner level: a textual deviation** from line 364. It covers 250 000 inner M2 fits.
   - **Apparent level: a gap in the plan**, not a textual deviation. It covers 2 500 apparent M2 fits, and the as-run choice has the same coupling property.
   - Through M2 λ, the deviation reaches:
     - C_app(M2), `opt_M2`, ΔC_app(k), ΔOpt(k) and ΔC_corr(k);
     - every interval, z₀, BC, kernel SE, between-run SD, EP-1, W3–W7 and figure.
   - Even the failure counts (H1, H2, H6, W1) depend on M2 fits, although all were 0 as run.
   - **Not affected:**
     - the outer and inner resample indices (stored `multiplicity_sha`, H5);
     - M1 (stored apparent λ = 0 in 2 500 / 2 500; the inner level is inferred from F13);
     - θ̂, which was computed by a different script. Its *comparability* is a separate question (Decision 4).

### 3.2 What is established

- **Reproducibility is not affected.** The driver's streams remain a pure function of `(BASE, run_id, k, b)`. All four consequences in lines 369–372 hold, and H5 confirms bit-reproducibility. **H5 PASS establishes reproducibility of the as-run procedure, not conformance** (completion report §3; draft v2 §9).
- **The statistical procedure is affected.** Under the plan, the fold assignment is pseudo-independent of the resample. In the driver it is a deterministic function of the same random bits (F14). The fold algorithm (`permutation` over row positions, F13) differs from the resample algorithm (`choice` with replacement), so the folds are not a copy of the resample. But nothing guarantees that the two are independent.
- **The driver's line-3 claim** ("implements v2.1 exactly") is false for this element.
- **The resume authorisation did not amend the plan.** `RESUME_AUTHORIZATION.json` waives the runtime ceiling only (scope: "runtime ceiling only"). It froze the code as run; it does not make the code's seeding the protocol.

### 3.3 What remains uncertain

- Whether the coupling changes the **distribution** of selected λ and therefore of ΔC_corr(k), or only which random fold split each fit receives.
- **One mechanism by which it *could* be systematic** (a hypothesis, untested):
  - the folds partition row positions, and the resample sets which positions hold duplicate patients;
  - a deterministic link between the two could change how often duplicates fall on both sides of a CV split;
  - that kind of leakage tends to favour weaker penalties.
- Nothing links the deviation to the observed z₀ shift (draft v2 §7), and nothing rules a link out.
- Whether Phase 9 primary and S-7 conform to their own protocol's seeding text (F24).

### 3.4 Options

| Option | Description | Advantages | Disadvantages |
|---|---|---|---|
| **1A** | Treat the driver as binding: amend §8.2 retroactively to match the as-run seeding | no computation | Changes a pre-specified protocol after results were seen. Contradicts the plan's status as the authoritative pre-specification. The resume record gives no basis for it. **Not defensible.** |
| **1B** | Plan binding; record the deviation; keep the as-run results as the reported results, with disclosure (draft v2 §12); no further computation | cheap and honest | The results rest on an untested assumption (that the coupling is distributionally inert) and cannot be called protocol-conforming |
| **1C** | 1B plus a **pre-specified paired sensitivity** on a fixed subset: the 20 H5 indices of run 0, or a larger fixed set, recomputed with `spawn(2)[0]` inner-CV seeding | Measures the effect replicate by replicate. The resamples are identical, so differences isolate the seeding. M1 values must reproduce bit-identically, which is a built-in check. Cost ≈ 10 min for 20 replicates at the observed throughput (F23). | Cannot recompute interval endpoints. Statistical power is limited by the unknown SD of the paired differences. Needs a decision rule fixed in advance. Still produces no protocol-conforming interval. |
| **1D** (recommended) | Plan binding; record the deviation; **full protocol-conforming rerun** (5 × 500 × 100) with a new, separately frozen driver into a new directory. **Declared primary in advance, whatever it shows.** The as-run results are kept as the historical, deviation-affected analysis, and the full paired comparison becomes the sensitivity analysis. | The only option that produces protocol-conforming results. The paired design isolates the deviation exactly, at full size, with no extra modelling assumption. Pre-declaration removes the risk of choosing whichever set looks better. | ≈ 19.5 h at the observed throughput (the plan's 7.4 h budget proved too low, F23). Needs new code (review, freezing, an H5-type check), an owner decision on the apparent-level CV seed (the plan is silent), and a plan clarification recorded **before** execution. |

### 3.5 Recommendation: **1D**, with 1C as the minimum acceptable alternative

**Basis.**
- The plan is the authoritative pre-specification. Its line 364 is unambiguous.
- Only a conforming computation yields protocol-conforming results.
- The pairing by identical resamples makes the comparison unusually clean.

**Conditions to settle before any run:**
1. **The apparent-level CV seed convention.**
   - It must leave the outer and inner resample streams bit-identical. For example, derive it from an explicit `spawn_key`, not from the order of `spawn()` calls: spawning a new child before the inner children would renumber the inner streams.
   - Whether to keep the as-run apparent convention or give it its own child stream is a genuine owner choice. The plan is silent.
2. **The primary designation, written before execution:** the conforming results are primary, whatever their direction.
3. **A verification gate for the conforming run:**
   - the multiplicities match the stored `multiplicity_sha` values;
   - the M1 quantities are bit-identical to the as-run values;
   - an H5-equivalent spot-check passes.

**What would change if 1D is chosen.**
- Every ΔC_corr-derived number is recomputed from the conforming checkpoints. Decisions 2–4 are then applied to that dataset.
- The as-run results move to a "historical / deviation-affected" section, beside the paired comparison.
- The PROVISIONAL label on the primary results would be removed only after the conforming run passes its own H1–H6 checks.

**If 1C is chosen instead:** the as-run results stay primary, but remain conditional on the deviation. The subset only says whether the coupling *appears* material, and the decision rule must be fixed before the subset is run.

---

## 4. Decision 2: EP-1 and Monte-Carlo error

### 4.1 Evidence

- **EP-1 (F3)** names "MC-SE(endpoint)". The only quantity the plan itself calls an "endpoint MC-SE" in §2.1 is row 2: the **within-run kernel quantile SE** (F4).
- **Row 3, the between-run SD,** is called "variability". Its precision is described as coarse (df = 4, ≈ 35 % relative SE; line 112). The W4 row (line 468) does call both "endpoint MC-SE", but only to compare them.
- **The plan's EP-1 sizing** (≈ 0.00073, ≈ 3.0 % of width, lines 107–108) is computed with the kernel formula.
- **The §2.4 scaling rule** `B_out′ = 500 × (ratio / 0.05)²` (line 118) assumes an SE ∝ 1/√B_out estimated from one run.
- **The bandwidth (F5):** the rule is "Silverman's rule". The ×{0.5, 1, 2} requirement is "recomputed … and all three reported", and the stated reason is that "a single number would overstate precision".

### 4.2 Answers to the five questions

1. **Which estimator does the plan prescribe?**
   - The **kernel quantile SE at Silverman's bandwidth** (lines 95, 105–108).
   - This is the best-supported reading, **but the plan never states explicitly** "EP-1 uses the ×1 kernel SE". `STATISTICAL_CONVENTIONS.md` #10 recorded this as unstated and deliberately evaluated every estimator instead of picking one.
2. **Are the bandwidth evaluations part of the acceptance criterion?**
   - The text makes them a **reporting requirement** ("recomputed … and all three reported").
   - Nothing says EP-1 must hold at all three. The stated purpose is to show how uncertain the SE estimate itself is.
   - **Supplementary**, on the plain wording.
3. **Does the between-run SD have a formal role in EP-1?**
   - **No.** It is quantity 3, reported beside the kernel SE and compared with it in W4. No sentence makes it an input to EP-1.
4. **Is the B_out top-up mandatory?**
   - **Only if EP-1 fails** (§2.4 heading, "If EP-1 fails — the exact rule").
   - Under the ×1 kernel reading, EP-1 **PASSES** in 10 / 10 evaluations on the as-run data (maximum ratio 0.0476, run 1 upper), so **§2.4 mandates no top-up**.
   - Under an "every bandwidth must pass" reading it fails (0.0535), and §2.4 requires B_out = 750 (required 573, rounded up to the next rung).
5. **What preserves the error-control objective without an unauthorised change?**
   - Govern EP-1 with the plan's named estimator. Report the sensitivity grid in full, and disclose its single exceedance.
   - That neither tightens nor relaxes a stated threshold, and it reports every number.

### 4.3 What remains uncertain

- **The margin under the recommended reading is small:** 0.0476 against 0.05 at run 1 upper. The ×0.5 bandwidth gives 0.0535 there. So **the conclusion "EP-1 PASS" is itself bandwidth-sensitive at one endpoint**, which is exactly what §2.2 warns about.
- **W5 is a separate signal.** Run 1 upper breaches (0.002986 against 0.002796, F22), and the W5 row says "raise B_out".
  - The §11.2 heading classes W5 as a warning ("investigate and report; never auto-validate or auto-invalidate").
  - The breach depends on the three-point grid of Convention 12.
  - Whether "raise B_out" is a directive or a suggested response is **ambiguous in the plan**.
- **Adjudicating the reading now is not blind.** The results are already known, so any adjudication happens after they were seen. The strict-reading outcome must be reported beside it.
- **If Decision 1 = 1D,** EP-1 is evaluated afresh on the conforming dataset. The as-run EP-1 result becomes historical.

### 4.4 Options

| Option | Description | Advantages | Disadvantages |
|---|---|---|---|
| **2A** (recommended) | EP-1 governed by the kernel SE at Silverman ×1. ×0.5, ×2 and the between-run SD reported as supplementary. The ×0.5 exceedance and the W5 breach disclosed. The W5 breach investigated from the **stored** fine trace (no new replicates). | Follows the plan's wording. No unauthorised change. No computation. | Adopted after the results were seen. The margin is thin. W5 is not acted on by a top-up. |
| **2B** | Strict: EP-1 must pass at every bandwidth, so top up every run to B_out = 750 per §2.4 | Most conservative error control. Also answers W5's "raise B_out". | Reads a reporting requirement as an acceptance criterion, which the text does not support. ≈ 1 250 new outer replicates (≈ 9.8 h at the observed throughput). Appends to checkpoints that the current rules protect. If Decision 1 is open, the top-up must use the **same** driver as the data it extends, which would propagate the deviation. |
| **2C** | 2A, plus a top-up triggered by W5 alone | Answers W5 | Treats a warning as a mandatory action, against the §11.2 heading. Same cost and checkpoint issues as 2B. |

### 4.5 Recommendation: **2A**

**Basis.** Plan lines 84–85, 95, 105–108 and 116–118.

**Record the following with the adjudication:**
- it was made after the 40-row grid was seen;
- the strict reading would require B_out = 750;
- the run-1 upper W5 breach and its dependence on the grid convention.

**What would change.**
- EP-1 moves from **BLOCKED** to **PASS** (PROVISIONAL through Decision 1). T09 and T10 in the ledger would be closed as "no top-up mandated".
- **H3 remains NOT TRIGGERED.**
- Under 1D, the same rule (×1 governs) is applied unchanged to the conforming dataset, and the top-up question is re-evaluated there and only there.

---

## 5. Decision 3: Interval interpretation

### 5.1 What A7 and A8 mean here

**Facts and derivation.** The basic interval is the reflection of the percentile interval through θ̂: `[2θ̂ − q.975, 2θ̂ − q.025]`. So at **each** endpoint,

`|percentile − basic| = |q.025 + q.975 − 2θ̂|`,

and **A7 = 2 · |midpoint − θ̂| / width**.

That explains why A7 is identical at both endpoints in every run (F20). Check on run 0 (*derived*):
- midpoint = (−0.008856 + 0.051630) / 2 = 0.021387;
- 2 × (0.021387 − 0.004316) / 0.060486 = **0.564**, matching the stored 56.4 %.

**Interpretation.**
- **A7 and A8 are not two independent confirmations.** Both measure how far the outer ΔC_corr(k) distribution sits from θ̂:
  - A7 measures the offset of the interval **midpoint**, which combines right skew (skewness 0.79 pooled) with the median offset;
  - A8 measures the **median** offset alone. Pooled, 30.4 % of the replicates fall below θ̂, so z₀ = −0.51, about 19 pooled SEs (*derived*: 0.512 / 0.026).
- **Under the bootstrap principle (§3.2),** an offset this large means one of the following:
  - (a) θ̂ is substantially biased for θ;
  - (b) the bootstrap distribution does not mimic the sampling distribution of this tuned, rank-based estimator at n = 452;
  - (c) θ̂ and the replicates are not realisations of the same functional (Decisions 1 and 4).
- **The evidence available does not separate (a), (b) and (c).**
  - **Supportive of (b), not decisive:** the shift is mainly in ΔC_app (+0.0083 of +0.0097, F18), and the median apparent M2 λ on D*ₖ (0.152) is one grid step below λ on D (0.351, F17). That fits weaker penalisation on resamples that contain duplicates, but it is not a demonstration.
  - **Bearing on (c):** θ̂'s own MC-SE cannot produce the offset (F19). The apparent-level CV randomness in θ̂ is unmeasured (Decision 4).
- **Plan §7.2 does not explain the offset.**
  - It predicted larger within-D*ₖ optimism and a distribution *below* θ̂. The observed optimism is *smaller* (−0.0014) and the distribution sits *above* θ̂.
  - §7.2 is also internally inconsistent: a distribution below θ̂ gives z₀ > 0 under the §7 definition, not "negative z₀" (final report §24).

### 5.2 Can an interval be designated primary?

- **Yes, at the level of method.**
  - The plan designates the **percentile** interval as primary (line 309).
  - A W6 / W7 breach does not change its label (F8, F9, F10). UNRELIABLE and IMPRECISE are reserved for H1–H3.
  - Line 314 ("neither interval should be trusted") is therefore a **warning that must be reported prominently**, not a reclassification.
  - The plan is in mild tension here, between line 314 and lines 330–331. The reading that satisfies both is: report it, with the warning, and do not validate it.
- **Basic and BC remain diagnostics** (lines 312–316; BC per audit line 98).
- **The percentile interval does not depend on θ̂.** Basic, BC, z₀, A7 and A8 do.

### 5.3 Five runs: separate or summarised?

**Facts.**
- The runs differ only in seed (lines 242–243).
- Between-run variability is "measured directly **at the size actually reported**" (line 244), which implies a reported interval at B_out = 500, the size of a single run.
- §2.1 row 5 requires the five endpoint values to be "reported in full" (line 98).
- **No rule picks one run or combines them.**

**Options.**

| Option | Description | Advantages | Disadvantages |
|---|---|---|---|
| **3A** (recommended) | Report all five per-run percentile intervals as the primary result set, with the between-run SD and range of each endpoint (quantities 3 and 5). No single combined headline. The pooled 2 500-draw interval only as **descriptive**, so labelled. | Uses only plan-defined quantities. Involves no selection. The five intervals agree closely (lower-endpoint range 0.0009, upper-endpoint range 0.0048). | No single headline number |
| **3B** | One run as the headline by a fixed rule (for example run 0, the lowest index) | One number | The rule is adopted after the results were seen. It discards four equally valid runs. |
| **3C** | Mean of the endpoints, or the pooled 2 500-draw quantiles, as the headline | Smaller MC error | A new estimator, not in the plan, and at a different size from the one the plan says is reported (line 244). Would need an amendment. |
| **3D** | Withhold any interval until the A7 / A8 investigation concludes | Maximally cautious | The plan forbids auto-invalidation on a warning (line 461). It is legitimate to withhold a *scientific conclusion*, but not to suppress the computed interval. |

### 5.4 How to communicate the uncertainty

Suggested wording, for PROVISIONAL use:

> *Conditional on the observed 136 recurrence events and 316 censored observations in the 452-patient A-19 population* (mandatory, lines 206–208): across five independent runs, the 95 % percentile intervals for ΔC_corrected ranged from about −0.010 to +0.052, and every interval, including the diagnostic basic and bias-corrected variants, contained zero. The analysis provides **no evidence that M2 improves optimism-corrected discrimination over M1**, and does not exclude improvements of up to about +0.05. The bootstrap distribution was **not centred on the point estimate** (+0.0043; z₀ ≈ −0.5 in every run), which the plan treats as a warning that neither the percentile nor the basic interval should be trusted for location. The interval method changes the endpoints materially: the basic interval's upper limit is about +0.018. These results were computed with an implementation that deviates from the pre-specified CV seeding, and they are provisional.

**Principles to keep:**
- Say "no evidence of improvement", never "no improvement".
- Quote no single magnitude as established.
- State the method dependence.
- Never present the zero-inclusion as resolving W6 / W7.
- Keep the §13 limitations attached to any reported interval, including the conservative-by-design penalisation (P6B-A1).

### 5.5 Is more analysis needed before a scientific conclusion?

- **For the qualitative "no evidence of improvement":**
  - it holds across all 15 intervals as run, because zero-inclusion survives every method;
  - it still depends on Decision 1. It becomes final only once the primary dataset is settled.
- **For any quantitative statement** (magnitude, or which interval locates the effect):
  - **yes**: the plan requires the W6 / W7 breaches to be **investigated** (lines 330–331), and so far only the descriptive shift decomposition has been done.
  - Candidate investigations, all from **stored** records, with no new replicates (§9, C5):
    - the λ distribution on D*ₖ against λ on D;
    - ΔC_app(k) against λ and against the number of distinct patients;
    - the contribution of skewness against the median offset to A7.
  - Plus the θ̂ sensitivity of Decision 4 (C4).

### 5.6 Recommendation: **3A**, with the wording in §5.4 and the C5 investigation before any final quantitative claim

**Basis.** Plan lines 98, 206–208, 244, 309, 314, 330–331 and 461.

**What would change.** The final report's headline becomes the five-run table with between-run variability. Any "mean across runs" phrasing (for example "averages [−0.0092, +0.0499]") would be explicitly **descriptive**. The W6 / W7 warnings move into the headline paragraph.

---

## 6. Decision 4: θ̂ reference value

### 6.1 Evidence

1. **The plan's definition.**
   - θ̂ = ΔC_corrected(D) = ΔC_app(D) − ΔOpt(D), stated numerically as +0.0043 (lines 6, 139, 165).
   - The plan asserts that the nested replicates apply "the identical functional" to D*ₖ (line 166).
   - θ̂ enters z₀ (line 316), the basic interval and A7. It is the stated centre of the analysis ("an interval around +0.0043", line 493).
2. **How the existing θ̂ was computed (F15, F16).**
   - A single apparent fit per model on D, with M2 λ = 0.351 chosen by one CV seeded `default_rng(20261102)`.
   - Then 5 × 500 single-layer optimism replicates, evaluated on D, using the advancing-generator CV convention.
   - θ̂ = mean over seeds of the corrected difference.
3. **Same functional?**
   - **At the level of definition, yes.** It uses the same `pipeline_fit`, the same event-stratified resampling algorithm (`boot_indices` ≡ `strat_boot`), the same C-index, and optimism evaluated on the analysis sample. That matches §3.2.
   - **At the level of the randomised tuning step, not demonstrably.** Three CV conventions are involved:
     - θ̂ uses folds drawn from bits *after* the resample, which is the convention closest to the plan's intent of independence;
     - the replicates use folds drawn from the *same* bits as the resample;
     - the plan specifies a separate child stream.
   - If the coupling is distributionally inert, θ̂ and the replicates estimate the same functional. If it is not, they differ, and the difference would show up exactly where z₀, A7 and A8 look: in the location.
   - **The difference in optimism replicates (2 500 against 100) affects Monte-Carlo noise only.** ΔOpt is a mean, and W3 already sizes that noise (SD inflation 0.34–0.36 %).
4. **Effect on the centre and its interpretation.**
   - **θ̂'s reported MC-SE (0.00019) leaves out the randomness of the apparent-fit CV on D.**
     - ΔC_app(D) is held fixed across the five seeds (`38:145`), so the between-seed SD measures optimism noise only.
     - The existing ±2 SE check on z₀ (F19) therefore does not cover the λ-selection randomness in ΔC_app(D).
     - Given that λ on D*ₖ spans the whole grid (F17), that randomness could be non-trivial. Its size is **unmeasured**.
   - **Unaffected by θ̂:** the percentile interval (the primary) and EP-1.
   - **Affected by θ̂:** basic, BC, z₀, A7 and A8.
5. **Is a same-functional reference needed?**
   - **Not as a replacement.**
     - θ̂ is the pre-specified estimate (line 6, "Estimand (unchanged)").
     - Re-deriving the centre after seeing that the diagnostics fire would be a result-driven re-centring, which runs against the plan's stated discipline ("never … re-centred", lines 129–130, in the EP-1 context).
     - There is also no plan-defined way to apply the *nested driver's* convention to D. D is not an outer replicate, so the apparent-level CV seed for D would have to be invented.
   - **Yes, as a sensitivity analysis:** to learn whether the W6 / W7 breaches are robust to θ̂'s unmeasured CV component.

### 6.2 What remains uncertain

- How much ΔC_app(D), and therefore θ̂, varies with the apparent-fit CV seed.
- Whether that variation is large enough to move z₀ across the A8 threshold. That would require θ̂ to rise to roughly the 40th percentile of the outer distribution, far above its current position at the 30th. This is plausible only if ΔC_app(D) is very sensitive to λ.

### 6.3 Options

| Option | Description | Advantages | Disadvantages |
|---|---|---|---|
| **4A** | Keep θ̂ = +0.0043 unchanged; disclose its seeding convention and the missing MC-SE component | Pre-specified. No computation. | Leaves the W6 / W7 attribution open |
| **4B** | Replace θ̂ with a re-derived value under the nested driver's (or the plan's) convention | Same seeding convention as the replicates | Post hoc re-centring. Needs a seed convention for D that the plan doesn't define. Changes the pre-specified estimate. **Not recommended.** |
| **4C** (recommended) | 4A, plus a labelled sensitivity analysis: refit the apparent M1 and M2 on D under a fixed, pre-declared set of CV seeds (for example 50 or 100); report the distribution of ΔC_app(D), the implied θ̂ (holding the stored ΔOpt(D)), and z₀ / A7 / A8 across that range. Optionally the same with plan-conformant seeding for the optimism replicates. | Keeps the pre-specified estimate. Directly measures the missing uncertainty component. Cheap: apparent fits only, minutes (C4). | Needs authorisation. Adds a sensitivity table. Must be pre-declared so its results cannot replace θ̂. |

### 6.4 Recommendation: **4C**

**Basis.** Plan lines 6, 129–130, 139, 165–166 and 316; `38:145`.

**What would change.**
- **θ̂, the percentile interval and EP-1:** nothing.
- **The basic interval, BC, z₀, A7 and A8:** each gains a sensitivity band. The final report would state that θ̂'s reported MC-SE excludes the apparent-CV component, and give the measured size of that component.

---

## 7. Cross-decision dependencies

| Decision | Depends on | Why | Can it be adjudicated now, independently? |
|---|---|---|---|
| **1** CV seed | none (the root) | It decides which dataset all later decisions are applied to | **Yes.** Binding text, deviation status, remediation option and primary designation are all decidable before any computation. |
| **2** EP-1 | 1 (application) | EP-1 must be evaluated on whichever dataset is primary. A top-up must use the same driver as the data it extends: topping up as-run data propagates the deviation, and topping up with a conforming driver mixes conventions within a run, which should not be done. | **The rule (which estimator governs): yes, now**, ideally before any conforming data are seen. **Its outcome: no.** |
| **3** Interpretation | 1 (dataset), 4 (θ̂ enters basic, BC, z₀, A7, A8), 2 (an IMPRECISE label or a top-up would change the reported interval) | A7 and A8 are functions of θ̂ (§5.1); the percentile interval is not | **The reporting rule (3A): yes, now. The final scientific wording: no.** It must wait for 1, the C5 investigation and 4C. |
| **4** θ̂ | 1 (whether the replicates' convention is the as-run or the plan's changes which comparison is relevant) | Under 1D the replicates follow the plan's convention, and θ̂'s advancing-generator convention is then closer to them in its independence property | **Retaining θ̂ (4A): yes, now.** The sensitivity (4C) can run independently of 1. |

**Answers to the four questions in the request.**
- **Does the CV-seed decision affect the interpretation of θ̂?** Yes. Whether θ̂ and the replicates realise the same functional depends on whether the coupling matters, which only 1C or 1D can show.
- **Does the choice of θ̂ affect z₀ and the interval diagnostics?** Yes: z₀, A7, A8, basic and BC all depend on it. The primary percentile interval does not.
- **Does the EP-1 decision affect whether more outer replicates are necessary?** Yes. Under 2A, none are mandated by §2.4. Under 2B, 1 250 are. W5 is a separate, warning-level pointer either way.
- **Can any conclusion be finalised independently of the others?**
  - Only the *rules*: the binding text for Decision 1, the EP-1 estimator, the reporting format and θ̂ retention.
  - No *scientific* conclusion can be finalised before Decision 1's remediation is complete.

---

## 8. Recommended decision sequence

1. **Adjudicate the rules now, before any new computation**, so that no rule is chosen on the basis of new results:
   - (a) Decision 1: plan binding, deviation recorded (closes "NOT YET ADJUDICATED" in draft v2); remediation option (1D recommended); primary designation fixed in writing.
   - (b) Decision 1D prerequisite: the apparent-level CV seed convention, as a recorded plan clarification.
   - (c) Decision 2: the EP-1 estimator rule (2A).
   - (d) Decision 3: the reporting rule (3A).
   - (e) Decision 4: θ̂ retained (4A), with the 4C sensitivity authorised or declined.
2. **If authorised, build and freeze the conforming driver** (a new file; the original stays untouched). Verify that the resample streams are unchanged *before* the full run, using the multiplicity hashes and bit-identical M1 values on a few replicates.
3. **Run the authorised computations** (§9) into new directories only.
4. **Apply the Decision 2 rule** to the primary dataset. A top-up follows only if that rule says so.
5. **Run C4 and C5** and complete the §7.1 investigation.
6. **Write the final interpretation** under the Decision 3 rule. Remove the PROVISIONAL labels only for results whose dependencies are all resolved.

---

## 9. Additional computation required, if any

**None of these has been performed. Each needs explicit authorisation.** Cost estimates use the observed throughput of ≈ 2.13 outer replicates/min on 8 cores (F23), not the plan's original budget, which proved too low.

| ID | Computation | For | Needed when | Writes to | Approximate cost |
|---|---|---|---|---|---|
| C1 | Paired CV-seed sensitivity on fixed indices (for example the 20 H5 indices of run 0), `spawn(2)[0]` inner CV | D1 | option 1C | a new directory | ≈ 10 min (20 replicates) |
| C2 | Full protocol-conforming rerun, 5 × 500 × 100, new frozen driver, plus an H5-equivalent check | D1 | option 1D | a new directory | ≈ 19.5 h |
| C3 | §2.4 top-up to B_out = 750, each run, with the driver that produced the data being extended | D2 | option 2B or 2C only | the run checkpoints (append), or new files; needs a decision either way | ≈ 9.8 h |
| C4 | θ̂ apparent-CV sensitivity: apparent M1 / M2 refits on D under a pre-declared seed set | D4 | option 4C | a new directory | minutes (apparent fits only) |
| C5 | A7 / A8 investigation from stored records: λ (D*ₖ vs D), ΔC_app(k) against λ and distinct-patient count, skew against median offset | D3 | before any quantitative claim | `final_project/` (new files) | seconds; no new replicates |
| C6 | W5 re-examination from the stored fine trace (step 10) | D2 | option 2A | `final_project/` (new files) | seconds; no new replicates |

---

## 10. Questions requiring project-owner authorisation

1. **Decision 1, binding text.** Is plan §8.2 (line 364) binding, so that the as-run seeding is a recorded protocol deviation? *(Recommended: yes.)*
2. **Decision 1, remediation.** Which option: 1B (disclose only), 1C (subset sensitivity, C1) or 1D (full conforming rerun, C2)? *(Recommended: 1D; minimum 1C.)*
3. **Decision 1, primary designation.** If 1D is chosen, do you confirm in writing, before the run, that the conforming results become primary whatever they show?
4. **Decision 1, apparent-level CV seed** (the plan is silent). Under 1D, should the apparent-fit CV keep the as-run convention or receive its own child stream, defined so that no resample stream changes?
5. **Decision 2.** Is EP-1 governed by the kernel quantile SE at Silverman ×1, with ×0.5, ×2 and the between-run SD supplementary (2A)? Or must every bandwidth pass, triggering the B_out = 750 top-up (2B)? *(Recommended: 2A, with disclosure that it was adjudicated after the grid was seen.)*
6. **Decision 2, W5.** Is W5's "raise B_out" a directive or a warning-level suggestion? *(Recommended: a warning; examine with C6.)*
7. **Decision 3.** Do you approve reporting all five per-run percentile intervals, with no single combined headline and the pooled interval descriptive only (3A), and the §5.4 wording? *(Recommended: yes.)*
8. **Decision 3, investigation.** Do you authorise C5 (stored records only) to complete the §7.1 investigation of the W6 / W7 breaches?
9. **Decision 4.** Do you confirm that θ̂ = +0.0043 is retained unchanged (4A), and authorise the C4 sensitivity (4C)? *(Recommended: yes to both.)*
10. **Scope.** Should Phase 9 primary and S-7 receive a code-level conformance audit against their own protocol's seeding text (F24)?

---

*This document changes nothing. No checkpoint, frozen file, plan, H5 artefact, existing result or completed report was modified, and nothing was committed, pushed or transmitted.*
