# Phase 9 — Methodological Clarification and Owner Adjudication Record

| | |
|---|---|
| **Project** | COAD transportability analysis, Phase 9: nested-bootstrap confidence interval for the optimism-corrected ΔC (`ΔC_corrected = C_corrected(M2 four-set) − C_corrected(M1)`) |
| **Purpose** | To record, **before any further computation**, explicit resolutions of the methodological ambiguities and deviations identified in `documentation/METHODOLOGICAL_DECISION_REVIEW.md` and `documentation/PRE_AUTHORIZATION_METHODOLOGY_AUDIT.md`, so that a protocol-correct rerun can be specified unambiguously |
| **Status** | **DRAFT — PENDING PROJECT-OWNER APPROVAL.** No decision in this record is approved. Every checkbox is unticked. |
| **Authoritative plan** | `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1, sha256/16 `d27c52a29930a338` (unchanged) |
| **Prepared** | 2026-11-04 |
| **Computation** | None was performed, and none is authorised by this record (§10). |

## Governing methodological principle

**The authoritative plan governs.**
- An ambiguity or silence in the plan is resolved only by an explicit, dated clarification recorded here, **before** any data produced under the clarified rule exist.
- A clarification is labelled as a clarification. It is never presented as something the original plan unambiguously required.
- No rule may be chosen, relaxed or replaced because of a result. "No threshold is relaxed because a result is inconvenient. B_out and B_in may only increase." (plan line 473)

Each decision below states: the **verified basis** (plan, code, stored evidence); the **proposed decision**, with its label (**plan requirement**, **clarification of ambiguous wording**, or **new methodological decision**); the **options**; and the **owner adjudication** boxes.

---

## Decision 1 — Authority of the CV-seed specification

**Verified basis.**
- **The plan:** §8.2, line 364, `cvseed = inner[b].spawn(2)  # [0] = inner 5-fold CV, [1] = reserve`.
- **The frozen driver** (`scripts/40_phase9_dc_interval.py`, sha256/16 `4b1738a4d0ee0e39`) seeds the inner CV with a fresh `np.random.default_rng(ss_inner[b])` (lines 141–142). That is the same initial state as the inner resample generator (lines 134–135). Equality was confirmed in 9 / 9 tested cases (`validation/cv_seed_state_evidence.json`).
- **No `spawn(2)` call** exists in the driver.
- **The driver's docstring** claims it "Implements docs/PHASE9_DELTAC_INTERVAL_PLAN.md v2.1 exactly" (line 3). **That claim is inaccurate for this CV-seed element.**
- **The full deviation record** is `documentation/CV_SEED_DEVIATION_DRAFT_v2.md` (DRAFT — NOT YET ADJUDICATED; sha256/16 `1e4200c100358161`). Adjudicating Decision 1 adjudicates that record. The draft itself is not edited.

**Proposed decision** (**plan requirement**):
- Plan §8.2 (line 364) is **binding** for inner cross-validation.
- The frozen driver's inner-CV seeding is recorded as a **methodological deviation** from that requirement.
- The existing 5 × 500 results are preserved unchanged as the **historical as-run analysis**. They are **not** relabelled as protocol-compliant.
- **H5 PASS establishes the reproducibility of the as-run procedure, not its conformance.**

**Owner adjudication:**
- [ ] Approve
- [ ] Reject
- [ ] Require further review

---

## Decision 2 — Apparent-fit CV random-number convention

**Verified basis.**
- **The plan specifies three streams:** the outer stream `outer[k]` (line 362), the inner stream `inner[b]` (line 363), and the inner-CV child (line 364).
- **It specifies no stream for the apparent-fit CV**, either on an outer resample D*ₖ or on the original dataset D. A search of the plan for `cvseed`, `spawn`, `SeedSequence`, `apparent` and `reserve` finds no other RNG text.
- **The frozen driver** seeds the apparent CV with a fresh `default_rng(ss_outer)` (lines 116–117). That is the same initial state as the outer resample (lines 97–100).
- **The plan anticipates growth of both indices.** B_out may grow through top-ups drawn "from the *same* `SeedSequence` spawn list at indices `B_out … B_out′−1`" (lines 127–128; ladder up to 3000, line 119). B_in "may only increase" (line 473).

**Why this is a new methodological decision.**
- The plan is silent at the apparent level, so **no stream is required by the plan**.
- The proposal extends, by analogy, the principle visible in line 364: the CV fold randomisation should be independent of the resampling draw that created the fitted sample. That principle supports "the identical functional" being applied to D*ₖ as to D (line 166).
- **It therefore requires owner approval** and is labelled a **new methodological decision**.

**Proposed construction** (**new methodological decision**):

| Element | Specification |
|---|---|
| Seed object | `numpy.random.SeedSequence(entropy=[BASE, run_id, T], spawn_key=(k,))`, with `BASE = 2026110300` (driver line 35) and `run_id ∈ {0, …, 4}` |
| Domain identifier `T` | A fixed, **non-zero** integer, **reserved exclusively** for this stream. Its exact value is frozen in the implementation specification **before execution**, and never chosen or changed after any rerun result is seen. |
| Generator | A **fresh** `numpy.random.default_rng(·)` on that SeedSequence **for each model** (M1, M2), mirroring the driver's per-model pattern. Only M2 consumes CV draws (`src/phase9_execution.py:196–201`). |
| Scope | Apparent-fit CV on each outer resample D*ₖ only. No other stream is changed: the outer resample, inner resample and inner-CV streams follow Decision 1 and plan lines 361–364. |

**Future collision protection (mandatory).**
- Separation must hold for the **full permitted future range**, not only for the replicates executed so far:
  - k ∈ [0, 3000), the top of the §2.4 ladder;
  - b ∈ [0, B_in,max), where **B_in,max must be declared** in the implementation specification, because the plan sets no upper bound on B_in (line 473).
- **A non-zero final entropy word alone is not accepted as proof of separation.** NumPy's handling of entropy words must itself be verified: no padding, truncation or other interpretation may make `[BASE, run_id, T]` produce the same pool as `[BASE, run_id]` or as any other project seed. That check is part of V2 (Decision 8).

**The original dataset D.**
- D is the single observed sample, **not an outer bootstrap replicate**. There is no `outer[k]` for D, and the bootstrap requires the same *procedure* on D*ₖ as on D, not the same random numbers.
- **Decision 3 retains the existing θ̂.** No primary-analysis refit on D occurs, and **no new seed convention for D is introduced.**
- The historical D convention (`default_rng(20261102)` per model, `scripts/38_phase9_execute.py:90–92`) is disclosed under Decision 3.

**Optional paired sensitivity comparison: as-run apparent-fit convention.**
- Within the same rerun, the apparent fits on each D*ₖ could additionally be computed under the as-run convention (`default_rng(SeedSequence([BASE, run_id], spawn_key=(k,)))`, i.e. driver lines 116–117).
- **Purpose:** to separate the effect of correcting the inner-CV deviation from the effect of the new apparent-fit choice.
- **It would be a sensitivity analysis only.** The independent-stream convention stays primary whatever it shows.
- If authorised, **V4 applies only to this comparison variant** (Decision 8). The independent-stream primary rerun is **not** required to reproduce the as-run ΔC_app.

**Owner adjudication — apparent-fit convention:**
- [ ] Approve independent apparent-fit stream (as specified above)
- [ ] Approve as-run convention
- [ ] Require further review

**Owner adjudication — optional paired sensitivity comparison:**
- [ ] Compute the as-run apparent-fit convention as a paired sensitivity analysis
- [ ] Do not compute this additional variant

**Approved T value:** ______________________  
**Declared B_in,max for the collision check:** ______________________

---

## Decision 3 — Treatment and disclosure of θ̂

**Verified basis.**
- **The plan's θ̂:** `Estimand (unchanged) … = +0.0043` (line 6); `theta_hat = DeltaC_corrected(D)  # +0.0043` (line 139).
- **The stored value:** `monte_carlo.M2_primary.corrected_C_mean − monte_carlo.M1.corrected_C_mean` from `reports/phase9/execution/phase9_execution_results.json` (sha256/16 `05d13b8b339fa7d4`), which is **0.004316342804030704** (`results/s1_interval_results.json#theta_hat`).
- **In `scripts/38_phase9_execute.py`:**
  - the corrected C for each seed is `apparent[name]["C"] - np.nanmean(o)` (line 145), so the apparent C on D is **the same value for all five seeds**;
  - the Monte-Carlo SE is the SD of the five seed-level values (line 175).

**Proposed decision** (**plan requirement**, with a disclosure added):
- Retain **θ̂ = +0.004316342804030704** as the pre-specified reference estimate for the primary analysis.
- Do **not** replace or re-centre it on the basis of observed bootstrap diagnostics.
- **Mandatory disclosure:** *The reported Monte Carlo standard error of approximately 0.00019 does not include the randomness of the cross-validation on the original dataset D that selected λ. The apparent ΔC on D was held fixed across all five seeds (`scripts/38_phase9_execute.py:145`).*
- **Further disclosure:**
  - θ̂ was produced under the Phase 9 primary seeding convention: CV draws are taken from one advancing generator after the resample (`38:117–129`), and the apparent CV on D is seeded `default_rng(20261102)` (`38:90–92`).
  - That apparent-CV generator starts from the same state as the generator of the first bootstrap resample of seed 20261102 (`38:91`, `38:118`). This cross-fit link has not been evaluated.
- **The percentile interval (the primary interval method) does not use θ̂.** The basic interval, the BC interval, z₀, A7 and A8 do.
- Any investigation of the full-data CV seed is labelled a **sensitivity analysis**, never a replacement estimate. It needs separate authorisation and is not authorised here.

**Owner adjudication:**
- [ ] Approve
- [ ] Reject
- [ ] Require further review

---

## Decision 4 — Formal EP-1 Monte-Carlo error estimator

**Verified basis.**
- **EP-1:** "`MC-SE(endpoint) ≤ 0.05 × (upper − lower)` for **both** endpoints" (lines 84–85).
- **The within-run endpoint MC-SE** is the kernel quantile SE (§2.1, line 95; §2.2, line 105).
- **The bandwidth:** "Bandwidth by **Silverman's rule**, with the estimate **recomputed at bandwidth × {0.5, 1.0, 2.0} and all three reported** — tail density estimates are bandwidth-sensitive and a single number would overstate precision" (lines 106–107).
- **§2.4** (lines 116–132): `B_out′ = 500 × (observed ratio / 0.05)²`, rounded **up** to the ladder {750, 1000, 1500, 2000, 3000}. "Each run is **topped up**". The interval is "never narrowed, re-centred or re-quantiled".
- **H3:** "EP-1 unmet after exhausting the B_out ladder … labelled **IMPRECISE**" (line 456).
- **What the plan does not state:** which of the three bandwidth-specific estimates governs acceptance (`PRE_AUTHORIZATION_METHODOLOGY_AUDIT.md` §4.2).

**Options.**
- **Option A:** the Silverman bandwidth h alone governs formal EP-1 acceptance. The 0.5h and 2h estimates are reported as sensitivity analyses.
- **Option B:** the **maximum** of the three bandwidth-specific estimates governs formal EP-1 acceptance.

**Recommendation: Option B,** as a **prospective conservative clarification**, subject to owner approval.
- It stops a bandwidth-sensitive precision estimate from passing merely because the nominal bandwidth passes. That is consistent with the plan's stated reason for the three bandwidths (lines 106–107).
- **Option B is not described as an unambiguous requirement of the original plan.** It is a clarification of ambiguous wording.

**Clarification text** (applies under Option B; items 1–4 and 6 also apply under Option A):
1. **Bandwidth rule:** Silverman (1986), eq. 3.31, h = 0.9 · min(s, IQR / 1.34) · n^(−1/5), with s computed with ddof = 1 and the IQR from type-7 quantiles, as in Statistical Convention 7 (`documentation/STATISTICAL_CONVENTIONS.md`). **The SciPy `silverman` variant (≈ 1.06 · s · n^(−1/5)) is not substituted.** It may appear only as a labelled cross-check.
2. **Estimates:** the kernel quantile SE, `√(p(1−p)/B_out) / f̂(q_p)`, with f̂ a Gaussian kernel density evaluated at the type-7 endpoint estimate. It is computed at bandwidths **0.5h, h and 2h**, with h recomputed on the run's full current B_out replicates.
3. **Interval:** the 95 % percentile interval of each run (p = 0.025 and 0.975), with width W = upper − lower at the run's full current B_out.
4. **Scope:** EP-1 is applied **separately in each run and at each endpoint**.
5. **Governing estimate (Option B):** for each run r and endpoint e, the ratio ρ(r, e) = max over the three bandwidths of [SE(r, e, bandwidth) / W(r)].
6. **The between-run SD** is a separate Monte-Carlo variability measure (§2.1, row 3). It is reported, and compared with the kernel estimate in W4. **It does not enter EP-1 acceptance.**
7. **The §2.4 "observed ratio"** is ρ_max = max over all runs r and endpoints e of ρ(r, e). That is the maximum of the three bandwidth-specific endpoint MC-SE-to-width ratios, taken over every run and endpoint.
8. **Acceptance threshold:** EP-1 **passes** at the current B_out if and only if ρ(r, e) ≤ 0.05 for **every** r ∈ {0, …, 4} and e ∈ {lower, upper}. **The 0.05 threshold is unchanged.**
9. **Sizing and stopping rule:**
   - If EP-1 fails at the current size B_cur, compute B_req = B_cur × (ρ_max / 0.05)². At B_cur = 500 this is exactly the plan's formula.
   - The new size is the smallest ladder rung ≥ B_req **and** > B_cur.
   - **All five runs are topped up to that same rung**, so that the runs remain identical in design apart from the seed (lines 242–243). The new outer replicates have indices B_cur … B_new − 1, drawn from the same spawn list (lines 127–128).
   - The **same rule is then re-evaluated in full at the new size**: h recomputed, all three bandwidths, every run and endpoint.
   - **Stop** at the first size at which EP-1 passes.
   - If EP-1 still fails at B_out = 3000, stop: the interval is reported with its measured endpoint MC-SE and labelled **IMPRECISE** (H3, line 456; lines 130–132).
   - B_out only increases. The interval is never narrowed, re-centred or re-quantiled.

**Historical note (disclosed).** Option B would **fail** the historical as-run data at B_out = 500.
- **Evidence:** run 1, upper endpoint, 0.5h: SE 0.002993 against W 0.055924, a ratio of **0.0535 > 0.05** (`results/s1_interval_results.json#EP1`). All 29 other kernel evaluations (of 30) pass.
- This was known when the option was proposed. It does not prevent prospective adoption for the corrected rerun, whose data do not yet exist.
- The historical analysis is not topped up (Decision 6).

**Owner adjudication:**
- [ ] Option A — Silverman's bandwidth alone
- [ ] Option B — Maximum across all three bandwidths
- [ ] Require further review

---

## Decision 5 — W5 definition and required action

**Verified basis.**
- **W5 row:** "MC-1 endpoint running trace over the final 40 % of B_out · spread > EP-1 bound · **raise B_out**" (line 469).
- **Section heading:** "Warning conditions — investigate and report; never auto-validate or auto-invalidate" (line 461). The column holding "raise B_out" is headed "Justification of the threshold" (line 463).
- **MC-1:** the "running-estimate trace … B = 50, 100, 200, 300, 400, 500, …"; "visibly stable over the final ~40 % of replicates" (`docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md:223`).
- **Current implementation:** spread = max − min over B ∈ {300, 400, 500} (`analysis/s1_interval_analysis.py:81–82`; `analysis/stats_core.py:141–142`). The bound is 0.05 × the full-B_out width (`s1_interval_analysis.py:138–139`).
- **The plan does not define:** "spread"; the grid for B_out ≠ 500; whether a breach mandates an increase; or how large an increase must be (`PRE_AUTHORIZATION_METHODOLOGY_AUDIT.md` §4.2, Q6–Q10).

**Proposed W5 definition** (**clarification**):

> For each run and each endpoint, **W5 spread** is the maximum minus minimum of the type-7 percentile endpoint, computed on the first B replicates in k order, for B in the trace grid intersected with [0.6 × B_out, B_out].
>
> The trace grid is `{50, 100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000}`.
>
> A breach occurs when `W5 spread > 0.05 × full-B_out percentile interval width`.

**Consistency with the current implementation and the 750 case.** The interval [0.6 × B_out, B_out] is closed.

| B_out | [0.6·B_out, B_out] | Grid points used |
|---|---|---|
| 500 | [300, 500] | **{300, 400, 500}**, identical to the current implementation |
| 750 | [450, 750] | **{500, 750}** |
| 1000 | [600, 1000] | {750, 1000} |
| 1500 | [900, 1500] | {1000, 1500} |
| 2000 | [1200, 2000] | {1500, 2000} |
| 3000 | [1800, 3000] | {2000, 3000} |

Above 500 the definition uses only two trace points, so W5 is a coarser stability check at those sizes. This is disclosed, not corrected: changing the grid would be a further decision.

**The interpretive point, stated plainly.**
- **Literal text:** W5's row says "raise B_out" without conditions (line 469).
- **Option W5-B departs from it.** Treating W5 as an investigation and reporting trigger, with EP-1 alone determining the mandatory top-up, is a **prospective clarification of that wording**. It is **not** a literal restatement of the plan.

**Options.**
- **Option W5-A: W5 independently mandates a top-up.**
  - Proposed corrective-action rule, which the owner must confirm or amend:
    - a W5 breach at size B_cur triggers a top-up of all five runs to the next ladder rung above B_cur;
    - W5 and EP-1 (Decision 4) are both re-evaluated at the new size;
    - the target size is the **larger** of the W5 target and the EP-1 target, whenever both trigger;
    - **stop** when EP-1 passes and W5 shows no breach;
    - if W5 still breaches at 3000, stop and report the persisting W5 breach as a warning. No plan label applies, and IMPRECISE remains reserved for EP-1 under H3.
- **Option W5-B: W5 triggers investigation and reporting; EP-1 governs mandatory top-up sizing.**
  - A W5 breach is investigated, using the running trace and the endpoint involved, and reported prominently. It does not by itself change B_out.
  - If a top-up occurs because EP-1 fails, W5 is re-evaluated and reported at the new size.

**Consistency with Decision 4.**
- **Under either option,** EP-1's rule is applied unchanged.
- **Under W5-A,** W5 can only *add* a reason to increase B_out. It can never reduce or override an EP-1-mandated size.
- **Under W5-B,** only EP-1 changes B_out.
- **Neither option relaxes a threshold.**

**Proposed recommendation:** Option W5-B, recorded explicitly as a clarification that departs from the literal wording. **The choice is the owner's.**

**Owner adjudication:**
- [ ] W5 independently mandates a top-up under a clarified sizing rule (Option W5-A)
- [ ] W5 triggers investigation and reporting; EP-1 determines mandatory top-up (Option W5-B)
- [ ] Require further review

**Approved W5 definition:** [ ] as proposed above · [ ] amended: ______________________  
**Approved W5 corrective-action rule:** ______________________

---

## Decision 6 — Existing checkpoints

**Verified basis.**
- **The checkpoints:** `reports/phase9/dc_interval/run{0..4}_outer.jsonl`, sha256/16 `885bbd1bd891b570`, `1bb27c56612a5f9d`, `0f12dce12bfce95e`, `4115bd4cd014d345`, `2202d559b0418392`. Each holds 500 records, and all five were verified unchanged by final QC (`audit/final_qc.json`).
- **The driver appends** to these files on any rerun (driver lines 179–195).

**Proposed decision:**
- The five existing 500-replicate checkpoints remain **frozen and unchanged**.
- **No top-up or append operation is authorised.**
- Any corrected rerun writes to a **separate output directory**. Its records are never merged or pooled with the historical checkpoints.
- Any increase in B_out applies only within the corrected rerun specification (Decisions 4 and 5).

**Owner adjudication:**
- [ ] Approve
- [ ] Reject
- [ ] Require further review

---

## Decision 7 — Corrected rerun becomes primary

**Proposed decision** (declared before execution):
- The protocol-correct rerun becomes the **primary analysis**, regardless of the direction, magnitude, statistical significance or favourability of its results.
- The existing analysis remains the **historical record** of the as-run implementation and its documented deviation.
- The corrected rerun **must not** be selectively accepted, rejected or replaced on the basis of whether its findings agree with the historical results.

**Failure-handling rule:**
- If the corrected rerun fails a verification gate (V1–V5, Decision 8) or a hard-failure criterion (H1–H6, plan lines 454–459), the implementation is **investigated and corrected**, and the rerun is repeated under an **explicitly versioned and documented procedure**: a new driver version, a new frozen hash and a dated record of the defect and its correction.
- **A failed gate or a hard failure never reinstates the historical as-run analysis as primary.**
- A correction may fix a defect in the implementation. It may not change any rule fixed in this record unless that change is itself adjudicated, dated and recorded before the repeated run.

**Owner adjudication:**
- [ ] Approve
- [ ] Reject
- [ ] Require further review

---

## Decision 8 — Verification gate and output location

**Canonical V1–V5** (verbatim from `PRE_AUTHORIZATION_METHODOLOGY_AUDIT.md` §3.4):
- **V1:** the outer and inner multiplicities reproduce the stored `multiplicity_sha` values.
- **V2:** every new CV state differs from its resample state and from every other stream in the tree, checked on the H5 indices.
- **V3:** all M1 quantities are bit-identical to the as-run records.
- **V4:** if A-R is chosen, ΔC_app(k) is bit-identical to the as-run records.
- **V5:** an H5-equivalent bit-identical recomputation of the new driver's own output.

**Explanatory gate summary**, as this record requires it to be operationalised:

| Check | Requirement | Notes |
|---|---|---|
| **V1** | Identical outer and inner resamples wherever paired comparison with the historical records requires them | **Operational gap, recorded and not resolved here.** The historical checkpoints store `multiplicity_sha` for the **outer** resample only (driver lines 103–107). No inner-resample hash is stored. The implementation specification must therefore state exactly how inner-resample identity is established, for example by regenerating the as-run inner indices from the unchanged `ss_inner[b]` streams with the frozen `strat_boot` and comparing them index for index. It must do so without weakening V1. |
| **V2** | Prescribed RNG stream assignment for every stream (Decisions 1 and 2); distinct states; deterministic reproducibility; **no collision with existing or future seed domains over the full permitted range** (k < 3000, b < declared B_in,max; all run_id; the Phase 9 primary seeds 20261102–20261106); verification of `T`, including that NumPy's entropy handling does not pad, truncate or otherwise interpret `[BASE, run_id, T]` so as to coincide with another domain | This **extends** the canonical V2 wording ("checked on the H5 indices") to the full future range, as Decision 2 requires. It strengthens the check and weakens nothing. |
| **V3** | M1 values identical where the methodology predicts invariance: the apparent C, the apparent λ (= 0), `opt_M1_mean`, and the M1 inner C values | Stored inner C values are rounded to 8 decimals (driver line 167). The comparison for those is therefore at the stored precision, and is stated as such. |
| **V4** | Applies **only** to the optional as-run apparent-fit comparison variant, if authorised (Decision 2): its ΔC_app(k) must be bit-identical to the historical records | **Never applied to the independent-stream primary rerun.** Its ΔC_app is expected to differ, by design. |
| **V5** | A reproducibility spot-check of pre-specified replicates: k = 0, 25, 50, …, 475 of run 0 (the plan §8.4 rule, lines 383–385), recomputed and bit-identical | — |
| **H4** (separate; not part of V1–V5) | Frozen-input integrity: every frozen-artefact hash is unchanged before and after execution (plan line 457; driver `check_frozen`, lines 38–58) | A hard failure in its own right, evaluated separately |

**Rules.**
- The exact canonical V1–V5 wording and thresholds, together with the operational specifications above, **must be copied into the implementation specification before execution**.
- Each check reports an explicit **PASS / FAIL**.
- **No criterion may be weakened, re-scoped or reinterpreted after any result is seen.**

**Approved output directory:** ______________________ *(owner to confirm; it must be new and separate from `reports/phase9/dc_interval/`)*

**Owner adjudication:**
- [ ] Approve V1–V5 after exact specification
- [ ] Require changes
- [ ] Require further review

---

## Decision 9 — A7 / A8 interpretation and interval reporting

**Verified basis.**
- **Percentile is primary** (plan line 309). The basic interval is diagnostic, with "a warning that **neither** interval should be trusted" on material divergence (lines 312–314). z₀ is defined at line 316.
- **A7 and A8 are warnings:** "Neither a breach nor a pass establishes that the interval is valid or invalid" (lines 330–331); warnings "never auto-validate or auto-invalidate" (line 461).
- **The A7 algebra.** With the basic interval reflecting the same two quantiles, A7 = 2·|m − θ̂| / W at **both** endpoints, where m is the percentile midpoint and W its width. That is why the two endpoint values are always identical. A7 splits into a tail-asymmetry term that does not depend on θ̂ and a median-offset term.
- **A8 is different.** A8 = |Φ⁻¹(#{ΔC_corr(k) < θ̂} / B)| depends only on θ̂'s rank among the replicates.
- **Hence they are neither equivalent nor proportional**, except under exact normality (`PRE_AUTHORIZATION_METHODOLOGY_AUDIT.md` §5).

**Proposed decision** (this resolves the earlier contradiction):
- **The five per-run percentile intervals are the primary interval results.**
- **No combination rule is needed, because no combined interval is constructed or designated as primary.**
- **Required reporting format:**
  1. The five individual run-level 95 % percentile intervals.
  2. Between-run variability: the SD (ddof = 1, df = 4) and the range of each endpoint (plan §2.1, quantities 3 and 5).
  3. **A7 and A8 reported separately**:
     - A7 as one statistic, because its two endpoint values are algebraically identical;
     - A8 as |z₀|.
  4. The interval asymmetry and displacement observations: A7's tail-asymmetry and median-offset components, and the shift of the outer distribution relative to θ̂.
  5. The pooled 2 500-draw interval, if included, **only as a descriptive summary**, labelled as not being a primary or separately established interval.
  6. An explicit statement that **intervals containing zero do not establish equivalence and do not prove the absence of an effect**.
- **Language:**
  - A7 and A8 are **not** described as mathematically equivalent.
  - Both are **warning diagnostics requiring investigation and qualification**. Their activation alone does not automatically invalidate an interval under the plan's warning framework. The labels UNRELIABLE and IMPRECISE remain reserved for H1–H3 (lines 454–456).
  - The mandatory conditional statement (lines 206–208) accompanies every reported interval.

**Owner adjudication:**
- [ ] Approve
- [ ] Reject
- [ ] Require further review

---

## 10. Authorisation boundary

**Approving any or all of the nine decisions does not authorise computation.** A separate, explicit rerun authorisation is required. It may be issued only after all of the following:

- every methodological choice in Decisions 1–9 has been adjudicated and dated;
- the RNG specification is documented, including the approved `T` and the declared B_in,max;
- the EP-1 and W5 rules are unambiguous and approved;
- a new driver has been written as a **new file**, reviewed and hash-frozen (the original driver `4b1738a4d0ee0e39` stays untouched);
- V1–V5 (including the V1 operational specification) and H4 are specified in the implementation specification;
- the output directory is confirmed;
- the runtime and resource commitment (§11) is accepted.

**Until then:**
- no bootstrap is started;
- no B_out top-up is performed;
- no frozen checkpoint is modified;
- no historical result is overwritten;
- no primary scientific conclusion is finalised.

---

## 11. Runtime and resource commitment

**These are planning estimates based on observed throughput, not guarantees.** The basis is ≈ 2.13 outer replicates per minute on 8 cores (`reports/phase9/dc_interval/HALT_RECORD.json`). The as-run runs 1–4 took 199–264 min each (`resume_driver.log`). The plan's original 7.4 h budget proved too low (HALT_RECORD diagnosis).

| Stage | New outer replicates | Estimated wall time (8 cores) |
|---|---|---|
| Corrected rerun, B_out = 500 (5 runs) | 2 500 | **≈ 19.5 h** |
| If EP-1 requires the next ladder rung, B_out = 750 | + 1 250 | **≈ 9.8 additional hours** |
| *For information:* further rungs, if required | 1000: + 1 250 (≈ 9.8 h); 1500: + 2 500 (≈ 19.6 h); 2000: + 2 500 (≈ 19.6 h); 3000: + 5 000 (≈ 39.1 h) | increments sum to ≈ 97.9 h; total including the base run ≈ 117 h |
| Optional as-run apparent-fit variant (Decision 2) | 2 extra apparent fits per outer replicate | small next to the 200 inner fits per replicate; not separately estimated |
| Gate checks V1–V5, H4 | — | not separately estimated |

**Commitment required.**
- **Before authorising execution, the owner must accept** the base estimate and the possibility of at least the additional ≈ 9.8 h for B_out = 750.
- If the owner adopts Option B (Decision 4) or Option W5-A (Decision 5), they should also accept that further rungs may follow under the approved rule.
- **No claim is made** about whether the corrected rerun will trigger any top-up, or about any change in its results.

---

## 12. Final owner approval

**Project-owner decision:**
- [ ] Approve all nine decisions as written.
- [ ] Approve with amendments.
- [ ] Defer authorisation pending further review.

**Approved amendments:**  
______________________________________________________________________

**Date of adjudication:** ______________________

**Approved T value:** ______________________

**Approved output directory:** ______________________

**RERUN AUTHORIZATION: NOT GRANTED BY THIS RECORD.**

**FINAL STATUS: PENDING PROJECT-OWNER ADJUDICATION.**

---

*This record is a draft for adjudication. It changes no plan, driver, historical implementation, checkpoint, result, report or audit record, and it authorises no computation.*
