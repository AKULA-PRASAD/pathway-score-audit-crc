# Pre-authorisation methodology audit: CV-seed specification, EP-1 and W5, and A7 / A8

| | |
|---|---|
| Status | **INVESTIGATION AND RECOMMENDATIONS ONLY. NOTHING HERE IS ADJUDICATED, AMENDED OR EXECUTED.** |
| Prepared | 2026-11-04, for the project owner, before any rerun authorisation |
| Relation to earlier work | Follows `documentation/METHODOLOGICAL_DECISION_REVIEW.md`, which is left unchanged. Corrections to it are listed in §1.2 and are made **here only**. |
| Computation | **None.** No bootstrap, sensitivity analysis or top-up was run, and no RNG was instantiated. The only arithmetic combines values already stored in `results/s1_interval_results.json` (percentile endpoints, medians, θ̂, z₀). Every such value is marked *derived*. |

Each issue separates **Verified source fact**, **Interpretation**, **Unresolved ambiguity**, **Recommendation**, and **Consequence of each option**. Line numbers refer to the files listed in §2.

---

## 1. Executive conclusion

### 1.1 Conclusions

**Issue A: CV-seed specification.**
- The plan explicitly prescribes the random-number stream for the **inner** fits only: resample from `inner[b]`, CV from `inner[b].spawn(2)[0]` (plan lines 363–364).
- It is **silent** on the apparent-fit CV on each outer resample D*ₖ.
- It has **no RNG requirement at all** for the original-data fit on D. θ̂ is taken as the fixed, already-computed +0.0043 (lines 6, 139) and is never recomputed by this plan.
- The driver is internally consistent: at both levels, CV restarts from the state of the resample stream. It conforms to the plan at the outer and inner resample level and deviates at the inner CV.
- **Giving the apparent CV its own stream is a new methodological decision, not compliance with a plan requirement.** There is a principled argument for it (§3.2), and a rerun cannot proceed until the owner has chosen. A precise, collision-free specification is proposed in §3.4.

**Issue B: EP-1 versus W5.**
- **EP-1 names an estimator family, not a single estimate.** It names "MC-SE(endpoint)", and §2.1 defines the within-run endpoint MC-SE as the kernel quantile SE. §2.2 makes Silverman's rule the bandwidth rule but requires ×0.5 and ×2 to be "recomputed … and all three reported", because "a single number would overstate precision" (lines 106–107).
- **Which of the three enters EP-1 is not stated.** The rationale clause can be read for either answer.
- **The two rules are of different kinds.** §2.4 is a mandatory corrective rule, triggered only when EP-1 fails. W5 is a warning row whose last cell reads "raise B_out" (line 469), inside a section headed "investigate and report; never auto-validate or auto-invalidate" (line 461). No sizing rule is attached to it.
- **The plan does not give an unambiguous decision rule for this situation.** The earlier review's recommendation (2A) is **withdrawn**: the protocol should be **clarified before any further computation**, and the as-run data should not be topped up.

**Issue C: A7 / A8.**
- With the basic interval built by reflecting the same two quantiles through θ̂, the percentile-minus-basic difference is **exactly equal at both endpoints**. So A7 = 2·|m − θ̂| / W, where m is the percentile midpoint and W its width. This is an algebraic identity, true for every dataset under the current construction (§5.1).
- **A8 is not equivalent and not proportional to A7.**
  - A8 is a rank (probability-scale) measure of where θ̂ falls relative to the outer distribution's median.
  - A7 is a data-scale measure that adds the 2.5 / 97.5 % quantile asymmetry.
  - Only under an exactly normal outer distribution does A7 = |z₀| / 1.96 hold.
- **In these data, the tail asymmetry term *alone* exceeds the 25 % threshold in runs 0, 1, 3 and 4** (*derived*, §5.3). A7 would breach there even if θ̂ were exactly the median.
- **The earlier review's statement "A7 and A8 measure the same thing here" is incorrect** and is corrected below.
- By the plan's rules, both remain warnings that require investigation and qualification. They do not invalidate the interval (lines 330–331, 461).

### 1.2 Corrections to `METHODOLOGICAL_DECISION_REVIEW.md`

| # | Earlier statement (location) | Correction |
|---|---|---|
| R1 | "A7 and A8 measure the same thing: how far the outer distribution sits from θ̂" (§1 table, row 3); "A7 and A8 are not two independent confirmations" (§5.1) | **Incorrect as stated.** They share a θ̂-dependent median-offset component, so they are *not independent*. But they are neither equivalent nor proportional: A7 also carries a tail-asymmetry component that does not depend on θ̂. That component is the larger share of A7 in runs 0, 1, 3 and 4, and alone exceeds 25 % in those four runs (§5.3). |
| R2 | Bandwidth evaluations "Supplementary, on the plain wording"; recommendation 2A (×1 governs, so no top-up) (§4.2–4.5) | **Overstated.** The plain wording requires the three bandwidths to be *reported*. It does not say which enters EP-1. The stated rationale ("a single number would overstate precision") is at least as compatible with requiring robustness across bandwidths. The plan is **genuinely ambiguous**, and 2A is withdrawn in favour of a pre-rerun clarification (§4.6). |
| R3 | The W5 breach "depends on the three-point grid of Convention 12"; suggestion to re-examine it from the stored fine trace (§4.3, C6) | **Misleading.** The grid {300, 400, 500} is the project's MC-1 grid (Phase 6B protocol line 223; `38:157–158`), not an arbitrary choice. The implemented "spread" is a range (`s1_interval_analysis.py:81–82`), and the range over a superset of trace points can only be ≥. **A finer grid cannot remove this breach.** What remains ambiguous is the definition of "spread", not the grid's coarseness. |
| R4 | Decision-1 condition 1: derive the apparent-level seed "from an explicit `spawn_key`, not from the order of `spawn()` calls" (§3.5) | **Incomplete.** Any key inside the root's child namespace or inside `outer[k]`'s child namespace can collide with a future replicate. The plan allows B_out to grow (lines 127–128) and B_in to grow (line 473). A disjoint namespace is required (§3.4). |
| R5 | 4C would show "whether the W6 / W7 breaches are robust to θ̂'s unmeasured CV component" (§6.1) | **Only partly true.** A θ̂ sensitivity bears on A8 and on A7's median-offset term. It cannot remove the A7 breach in the four runs where the θ̂-independent asymmetry term alone exceeds 25 %. |

### 1.3 Readiness

**Not yet ready for rerun authorisation.** The following must be decided and recorded first (§7, §10):
- the apparent-fit CV convention;
- the EP-1 governing-estimator rule, and whether W5 obliges an increase in B_out (and by how much);
- the written pre-designation of primary results;
- the verification gate.

None of these needs computation.

---

## 2. Source inventory and exact versions inspected

All read-only. sha256 prefixes (16 hex) were computed on 2026-11-04.

| File | sha256/16 | Inspected |
|---|---|---|
| `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` (v2.1, authoritative) | `d27c52a29930a338` | all 534 lines; cited lines re-read for this audit |
| `scripts/40_phase9_dc_interval.py` (frozen driver) | `4b1738a4d0ee0e39` | lines 1–206 |
| `src/phase9_execution.py` (frozen shared code) | `571d0c233f367be6` | lines 119–206 |
| `scripts/38_phase9_execute.py` (Phase 9 primary) | `82afbcf05f3b7514` | lines 1–200 |
| `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` (MC-1 definition) | `755e44d3c00d4a59` | lines 223, 228 |
| `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` | `ab56708919f3746f` | lines 543, 694 (keyword search) |
| `reports/phase9/execution/phase9_execution_results.json` | `05d13b8b339fa7d4` | `apparent` block |
| `reports/phase9/final_project/analysis/stats_core.py` | `26d60618fe05386c` | `quantile`, `percentile_ci`, `basic_ci`, `z0`, `running_trace`, `MC1_GRID` |
| `reports/phase9/final_project/analysis/s1_interval_analysis.py` | `abb0cacc40fe9f8a` | lines 57–58, 81–88, 120–142 |
| `reports/phase9/final_project/results/s1_interval_results.json` | `99ad5e20c1485425` | per-run endpoints, medians, z₀, A7, W5, kernel SE |
| `reports/phase9/final_project/validation/cv_seed_state_evidence.json` | `976ed0dc3a97672d` | 9 / 9 state-equality result |
| `reports/phase9/final_project/documentation/CV_SEED_DEVIATION_DRAFT_v2.md` | `1e4200c100358161` | whole |
| `reports/phase9/final_project/documentation/STATISTICAL_CONVENTIONS.md` | `fe0845dbce0c563c` | #10, #12 |
| `reports/phase9/final_project/documentation/METHODOLOGICAL_DECISION_REVIEW.md` | `ade6627a786c113f` | whole (not modified) |

---

## 3. Issue A: CV-seed specification

### 3.1 Verified source facts

**Plan §8.2 (lines 361–364), verbatim:**
```
root_k = numpy.random.SeedSequence(entropy=[BASE, run_id])
outer  = root_k.spawn(B_out)                 # outer[k]  -> outer replicate k
inner  = outer[k].spawn(B_in)                # inner[b]  -> inner replicate b of outer k
cvseed = inner[b].spawn(2)                   # [0] = inner 5-fold CV, [1] = reserve
```

**Other relevant plan text:**
- Line 367: "Each replicate's stream is a pure function of `(BASE, run_id, k, b)`."
- Lines 369–372: independence between runs and between layers, scheduling invariance, resumability.
- Line 395: each worker computes "apparent fit + 100 inner replicates, both models".
- Line 139: `theta_hat = DeltaC_corrected(D)  # +0.0043`.
- Line 6: "Estimand (unchanged) … = **+0.0043**".
- Lines 127–128: the top-up draws "from the *same* `SeedSequence` spawn list at indices `B_out … B_out′−1`".
- Line 473: "B_out and B_in may only increase."

A search of the plan for `cvseed`, `spawn`, `SeedSequence`, `apparent` and `reserve` finds **no other RNG text**: line 303 mentions the apparent fit only as a reliability criterion, and lines 170 and 422 only as a stored quantity.

**Implementations (code):**

| Fit | Frozen driver | Phase 9 primary |
|---|---|---|
| Outer resample on D | `rng_o = default_rng(ss_outer)`, `ss_outer = root.spawn(k+1)[k]` (lines 96–100) | n/a |
| Apparent CV on D*ₖ | fresh `default_rng(ss_outer)` per model (lines 116–117): **same initial state as the outer resample** | n/a |
| Inner resample | `default_rng(ss_inner[b])`, `ss_inner = ss_outer.spawn(B_in)` (lines 130, 134–135) | `boot_indices(e, rng)` from one advancing `rng = default_rng(seed)` (`38:117–122`) |
| Inner CV | fresh `default_rng(ss_inner[b])` per model (lines 141–142): **same initial state as the inner resample** | the **same advancing** `rng`, after the resample and after earlier models' CV (`38:126–129`) |
| Apparent CV on D | **not computed** | fresh `default_rng(SEEDS[0]=20261102)` per model (`38:90–92`) |

**Further facts:**
- **CV consumption:** only `tune_lambda` consumes the generator. It is called only when `p > n_clin` (`src/phase9_execution.py:196–201`), so only M2 consumes CV draws. The protocol confirms that "M1 … has no λ to tune" (`PHASE9_MODELING_PROTOCOL_DRAFT.md:543`).
- **State equality:** the executed CV state equals the resample state in 9 / 9 tested cases; planned `spawn(2)[0]` equals it in 0 / 9 (`cv_seed_state_evidence.json`).
- **An additional observation from the code.** In Phase 9 primary, the generator for the apparent CV on D (`default_rng(20261102)`, `38:91`) and the generator for the first bootstrap resample of seed 20261102 (`38:118`) start from the **same state**. The fold permutation on D and that first resample both read the event-index set of D.
  - This links the apparent fit on D to one of the 2 500 optimism replicates. It is a dependence *across* fits, not within one.
  - Its effect has not been evaluated. It is recorded for completeness only.

### 3.2 Interpretation

**Answers to the eight questions.**

1. **Inner fits.**
   - The plan explicitly prescribes `inner[b].spawn(2)[0]` for the inner 5-fold CV (line 364).
   - It implies, without stating it, that the inner resample is drawn from `inner[b]` itself ("inner[b] → inner replicate b", line 363). The driver's lines 134–135 implement that reading.
2. **Apparent fit on D*ₖ: not specified.**
   - Line 362 maps `outer[k]` to "outer replicate k". One could read that as covering all randomness of outer replicate k, apparent CV included, which is what the driver does.
   - But read the same way, line 363 would cover the inner CV, and the plan instead carved the inner CV out into its own child. Line 364 is evidence that the author intended CV to be separate from resampling, but only at the inner level.
   - **Neither reading is explicit.**
3. **Apparent fit on D: not specified, and not in this plan's scope.** The plan treats θ̂ as an existing number (lines 6, 139) and never instructs refitting on D.
4. **Is the driver consistent?** Yes, internally. At both levels it seeds the CV with a fresh generator on the SeedSequence that seeded that layer's resample. It is consistent with the plan at the apparent level (where the plan is silent) and inconsistent at the inner level (line 364).
5. **The scientifically defensible convention.**
   - **Principle.** The bootstrap rests on applying "the identical functional" to D*ₖ as to D (line 166). Fitting M2 includes a **randomised** step, the CV fold assignment. On D, that randomisation is independent of how D was sampled: nature drew the patients, and the analyst drew the folds.
   - **Consequence.** To mimic this, the fold randomisation inside each bootstrap fit should be independent of the bootstrap draw that created that fit's sample. That is what line 364 secures at the inner level. Applying the same principle to the apparent fit on D*ₖ calls for an apparent-CV stream independent of the outer resample stream.
   - **Why the as-run convention fails this.** It makes the fold assignment a deterministic function of the bits that drew the resample. By analogy, that is as if the folds on D depended on how nature sampled the patients.
   - **On D itself,** no resample stream exists, so this coupling cannot arise. Any fixed, pre-declared seed independent of the data satisfies the principle. The historical `default_rng(20261102)` does, apart from the cross-fit link noted in §3.1.
   - **Why D is not an outer replicate.** D is the single observed sample, and θ̂ is one realisation of the estimator on it. The outer replicates are simulated samples. They need the same *procedure*, not the same *random numbers*. So D's CV seed need not, and cannot meaningfully, follow the outer-replicate convention: there is no `outer[k]` for D.
6. **Is a separate apparent-CV stream compliance or a new decision?** It is a **new methodological decision**. It extends line 364's principle by analogy, but no plan text requires it. Keeping the as-run apparent convention is equally compliant *textually*. It is weaker on the principle above.
7. **The exact convention to freeze:** see §3.4.
8. **What the apparent-CV choice affects (from the code):**
   - **Affected:**
     - apparent M2 λ;
     - C_app(M2, k), ΔC_app(k), ΔC_corr(k);
     - every interval, z₀, A7, A8, kernel SE, between-run SD, EP-1, W4–W7 and figure;
     - the apparent-stage failure accounting (H1 via §6.5).
   - **Not affected:**
     - inner optimism. Inner fits tune their own λ (driver lines 141–142) and do not use the apparent λ, so `opt_M1`, `opt_M2` and ΔOpt(k) are unchanged;
     - the outer and inner resamples;
     - all of M1;
     - θ̂.
   - **Effect on the as-run comparison.** If the rerun keeps the as-run apparent convention, ΔC_app(k) is bit-identical to the as-run value. Paired differences then isolate the **inner-CV deviation alone**, which is exactly the textual deviation. If the apparent convention changes, the paired differences mix both changes.

**Does the choice change the estimand?** No. The estimand ΔC_corrected (line 6) is defined at the level of the functional. The CV convention is part of the **finite-sample implementation of a randomised estimator**. It can change the estimator's realisations, and, if the randomisation is not independent of the resample, possibly the distribution the bootstrap targets. It does not change what is being estimated.

### 3.3 Unresolved ambiguity

- Whether line 362's "outer[k] → outer replicate k" was meant to cover the apparent CV. The text cannot settle this.
- Whether the as-run coupling has any distributional effect. This is unmeasured, and only a computation can answer it.

### 3.4 Proposed RNG specification (for owner approval; not executed)

Notation: `SS(entropy, spawn_key)` = `numpy.random.SeedSequence(entropy=…, spawn_key=…)`; `G(s)` = `numpy.random.default_rng(s)`.

| Element | Proposed specification | Label |
|---|---|---|
| **Outer resample** `I*ₖ` | `G(SS([BASE, run_id], (k,)))`. Identical to `root_k.spawn(B_out)[k]` and to the driver's `spawn(k+1)[k]` (equivalence verified, `cv_seed_state_evidence.json`). | **Explicitly specified** (line 362) |
| **Inner resample** `J*ₖ,b` | `G(SS([BASE, run_id], (k, b)))`, identical to `outer[k].spawn(B_in)[b]` | **Specified** for the stream (line 363). **Inferred** that the resample uses `inner[b]` directly (driver 134–135). |
| **Inner-fit CV** (M2 only) | `G(SS([BASE, run_id], (k, b, 0)))`, identical to `inner[b].spawn(2)[0]`. A fresh generator per model; M1 never consumes it. | **Explicitly specified** (line 364). "Fresh per model" is **inferred** from the driver pattern. |
| **Reserve** | `SS([BASE, run_id], (k, b, 1))`, unused | **Explicitly specified** as reserve (line 364) |
| **Apparent-fit CV on D*ₖ**: option **A-I** (independent) | `G(SS([BASE, run_id, T], (k,)))` with a fixed, **non-zero** domain tag, for example `T = 1`. A fresh generator per model. | **New methodological decision** |
| **Apparent-fit CV on D*ₖ**: option **A-R** (as run) | `G(SS([BASE, run_id], (k,)))`, the same state as the outer resample | **Inferred from the implementation.** Textually compliant through silence; fails the principle of §3.2 Q5. |
| **Apparent-fit CV on D** | **No new specification.** θ̂ is retained (lines 6, 139). Historical convention: `G(20261102)` per model (`38:90–92`). A specification would be needed only if θ̂ were re-derived or the earlier review's C4 sensitivity were authorised, and that would be a separate decision. | **Inferred from the implementation (historical).** Any new D convention is a **new decision**. |

**Why A-I uses a separate entropy domain instead of a spawn key.**
- Every node in the `[BASE, run_id]` tree can be reached by a future, permitted extension:
  - `(k,)` for k ≥ 500 by a B_out top-up (lines 127–128);
  - `(k, b)` and `(k, b, ·)` for b ≥ 100 by a B_in increase (line 473).
- A disjoint entropy (`[BASE, run_id, T]`) cannot collide under any B_out or B_in.
- The tag must be **non-zero**. I understand that NumPy's SeedSequence pads short entropy with zero words, which would make `[BASE, run_id, 0]` coincide with `[BASE, run_id]`. That is an implementation detail I have **not verified in this audit**. Gate V2 below checks state-distinctness directly.

**Verification gate, to run before any full rerun (none of it run now).**
- **V1:** the outer and inner multiplicities reproduce the stored `multiplicity_sha` values.
- **V2:** every new CV state differs from its resample state and from every other stream in the tree, checked on the H5 indices.
- **V3:** all M1 quantities are bit-identical to the as-run records.
- **V4:** if A-R is chosen, ΔC_app(k) is bit-identical to the as-run records.
- **V5:** an H5-equivalent bit-identical recomputation of the new driver's own output.

### 3.5 Recommendation

**Adopt A-I as primary.** It applies line 364's evident principle consistently to both fit levels, and it is the only option under which "the identical functional" (line 166) holds in distribution by construction. **Record it as a plan clarification**, explicitly labelled a new methodological decision.

**Also pre-declare a secondary A-R apparent computation inside the same rerun.**
- Cost: two extra apparent fits per outer replicate, small next to the 200 inner fits.
- This lets the inner-CV deviation and the apparent-CV choice be **separated** in the paired comparison with the historical records.
- This is a recommendation, not a requirement. The owner may decline it.

### 3.6 Consequence of each option

| Option | Consequence |
|---|---|
| **A-I primary** (with or without the A-R secondary) | Both CV levels independent of their resample. A new decision is recorded. Paired differences against the historical records combine both changes, unless the A-R secondary is also computed. |
| **A-R primary** | The fewest departures from the run as executed. Only the textual deviation is corrected. Paired differences isolate the inner-CV deviation exactly. The apparent-level coupling remains and is disclosed as such. |
| **Undecided** | A rerun cannot be specified. **BLOCKED.** |

---

## 4. Issue B: EP-1 and W5

### 4.1 Verified source facts (minimum wording)

| # | Text | Location |
|---|---|---|
| B1 | EP-1: "`MC-SE(endpoint) ≤ 0.05 × (upper − lower)` for **both** endpoints, with the **absolute** endpoint MC-SE reported **separately from sampling uncertainty**." | §2, lines 84–85 |
| B2 | Row 2, "**Within-run endpoint MC-SE** · noise in each quantile from finite B_out · kernel quantile SE (§2.2)". Row 3, "**Between-run endpoint variability** · SD of each endpoint across the five independent runs". "2, 3 and 5 are **Monte-Carlo**." | §2.1, lines 95–96, 100 |
| B3 | "`SE(q_p) = √(p(1−p)/B_out) / f̂(q_p)` … Bandwidth by **Silverman's rule**, with the estimate **recomputed at bandwidth × {0.5, 1.0, 2.0} and all three reported** — tail density estimates are bandwidth-sensitive and **a single number would overstate precision**." | §2.2, lines 105–107 |
| B4 | "**If EP-1 fails — the exact rule.** … `B_out′ = 500 × (observed ratio / 0.05)²`, rounded **up** to the next rung … {750, 1000, 1500, 2000, 3000}"; "Each run is **topped up**"; "B_out may only increase." | §2.4, lines 116–119, 127, 130 |
| B5 | H3: "EP-1 unmet after exhausting the B_out ladder … reported and labelled **IMPRECISE**" | §11.1, line 456 |
| B6 | §11.2 heading: "Warning conditions — **investigate and report; never auto-validate or auto-invalidate**". Column headers: "Condition · Threshold · **Justification of the threshold**". | lines 461, 463 |
| B7 | W5: "MC-1 endpoint running trace over the final 40 % of B_out · spread > EP-1 bound · **raise B_out**" | line 469 |
| B8 | W3: "Inner-noise SD inflation · > 5 % · **raise B_in**; never deconvolve the variance". After the table: "**B_out and B_in may only increase.**" | lines 467, 473 |
| B9 | Singular references to "the endpoint MC-SE" (§1.3, §6.4, §6.5); W4 calls both the kernel estimate and the between-run SD an "endpoint MC-SE". | lines 77, 290, 304, 468 |
| B10 | MC-1: "Running-estimate trace … B = 50, 100, 200, 300, 400, 500 … the estimate must be visibly stable over the final ~40 % of replicates" | Phase 6B protocol line 223 |
| B11 | Implementation: W5 "spread" is the **range** (max − min) over the trace at B ∈ {300, 400, 500} (`s1_interval_analysis.py:81–82`; `stats_core.MC1_FINAL`). Phase 9 primary's MC-1 trace uses the same grid (`38:157–158`). | code |

**Stored values, run 1, upper endpoint** (width 0.055924, EP-1 bound 0.002796):
- kernel SE ×0.5: 0.002993 (ratio 0.0535);
- kernel SE ×1: 0.002663 (0.0476);
- kernel SE ×2: 0.001969 (0.0352);
- between-run SD: 0.002186 (0.0391);
- W5 range: 0.002986.

Every other run and endpoint passes under every estimator, and W5 holds everywhere else (`s1_interval_results.json#EP1`, `#runs[].warnings`).

### 4.2 Answers to the ten questions

1. **Operative EP-1 wording:** B1, verbatim.
2. **The §2.1 definition:** the within-run endpoint MC-SE is the **kernel quantile SE** (row 2). The between-run SD is "endpoint variability" (row 3). Both are classed as Monte-Carlo quantities (line 100), and W4 calls both "endpoint MC-SE" (line 468).
3. **Silverman's rule** is the named bandwidth rule for the kernel estimator (line 105). The variant is unstated (Convention 7).
4. **Half and double bandwidth: sensitivity analyses that must be reported. Whether they are also acceptance inputs is not stated.**
   - The rationale clause, "a single number would overstate precision", is compatible with two readings:
     - (i) report the spread so readers see how uncertain the SE estimate is, with EP-1 judged at Silverman's bandwidth;
     - (ii) judge EP-1 robustly, so that a pass at one bandwidth does not count as precision.
   - **The text does not choose between them.**
5. **§2.4 when EP-1 fails** (mandatory):
   - compute B_out′ from the "observed ratio";
   - round up to the ladder;
   - top up **each run** from the same spawn list;
   - never narrow, re-centre or re-quantile;
   - if the ladder is exhausted, H3 applies (IMPRECISE).
6. **W5 when breached:** the row's last cell says "raise B_out". The section says "investigate and report".
7. **Instruction, warning or recommendation?** **Ambiguous.**
   - W5 is classed as a warning (line 461).
   - "raise B_out" sits in a column headed "Justification of the threshold" (line 463), which it is not. It is a remedial action, like W3's "raise B_in" (line 467).
   - Line 473, "B_out and B_in may only increase", presupposes that increases are an anticipated response to warnings.
   - The most faithful description: **a warning with a stated remedial direction**. It is not marked mandatory and it has no sizing rule.
8. **Does W5 independently trigger a top-up?**
   - **Not by any explicit rule.** §2.4 is triggered only by EP-1 failure, and no text links W5 to the §2.4 ladder or to any amount.
   - It is not excluded either: "raise B_out" is an unconditional phrase in its row.
   - **The plan does not settle this.**
9. **Can the top-up rule be applied to the half-bandwidth result?**
   - **Arithmetically, yes:** 500 × (0.0535 / 0.05)² = 573, which rounds up to **750** (stored `required_B_out` 572.8).
   - **Procedurally, only if ×0.5 is an EP-1 input**, which is the unresolved point in (4).
   - **The scaling assumption** `SE ∝ 1/√B_out` (line 118) holds for the kernel SE only if f̂(q) is stable as B_out grows. Silverman's h itself shrinks as B_out^(−1/5), and the ×0.5 estimate has the highest density-estimation variance, so 750 is a rule-of-thumb projection, not a guarantee.
   - **For W5,** applying the §2.4 formula to a range statistic would be a new rule.
10. **Is there an unambiguous decision rule?** **No.** The plan leaves open:
    - which kernel bandwidth (or which Monte-Carlo estimator) EP-1 uses;
    - whether W5 is mandatory;
    - how far B_out must rise under W5;
    - how "spread" is defined (range is implemented; MC-1 itself is a visual criterion, line 223).

### 4.3 Interpretation

- **The two signals agree with each other.** The marginal endpoint is the same one in both: run 1 upper, ratio 0.0535 under ×0.5 and 0.0534 under W5's range (*derived*: 0.002986 / 0.055924).
- **Together they indicate that this upper quantile is the least stable endpoint at B_out = 500.** That is a Monte-Carlo-precision statement. It says nothing about sampling validity.
- **This coincidence of two independent signals argues against dismissing either one.** It does not convert either into a mandatory rule.

### 4.4 Unresolved ambiguity

These are Q4, Q7, Q8 and Q10 above. They are matters of **protocol interpretation**. They should be resolved by an explicit clarification record, made **before** any data under the clarified rule exist.

### 4.5 Options and consequences

| Option | Rule | Consequence for the as-run data | Consequence for a corrected rerun |
|---|---|---|---|
| **B-1** No top-up | EP-1 at Silverman ×1; ×0.5 / ×2 reported; W5 a warning, investigated and reported | EP-1 PASS (10 / 10); no new replicates | The same rule applied to the rerun at B_out = 500 |
| **B-2** Top-up to the ladder rung | EP-1 must hold at all three bandwidths (the maximum governs), **or** W5 treated as mandatory and sized by §2.4 analogy | EP-1 FAIL, so 750 per run (≈ 1 250 replicates). Appends to protected checkpoints. Would extend the CV-seed deviation. | Evaluated at 500 under the rule; topped up only if it triggers |
| **B-3** Clarify the protocol first (**recommended**) | The owner fixes, in a dated clarification record: (a) the EP-1 estimator rule (B-1's or B-2's); (b) whether W5 obliges an increase, and the sizing; (c) the definition of "spread" | **No action on the as-run data.** They become historical once a corrected rerun is primary. Their EP-1 status is reported under both readings. | The rule is fixed before the rerun's data exist, then applied mechanically. Because streams are index-derived (lines 127–128, 371–372), running 500 first and topping up later is bit-identical to running the larger size directly, so nothing is lost by evaluating at 500 first. |

### 4.6 Recommendation: **B-3**

**Do not top up the as-run checkpoints.** Under any reading:
- topping up would propagate the CV-seed deviation;
- it would modify protected files;
- it would spend computation on data that become historical once the corrected rerun is primary.

**The substance of the clarification is the owner's decision.**
- **Stated neutrally:** B-2's estimator rule matches the rationale clause (lines 106–107) and W5's remedial direction (line 469). B-1's matches the naming of Silverman's rule (line 105) and the singular "MC-SE" in EP-1.
- **What the clarification must not do:** be chosen on the basis of the corrected rerun's results. Hence it must be made first.

---

## 5. Issue C: A7 and A8, mathematical verification

### 5.1 Definitions and derivation

**Notation.**
- x₁…x_B are the outer replicates ΔC_corr(k), with B = B_out.
- q(p) is the type-7 empirical quantile, and q_a = q(α/2), q_b = q(1 − α/2), with α = 0.05.
- θ̂ is the fixed reference value.
- F̂(t) = #{xₖ < t} / B, with strict inequality.

**Interval definitions (plan §7, lines 309–316; `stats_core.py:26–35`).**
- Percentile: [L_p, U_p] = [q_a, q_b], with width W = q_b − q_a.
- Basic: [L_b, U_b] = [2θ̂ − q_b, 2θ̂ − q_a], built from **the same** q_a, q_b (`basic_ci` calls `percentile_ci`).

**A7 (W6).**
- Plan: "|percentile − basic| at either endpoint", threshold "25 % of width" (lines 327, 470).
- Implementation: A7_lo = |L_p − L_b| / W and A7_hi = |U_p − U_b| / W (`s1_interval_analysis.py:88`).

**A8 (W7).**
- z₀ = Φ⁻¹(F̂(θ̂)) (line 316), and A8 = |z₀| with threshold 0.25 (lines 328, 471; `stats_core.z0`).

**Derivation of A7.**
- Step 1: L_p − L_b = q_a − (2θ̂ − q_b) = q_a + q_b − 2θ̂.
- Step 2: U_p − U_b = q_b − (2θ̂ − q_a) = q_a + q_b − 2θ̂.
- Step 3: steps 1 and 2 are **identical as signed quantities**, so A7_lo = A7_hi.
- Step 4: with m = (q_a + q_b) / 2, the difference is q_a + q_b − 2θ̂ = 2(m − θ̂). Hence **A7 = 2·|m − θ̂| / W**.
- Step 5: the basic width is (2θ̂ − q_a) − (2θ̂ − q_b) = W, so normalising by either interval's width gives the same value.

**Assumptions** (all hold in the current construction):
- (i) the basic interval reflects the **same** two quantiles used for the percentile interval;
- (ii) both use the same α and quantile estimator;
- (iii) A7 is normalised by W.

The identity is **algebraic**. It holds for every dataset and every θ̂ under these assumptions, so it is not an empirical finding. It would fail if, for example:
- the basic interval used different quantile levels or a different quantile type;
- the basic interval were studentised;
- the two endpoints were normalised differently.

**Decomposition of A7.**

Let med = q(0.5). Then m − θ̂ = (m − med) + (med − θ̂), so

A7 = |T_asym + T_med|, with T_asym = 2(m − med) / W and T_med = 2(med − θ̂) / W.

- **T_asym** is the asymmetry of the 2.5 % and 97.5 % quantiles about the median. It does **not** involve θ̂.
- **T_med** is the data-scale median offset from θ̂.

**A8 by contrast.**
- z₀ depends on the data only through the **rank** of θ̂ among the xₖ. It is invariant to monotone transformations of the replicates.
- z₀ ≈ 0 if and only if θ̂ ≈ med; its sign is that of (θ̂ − med).
- It does not use q_a, q_b or W.

**Relationship between the two.**
- **In general:** A7 and A8 share only the θ̂ position. A8 contains no tail-asymmetry information, and A7 is in data units scaled by the interval width. **They are neither equivalent nor proportional.**
- **Special case:** if the xₖ were exactly normal N(μ, σ²) (and ignoring quantile-estimation noise), then m = med = μ, W = 2·1.96σ, and z₀ = (θ̂ − μ)/σ, giving A7 = |z₀| / 1.96. Proportionality holds **only** under this symmetric special case.

### 5.2 Answers to the ten questions

1. **A7:** |L_p − L_b| / W at each endpoint. Under the current construction this equals 2·|m − θ̂| / W.
2. **A8:** |Φ⁻¹(#{xₖ < θ̂} / B)|.
3. **Relationship:** neither equivalent nor proportional. They are related through θ̂'s position, and proportional only under exact normality. Empirically they are correlated, through T_med.
4. **The identity A7 = 2·|m − θ̂| / W:** yes, it holds exactly at **both** endpoints (steps 1–4).
5. **Same quantities and normalisation?** No. A8 uses rank, not the two quantiles or the width, and is scaled in standard-normal units.
6. **Why the recorded values are identical at both endpoints:** by algebraic necessity (step 3), not because of any property of these data. **The per-endpoint A7 is therefore one statistic, not two checks.**
7. **Generality:** it holds generally under assumptions (i)–(iii), and only under the current interval construction.
8. **Independent information?** Partly. They share the θ̂-dependent median-offset component, so they are not independent confirmations of a centring problem. A7 also carries the tail asymmetry T_asym, which A8 cannot detect. A8 is the cleaner median-bias detector, and the plan grounds it partially (≈ 4.5 SE, line 328). A7's threshold has "**NONE**" (no statistical justification) (line 327).
9. **What the plan instructs:** "A breach of either triggers investigation and explicit reporting" (line 330). The basic-versus-percentile divergence is "a warning that **neither** interval should be trusted" (line 314).
10. **Invalidation or qualification?** **Investigation and qualification.** "Neither a breach nor a pass establishes that the interval is valid or invalid" (lines 330–331). Warnings "never auto-validate or auto-invalidate" (line 461). The UNRELIABLE and IMPRECISE labels are reserved for H1–H3 (lines 454–456).

### 5.3 Project-specific values (*derived* from stored endpoints, medians, θ̂ and z₀; θ̂ = 0.0043163)

| Run | A7 (stored) | A7 recomputed, 2·\|m − θ̂\|/W | T_asym (θ̂-free) | T_med | \|z₀\| | \|z₀\|/1.96 (normal-case value) |
|---|---|---|---|---|---|---|
| 0 | 0.5644 | 0.5644 | **0.348** | 0.216 | 0.468 | 0.239 |
| 1 | 0.5217 | 0.5217 | **0.270** | 0.252 | 0.519 | 0.265 |
| 2 | 0.5286 | 0.5286 | 0.230 | 0.299 | 0.571 | 0.291 |
| 3 | 0.5393 | 0.5393 | **0.294** | 0.245 | 0.429 | 0.219 |
| 4 | 0.5546 | 0.5546 | **0.299** | 0.256 | 0.577 | 0.294 |
| pooled (descriptive) | — | 0.5526 | **0.300** | 0.253 | 0.512 | 0.261 |

**Interpretation.**
- **The identity reproduces every stored A7** to 4 decimal places.
- **Tail asymmetry is large.** T_asym (bold where it alone exceeds 0.25) is between 0.23 and 0.35, consistent with the right skew (pooled skewness 0.79). In four of five runs it would breach the 25 % threshold even if θ̂ were exactly the median, with z₀ = 0.
- **A7 is roughly double its normal-case value.** Observed A7 is about twice |z₀|/1.96, which quantifies how far from the symmetric case these data are.
- **Consequence:** A7's breach is driven **jointly** by skewness (θ̂-independent) and median offset (θ̂-dependent). A8's breach reflects median offset only.

### 5.4 Unresolved ambiguity

- **Normalising width.** The plan says "25 % of width" without naming which interval's width. Step 5 shows the choice does not matter here.
- **What "investigation" must contain.** The plan does not say, and no investigation beyond the descriptive decomposition has yet been performed.

### 5.5 Recommendation and option consequences

**Recommendation.**
- Report A7 and A8 as **distinct, partly overlapping** diagnostics.
- Report A7 as **one statistic**, since its two endpoint values are algebraically identical, with its decomposition into T_asym and T_med.
- Report A8 as the median-offset detector.
- Keep both as **warnings requiring qualification**. Neither invalidates the interval.

**Consequences of how they are reported.**
- **Treating them as one signal** (the earlier review's wording) would understate A7's skewness content.
- **Treating them as two independent confirmations** would overstate the evidence for a centring problem.

---

## 6. Cross-decision dependencies

| Question | Answer |
|---|---|
| **Which decisions must be frozen before the corrected rerun?** | (1) Plan §8.2 binding, and the as-run seeding recorded as a deviation. (2) The apparent-CV convention (A-I or A-R; §3.4), and whether the A-R secondary is computed. (3) The written pre-designation that the corrected results are primary, whatever they show. (4) The EP-1 / W5 clarification (B-3 content), because it governs whether the rerun is topped up. (5) The verification gate V1–V5. (6) A new output directory and a separately frozen driver hash; the original driver is untouched. |
| **Which can be adjudicated independently?** | The A7 / A8 reporting form (Issue C). The five-run reporting format. Retention of θ̂. None of these depends on rerun data or on each other. |
| **Does the apparent-CV convention change the estimand?** | **No.** It changes the finite-sample implementation of a randomised estimator. It can change realisations, and possibly the distribution the bootstrap targets if the randomisation is coupled to the resample. The estimand of line 6 is unchanged (§3.2). |
| **Should θ̂ remain fixed for the primary analysis?** | **Yes.** It is pre-specified (lines 6, 139) and its fold randomisation on D is independent of D's sampling. That matches the independence principle behind A-I. Any θ̂ sensitivity is supplementary and must not replace it. |
| **Is the five-run reporting recommendation still defensible?** | **Yes**, and unaffected by Issues A–C. It uses only plan-defined quantities (lines 98, 244, 309). |
| **Can the existing 2 500 replicates remain as historical evidence?** | **Yes**, labelled *as-run / CV-seed deviation*. They are reproducible (H5) but not protocol-conforming. They should not be topped up, and they should not be pooled with corrected replicates, because the two use different procedures. |
| **What must be documented before the next computational stage?** | The checklist in §10. |
| **What the rerun will show** | **No claim is made** that the corrected rerun will widen, narrow, move or preserve the interval, or change any warning. |

---

## 7. Recommended pre-rerun decisions

1. **Plan §8.2 is binding.** The as-run inner-CV seeding is recorded as a protocol deviation, and the v2 draft is adjudicated accordingly.
2. **Apparent-CV convention.** A-I as primary (`[BASE, run_id, T]` with non-zero T, a fresh generator per model), recorded as a **new methodological decision / plan clarification**. Optionally, A-R as a pre-declared secondary.
3. **D.** No new RNG specification. θ̂ is retained unchanged.
4. **EP-1 / W5.** The owner fixes, in writing and before the rerun:
   - the governing estimator (×1, or the maximum over the three bandwidths);
   - whether a W5 breach obliges an increase in B_out, and the sizing;
   - the definition of "spread" (range, as implemented, or another).
5. **No top-up of the as-run checkpoints.**
6. **Pre-designation.** The corrected results are primary, whatever their direction. The as-run results are historical.
7. **Gate.** V1–V5 must pass before the full rerun's results are used.
8. **Reporting.** Five per-run percentile intervals. A7 as one statistic with its T_asym / T_med decomposition. A8 separately. Both reported as warnings with qualification.

---

## 8. Decisions requiring explicit project-owner authorisation

1. Confirm that plan §8.2 is binding and that the as-run seeding is a recorded deviation.
2. Choose the apparent-fit CV convention, A-I or A-R. If A-I, approve the domain tag T. Decide whether the A-R secondary is computed.
3. Confirm that θ̂ = +0.0043 is retained, with no RNG specification for D.
4. Fix the EP-1 governing-estimator rule.
5. Rule on W5: is a breach an obligation to raise B_out? If so, by what rule? Define "spread".
6. Confirm that the as-run checkpoints will not be topped up.
7. Pre-designate the corrected results as primary, in writing.
8. Approve the verification gate V1–V5 and the new output location.
9. Approve the A7 / A8 reporting form (§5.5) and the five-run format.
10. Authorise the corrected rerun itself (≈ 19.5 h at the observed throughput), as a separate act after decisions 1–9 are recorded.

---

## 9. What must remain frozen

- The **authoritative plan** v2.1 (`d27c52a29930a338`). Clarifications are recorded as separate, dated documents, never as edits.
- The **frozen driver** (`4b1738a4d0ee0e39`) and `src/phase9_execution.py` (`571d0c233f367be6`). A corrected driver is a **new file** with its own frozen hash.
- The **Phase 9 primary script** and its results (`82afbcf05f3b7514`, `05d13b8b339fa7d4`), and therefore θ̂.
- The five **as-run checkpoints**, the halt and resume records, and `resume_driver.log`.
- The **H5 evidence** and deviation draft v1.
- The **seed definitions**: `BASE = 2026110300`, run_id 0–4, and every existing stream in the `[BASE, run_id]` tree.
- The **statistical methodology**: percentile primary; basic and z₀ diagnostic; BCa excluded; EP-1 threshold 0.05; the §2.4 ladder; the H and W thresholds.
- Every existing result, table, figure, report and audit record in `final_project/`, including `METHODOLOGICAL_DECISION_REVIEW.md`.

---

## 10. Final pre-authorisation checklist

| # | Item | Status |
|---|---|---|
| 1 | §8.2 binding; deviation adjudicated | **NOT YET ADJUDICATED** |
| 2 | Apparent-CV convention chosen and recorded as a clarification | **NOT YET ADJUDICATED** |
| 3 | A-R secondary: yes or no | **NOT YET ADJUDICATED** |
| 4 | θ̂ retained; no D specification | **NOT YET ADJUDICATED** (recommended: retain) |
| 5 | EP-1 governing-estimator rule recorded | **NOT YET ADJUDICATED** |
| 6 | W5 obligation, sizing and "spread" definition recorded | **NOT YET ADJUDICATED** |
| 7 | No top-up of as-run checkpoints confirmed | **NOT YET ADJUDICATED** |
| 8 | Corrected results pre-designated primary, in writing | **NOT YET ADJUDICATED** |
| 9 | Corrected driver written, reviewed and hash-frozen as a new file | **BLOCKED** (by 1–8) |
| 10 | Gate V1–V5 specified, passed before use | **BLOCKED** (by 9) |
| 11 | A7 / A8 reporting form and five-run format approved | **NOT YET ADJUDICATED** |
| 12 | Rerun authorised | **BLOCKED** (by 1–10) |

**The project is not ready for rerun authorisation until items 1–8 are recorded.**

---

*This audit created one file only. No computation was initiated. No checkpoint, frozen file, plan, driver, historical implementation, H5 artefact, existing result, report or audit record was modified, and nothing was committed, pushed or transmitted.*
