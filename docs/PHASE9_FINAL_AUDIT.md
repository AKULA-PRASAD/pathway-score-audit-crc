# Phase 9 — Final Protocol Audit

**Audited version:** 2.0 · **Issued:** **v2.1** · **Date:** 2026-11-03
**Script:** `scripts/37_phase9_final_audit.py` · **Tests:** `tests/test_phase9_audit.py`
**Artifact:** `reports/phase9/final_audit_summary.json`

> **0 models fitted · 0 decisions applied · 0 endpoint access · 0 commits.**
> No event count was read, computed or inferred. No frozen record was modified.

---

## Finding 1 — the 453 vs 452 batch counts: **no error of fact; the premise conflates two sets**

**The premise does not hold.** The audit request states that *"the A-19 complete-case cohort
contains 453 samples."* It does not. Two distinct sets are in play:

| Set | n | Defined by |
|---|---:|---|
| **S3 cohort** | **453** | Phase 1 (A-01) + Phase 7C D2 — **no exclusions** |
| **A-19 complete-case model-fitting set** | **452** | A-19 drops **GSM972262**, missing `chemotherapy.adjuvant` |

The nine counts in Phase 9 v2.0 sum to 452 because they are the **complete-case** counts, and v2.0
labelled them *"(complete case, n = 452)"*.

### Verification (outcome-blind, explicit sample IDs)

| Check | Result |
|---|---|
| Complete-case set derived independently over the five M1 covariates | **452**, and the derived drop **equals** the A-19 record exactly |
| S3 IDs absent from batch metadata | **0** |
| Duplicate `sample_id` in batch metadata | **0** |
| S3 IDs with blank / `N/A` batch | **0** |
| S3 IDs mapping to more than one batch row | **0** |
| Distinct batch labels over S3 | **9** |
| Bundle batch == ArrayExpress SDRF batch, all 453 | **True** |

| batch | S3 cohort (453) | complete-case (452) |
|---|---:|---:|
| batch01 | 33 | 33 |
| batch02 | 14 | 14 |
| batch03 | 102 | 102 |
| batch04 | 103 | 103 |
| **batch05** | **49** | **48** |
| batch07 | 81 | 81 |
| batch08 | 41 | 41 |
| batch09 | 15 | 15 |
| batch10 | 15 | 15 |
| **total** | **453** | **452** |

The single difference is **GSM972262**, whose batch is **batch05** in both the Phase 3 bundle and
the SDRF, and whose `chemotherapy.adjuvant` is `N/A` in the original metadata and `""` in the
bundle — exactly the A-19 drop.

**Verdict: not a data error, not a missing or ambiguous batch assignment, and not an error of
fact in v2.0.** Nothing was invented and no batch was assigned. **All 453 S3 samples remain in the
cohort**; A-19 affects only which rows enter a model fit. The event status of GSM972262 was **not
inspected**.

**Correction issued — clarity, not fact.** That this generated an audit request is itself evidence
the v2.0 wording was easy to misread. **v2.1** prints **both** count columns side by side
(`phase9_summary.json` → `batch_counts_both_sets`) and states the cohort / complete-case
distinction explicitly wherever a batch count appears.

---

## Finding 2 — the D8D-4 "50 % of events" rule: **WITHDRAWN**

**Verdict: the rule is both logically ill-defined and statistically unfounded. It is withdrawn and
replaced.**

### Defect 1 — logical: the criterion is near-vacuous

As worded, it fires when *"fewer than 50 % of events lie in strata that contribute to the partial
likelihood."* **Every event lies, by construction, in a stratum containing at least one event —
itself.** If "contributes" means "has ≥ 1 event", the measured quantity is **100 % always** and the
threshold can never fire. The criterion does not do what it appears to do.

### Defect 2 — statistical: 50 % has no basis

It is not derived from power, from information, or from any property of the partial likelihood. A
stratified fit with 40 % of events in comparable risk sets can be perfectly well estimated; one
with 60 % can still be unstable. Worse, it applies the word **"estimable"** to what is really a
statement about **precision**.

### Three concepts the rule conflated

| Concept | What it actually is |
|---|---|
| **Mathematical estimability** | β is identified iff **at least one** stratum contains an event whose within-stratum risk set holds ≥ 2 subjects **and** those subjects differ on a covariate. An extremely weak condition, essentially always met at this cohort size. A property of design + data, not of a percentage |
| **Numerical convergence** | separate: Newton–Raphson can fail through monotone likelihood, separation, or a rank-deficient design even when β *is* identified |
| **Useful information** | **continuous**. Stratification only removes between-stratum comparisons, never adds any, so the stratified fit carries **less** information about β. *How much less* is measurable — and is the thing worth reporting |

### How strata actually contribute

- **Zero-event stratum** → contributes a likelihood factor of **1**. It drops out entirely: no term, no score, no information. Its patients are not "lost" in any estimability sense; they simply never enter a comparison.
- **Sparse-event stratum** → contributes its own event terms, informative only insofar as the within-stratum risk set at that event time holds other subjects who **differ on the covariates**.
- **Event with a degenerate risk set** (the subject alone at risk) → contributes **log(1) = 0**: zero score, zero information — *while still being counted* in any naive events-per-stratum tally. **This is precisely the trap the withdrawn rule falls into.**

The informative quantity is therefore **not** "events per stratum" but **the number of events whose
within-stratum risk set yields a non-degenerate, covariate-varying comparison**.

### ⛔ Proposed replacement **D8D-4-R1** — **PENDING USER APPROVAL**

> **Principle.** Declare non-estimability on **mathematical and numerical grounds only** — criteria
> whose *definitions* are outcome-blind — and report information loss as a **continuous
> diagnostic** rather than gating on an invented threshold.
>
> **The batch-stratified Cox sensitivity model is declared NOT ESTIMABLE if and only if:**
> 1. the fit fails to converge under the declared criteria (log-likelihood change < 1e-9, max 100 iterations), **or**
> 2. the observed information matrix for β is singular or not positive definite (numerical rank < number of coefficients), **or**
> 3. the separation sentinel trips (any \|β\| > 10 or SE > 10).
>
> **Reported always, never gated on:**
> - per-stratum **sample** counts (available now, outcome-blind — Finding 1's table);
> - per-stratum **event** counts (at authorised execution);
> - per-stratum count of events with a **non-degenerate within-stratum risk set** — the events that actually generate a comparison;
> - the number of strata contributing ≥ 1 such event;
> - convergence status, iteration count, final gradient norm;
> - the **information-retention ratio**: median over coefficients of (SE<sub>unstratified</sub> / SE<sub>stratified</sub>)², a direct and interpretable measure of what stratification cost;
> - an explicit statement of the limitation beside the estimate.
>
> **No threshold is placed on the information-retention ratio.** Any cut-point would reintroduce
> the arbitrariness just removed.
>
> **If declared not estimable:** that fact is **reported**, and the leave-one-batch-out diagnostic
> stands alone. **We do not fall back to pooling** — a pooling choice made after seeing
> events-per-stratum would be an outcome-informed design change.

**Unchanged:** all nine ComBat batches remain individual strata · **no pooling, by count or
otherwise** · no post-outcome pooling and no redefinition of strata based on the observed event
distribution. The three estimability criteria are defined without reference to any outcome
quantity; their *evaluation* requires a fit, which happens only at authorised execution — exactly
as event-stratified fold construction does.

---

## Finding 3 — D9-1 penalisation scope: **SOUND, retain**

**Verdict: sound. Retain, with two interpretive statements made explicit.**

| Question | Finding |
|---|---|
| Is the M1/M2 comparison still interpretable? | **Yes — and improved.** With the clinical block unpenalised in **both** models, M1 is exactly the restriction of M2 obtained by setting the pathway coefficients to zero |
| Is the nested comparison coherent? | **Yes.** Penalising only the molecular block means λ → ∞ recovers M1 **exactly**. Had the clinical block been penalised with a model-specific λ, M1 would **not** be the limit of M2 and the nesting would be broken — so this choice is what **makes** the comparison coherent |
| Are clinical coefficients comparable across models? | **Improved.** Unpenalised in both, they differ only through genuine **adjustment** for the pathway block. Under a penalised clinical block with model-specific λ, part of any difference would be a **shrinkage artefact** and the two sets would not be comparable at all |
| How do scaling and functional form interact with the penalty? | D8D-1 makes the penalty scale-fair **across** the pathway scores; D9-1 removes the clinical scaling question **entirely**, since an unpenalised coefficient is **scale-equivariant**. Leaving the age basis unpenalised keeps the spline at its declared **2 df**; penalising it would shrink the *nonlinearity* — implicit smoothing with an effective df nobody specified. **The two decisions are complementary and neither alone is sufficient** |
| Does it change the molecular interpretation? | **Yes, and this must be stated.** ΔC(M2 − M1) estimates the gain from adding a **shrunk** pathway block to an **unshrunk** clinical model. The penalty pulls the pathway coefficients toward zero — toward M1 — so the comparison is **conservative for the molecular block**: biased *against* finding a molecular contribution, not toward one. That is the right direction of error for this study, but it must be stated, not relied on silently |
| Does Phase 6B permit it? | **Not as written.** Penalisation **scope** is not a Register A item (A-11 fixes the predictor *list*, A-15 the *roles*, A-17 the resampling). It is protocol text, so it requires an explicit amendment — **not** a Register A change |

### ⛔ Proposed amendment **P6B-A1** — **PENDING USER APPROVAL** · **NOT APPLIED**

| | |
|---|---|
| **Target** | `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` §3.1, the *Penalisation* row of the M1 table |
| **From** | `ridge, λ tuned inside resampling` |
| **To** | `the clinical block (four binary covariates and the age spline basis) is FORCE-INCLUDED and UNPENALISED in M1, M2 and M3. The ridge penalty applies only to the M2 pathway scores and the elastic-net penalty only to the M3 candidate genes; λ (and α for M3) are tuned by inner cross-validation inside each training fold. M1 therefore has no penalised term and no λ to tune.` |
| **Register A item?** | **No** |
| **Consequence to state** | M1 ceases to be a penalised model; any comparison with the 6B draft text must say so |

**Residual caveats.** M1's apparent performance is unshrunk and carries the usual optimism — which
the bootstrap correction estimates, with the uniform shrinkage factor reported regardless. If a
future protocol wished to penalise the clinical block, the scaling question returns.

---

## Finding 4 — internal consistency: **33 / 33 checks pass**

Machine-readable: `reports/phase9/final_audit_summary.json` → `finding_4_consistency`.

| Area | Checks | Result |
|---|---:|---|
| Register A (A-01, A-04/05, A-13, A-14/M2-A1, A-19) | 5 | ✅ |
| M1/M2/M3 definitions and predictor roles | 2 | ✅ |
| Training-fold-only preprocessing; nested tuning | 3 | ✅ |
| M2-A1 vs proposed M2-A2 | 3 | ✅ |
| Age spline knots, reference age, fallbacks | 3 | ✅ |
| Batch sensitivity specification | 3 | ✅ |
| Bootstrap and Monte-Carlo stopping rules | 3 | ✅ |
| Non-estimable replicate handling | 1 | ✅ |
| U1 / U2 endpoint restrictions | 2 | ✅ |
| Competing-risk and composite-DFS estimands | 2 | ✅ |
| External-validation limitations | 2 | ✅ |
| Complete-case sample and event-count statements | 3 | ✅ |
| Execution boundary | 1 | ✅ |

Notable confirmations: **M2-A1 is the only amendment on disk**; `score_standardisation` still reads
`NONE`; the complete-case **event** count is left undetermined and event counts are **cited, never
recomputed**; no fixed-time quantity appears in the estimated set while all four U1-blocked
quantities appear in the *not-estimated* list; S-1 is marked "NOT comparable" with the primary and
S-3 "CONDITIONAL on U1".

---

## Amendments now awaiting approval

| ID | What | Status |
|---|---|---|
| **M2-A2** | score-level standardisation, training-fold only (unchanged from v2.0) | **PENDING** · not applied |
| **P6B-A1** | Phase 6B §3.1 penalisation-scope row (**new**) | **PENDING** · not applied |
| **D8D-4-R1** | replaces the withdrawn 50 % rule (**new**) | **PENDING** · not applied |

Plus the five standing decisions **D8D-1 … D8D-4, D9-1**, and approval of the protocol as a whole.

**Frozen records modified: NONE.**
