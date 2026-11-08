# Phase 6A — Endpoint Resolution (U1, U2) — rev. 2

**Date:** 2026-11-02 · **Revision 2** — corrections applied.
**Status:** Investigation complete. **No source value changed, no exclusion changed, no model built, Phase 3 bundle untouched.**

---

## 0. Outcome in one line

| Item | Verdict |
|---|---|
| **U1 — absolute time unit of `rfs.delay`** | **REMAINS UNRESOLVED.** Months is strongly supported but **never documented**. Status stays `[PROVISIONAL]` |
| **U2 — death without recorded recurrence** | **The published rule is UNRESOLVED.** The *deposited encoding* is censoring-at-death — and that censoring is **demonstrably informative** (new finding, §2.3) |

---

## 1. U1 — absolute time unit

### 1.1 Sources searched (eight; two new this phase)

| Source | Result | Tier |
|---|---|---|
| GEO SOFT `!Sample_characteristics_ch1` | `rfs.delay` **unlabelled**; `os.delay (months)` **labelled** | evidence |
| Marisa et al. 2013 main text | no unit statement for the column | evidence (absence) |
| Text S1 supplementary methods | "five years/60 months" — the **analysis horizon**, not the column | evidence |
| Table S1 per-patient annotation | single sheet, no legend, `RFS.delay` unlabelled | evidence (absence) |
| ArrayExpress E-GEOD-39582 SDRF | `Characteristics [rfs.delay]`, **no `Unit[]` column anywhere** | evidence (absence) |
| GEO SuperSeries GSE40967 | no RFS metadata at all | evidence (absence) |
| **NEW — published median follow-up** | §1.2 | **inference** |
| **NEW — downstream reanalyses of GSE39582** | §1.3 | **inference** |

### 1.2 New evidence — the published median follow-up

Marisa et al. report *"median follow-up of 51.5 mo"*. Computed from the deposit:

| Subset | median `os.delay` | median `rfs.delay` |
|---|---:|---:|
| all 585 deposited | 52.00 | 44.00 |
| **566 primary tumours** | **51.00** | 43.00 |
| 443 discovery | 50.00 | 41.00 |

**51.5 sits between the 566-sample (51.00) and 585-sample (52.00) medians of `os.delay`.** I could **not** reproduce 51.5 exactly on any subset tried, so the precise definition the authors used is not pinned down.

**What this establishes:** the paper's follow-up summary is consistent with `os.delay` in months — a column already labelled months. **It is corroboration for `os.delay`, not direct evidence about `rfs.delay`.** It is recorded as inference, not evidence.

### 1.3 New evidence — downstream reanalyses

A literature search surfaced a secondary analysis reporting *"The median follow-up time was **51.0 months** in the GSE39582 cohort"* — exactly the median `os.delay` computed above for the 566 primary tumours.

**What this establishes:** at least one independent group read the deposit's time column as months and obtained the same number. **It is third-party corroboration, not documentation** — they were reading the same unlabelled column and may have made the same assumption.

### 1.4 The internal-consistency argument

Across the 560 tumours with both times:

| Constraint | Result |
|---|---:|
| `rfs.delay > os.delay` (impossible if same unit) | **0 violations** |
| ranges | **identical**: both `[0, 201]` |
| both integer-valued | yes |
| `rfs.delay == os.delay` exactly | **387 (69.1%)** |

**The 387 exact equalities are the strongest constraint.** If the two columns were in different units, those 387 patients would each have a recurrence/last-contact time differing from their death/censoring time by exactly the unit ratio — e.g. recurring at precisely 1/30th of their OS time, 387 times over. That is not a credible data pattern.

### 1.5 Verdict

> **U1 REMAINS UNRESOLVED.** Months is supported by: a labelled sibling column in the same units; zero ordering violations across 560 records; identical ranges; 387 exact equalities; a consistent published median follow-up; and an independent group's reading. **None of these is documentation.** No source states the unit of `rfs.delay`.
>
> **Status: `[PROVISIONAL: months]`.** It must be labelled as an assumption in every artifact and may only be upgraded by an explicit statement from the CIT programme or the corresponding author.

The one untried avenue remains **direct correspondence with the authors**, which is outside what this investigation can do unilaterally.

---

## 2. U2 — death without recorded recurrence

### 2.1 The published rule — unresolved

The main text, Text S1 and Table S1 were all searched. **"censor" appears zero times in Text S1.** No source states whether death without recurrence was treated as an event, as censoring, or as a competing risk.

**Explicitly acknowledged:** matching RFS and OS follow-up values **cannot establish a published rule**. The 49/49 pattern (46/453 within S3) shows what the *deposited data encodes*; it says nothing about what the authors intended or wrote.

### 2.2 What the deposited encoding does

| Observation | Value |
|---|---|
| Deaths without recorded recurrence, S3 | **46 / 453 (10.2%)** |
| Of those, `rfs_event = 0` | **46/46** |
| Of those, `rfs_time == os_time` | **46/46 (100%)** |

**Operationally this is settled regardless of intent:** any analysis using `rfs.event`/`rfs.delay` as deposited treats death without recurrence as **censoring at the date of death**. There is no ambiguity about what the fields *do*; the ambiguity is about what was *meant*.

### 2.3 NEW FINDING — that censoring is informative

Censoring at death is only valid if those patients would have carried the same future recurrence hazard as those remaining at risk. They do not:

| Group | n | median age | stage III | adjuvant chemo |
|---|---:|---:|---:|---:|
| **died without recurrence** (censored) | **46** | **75.0** | **32.6%** | **23.9%** |
| alive, no recurrence (true censoring) | 271 | 68.0 | 38.7% | 42.1% |
| recurrence | 136 | 68.5 | 58.8% | 55.1% |

**The censored-at-death group is ~7 years older, less often stage III, and roughly half as likely to have received chemotherapy.** That is the profile of competing non-cancer mortality, not of a random sample of those at risk.

> **Conclusion: censoring at death in this cohort is INFORMATIVE.** Death is a **competing event**, not a nuisance censoring mechanism. This is a methodological finding about the data, independent of what the authors intended.

### 2.4 Verdict

> **U2 REMAINS UNRESOLVED as a published rule**, and is now **resolved as a methodological problem**: 10.2% of S3 experiences a competing event whose occurrence is associated with age, stage and treatment. The analysis must address this explicitly rather than inherit it silently.

---

## 3. Consequences

### 3.1 U1 — what it does and does not block

**Demonstrated** (synthetic data, `scripts/`-free demonstration): the Cox partial likelihood depends only on the **ordering** of times. Rescaling by any positive constant leaves β, SE and C **identical to 10 decimal places**:

| Rescaling | β | SE | C |
|---|---|---|---|
| as-is | 0.9252625712 | 0.0579913478 | 0.666870 |
| × 30.44 (months→days) | 0.9252625712 | 0.0579913478 | 0.666870 |
| × 1/12 (months→years) | 0.9252625712 | 0.0579913478 | 0.666870 |

| Analysis | Affected by U1? |
|---|---|
| **Cox relative-risk (hazard ratios)** | **NO — fully invariant** |
| **C-index / discrimination** | **NO — rank-based, invariant** |
| Calibration slope on the log-hazard scale | **NO — invariant** |
| **Absolute survival / recurrence-risk predictions** | **YES — baseline hazard is expressed in the unknown unit** |
| **Fixed-time predictions (t-year risk)** | **YES — blocked. "60 units" is 5 y, 1.15 y or 0.16 y depending on the unit** |
| **cloglog calibration at a horizon (Riley's estimator)** | **YES — blocked; needs a horizon** |
| **Time-dependent AUC / Brier at t** | **YES — blocked** |
| The conditional 60-unit endpoint | **YES — conditional** |
| **Event counts** | **NO — counts are unit-free** |

### 3.1.1 The operative distinction — relative risk vs fixed-time absolute risk

> **U1 status: `[PROVISIONAL: months]` — supported but NOT documented.** Eight sources searched; none labels the column (§1). This status is unchanged by anything in this section.

**Two classes of analysis, and the unit matters only to one:**

| Class | Examples | Depends on the unit? | Permitted under `[PROVISIONAL]`? |
|---|---|---|---|
| **RELATIVE-RISK** — rank- and ordering-based | hazard ratios · C-index · calibration slope on the log-hazard scale · likelihood-ratio tests · ΔC between models | **NO.** The Cox partial likelihood uses only the ordering of times, and any unit change is a strictly monotone rescaling | **YES — fully permitted.** No assumption about the unit is invoked, so none can be violated |
| **FIXED-TIME ABSOLUTE-RISK** — referencing a point on the time axis | t-year recurrence risk · predicted survival/CIF curves · cloglog calibration at a horizon · time-dependent AUC/Brier · the 60-unit endpoint | **YES.** "60 units" is 5 years, 1.15 years or 0.16 years depending on whether the unit is months, weeks or days | **NO — prohibited** while U1 is unresolved, except as explicitly labelled conditional analyses |

**Demonstrated above:** rescaling by ×30.44, ×1/12 or ×0.0329 leaves β, SE and C identical to 10 decimal places. The invariance is exact, not approximate.

**Consequence for the protocol:** the primary analysis is confined to the relative-risk class, which is why it can proceed with U1 unresolved. Every fixed-time analysis is either prohibited or carries an explicit conditional label. **This is a restriction on the claims, not a resolution of the unit.**

### 3.2 U2 — censoring assumptions and competing risks

#### 3.2.1 Two different censoring concepts — the distinction that matters

| | **Independent (non-informative) censoring** | **Censoring at a competing death** |
|---|---|---|
| **What it assumes** | a censored patient's future recurrence hazard is the same as that of comparable patients remaining under observation | the same assumption is *made*, but it is **counterfactual**: the patient cannot recur, ever |
| **Typical cause** | administrative end of follow-up, loss to follow-up | death from another cause |
| **Is the assumption testable?** | not directly, but plausible when censoring is by study design | **demonstrably false as a statement about the real world** — a dead patient has zero future recurrence hazard |
| **Does it invalidate the analysis?** | no, if it holds | **it does not invalidate the cause-specific hazard, but it changes what that hazard means** (§3.2.2) |

**In this cohort, censoring is of two distinct kinds:**

| Kind | n (S3) | Nature |
|---|---:|---|
| alive and recurrence-free at last contact | **271** | genuine administrative/follow-up censoring |
| **died without recurrence** | **46** | **competing event, censored as if it were follow-up censoring** |

And the competing-death group is **not exchangeable** with those remaining at risk: ~7 years older, less often stage III, roughly half as likely to have received chemotherapy (§2.3).

#### 3.2.2 What each estimator actually estimates

| | **Cause-specific Cox** (death censored) | **CIF / Fine–Gray** (death as competing) | **Composite DFS** (death counted as event) |
|---|---|---|---|
| **Quantity** | **cause-specific hazard** of recurrence | **cumulative incidence** of recurrence | hazard / risk of *recurrence-or-death* |
| **Plain statement** | *"Among patients still alive and recurrence-free, what is the instantaneous rate of recurrence?"* | *"What proportion of patients will actually experience recurrence, given that some will die first and can never recur?"* | *"What proportion will experience either event?"* |
| **Risk set** | death **removes the patient from the risk set** — they stop contributing, and are treated as though they *could* have recurred later | the patient **remains in the denominator** of the subdistribution; their death is recorded as a competing exit | death is an event; the patient exits by having the event |
| **Appropriate use** | **aetiology** — which features drive the recurrence mechanism | **prediction** — absolute probability a patient actually recurs | prognosis for a combined clinical outcome |
| **Relationship to 1−KM** | 1−KM **overstates** recurrence probability here, because it credits the 46 deaths with hypothetical future recurrences | CIF is the correct absolute-risk estimator | n/a — different endpoint |

**The precise point:** for the **cause-specific hazard**, death legitimately removes a patient from the recurrence risk set — that is exactly the quantity being defined, and the estimator is valid. But the resulting hazard **cannot be converted into an absolute recurrence probability via 1−KM**, because doing so implicitly assumes those 46 patients would have remained at risk. **The CIF estimates the actual recurrence probability in the presence of competing death; the cause-specific hazard does not.**

#### 3.2.3 Consequences by analysis type

| Analysis | Consequence of U2 |
|---|---|
| Cox relative risk | **valid as a cause-specific hazard ratio.** Not an absolute-risk statement |
| Absolute recurrence risk | **1−KM is biased upward.** Use CIF |
| Fixed-time predictions | require CIF, and are separately blocked by U1 |
| Calibration | predicted risks from a cause-specific model are not the quantity observed; calibrating them against observed CIF is an estimand mismatch unless stated |
| Event counts | 136 recurrences unchanged; **46 competing events** that a censoring framing discards as uninformative |
| Sensitivity analyses | the composite gives **182** events — a **different estimand**, not a larger version of the same one |

---

---

## 4. Primary endpoint recommendation

**No source value is changed and no exclusion is changed.** All variants are derived columns on the existing Phase 3 bundle.

### 4.1 Primary

> **Cause-specific recurrence, full follow-up, cohort S3.**
> `rfs_event_full` / `rfs_time_full`, death without recurrence treated as censoring (as deposited), analysed as a **cause-specific hazard model**.
> **Estimand: relative risk (hazard ratios) and discrimination only.**

**Why this one:**
- It is **unit-free**, so U1 cannot invalidate it (§3.1).
- It uses the endpoint **exactly as deposited** — no recoding.
- It is the largest unconditional event set (136).
- Its estimand is honestly statable despite U2: a cause-specific hazard, explicitly not an absolute-risk claim.

**Mandatory accompanying statement:** *"Hazard ratios are cause-specific. Death without recurrence (10.2% of the cohort) is treated as censoring; that censoring is associated with age, stage and treatment, so these estimates must not be read as absolute recurrence probabilities."*

### 4.2 Verified event accounting for the composite endpoint

**Checked, not assumed.** Full partition of S3 (n = 453), from `endpoint.csv`:

| Category | Definition | n |
|---|---|---:|
| **A** | recurrence recorded (`rfs_event = 1`) | **136** |
| **B** | died, no recorded recurrence (`rfs_event = 0`, `os_event = 1`) | **46** |
| **C** | alive, no recurrence (`rfs_event = 0`, `os_event = 0`) | **271** |
| **D** | `rfs_event = 0`, `os_event` missing | **0** |
| | **TOTAL** | **453** ✅ exhaustive |

| Check | Result |
|---|---|
| Missingness in `rfs_event`, `rfs_time`, `os_event`, `os_time` | **0, 0, 0, 0** |
| Categories mutually exclusive | **yes — no overlapping pairs** |
| Additional event categories beyond A–D | **none** |
| Patients with **both** recurrence and death | **89** |
| Of those, `rfs_time ≤ os_time` | **89 / 89** |
| **Double counting in the composite?** | **No** — the composite takes the *first* event, so these 89 count **once**, at recurrence |
| **Composite event count, computed** | **182** |
| 136 + 46 | **182** ✅ **confirmed** |
| Composite censored | **271** (= category C) |
| `rfs_time ≤ os_time` for all 453 | **yes**, so composite time = `min(rfs_time, os_time)` = `rfs_time_source` |
| Records where `rfs_time ≠ os_time` | 132 |

> **Confirmed: the 182 composite DFS events are exactly the 136 category-A recurrences plus the 46 category-B deaths without recorded recurrence.** No overlap, no missingness, no further event categories.

### 4.3 Sensitivity analyses — pre-specified

| # | Variant | Events | Estimand | Addresses |
|---|---|---:|---|---|
| **S-1** | CIF + Fine–Gray subdistribution | 136 rec / 46 competing | **absolute recurrence probability** | U2 |
| **S-2** | Composite DFS | **182** (136 + 46, verified) | **recurrence-or-death** | U2 |
| **S-3** | 60-unit administrative horizon | 130 | same as primary, truncated | **conditional on U1** |
| **S-4 … S-6** | Cohorts S0 / S1 / S2 | 139 / 136 / 139 | same as primary | zero-time rule |
| **S-7** | Representative-IQR representation | 136 | same as primary | aggregation choice |

### 4.4 How the three framings will be compared — and what comparison can and cannot show

> **The primary (cause-specific), S-1 (subdistribution) and S-2 (composite) estimate THREE DIFFERENT QUANTITIES.** They are not three measurements of one thing.
>
> **Agreement between them does NOT prove that death handling has no effect**, and will never be reported as if it did. Two estimands can coincide numerically in a given sample while remaining conceptually distinct, and the coincidence carries no inferential weight.

**What will be compared, and how:**

| Dimension | Comparison | Interpretation rule |
|---|---|---|
| **Effect direction** | sign of each coefficient across the three framings | **Concordant signs** = the predictor's association is not an artifact of the death-handling choice. **Discordant signs** are a substantive finding and are reported as such, not averaged away |
| **Effect magnitude** | cause-specific HR vs subdistribution HR vs composite HR, side by side with CIs | **Magnitudes are NOT expected to match** and no agreement threshold is defined. A predictor acting on both recurrence and mortality will have a subdistribution HR shrunk toward the null relative to its cause-specific HR — that is the expected behaviour, not a discrepancy |
| **Uncertainty** | CI width per framing | the composite has more events (182 vs 136) and will generally be **more precise**. **Greater precision is not greater validity** — it is a different question answered with more events |
| **Discrimination** | C-index per framing, each with CI | compared descriptively only; the outcomes being discriminated differ |
| **Interpretation** | each estimate reported with its estimand stated in words | no single "headline" number across framings |

**Reporting commitment:** all three are reported together, each labelled with its estimand in plain language, and **no framing is designated the "true" answer**. The primary is designated in advance (§4.1) on grounds of being unit-free and using the endpoint as deposited — **not** on grounds of being more correct.

**What would change the conclusion:** if the cause-specific and subdistribution analyses disagree in *direction* for a predictor, that predictor's apparent recurrence association is plausibly mediated by competing mortality, and it must be reported with that caveat rather than as a recurrence predictor.

### 4.3 Explicitly prohibited until U1 resolves

Absolute t-year risk · predicted survival curves presented as calibrated risk · time-dependent AUC/Brier at a named time · a faithful Riley cloglog calibration · **any statement of the form "5-year recurrence risk"**.

---

## 5. What stays unresolved

| # | Item | Status | Only route to resolution |
|---|---|---|---|
| U1 | Absolute unit of `rfs.delay` | `[PROVISIONAL: months]` — eight sources, none documents it | author / CIT correspondence |
| U2 | Published death-handling rule | **UNRESOLVED** | author correspondence |
| U2′ | *Methodological* handling | **RESOLVED by design** — competing event, addressed by S-1/S-2 | — |
| U4 | CIT446 stage/M conflict | unresolved; excluded by Criterion 1 | — |
| U5 | `os.delay = 0` sentinel vs value | unresolved; those records excluded by Criterion 2 | — |
