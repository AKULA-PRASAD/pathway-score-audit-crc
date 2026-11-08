# Phase 1 — Cohort and Endpoint Protocol (rev. 4, APPROVED)

**Date:** 2026-11-02 · **Revision 4** — approved with qualifications
**Status:** Approved. Decisions recorded below; **no source value altered, no dataset filtered.**
**Reproduce:** `.venv/bin/python scripts/15_investigate_cohort_endpoint.py` → `reports/phase1/`

## STATUS SUMMARY — read this first

Three categories, deliberately never merged. **A decision being approved does
not make an assumption underneath it true.**

### ✅ APPROVED DECISIONS

| # | Decision | Basis |
|---|---|---|
| A1 | **RFS = time from surgery to first recurrence** | `[DOCUMENTED]` — Text S1, quoted §1.1 |
| A2 | **Primary cohort = stage II–III primary tumours** | `[DOCUMENTED]` — authors' rationale, §4.1 |
| A3 | **Scenario S3 = proposed primary cohort** (Criterion 1 + Criterion 2 applied, kept separate) | §2.6; counts §7 |
| A4 | **S0, S1, S2 retained as sensitivity scenarios** | §2.5, §7 |
| A5 | **Full-follow-up analysis retained as sensitivity** | §1.4; originals never overwritten |
| A6 | **CIT subset = internal validation only, never external** | §4.2 |
| A7 | **Chemotherapy = prognostic covariate candidate, never a causal exposure** | §4.5 |
| A8 | **CIT446 flagged; stage and M status NOT recoded** | §3 |
| A9 | **Event count is a modelling constraint; no rigid 10-EPV rule** | §4.3 |
| A10 | **Nine-batch ArrayExpress assignments carried to Phase 4; no batch/normalization decision made** | §5.1 |

### ⚠️ PROVISIONAL ASSUMPTIONS — labelled, not facts

| # | Assumption | Evidence | If wrong |
|---|---|---|---|
| P1 | **Death without recurrence = censoring** | `[OBSERVED]` 49/49 pattern. **NOT a published rule** — "censor" appears zero times in Text S1 | Endpoint becomes DFS; events 136 → ~185 |
| P2 | **`rfs.delay` is in months** | `[OBSERVED]` 4 converging patterns. **No source labels the column** (6 searched, §1.3.1) | Every time-based result invalid |
| P3 | **60-unit horizon = 5 years** | **Conditional on P2**, not independent | Truncation at the wrong time; silent corruption |

**P2 remains the one blocking item.** P3 inherits its status entirely: the
60-month cutoff is used as a *working rule under a declared assumption*, never
cited as verified. Any artifact derived under it carries the label.

### ❓ UNRESOLVED LIMITATIONS — carried, not closed

| # | Limitation | Why it cannot be resolved |
|---|---|---|
| U1 | CIT446 stage II vs M1 | Table S1 has no M column; FOLFIRI circumstantial only (§3) |
| U2 | CIT106, CIT250 zero times | No positive evidence of baseline disease |
| U3 | `os.delay = 0`: true value or sentinel | No spike at zero; counts too small to exclude |
| U4 | R0 / curative resection status | Not stated in any source, not derivable |
| U5 | GEO vs Table S1 stage disagreement (4 records) | Outside the cohort; proves sources differ |
| U6 | Authors' handling of t = 0 records in fitted models | Denominator includes them; fitting unreported |
| U7 | Era effects, 1987–2007, seven centres | No per-sample surgery date exists |

---

### Evidence tiers used throughout

| Tier | Meaning |
|---|---|
| **[DOCUMENTED]** | Stated verbatim in a cited primary source, quoted below |
| **[OBSERVED]** | A data pattern measured in the deposited files, shown with counts |
| **[PROVISIONAL]** | Working assumption supported by [OBSERVED] evidence but **not** documented |
| **[UNRESOLVED]** | Cannot be determined from any source available to this investigation |

**A [PROVISIONAL] item is not a fact.** It is an assumption carried forward with
its evidence attached, and it must be re-labelled if documentation is found or
contradicted.

---

## 0. Sources

**Primary publication**
Marisa L, de Reyniès A, Duval A, Selves J, Gaub MP, Vescovo L, Etienne-Grimaldi MC,
Schiappa R, Guenot D, Ayadi M, Kirzin S, Chazal M, Fléjou JF, Benchimol D, Berger A,
Lagarde A, Pencreach E, Piard F, Elias D, Parc Y, Olschwang S, Milano G,
Laurent-Puig P, Boige V. *Gene expression classification of colon cancer into
molecular subtypes: characterization, validation, and prognostic value.*
PLoS Medicine 2013;10(5):e1001453.
DOI [10.1371/journal.pmed.1001453](https://doi.org/10.1371/journal.pmed.1001453)
· PMID 23700391 · PMCID PMC3660251

**Text S1** (Supplementary methods) — <https://doi.org/10.1371/journal.pmed.1001453.s016>
**Table S1** (Patient and Tumor Characteristics, **per-patient**, 566 CIT ids) —
<https://doi.org/10.1371/journal.pmed.1001453.s009>
Fetched to `data/external/` (gitignored), SHA-256 `38b6acb318022c27…`

**GEO** GSE39582 (SubSeries of GSE40967), GPL570 — local SOFT file.

**Methodological references** (for §4.3 only)
- Riley RD, Snell KIE, Ensor J, Burke DL, Harrell FE Jr, Moons KGM, Collins GS.
  *Minimum sample size for developing a multivariable prediction model: PART II —
  binary and time-to-event outcomes.* Stat Med 2019;38(7):1276–1296.
  DOI [10.1002/sim.7992](https://doi.org/10.1002/sim.7992)
- van Smeden M, de Groot JAH, Moons KGM, Collins GS, Altman DG, Eijkemans MJC,
  Reitsma JB. *No rationale for 1 variable per 10 events criterion for binary
  logistic regression analysis.* BMC Med Res Methodol 2016;16:163.

---

## 1. RFS endpoint semantics

### 1.1 Event definition — **[DOCUMENTED]**

> "Relapse-Free Survival was used and defined as the time from surgery to the
> first recurrence."
> — Text S1, *Survival Analyses*

**Retained as the endpoint definition.** The event is first recurrence; the
time origin is **date of surgery**. Both are documented.

### 1.2 Death without documented recurrence

**Two separate things, deliberately not merged:**

**(a) Published evidence — [UNRESOLVED].**
**The publication does not specify this rule.** I searched the main text twice
with different queries and extracted and searched the full Text S1. The string
"censor" appears **zero times in Text S1**. No sentence anywhere in either
source states how death without recurrence was handled.

**(b) Observed data pattern — [OBSERVED].**
Among stage II–III records with complete RFS and OS (n = 461):

| rfs.event | os.event | n | `rfs.delay == os.delay` |
|---:|---:|---:|---:|
| 0 | 0 | 273 | 273 (100%) |
| 0 | 1 | **49** | **49 (100%)** |
| 1 | 0 | 48 | 1 (2%) |
| 1 | 1 | 91 | 3 (3%) |

All 49 patients who died without recorded recurrence carry `rfs.event = 0`
with RFS time set exactly to their OS time. In the recurrence rows, by
contrast, 135 of 139 have `rfs.delay < os.delay`, so `rfs.delay` is a genuine
recurrence time and not a copy of the OS column.

**Conclusion — [PROVISIONAL]:** death without recurrence is treated as
**censoring at the date of death**. This is an inference from the data pattern
above, carried forward as a documented assumption. **It is not a published
rule and must not be cited as one.**

**Why it matters:** counting death as an event would take the event count from
**139 → 188 (+35%)** and change the estimand from RFS to DFS.

### 1.3 Time unit — **[PROVISIONAL]**, and a blocking dependency

**Direct documentary evidence (about the analysis, not the column):**

> "...assessed by Log-rank test using an endpoint of **five years/60 months**."
> — Text S1

This documents that the authors' *analysis horizon* was expressed in months.
It does **not** document the unit of the deposited `rfs.delay` column.

**Inferential evidence — [OBSERVED]**, 560 tumours with both times:

| Test | Result | Would fail under a unit mismatch? |
|---|---|---|
| `rfs.delay > os.delay` | **0 violations** | yes |
| range vs `os.delay` range | identical, 0–201 | yes |
| both integer-valued | true | likely |
| `rfs.delay == os.delay` | 387 (69.1%) | yes |

### 1.3.1 Documentary search — sources now exhausted

| Source | Result |
|---|---|
| GEO SOFT `!Sample_characteristics_ch1` | `rfs.delay` has **no** unit label; `os.delay` **is** labelled `(months)` |
| Publication main text | no unit statement for the column |
| Text S1 supplementary methods | "five years/60 months" — the **analysis horizon**, not the column |
| **Table S1** per-patient annotation | single sheet, no legend, no units; `RFS.delay` unlabelled |
| **ArrayExpress E-GEOD-39582 SDRF** | `Characteristics [rfs.delay]` present; **no `Unit[]` column anywhere** in the SDRF; no overall-survival column at all |
| **GEO SuperSeries GSE40967** | no RFS metadata whatsoever |

**Outcome: no source labels the unit.** Status remains **[PROVISIONAL]**.
The one avenue not attempted is direct correspondence with the CIT programme
or the corresponding author, which is outside what this investigation can do
unilaterally.

### 1.3.2 The 60-month cutoff is *conditional* on this unit — blocking

**A "60" cutoff is only a 5-year cutoff if `rfs.delay` is in months.** The
censoring rule in §1.4 and the endpoint definition are therefore **not
independent decisions**: the cutoff inherits the provisional status of the
unit.

Concretely, if the column were in weeks, "60" would truncate at ~14 months and
reclassify a large fraction of events; if in days, at ~2 months. Either would
silently destroy the analysis while producing plausible-looking output.

**Therefore: the 60-month cutoff cannot be confirmed as the primary analysis
rule until the time unit is established.** It is carried forward as a
*conditional* proposal — "IF months, THEN truncate at 60" — and the condition
must be discharged, not assumed. Two acceptable ways to discharge it:

1. **Documentary** — obtain an explicit statement (author correspondence, CIT
   data dictionary). Upgrades to `[DOCUMENTED]`.
2. **Declared assumption** — proceed with months as a stated, pre-registered
   assumption, and run the full pipeline under an alternative-unit sensitivity
   analysis to demonstrate the conclusions do not hinge on it.

Until one of these is done, §1.4 is a **proposal with an unmet precondition**,
not an approved rule.

### 1.4 Five-year administrative censoring — **conditional proposal**

> **PRECONDITION (§1.3.2): this rule presumes `rfs.delay` is in months, which
> is [PROVISIONAL]. It cannot be confirmed as the primary analysis rule until
> that unit is established.**

> "...using an endpoint of five years/60 months." — Text S1

**The deposited data is NOT censored at 5 years** — [OBSERVED]. 180 of 464
stage II–III records (38.8%) exceed 60 months, up to 201 months (16.8 years).
The authors applied the horizon at analysis time; GEO holds raw follow-up.

**Exact mechanics of the proposed primary analysis.** For each record, derive
new columns — originals are never overwritten:

```
if rfs.delay > 60:
    rfs_time_5y  := 60        # time truncated to the horizon
    rfs_event_5y := 0         # RECLASSIFIED AS CENSORED, whatever it was
else:
    rfs_time_5y  := rfs.delay # unchanged
    rfs_event_5y := rfs.event # unchanged
```

A recurrence observed at month 75 becomes a **censored observation at month
60**. This is standard administrative censoring and it is *not* data deletion:
the patient still contributes 60 months of event-free person-time.

**Exact numerical consequences — [OBSERVED]:**

| subset | n | events before | events after | events reclassified as censored | times truncated to 60 |
|---|---:|---:|---:|---:|---:|
| discovery | 359 | 109 | **103** | 6 | 134 |
| validation | 102 | 30 | **30** | 0 | 46 |

### 1.4.1 Record-level verification of the six reclassified events

All six are **discovery**, all carry `rfs.event = 1`, all have
`rfs.delay > 60`, and all satisfy `rfs.delay <= os.delay`:

| sample_id | title | stage | rfs.event | rfs.delay | os.event | os.delay | chemo |
|---|---|---|---:|---:|---:|---:|---|
| GSM972136 | CIT180 | III | 1 | 61 | 0 | 130 | Y |
| GSM972322 | CIT366 | II | 1 | 66 | 1 | 86 | Y |
| GSM972352 | CIT396 | III | 1 | 74 | 0 | 139 | Y |
| GSM972328 | CIT372 | III | 1 | 78 | 0 | 82 | Y |
| GSM971968 | CIT012 | II | 1 | 93 | 0 | 112 | N |
| GSM972082 | CIT126 | II | 1 | 119 | 1 | 132 | Y |

Each would become a **censored observation at 60**, retaining 60 units of
event-free person-time. Note GSM972136 recurs at 61 — one unit past the
boundary. Full records: `reports/phase1/five_year_reclassified.csv`.

### 1.4.2 Why validation reclassifies ZERO despite 46 records beyond 60

This looked anomalous and is not. Truncation changes **time** for every record
beyond the horizon, but changes **event status** only for records that were
*events*.

| subset | records >60 | of those, events | of those, censored | latest event time | latest censored time |
|---|---:|---:|---:|---:|---:|
| discovery | 134 | **6** | 128 | 119 | 201 |
| validation | 46 | **0** | 46 | **45** | 146 |

**Validation's latest recurrence occurs at month 45** — comfortably inside the
horizon. Every one of its 46 records beyond 60 months is therefore a *censored*
observation whose event status was already 0. Truncating their time to 60 is a
no-op for the event count, though it does shorten their contributed
person-time. Discovery's events extend to month 119, so it has six beyond the
boundary.

A consequence worth carrying forward: **the two subsets have materially
different follow-up profiles** — late recurrences exist in discovery and not
in validation — which is one more reason the CIT subset is not exchangeable
with discovery (§4.2).

**Sensitivity analysis:** the identical pipeline on full uncensored follow-up
(139 events). If a finding holds under the 5-year horizon but vanishes with
full follow-up — or the reverse — that is itself a reportable result about
the time-dependence of the signal, not a nuisance.

### 1.5 Cohort inclusion criteria — **[DOCUMENTED]**

> "750 patients with stage I to IV CC who underwent surgery between 1987 and
> 2007 in seven centers. Patients who received preoperative chemotherapy
> and/or radiation therapy and those with primary rectal cancer were excluded."

Median follow-up **51.5 months**.

**[UNRESOLVED]:** no statement on R0 resection or curative intent for
stage II–III. The phrase "curative surgery" appears only in the rationale for
*excluding* stage I. Whether R1/R2 resections are present cannot be determined.

---

## 2. Zero-time records — re-investigation

**No observation was shifted, recoded or removed. The 0.5-month shift
previously proposed has been withdrawn.** All original values preserved.

### 2.1 Source tracing — **[OBSERVED]**

I traced all eight records back to the raw SOFT `!Sample_characteristics_ch1`
lines and forward to the authors' per-patient Table S1.

**Result: all three sources agree exactly.** The GEO CSV faithfully reflects
the SOFT deposit (no parsing artefact), and across all 566 tumours there are
**0 disagreements** between GEO and Table S1 on `rfs.event` and `rfs.delay`.

| sample_id | title | subset | stage | T | N | M | rfs.event | rfs.delay | os.event | os.delay | chemo | chemo type |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSM972062 | CIT106 | discovery | II | T4 | N0 | M0 | **1** | 0 | 1 | 23 | Y | FOLFOX |
| GSM972206 | CIT250 | discovery | II | N/A | N/A | N/A | **1** | 0 | 0 | 42 | N | — |
| GSM972402 | CIT446 | validation | II | T3 | N0 | **M1** | **1** | 0 | 1 | 31 | Y | **FOLFIRI** |
| GSM972229 | CIT273 | discovery | II | T3 | N0 | M0 | 0 | 0 | 0 | **0** | N | — |
| GSM972246 | CIT290 | discovery | III | T3 | N1 | M0 | 0 | 0 | 0 | **0** | N | — |
| GSM972258 | CIT302 | discovery | II | T4 | N0 | M0 | 0 | 0 | 1 | **0** | N | — |
| GSM972401 | CIT445 | validation | II | T3 | N0 | M0 | 0 | 0 | 1 | **0** | N | — |
| GSM972456 | CIT500 | validation | II | T3 | N0 | M0 | 0 | 0 | 1 | **0** | N | — |

### 2.2 Pattern A — `rfs.event = 1, rfs.delay = 0` (n = 3)

**New evidence found.** This pattern is not an anomaly in this dataset — it is
the dominant coding for **stage IV**:

| | n |
|---|---:|
| stage IV tumours | 60 |
| stage IV with `rfs.delay = 0` | 30 (50%) |
| of those, `rfs.event = 1` | **29 / 30** |

Stage IV patients are metastatic at baseline — there is nothing to relapse
*from*. The CIT coders evidently recorded "disease present at time zero" as
`rfs.delay = 0, rfs.event = 1`. **[OBSERVED]**

**Therefore the supported interpretation of Pattern A is "disease present at
baseline", not "a recurrence observed during month 0."** My earlier
"synchronous metastasis / rounding / data-entry error" list was a set of
untested hypotheses; this convention is the one with actual evidence behind it.

**But the evidence is per-record, and it is not uniform:**

| record | supports baseline disease | contradicts | strength |
|---|---|---|---|
| **CIT446** | `tnm.m = M1`; FOLFIRI (a metastatic regimen, not standard stage II adjuvant therapy) | `tnm.stage = 2` in **both** GEO and Table S1 | **moderate–strong** |
| **CIT106** | fits the t=0 convention | `tnm.m = M0`; FOLFOX is plausible for high-risk T4 stage II | **weak — unexplained** |
| **CIT250** | fits the t=0 convention | T/N/M all `N/A`; no chemo; alive at 42 months | **weak — unexplained** |

For CIT106 and CIT250 there is **no positive evidence** of baseline metastatic
disease. Their `rfs.delay = 0` is **[UNRESOLVED]**.

### 2.3 Pattern B — `rfs.event = 0, rfs.delay = 0, os.delay = 0` (n = 5)

These have **zero recorded follow-up on both endpoints**. Three have
`os.event = 1` (death recorded at month 0); two have `os.event = 0`.

**Is `0` a true value or a missing-data sentinel? — [UNRESOLVED].** I tested
it: across all 566 tumours, `os.delay` takes the value 0 in 6 records, 1 in 6,
2 in 6, 3 in 1. **There is no spike at zero**, so the data do *not* support a
sentinel reading — but with counts this small they do not exclude one either.
I could not resolve this, and I am not going to assert it either way.

What *is* certain: **these five records contribute zero person-time**. Under
any survival model they are uninformative, whichever interpretation is true.

### 2.4 How the implementation actually behaves at t = 0 — **[OBSERVED]**

Before choosing a rule, I probed what `lifelines 0.30.3` does with a
zero-duration row, using **synthetic data only** (no project data touched):

| Scenario | `CoxPHFitter` | `KaplanMeierFitter` |
|---|---|---|
| event at t = 0 | **fits, n_used = all rows, NO warning** | fits; timeline starts at 0; **S(0) = 0.995 < 1** |
| censored at t = 0 | **fits, n_used = all rows, NO warning** | fits; S(0) = 1.0 |

**The record is not dropped and not flagged.** It enters the partial
likelihood. In a probe where the t = 0 event carried an extreme covariate, the
fitted coefficient moved from −0.0652 to −0.0108 — a single row shifting the
estimate by 83%.

**Why that happens:** an event at the time origin has the **entire cohort** in
its risk set, which is the largest possible risk set, so it is maximally
influential. And a survival curve that starts below 1 is incoherent — every
patient is by construction recurrence-free at the moment of surgery.

**This is the dangerous case precisely because it is silent.** "The software
will catch it" is false.

### 2.5 The four scenarios — **[OBSERVED]** counts, nothing applied

| Scenario | n | events | censored | disc n / events | val n / events |
|---|---:|---:|---:|---:|---:|
| **S0** retain all 8 | 461 | 139 | 322 | 359 / 109 | 102 / 30 |
| **S1** exclude Group A only (3) | 458 | 136 | 322 | 357 / 107 | 101 / 29 |
| **S2** exclude Group B only (5) | 456 | 139 | 317 | 356 / 109 | 100 / 30 |
| **S3** exclude both (8) | 453 | 136 | 317 | 354 / 107 | 99 / 29 |

### 2.6 Recommended primary approach — **two explicit criteria, not one**

My revision-2 recommendation ("exclude all 8 under a single `rfs.delay > 0`
rule") was convenient but poorly justified: it bundled two unrelated problems
under one arbitrary-looking threshold. **S3 remains the recommendation, but it
should be stated as two separate, independently defensible eligibility
criteria** — because the groups fail for different reasons, and an external
cohort may contain one without the other.

**Criterion 1 — "recurrence-free at the time origin" (removes Group A, n = 3).**
Not a technical fix. The endpoint is *recurrence-free* survival; a patient
recorded as having recurred at the moment of surgery was **never
recurrence-free**, so was never at risk and does not belong in the risk set.
The stage IV coding evidence (29/30) supports reading these as *disease
present at baseline* rather than month-0 recurrences. This criterion is
conceptual and transfers to any cohort.

**Criterion 2 — "> 0 units of recorded follow-up" (removes Group B, n = 5).**
A standard follow-up eligibility requirement. These five have zero recorded
follow-up on **both** endpoints and contribute zero person-time under every
interpretation, so the unresolved sentinel question (§2.3) never has to be
settled. They cannot inform any model.

**Combined cost: 8 of 461 records (1.7%), 3 of 139 events (2.2%).**

**Why not the alternatives:**

- **S0 (retain all)** — rejected as *primary* because of §2.4: the pathology is
  silent, the three t = 0 events are maximally influential, and S(0) < 1 is
  incoherent. It remains a **required sensitivity analysis** precisely so the
  influence is measured rather than assumed.
- **S1 (Group A only)** — leaves five zero-follow-up records that cannot
  contribute, inflating the apparent denominator.
- **S2 (Group B only)** — leaves the three t = 0 events, i.e. keeps the actual
  statistical hazard while removing the harmless rows.

**If retention is nonetheless chosen (S0 or S2), the handling is specified as:**
keep the source value `0` **unchanged** — no offset, no imputation, no
recoding — and in the modelling layer either (a) use a left-truncation /
delayed-entry specification with entry time 0 and explicit tie handling, or
(b) report per-record influence diagnostics (dfbeta) for the three t = 0
events so their leverage is visible rather than silent. **No arbitrary offset
such as 0.5 months will be introduced under any scenario.** The 0.5-month
shift proposed in revision 1 is withdrawn permanently.

**All four scenarios must be pre-specified now**, before any model is fitted.

---

## 3. GSM972402 / CIT446 — stage II vs M1

**Investigated as requested. Not resolvable. Not recoded.**

**Original source fields traced.** From the raw SOFT record:

```
TITLE: CIT446
tnm.stage: 2          tnm.t: T3      tnm.n: N0      tnm.m: M1
chemotherapy.adjuvant: Y            chemotherapy.adjuvant.type: FOLFIRI
rfs.event: 1          rfs.delay: 0
os.event: 1           os.delay (months): 31
```

The GEO CSV is a faithful transcription — **the conflict exists in the deposit
itself**, not in our processing.

**Can Table S1 arbitrate? No.** The authors' per-patient table carries
`TNM.Stage` but **no `TNM.T` / `TNM.N` / `TNM.M` columns and no OS columns**.
It records `TNM.Stage = 2` for CIT446, agreeing with GEO's stage field — which
is the side of the conflict that was never in doubt. **The M field has no
independent corroboration anywhere.**

| Evidence for **M1 / stage IV** | Evidence for **stage II** |
|---|---|
| `tnm.m = M1` in the GEO deposit | `tnm.stage = 2` in GEO |
| `rfs.delay = 0, rfs.event = 1` matches the stage IV convention (29/30) | `TNM.Stage = 2` in the authors' Table S1 |
| FOLFIRI — *weak, circumstantial only* (see below) | T3 N0 is consistent with stage II if M0 |

**On FOLFIRI — downgraded, and explicitly not treated as definitive.** My
revision-2 text called FOLFIRI "a metastatic regimen", implying it settles the
question. It does not, and that inference was overstated. Reasons it is weak:

- Chemotherapy regimen is **not a staging instrument**. Regimen choice depends
  on clinician preference, trial enrolment, centre protocol, comorbidity and
  era — this cohort spans surgeries from 1987 to 2007 across seven centres.
- The GEO field is `chemotherapy.adjuvant.type` with
  `chemotherapy.adjuvant = Y`, i.e. it is recorded as **adjuvant** therapy.
  Adjuvant intent is inconsistent with metastatic disease, so this field
  arguably argues *against* M1 as much as the regimen argues for it.
- Irinotecan-containing regimens were investigated in the adjuvant setting
  during this period, so adjuvant FOLFIRI is not *prima facie* impossible.

FOLFIRI is therefore recorded as a **circumstantial observation, carrying no
decisive weight**, not as evidence of stage.

**Verdict: [UNRESOLVED].** Two primary sources say stage II; two independent
clinical signals (M1, FOLFIRI) say metastatic. I cannot determine which field
is in error.

**Handling: flag, never recode.** The source values stay exactly as deposited.

CIT446 is removed from the primary analysis under **Criterion 1** of §2.6
(recurrence recorded at the time origin) — *as a consequence of the zero-time
rule, not as a staging judgement*. This keeps the two issues strictly
separate: we are not deciding its stage, and we are not using the zero-time
rule as a backdoor for a staging decision.

**Explicit CIT446 sensitivity analysis — [OBSERVED] counts:**

| | n | events |
|---|---:|---:|
| cohort **with** CIT446 | 461 | 139 |
| cohort **without** CIT446 | 460 | 138 |

Because CIT446 is already removed by S1 and S3, the sensitivity analysis that
matters is **within S0 and S2** (the retention scenarios), run both with and
without it. If any conclusion changes on the strength of this single
unresolvable record, that fact is itself the finding and must be reported.

The stage/M conflict goes in the limitations section regardless of which
scenario is chosen.

### 3.1 New finding — GEO and Table S1 disagree on four stage values

| title | sample_id | GEO `tnm.stage` | Table S1 `TNM.Stage` |
|---|---|---:|---:|
| CIT029 | GSM971985 | 0 | 1 |
| CIT317 | GSM972273 | 0 | 1 |
| CIT340 | GSM972296 | 0 | 1 |
| CIT344 | GSM972300 | 0 | 1 |

These are the four records GEO codes as "stage 0". **None is in the stage II–III
cohort, so the primary analysis is unaffected** — but it establishes that the
two sources are not identical, which is why the M1 conflict could not simply be
settled by consulting Table S1.

---

## 4. Cohort, treatment and covariates

### 4.1 Stage II–III as primary cohort — **retained** — **[DOCUMENTED]** rationale

> "Survival analyses were restricted to the subgroup of patients with stage
> II-III tumors. Additional prognostic biomarkers are most needed for these
> patients. This is because the vast majority of stage I CRC patients will
> never relapse after curative surgery and will not derive benefit from
> adjuvant chemotherapy because the prognosis is excellent. Also, most stage
> IV CRC patients are already metastatic and will die from their disease."
> — Text S1

This is the authors' own clinical-utility argument and it applies directly to
our question: stage II–III is where the adjuvant decision is genuinely
uncertain, so incremental predictive value has real decision relevance. It
also keeps us comparable to the source study and to external cohorts using the
same convention.

**Cost:** 97 tumours excluded; 139 events available.

```
deposited samples                          585
primary tumours                            566   (-19 non-tumoural mucosa)
stage II-III primary tumours               469   (-97)
with documented RFS (event and time)       461   (-8)
  of those, rfs.delay > 0                  453   (-8)
```

Tumour selection is unambiguous — exactly three `source_name` values exist
(443 discovery / 123 validation / 19 non-tumoural mucosa).

**Paired normals:** 19 non-tumoural samples, **17 carrying a `dependancy_sample`
link that resolves to a tumour in the series** — i.e. 17 same-patient
tumour/normal pairs. A tumour-only analysis is unaffected, but the normals are
not independent observations and must never be pooled with tumours.

### 4.2 CIT "validation" is **internal**, not external — **retained**

1. **Same study, same consortium.** `!Series_overall_design`: "The 585 tumors
   was divided into a discovery dataset of 443 CC and a validation dataset of
   123 CC" — one collection, split.
2. **Jointly preprocessed.** `!Sample_data_processing`:
   `RMA normalisation (Bioconductor R affy package) + ComBat (R sva package)`
   — fitted across all 585 arrays together, so the subsets are statistically
   entangled.
3. **Technically confounded, with the evidence removed.** RNA extraction method
   differs sharply (Trizol 22.7% vs 4.7%; RNA NOW 0.3% vs 15.0%; Rneasy Micro
   12.2% vs 37.4%) — exactly what ComBat was applied to erase.
4. **Underpowered regardless.** 30 events.

Note the publication's own "validation set (n = 1,029)" is the 123 CIT samples
**plus 906 public samples** — so the GEO `validation` column is not the paper's
validation set either.

**The CIT subset may serve as an internal holdout at most. It must never be
described as external validation and cannot support the transportability
claim.** Genuine external cohorts remain a Phase 4 requirement.

### 4.3 Sample size as a modelling constraint — not a 10-EPV rule

My previous framing ("139 events → ~13 parameters at 10 EPV") overstated a
rule of thumb as a limit. Revised:

**The 10-events-per-variable rule has been directly challenged.** van Smeden
et al. (2016) found EPV has only a weak relation to prediction error and
concluded it is not an appropriate sample-size criterion. Riley et al. (2019)
replaced the single ratio with targets tied to expected model performance —
anticipated shrinkage, optimism in apparent fit, and precision of the
predicted risks.

**What this means operationally:**

- **139 events is the binding resource** (103 under the 5-year horizon). That
  is a real constraint on model complexity, but the number of admissible
  parameters is **not** a fixed 13.
- The relevant question is **expected optimism and shrinkage**, not a ratio.
  With penalisation (ridge/elastic net), effective degrees of freedom are far
  below the nominal parameter count, so more candidate predictors can be
  entertained than an EPV count suggests.
- A formal Riley-style sample-size calculation should be run in **Phase 3**,
  once the candidate predictor set and a target C-statistic are defined. It
  cannot be run now, because it needs both.
- Regardless of the rule, **feature selection must be nested inside
  cross-validation**. Selecting predictors on the full data before CV inflates
  apparent performance no matter how favourable the event count.

### 4.4 Covariate availability (stage II–III, documented RFS, n = 461)

**No variable is proposed for inclusion here.** Availability audit only.

| covariate | complete | missing | missing % | levels | events w/ data |
|---|---:|---:|---:|---:|---:|
| age at diagnosis | 461 | 0 | 0.0% | cont. | 139 |
| sex | 461 | 0 | 0.0% | 2 | 139 |
| tnm.stage | 461 | 0 | 0.0% | 2 | 139 |
| tumour location | 461 | 0 | 0.0% | 2 | 139 |
| chemotherapy.adjuvant | 460 | 1 | 0.2% | 2 | 139 |
| tnm.t / tnm.n / tnm.m | 441 | 20 | 4.3% | 4/5/3 | 119 |
| kras.mutation | 440 | 21 | 4.6% | 2 | 131 |
| mmr.status | 417 | 44 | 9.5% | 2 | 131 |
| braf.mutation | 407 | 54 | 11.7% | 2 | 112 |
| cimp.status | 394 | 67 | 14.5% | 2 | 105 |
| cin.status | 372 | 89 | 19.3% | 2 | 94 |
| tp53.mutation | 288 | 173 | 37.5% | 2 | 101 |
| chemotherapy.adjuvant.type | 143 | 318 | **69.0%** | 5 | 51 |
| cit.molecularsubtype | 461 | 0 | 0.0% | 6 | 139 |
| rna.extraction.method | 461 | 0 | 0.0% | 6 | 139 |

**Must not enter a clinical baseline:**
- **`cit.molecularsubtype`** — derived *from the expression data*. Including it
  would absorb the very signal we measure incremental value against.
- **`rna.extraction.method`** — technical. A batch *diagnostic*, never a predictor.
- **`chemotherapy.adjuvant.type`** — 69% missing, 5 levels, 51 events.

**Stage-specific risk:** stage II 59/260 (22.7%), stage III 80/201 (39.8%) —
nearly two-fold, so stage belongs in any baseline model.

### 4.5 Adjuvant chemotherapy — an association, explicitly not a treatment effect

| chemo | stage | n | events | event rate |
|---|---|---:|---:|---:|
| N | II | 203 | 40 | 19.7% |
| Y | II | 56 | 19 | 33.9% |
| N | III | 55 | 22 | 40.0% |
| Y | III | 146 | 58 | 39.7% |

**This table does not estimate a treatment effect and must not be read as one.**

Within stage II, treated patients recur more often than untreated. That ordering
is the signature of **confounding by indication**: in this era, stage II
patients were selected for adjuvant chemotherapy *because* they carried
high-risk features (T4, obstruction, perforation, inadequate node sampling,
poor differentiation). Treatment assignment is a function of prognosis, so the
crude contrast compares sicker treated patients against healthier untreated
ones.

**Three reasons no causal reading is available here**, and none is fixable
within this dataset:

1. **Observational, non-randomised** allocation.
2. **Unmeasured confounding by the very features that drove the decision** —
   obstruction, perforation, lymph-node yield and differentiation grade are
   **not in the GEO metadata at all**, so they cannot be adjusted for.
3. **Era effects** — surgery spans 1987–2007 across seven centres, and adjuvant
   practice changed substantially over that period. No per-sample surgery date
   is available, so era cannot be modelled.

**What chemotherapy is for in this project:** it is a **prognostic covariate in
a prediction model**, not an exposure in a causal model. We are asking whether
transcriptomic features add predictive value beyond routinely available
clinical information, and treatment received is part of that information. Its
coefficient must never be interpreted as the effect of giving chemotherapy.

Stage stratification is the minimum response; whether to stratify, interact, or
restrict is a **Phase 3 decision** and is not made here.

---

## 5. Unresolved questions

| # | Question | Status | Blocking? |
|---|---|---|---|
| 1 | `rfs.delay` units | `[PROVISIONAL: months]` — six documentary sources exhausted (§1.3.1) | **YES — blocks confirming the 60-month cutoff (§1.3.2)** |
| 2 | Death without recurrence | `[PROVISIONAL: censoring]` — 49/49 observed; no published rule | No — assumption is stated and its effect quantified |
| 3 | CIT446 stage II vs M1 | `[UNRESOLVED]` — Table S1 has no M column; FOLFIRI is circumstantial only | No — removed by Criterion 1, with explicit sensitivity |
| 4 | CIT106, CIT250 zero-time | `[UNRESOLVED]` — no positive evidence of baseline disease | No — removed by Criterion 1 |
| 5 | `os.delay = 0`: true value or sentinel? | `[UNRESOLVED]` — no spike at zero, counts too small to exclude | No — Criterion 2 removes them either way |
| 6 | R0 / curative resection status | `[UNRESOLVED]` — not stated, not derivable | No — limitation |
| 7 | GEO vs Table S1 stage disagreement (4 records) | `[OBSERVED]` — stage 0 vs 1; outside the cohort | No — but proves sources are not identical |
| 8 | Did the authors fit models including t = 0 records? | `[UNRESOLVED]` — denominator 359 includes them; fitted handling unreported | No |
| 9 | 1987–2007 surgery window, 7 centres, no per-sample date | `[UNRESOLVED]` — era effects unmodellable | No — limitation |

**Only item 1 is blocking.** Everything else is either resolved by a
pre-specified rule or carried as a documented limitation.

### 5.1 Carried forward to Phase 4

The ArrayExpress SDRF for E-GEOD-39582 exposes **the actual ComBat batch
assignments**, which the GEO deposit does not:

| ComBat batch | discovery | validation |
|---|---:|---:|
| batch01 | 60 | 1 |
| batch02 | 14 | 0 |
| batch03 | 102 | 10 |
| batch04 | 87 | 32 |
| batch05 | 81 | 5 |
| batch07 | 49 | 47 |
| batch08 | 39 | 5 |
| batch09 | 1 | 18 |
| batch10 | 10 | 5 |

The batch structure is sharply unbalanced across the two subsets. This is
direct documentary evidence for the entanglement argued in §4.2, and it means
Phase 4 can **model** the batch structure rather than infer it. Recorded here;
no action taken in Phase 1.

## 6. Reproducibility record

### 6.1 Installed dependencies

| Package | Version | Phase | Purpose |
|---|---|---|---|
| `pyarrow` | 26.0.0 | 0 | Parquet backend for the probe matrix cache |
| `pytest` | 9.2.0 | 0 | Phase 0 verification suite |
| `iniconfig` | 2.3.1 | 0 | transitive (pytest) |
| `pluggy` | 1.7.0 | 0 | transitive (pytest) |
| **`xlrd`** | **2.0.2** | **1** | **Read-only parser for the legacy `.xls` Table S1 (per-patient annotation, Marisa et al. 2013). Pure Python, no native deps. Used only by `scripts/15_investigate_cohort_endpoint.py` to cross-check GEO against the published table. Not imported by the Phase 0 suite.** |

Pre-existing in the environment and unchanged: numpy 2.5.4, pandas 2.3.3,
scipy 1.18.1, scikit-learn 1.9.1, lifelines 0.30.3, GEOparse 2.0.4
(Python 3.13.13).

**Deferred, not installed:** `statsmodels` 0.15.0 and `scikit-survival` 0.28.0
(Phase 5+). Both verified installable by `pip install --dry-run`;
scikit-survival has no prebuilt wheel for this Python/numpy pair and compiles
from source.

`requirements.txt` is the original environment capture and remains unmodified.
Additions are recorded in `requirements-analysis.txt`.

### 6.2 Integrity

| Check | Result |
|---|---|
| Phase 0 test suite | **21/21 passed** |
| Raw data SHA-256 (`data/raw/CHECKSUMS.sha256`) | both files **OK** |
| Probe matrix cache | **current** (SHA-256 matches source) |
| Scripts 01–14 | **unmodified** |
| `metadata/` source values | **byte-identical** (asserted in `scripts/16_build_cohort_manifest.py`) |
| `config/analysis_config.yaml` | **no scientific parameter present** |

### 6.3 Reproduce

```
.venv/bin/python scripts/15_investigate_cohort_endpoint.py   # evidence
.venv/bin/python scripts/16_build_cohort_manifest.py         # cohort labels
```

Outputs land in `reports/phase1/`. Neither script writes to `metadata/`,
`data/raw/`, `data/processed/`, or `config/`.

---

## 7. Final cohort flow and counts

### 7.1 Flow to the approved primary cohort (S3)

```
deposited samples                             585
primary tumours                               566   (-19  non-tumoural mucosa)
stage II-III primary tumours                  469   (-97  stage 0/I/IV, missing)
with documented RFS (event and time)          461   (-8   no rfs.event/delay)
minus Criterion 1 failures (Group A, n=3)     458   (-3   not recurrence-free at origin)
minus Criterion 2 failures (Group B, n=5)     453   (-5   zero recorded follow-up)
                                              ===
APPROVED PRIMARY COHORT (S3)                  453
```

The two criteria are applied and reported **separately**, never collapsed into
one threshold.

### 7.2 Exact counts, all scenarios

`events_5y` is **conditional on provisional assumption P2** (months).

| scenario | description | n | events (full) | censored | events (5y)* | median FU |
|---|---|---:|---:|---:|---:|---:|
| **S3** | **PRIMARY** — C1 + C2 | **453** | **136** | **317** | **130** | 49 |
| S0 | sensitivity: retain all 8 | 461 | 139 | 322 | 133 | 48 |
| S1 | sensitivity: C1 only | 458 | 136 | 322 | 130 | 48 |
| S2 | sensitivity: C2 only | 456 | 139 | 317 | 133 | 48.5 |

By subset:

| scenario | disc n | disc ev (full) | disc ev (5y)* | val n | val ev (full) | val ev (5y)* |
|---|---:|---:|---:|---:|---:|---:|
| **S3** | **354** | **107** | **101** | **99** | **29** | **29** |
| S0 | 359 | 109 | 103 | 102 | 30 | 30 |
| S1 | 357 | 107 | 101 | 101 | 29 | 29 |
| S2 | 356 | 109 | 103 | 100 | 30 | 30 |

\* conditional on P2

### 7.3 The eight excluded records — values verbatim, nothing recoded

| sample_id | title | subset | stage | M | rfs.event | rfs.delay | os.event | os.delay | criterion |
|---|---|---|---|---|---:|---:|---:|---:|---|
| GSM972062 | CIT106 | discovery | II | M0 | 1 | 0 | 1 | 23 | **C1** (Group A) |
| GSM972206 | CIT250 | discovery | II | — | 1 | 0 | 0 | 42 | **C1** (Group A) |
| GSM972402 | CIT446 | validation | II | **M1** | 1 | 0 | 1 | 31 | **C1** (Group A) |
| GSM972229 | CIT273 | discovery | II | M0 | 0 | 0 | 0 | 0 | **C2** (Group B) |
| GSM972246 | CIT290 | discovery | III | M0 | 0 | 0 | 0 | 0 | **C2** (Group B) |
| GSM972258 | CIT302 | discovery | II | M0 | 0 | 0 | 1 | 0 | **C2** (Group B) |
| GSM972401 | CIT445 | validation | II | M0 | 0 | 0 | 1 | 0 | **C2** (Group B) |
| GSM972456 | CIT500 | validation | II | M0 | 0 | 0 | 1 | 0 | **C2** (Group B) |

**CIT446 is excluded by Criterion 1 only.** Its stage/M conflict is neither
recoded nor used as a reason for exclusion — the two issues stay separate.

### 7.4 Artifacts

| File | Contents |
|---|---|
| `reports/phase1/cohort_manifest.csv` | per-sample eligibility flags and scenario membership for all 585 samples; source values copied verbatim; derived `rfs_time_5y` / `rfs_event_5y` as **separate** columns |
| `reports/phase1/cohort_counts.csv` | the tables in §7.2 |
| `reports/phase1/excluded_records.csv` | the eight records in §7.3 |
| `reports/phase1/cohort_summary.json` | machine-readable summary |

The manifest is a **label table**, not a filtered dataset. No analysis dataset
has been materialised.

---

## 8. Status

Phase 1 is complete and approved with qualifications. Approved decisions,
provisional assumptions and unresolved limitations are separated in the
**Status Summary** above.

`config/analysis_config.yaml` still contains **no** endpoint, eligibility,
cohort, filtering, aggregation or normalization parameter. Scripts 01–14, raw
data, the cached matrix and `metadata/` are unchanged.

**Not started:** survival modelling, feature selection, probe aggregation,
expression filtering, normalization strategy, external validation.
