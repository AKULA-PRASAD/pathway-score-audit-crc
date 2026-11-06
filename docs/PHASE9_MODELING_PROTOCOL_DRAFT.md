# Phase 9 — Modeling Protocol Draft

**Version:** 3.1 · **Status:** **ALL FIVE DECISIONS APPROVED — awaiting execution authorisation** · **Date:** 2026-11-03
**Script:** `scripts/36_phase9_protocol_plan.py` · **Tests:** `tests/test_phase9_protocol.py`
**Artifacts:** `reports/phase9/decision_register.json` ·
`reports/phase9/leakage_control_matrix.csv` · `reports/phase9/model_comparison_plan.csv` ·
`reports/phase9/validation_specification.json` · `reports/phase9/estimand_separation.csv` ·
`reports/phase9/phase9_readiness_checklist.csv` · `reports/phase9/age_spline_specification.json` ·
`reports/phase9/protocol_change_log.json` · `reports/phase9/phase9_summary.json`

### Revision history

| Version | Date | Status | Summary |
|---|---|---|---|
| **1.0** | 2026-11-03 | superseded | Initial planning draft: recommendations for D8D-1…4, leakage matrix, validation specification, estimand separation. |
| **2.0** | 2026-11-03 | superseded | D8D-1 expanded (scaling guards, clinical-scaling answer, interpretation) · **D8D-4 pooling rule REVERSED** · **D9-1 added** (penalisation scope) · age spline finalised · optimism correction, convergence and non-estimable-replicate policy written out · readiness checklist renamed. **No frozen record modified.** |

| **2.1** | 2026-11-03 | superseded | Issued by the Phase 9 Final Protocol Audit (`docs/PHASE9_FINAL_AUDIT.md`). Batch counts clarified (both sets shown; **no error of fact**) · **D8D-4 50 % rule WITHDRAWN**, replaced by **D8D-4-R1** · D9-1 audited **sound**, with amendment **P6B-A1** drafted. **No frozen record modified.** |

| **3.0** | 2026-11-03 | superseded | **Three approved amendments APPLIED: M2-A2** (training-fold-only pathway-score standardisation; panels and freeze record → **specification revision 3**) · **P6B-A1** (Phase 6B §3.1 penalty scope; original wording preserved) · **D8D-4-R1** (batch estimability; the squared-SE ratio **rejected** and replaced). **D8D-2 and D8D-3 remain PENDING.** |

| **3.1** | 2026-11-03 | **current** | **D8D-2 and D8D-3 APPROVED** and their reporting protocols formalised. Four-set panel = pre-specified **primary**; five-set **sensitivity-only**. Median representation **primary**; representative-IQR **sensitivity-only**. **No specification amendment required** — reporting protocol only; membership, mappings, scoring and the amendment chain unchanged. **No decision remains pending.** |

Machine-readable: `reports/phase9/protocol_change_log.json` · audit:
`reports/phase9/final_audit_summary.json` · Phase 6B amendment record:
`config/analysis/PHASE6B_AMENDMENT_RECORD.json`.

---

## ✅ Applied amendments (v3.0)

| Amendment | Decision | What it settles | Applied to |
|---|---|---|---|
| **M2-A2** | D8D-1 (Option B) | pathway scores standardised **training-fold only** | `m2_panel_4set.json`, `m2_panel_5set.json`, `M2_FREEZE_RECORD.json` → **rev 3**; `src/score_standardisation.py` |
| **P6B-A1** | D9-1 (Option B) | clinical block **unpenalised**, molecular block penalised | `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` §3.1 + new §3.1.1; `config/analysis/PHASE6B_AMENDMENT_RECORD.json` |
| **D8D-4-R1** | D8D-4 | batch-stratification **estimability** rule | this protocol §3.4 |

**Approved 2026-11-03 — reporting protocols, no specification amendment required:**

| Decision | What it formalises |
|---|---|
| **D8D-2** | four-set panel = **pre-specified primary**; five-set = **sensitivity-only** (§3.2) |
| **D8D-3** | median aggregation = **primary** representation; representative-IQR = **sensitivity-only** (§3.3) |

**No Phase 9 decision remains pending.** Execution still requires separate authorisation.

### Frozen and unchanged by v3.0

M2 panel **membership** · panel **roles** (4-set PRIMARY, 5-set SENSITIVITY-only) · gene-level
`sigma_ddof = 1` (M2-A1) · `src/pathway_scoring.py` (byte-identical) · **A-19** complete-case ·
**Phase 7 QC** decisions and records · **Phase 2/3** artifacts · cohort **S3 = 453** · endpoint
definitions · median representation primary, representative-IQR sensitivity-only · **Phase 8A–8D**
findings and artifacts · Phase 6B **Register A**.

> **No model was fitted. No outcome column was read.** `endpoint.csv` is blocked at the I/O
> layer. Every documented event count below is **cited by regex from its governing report**, never
> recomputed — see `phase9_summary.json` → `documented_figures_cited`.
>
> **As of v3.1 all five Phase 9 decisions are APPROVED.** Three carried specification amendments —
> **M2-A2, P6B-A1, D8D-4-R1** — which are applied. **D8D-2 and D8D-3 formalise reporting protocol
> and required no amendment.** The M2
> panel specifications and freeze record are at **specification revision 3** under M2-A2, with
> every superseded hash retained in the amendment chain. **D8D-2 and D8D-3 remain
> `PENDING_USER_APPROVAL`.** A-19, the QC records, Phase 2/3 artifacts, panel membership, cohort,
> endpoint definitions and the Phase 8 findings are **unchanged**.

---

## 0. Governing records — one correction

The brief named `docs/PHASE6C_HALLMARK_PANEL_REVIEW.md`. **That file does not exist.** The actual
Phase 6C report is **`docs/PHASE6C_PATHWAY_PANEL_REVIEW.md`**, and it is the record used here.
All twelve other named records were verified present before any recommendation was formed.

**Phase 6B is further advanced than the brief assumes.** Register A is **fully frozen** (A-01…A-20),
and §5/§5.1 already specify bootstrap optimism correction and the MC-1…MC-4 stability rules. Phase
9 therefore *completes* that protocol rather than replacing it: what was genuinely open is the
resampling unit, stratification, seed policy, M1 functional form, the M3 nested-selection detail,
and the four D8D decisions.

---

## 1. What is frozen and carried in unchanged

| Item | Value | Authority |
|---|---|---|
| Cohort | **S3, 453 samples, 0 exclusions** | A-01, Phase 7C D2 |
| Primary endpoint | cause-specific recurrence, as deposited | A-07 |
| Estimand | cause-specific hazard ratio + discrimination | A-09 |
| M1 | clinical-only Cox baseline | A-11 |
| M2 primary | clinical + frozen **4-set** panel | A-13 |
| M2 sensitivity | clinical + **5-set** panel (+PI3K), sensitivity-only | A-13, M2-A1 |
| M3 | transcriptome-wide, **exploratory only** | A-15 |
| Representation | median aggregation primary; representative-IQR sensitivity | A-04/A-05 |
| Expression matrix | deposited, as-is | B-03 |
| Scoring | `src/pathway_scoring.py`, frozen membership, **`sigma_ddof = 1`** | A-14, M2-A1 |
| Missing data | **complete-case** | A-19 |
| External validation | **none — not feasible** | Phase 4 G4 |

---

## 2. Cohort and parameter budget

| Model | p | complete-case n | lost |
|---|---:|---:|---:|
| M1 (clinical only) | 5 | **452** | 1 |
| M1 + optional TNM/MMR | 9 | 390 | 63 |
| M2 primary (M1 + 4 scores) | 9 | **452** | 1 |
| M2 sensitivity (M1 + 5 scores) | 10 | **452** | 1 |

Source: `reports/phase6d/definitive_audit.json`. The single A-19 drop is **GSM972262** (CIT306),
missing `chemotherapy.adjuvant`.

### 2.1 The event count in the complete-case set is not determined here

S3 carries **136 events** (cited from `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` Register A-01).
A-19 drops one sample. **Whether that sample is an event is outcome data and was not read.**
Complete-case events therefore lie in **[135, 136]**, and M2-primary EPV lies in
**[15.0, 15.1]** — comfortably above the conventional floor of 10 either way, so the budget
decision does not depend on the unknown. The exact count is confirmed at execution time.

### 2.2 Recommendation on the optional covariates

Adding T, N and MMR raises p from 5 to 9 but costs **63 observations** (452 → 390), and MMR is
the variable most entangled with batch and with `combat_mod` (Phase 8B). **M1 primary therefore
uses the five pre-specified covariates only.** The 9-predictor variant is retained as a
pre-specified secondary fit, reported with its reduced n, exactly as A-11 permits.

---

## 3. Decision resolutions — recommendations **PENDING APPROVAL**

### 3.1 D8D-1 — Score-level standardisation · **RECOMMEND OPTION B**

**Currently frozen:** `"NONE before the Cox model (penalisation handles scale). If added later it
is training-fold-only."`

**Recommendation: standardise the pathway scores to unit variance using training-fold parameters
only, inside every resampling boundary, retaining the unstandardised scores as the auditable
source features.**

**Why.** The frozen justification — *"penalisation handles scale"* — **is not correct as stated.**
Ridge and elastic-net penalties are **not scale-invariant**: the penalty λ·Σβ² is applied per unit
coefficient, so a predictor with smaller dispersion needs a *larger* coefficient to express the
same effect and is therefore shrunk *harder*. Phase 8C measured a **3× dispersion range** across
the five scores (SD 0.554 E2F → 0.185 PI3K). Under the status quo the penalty would treat them
unequally **for a reason unrelated to their content**. Standardised coefficients are also
interpretable on a common footing — log-hazard **per 1 SD of score** — which is precisely what a
cross-pathway comparison requires.

**The honest cost.** Standardising **deliberately erases the dispersion difference from the
coefficient scale**, and that difference is real information about set coherence (implied
inter-gene r̄: E2F 0.304, EMT 0.241, WNT 0.041, P53 0.033, PI3K 0.025). The recommendation is
only defensible because that information is **preserved and reported separately** from
`reports/phase8c/score_distributions.csv`, not because it is unimportant.

**Low dispersion is not low biological relevance.** Standardisation must not be read as a claim
that the five scores are equally informative. **No outcome association has been tested.**

Further limitations: it adds a quantity estimated inside each fold (more Monte-Carlo variance,
covered by MC-1…MC-4); and any comparison with an unstandardised fit is **not like-for-like**.

**Leakage constraints if approved:** μ and σ from **training samples only** · applied unchanged to
held-out data · recomputed independently inside **every resampling, bootstrap and inner-CV tuning
operation** — the inner loop, not only the outer split · held-out samples never influence the
parameters · unstandardised scores preserved **unchanged**, written to disk and reported.

#### Scaling guards — defined, not assumed away

| Condition | Rule |
|---|---|
| Training-fold score SD = 0 (or < 1e-12) | the score is **dropped for that fold** and the drop recorded — the same convention the frozen gene-level rule uses. **Never** silently set to 0 or 1 |
| Non-finite training-fold mean or SD | the replicate is **aborted and recorded as non-estimable** (§6.3), never imputed away |

Neither is expected — Phase 8C observed 0 zero-variance genes and 0 non-finite values across 453
samples. The rules exist so the behaviour is **defined**, not because the conditions are
anticipated.

#### Does this settle the clinical predictors? **No.**

Standardising the pathway scores settles **only the pathway block**. The clinical block is four
binaries plus a 2-column age spline basis **on the year scale**; if penalised, those columns
inherit exactly the scale-dependence D8D-1 is fixing — in the *opposite* direction, since a
large-valued basis column carries a small coefficient and is barely penalised. **This is not
resolved inside D8D-1.** It is raised as **D9-1** (§3.5), where the clean answer is to leave the
clinical block *unpenalised* rather than to standardise it.

#### Coefficient interpretation after standardisation

The scale becomes **log cause-specific hazard ratio per 1 training-fold SD of the score**, which
makes coefficients comparable **across** pathways.

**What it does not mean.** A larger standardised coefficient does **not** mean a pathway is more
biologically important. Standardisation removes the dispersion difference from the coefficient
scale, so the coefficient no longer carries information about how coordinated the gene set is.
Score dispersion is a joint property of set coherence and the averaging transform (SD = 1/√k if
members were independent); a low-dispersion score such as P53 (implied inter-gene r̄ = 0.033) is a
**weakly coordinated average** — a statement about co-expression in this tissue under this
preprocessing, **not** about biological relevance and **not** about recurrence, which has not been
tested. Phase 8C's coherence table is reported alongside the coefficients precisely so the two are
not conflated.

#### Amendment **M2-A2** — ✅ **APPLIED 2026-11-03** (specification revision 2 → 3)

| | |
|---|---|
| **Targets** | `m2_panel_4set.json`, `m2_panel_5set.json`, `M2_FREEZE_RECORD.json` |
| **Field** | `scoring.score_standardisation` |
| **From (superseded, retained in the amendment history)** | `NONE before the Cox model (penalisation handles scale). If added later it is training-fold-only.` |
| **To** | `TRAINING-FOLD ONLY. Each pathway score is centred and scaled to unit variance using the mean and sample SD (ddof=1) of the TRAINING samples of that fold, applied unchanged to held-out samples. Recomputed inside every resampling replicate. The unstandardised scores are preserved as the auditable source features. This is a SECOND standardisation, applied to the scores, and does not alter the frozen gene-level standardisation (sigma_ddof=1) that produces them.` |
| Membership / roles changed | **No / No** |
| Outcome data used | **No** |
| Permitted under | `modification_policy` — non-outcome correction, recorded with reason and timestamp, prior version retained, re-freeze |
| **Implementation** | `src/score_standardisation.py` :: `standardise_scores_training_fold` |
| **Canonical scorer** | `src/pathway_scoring.py` **unmodified**; gene-level `sigma_ddof = 1` untouched |
| **Coefficient interpretation** | log cause-specific hazard ratio **per 1 training-fold SD** of the score |
| **Clinical scaling** | **out of scope** here — governed by P6B-A1 |

### 3.2 D8D-2 — Four-set vs five-set · ✅ **APPROVED 2026-11-03** — reporting framework formalised, **no role change**

**Comparisons reported:** M1 (benchmark) · **M2 primary vs M1 — the primary comparison**
(optimism-corrected ΔC with CI and MC-SE) · M2 sensitivity vs M1 (same quantities) · M2
sensitivity vs M2 primary **as a stability check** · M3 vs M1, exploratory-labelled.

**Rules fixed now, before any outcome is seen:**

1. The **4-set result is primary always**, whatever the 5-set shows.
2. The 5-set is reported **alongside, never instead of**, the 4-set.
3. A ΔC between panels **within its Monte-Carlo SE** is declared **"no evidence of a difference"** — not a tie to be broken by a second metric.
4. **A metric improvement in the 5-set is not grounds to switch panels.** An extra parameter and added flexibility can improve apparent fit; that is a property of the comparison, not evidence about PI3K.
5. P53 and PI3K coefficients are reported **together** with CIs and variance inflation, because **r = 0.631** makes individual coefficient instability **foreseeable**.
6. PI3K's descriptive redundancy (**R² = 0.540** on the primary four) is stated whenever its coefficient is discussed.
7. Divergence is evidence about the **stability of the primary estimate**, not a reason to prefer a panel.

**Pre-specified divergence metric:** optimism-corrected ΔC (5-set − 4-set) with bootstrap CI and
MC-SE, fixed before unblinding. **M2-A1's immutability clause stands.**

#### Formalised decision (approved 2026-11-03)

| | |
|---|---|
| **Primary — pre-specified** | the frozen **4-set** panel: **EMT · E2F_TARGETS · WNT_BETA_CATENIN · P53_PATHWAY** |
| **Sensitivity-only** | the frozen **5-set** panel, adding **PI3K_AKT_MTOR_SIGNALING** |
| **Unchanged** | pathway membership · gene mappings · scoring definitions (`sigma_ddof = 1`, M2-A1; score standardisation, M2-A2) · the frozen specification files |

**Prohibitions, binding:** the 5-set panel **must not** be promoted to primary on model performance
or statistical significance · the 4-set analysis is reported as the pre-specified primary result
**regardless of which panel performs better** · the 5-set is reported **as a sensitivity analysis**
· pathway-score dispersion or correlation **must not** be interpreted as proof of biological
activity or clinical relevance.

**Mandatory with any five-set result** (source: `reports/phase8c/score_correlation_pearson_5set_sensitivity.csv`,
`reports/phase8c/phase8c_summary.json`):

| Diagnostic | 4-set | 5-set |
|---|---:|---:|
| **P53 ↔ PI3K Pearson r** | — | **0.6314** |
| Condition number | 2.949 | **7.032** |
| Minimum eigenvalue | 0.4811 | **0.2777** |
| Pairs with \|r\| > 0.5 | 0 | **1** |
| PI3K R² on the primary four | — | 0.5398 |

Variance inflation factors and the coefficient CIs for **P53 and PI3K are reported together**.

> **Why PI3K is sensitivity-only — and what that does *not* mean.** Its correlation with P53
> (r = 0.6314) is the only \|r\| > 0.5 pair in the panel and raises the condition number from 2.949
> to 7.032, making individual coefficient instability **foreseeable**. **This is not evidence that
> PI3K/AKT/mTOR is biologically irrelevant.** It is a statement about the **conditioning of the
> score matrix** in this cohort under this preprocessing, and nothing more. **No outcome
> association has been tested for any pathway.**

**Hierarchy vs performance.** The reporting protocol distinguishes the **pre-specified analysis
hierarchy** — fixed before any outcome was seen, and immovable — from **observed model
performance**, which is a result reported *within* that hierarchy. A ΔC between panels within its
Monte-Carlo SE is reported as **"no evidence of a difference"**, not as a tie to be broken by
another metric.

### 3.3 D8D-3 — Gene representation · ✅ **APPROVED 2026-11-03** — median primary, representative-IQR sensitivity-only

**Do not create a new policy.** Phase 6A **already registers estimand S-7** — a representative-IQR
re-fit (A-08). The recommendation is to **attach a prespecified, WNT-focused comparison to S-7**,
not to invent a parallel mechanism.

**Observed (`reports/phase8c/representation_sensitivity.csv`):** WNT **0.858** against E2F 0.998,
EMT 0.992, P53 0.972, PI3K 0.910. WNT is the only score whose two representations correlate below
0.9.

**Comparison:** a complete re-fit of M2 primary under representative-IQR — **same pipeline, same
folds, same seeds**. Quantities: per-pathway standardised coefficient and CI, pathway by pathway ·
optimism-corrected C and ΔC vs M1 · rank correlation of the linear predictor between
representations · **the WNT coefficient specifically**.

**Interpretation.** Agreement means the result does not depend on the probe→gene rule. It is
**not independent confirmation** — same samples, same genes, same deposited matrix, summarised
differently. Disagreement confined to WNT is **representation sensitivity of a 41-gene,
low-coherence set (r̄ = 0.041)**, not biology. Disagreement across several pathways would call the
representation choice itself into question.

**Selection guard:** the **median** representation is reported as primary **regardless of which is
more favourable**, and the representative-IQR re-fit is reported in full whatever it shows.

#### Formalised decision (approved 2026-11-03)

| | |
|---|---|
| **Primary** | **median probe aggregation** (frozen Phase 2B, A-04) — `expression_median_S0.parquet` |
| **Sensitivity-only** | **representative-IQR**, discovery_443 selection (frozen Phase 2B, A-05; registered estimand **S-7**) — `expression_repiqr_disc443_S0.parquet` |

**Prohibitions, binding:** the primary representation **must not** be replaced on observed model
performance · the correlations below **must not** be read as evidence that one representation is
universally superior · neither may be chosen after seeing outcomes.

**Reported when the analyses are authorised:** whether key conclusions are **directionally** and
**quantitatively** consistent across the two representations · any differences in **pathway
scores, coefficients, uncertainty, discrimination** or other relevant model outputs · the
median-representation result as primary, in full.

**Representation agreement** (Pearson, 453 S3 samples; `reports/phase8c/representation_sensitivity.csv`):

| Pathway | r |
|---|---:|
| E2F_TARGETS | 0.9975 |
| EMT | 0.9920 |
| P53_PATHWAY | 0.9720 |
| PI3K_AKT_MTOR | 0.9096 |
| **WNT_BETA_CATENIN** | **0.8577** |

> **Acknowledged explicitly: WNT has the weakest representation agreement of the five scores
> (0.858) and is the only one below 0.9.** It is the smallest set (41 genes) and the least
> internally coherent of the primary four (implied inter-gene r̄ = 0.041), so each gene's probe
> choice matters more. This is carried into the reporting of any WNT result.
>
> These correlations say how much the two summarisation rules **agree**. They **do not rank the
> representations**, and a higher correlation does not make either one more correct.

> **Robustness ≠ performance.** **Representation robustness** — whether a conclusion survives
> changing the probe-to-gene rule — is a stability property of the pipeline. **Predictive
> performance** is an outcome-dependent quantity, and **none has been computed.** The two are
> reported separately and never substituted for one another. And agreement between representations
> is **not independent replication**: same samples, same genes, same deposited matrix, summarised
> differently.

### 3.4 D8D-4 — Batch sensitivity · **narrow diagnostic + one stratified model**

**Justified, but narrowly.** No batch covariate enters M1 or M2, and batch adjustment never
becomes primary.

**Eligible variables:**

| Variable | Eligibility |
|---|---|
| `combat_batch` (9 levels) | **ELIGIBLE** — the only variable denoting processing group |
| `combat_mod` | **NOT ELIGIBLE** — perfectly concordant with MMR; adjusting for it would adjust away the cohort's dominant biological axis |
| `rna.extraction.method` | **NOT ELIGIBLE** as adjustment — almost completely confounded with batch (χ² = 2265, df = 40), min cell 14; descriptive only |
| `dataset` (discovery/validation) | **NOT ELIGIBLE** — not independent cohorts; expression separation negligible (ε² ≤ 0.005) |

**Design.**
- **Diagnostic 1 — leave-one-batch-out stability.** Refit M2 primary nine times, each omitting one batch; report the spread of standardised coefficients and of ΔC. A stability diagnostic that **estimates no batch effect**. **Pre-registered caveat:** the batches differ ~7-fold in size (14 → 103), so omitting batch03 or batch04 removes ≈ 23 % of the cohort while omitting batch02 removes 3 %. The spread therefore **conflates batch identity with sample size**, and the omitted *n* is reported beside every refit.
- **Sensitivity model 1 — batch-stratified Cox.** M2 primary with the baseline hazard **stratified** by ComBat batch, estimating **no batch coefficients**. Stratification permits different baseline hazards without forcing batch to absorb biological variation through a coefficient.
- Neither runs before the primary analysis; neither can replace it.

#### Sparse strata — the v1.0 pooling rule is **WITHDRAWN**

v1.0 proposed pooling strata below 20 samples. **That was wrong, for a technical reason.**

**A stratified Cox model estimates no parameters for its strata.** The partial likelihood
factorises across strata, so a stratum consumes **no degrees of freedom**. The familiar "too few
events per parameter" concern — which motivates pooling for a *covariate* — simply does not apply
to *stratification*. A stratum with few events contributes little information; one with zero
events contributes a likelihood factor of 1 and drops out. **Neither is an estimability failure.**
The 20-sample threshold was a reflex imported from the covariate case.

**What pooling would have cost.** It forces batch02, batch09 and batch10 — which have nothing in
common except being small — to **share one baseline hazard**. That is precisely the assumption the
sensitivity analysis exists to test, so pooling would **conceal the between-batch heterogeneity
being probed**, and would do so in the three batches where heterogeneity is least observable. The
pooled stratum would also be an artefact with no technical provenance: *"small"* is not a
processing group.

| Option | Verdict |
|---|---|
| **1. Keep all nine strata individually** | ✅ **RECOMMENDED** — no parameters consumed, nothing non-estimable, no heterogeneity concealed. The cost is **efficiency, not validity** |
| 2. Prespecified pooling on counts + provenance | ❌ rejected — no technical provenance links the three small batches, so pooling would rest on sample count alone. Provenance-based pooling would be defensible *if* the SDRF documented shared processing; it does not |
| 3. Omit the stratified model | held in reserve as the **pre-registered failure response**, not as the design |

**Rule:** all **nine** ComBat batches enter as individual strata. **No pooling, by count or
otherwise.**

> **Two different sets — do not conflate them** *(clarified in v2.1)*. The **S3 cohort is 453**
> (Phase 7C D2, no exclusions). The **A-19 complete-case model-fitting set is 452**, because A-19
> drops **GSM972262** (missing `chemotherapy.adjuvant`). Stratum sizes differ by exactly that one
> sample, in **batch05**.

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

Verified against the Phase 3 bundle **and** the ArrayExpress SDRF: **0** missing, **0** duplicated,
**0** blank, **0** ambiguous batch assignments; bundle and SDRF agree for all 453
(`docs/PHASE9_FINAL_AUDIT.md` Finding 1). **All 453 remain in the cohort** — A-19 affects only
which rows enter a model fit.

> **Estimability depends on outcome information.** How much information the nine strata actually
> carry depends on **events per stratum** — outcome data, **not read here**. The *design* is fixed
> outcome-blind; its *information content* can only be measured at execution. Stated rather than
> papered over.

#### Estimability — **D8D-4-R1** ✅ **APPLIED**, replacing the withdrawn 50 % rule

> The v2.0 rule — *"not estimable if fewer than 50 % of events lie in strata that contribute to the
> partial likelihood"* — is **WITHDRAWN**. It was **near-vacuous** (every event lies in a stratum
> containing an event, itself) and **statistically unfounded** (50 % has no basis, and labelled a
> *precision* statement as *estimability*). See `docs/PHASE9_FINAL_AUDIT.md` Finding 2.

**Declared NOT ESTIMABLE if and only if** — criteria whose definitions are outcome-blind:

1. the fit **fails to converge** (ΔlogL < 1e-9, max 100 iterations), **or**
2. the **observed information matrix for β is singular / not positive definite** (rank < p), **or**
3. the **separation sentinel** trips (\|β\| > 10 or SE > 10).

**Reported always, never gated on:** per-stratum **sample** counts (table above) · per-stratum
**event** counts (at execution) · per-stratum count of events with a **non-degenerate
within-stratum risk set** · number of strata contributing ≥ 1 such event · convergence status,
iterations, final gradient norm · the **risk-set comparison retention** (below) · the stratified
and unstratified coefficient estimates with bootstrap intervals, **side by side and never reduced
to a ratio** · an explicit limitation statement beside the estimate.

#### The squared-SE ratio proposed in v2.1 is **REJECTED**

| Reason | Detail |
|---|---|
| **Estimand mismatch** | in a stratified Cox, β is a **within-stratum** log HR; the unstratified β is a different functional, equal only under effect homogeneity **and** proportional baselines. Comparing precisions of two different parameters is not an information statement about either |
| **Penalised SEs are not nominal** | under P6B-A1 the pathway coefficients are ridge-penalised; a penalised SE lacks the usual Wald meaning, and λ is tuned separately per fit — the ratio would partly measure a **difference in tuned penalty** |
| **No guaranteed direction** | stratification removes between-stratum comparisons (less information) but also removes stratum confounding (less residual variation). The ratio is **not bounded by 1**, and a value > 1 is uninterpretable as "retention" |
| **Unquantified** | it is a ratio of estimates with its own sampling variability, and no uncertainty quantification was specified |

**It is not forced into the protocol.**

#### Replacement — **risk-set comparison retention** (descriptive, no threshold)

> **N**<sub>strat</sub> **/ N**<sub>unstrat</sub>, where
> N<sub>unstrat</sub> = Σ<sub>events *i*</sub> (\|R(t<sub>i</sub>)\| − 1) and
> N<sub>strat</sub> = Σ<sub>events *i*</sub> (\|R<sub>s(i)</sub>(t<sub>i</sub>)\| − 1).

**Why defensible:** **model-free** — a combinatorial property of the design and event times,
needing no fit, no SE, no λ · **bounded by construction**, since R<sub>s</sub> ⊆ R, so the ratio
lies in [0, 1] · **directly interpretable** — "stratification retained X % of the risk-set
comparisons the partial likelihood could have used" · it **subsumes** the one-event-stratum
question automatically, because an event alone in its stratum contributes 0 comparisons, matching
its log(1) = 0 likelihood contribution.

**Limitations, stated not hidden:** it counts **comparisons, not Fisher information** — a retained
comparison between two subjects with identical covariates is counted yet contributes nothing · it
is reported **alongside**, not instead of, the count of events with a non-degenerate risk set · it
says nothing about any coefficient's precision and is **not** a substitute for the confidence
intervals · it needs event times, so it is computed at authorised execution, though its
**definition is fixed here, outcome-blind**.

**No pass/fail cut-point is placed on it, now or later.**

#### The diagnostics are **necessary, not sufficient**

- A **positive-definite information matrix does not prove adequate information or precision** — only that the parameters are locally identified and the Hessian invertible. A model can pass all three checks and still have intervals too wide to be useful.
- A **numerical failure is not a wide confidence interval.** A wide interval is a **result**, reported as such; a numerical failure means there is **no** result. The two are never reported as if they were the same.
- A **zero-event stratum** contributes no event information and **must not** be treated as a numerical failure.
- A **one-event stratum is not automatically informative**: whether its event has a non-degenerate risk set is **reported**, not assumed.

**Why events-per-stratum is the wrong quantity.** A zero-event stratum contributes a likelihood
factor of **1** and drops out. An event whose within-stratum risk set is itself alone contributes
**log(1) = 0** — zero score, zero information — *while still being counted* in a naive tally.

If declared not estimable, that is **reported**, and the leave-one-batch-out diagnostic stands
alone. **We do not fall back to pooling** — a pooling choice made after seeing events-per-stratum
would be an outcome-informed design change.

**Interpretation.** Agreement means the result is not driven by processing group. **It does not
mean batch effects are absent** — ComBat was applied before deposit, so between-batch variance was
removed by construction. Disagreement is reported as an unresolved limitation, **not** grounds to
adopt the adjusted model. Note that the two surviving Phase 8B signals are a **control-probe**
effect and a **dispersion** effect; neither is a location effect, so **neither would be addressed
by a batch coefficient** — a further reason the adjustment is not promoted.

**Why never primary:** batch, extraction method and MMR are mutually confounded · a batch covariate
would absorb biology (batch05 is 63.3 % dMMR vs 14.3 %) · the matrix is already ComBat-corrected,
so adjusting again risks removing signal twice · the primary analysis was pre-specified without
batch, and changing that after seeing results forfeits its pre-specified status.

### 3.5 D9-1 — Penalisation scope · ✅ **APPROVED (Option B) — P6B-A1 APPLIED**

**Raised by D8D-1**, which exposes the question for the clinical block without answering it.

Phase 6B §3.1 says *"ridge, λ tuned inside resampling"* **without saying which predictors the
penalty applies to**. Penalisation scope is **not** a Register A item (A-11 fixes the predictor
*list*, A-15 the model *roles*, A-17 the resampling), so Phase 9 may finalise it — but it changes
the 6B draft text, so it is **raised for approval, not settled silently**.

**The problem.** Ridge is not scale-invariant. If the four binaries and the year-scale age spline
basis are penalised alongside the scores, the basis columns carry small coefficients and are barely
shrunk while the binaries are shrunk comparatively hard — an ordering driven by **units**.

| Option | Verdict |
|---|---|
| A. Penalise everything, standardising the clinical block too | unusual for binaries, and makes a 0/1 contrast depend on the training fold's prevalence |
| **B. Force-include the clinical block UNPENALISED; penalise only the pathway scores (M2) and selected genes (M3)** | ✅ **RECOMMENDED** |
| C. Leave the scope unspecified | not viable — the implementation must do something |

**Why B.** An unpenalised coefficient is **scale-equivariant**: the fit is invariant to any linear
rescaling of that predictor, so option B makes the clinical scaling question **disappear** rather
than answering it with a convention. M1 is the benchmark, and with 6 df against 135–136 events
(EPV ≈ 22) a classical unpenalised Cox fit is appropriate. It makes **ΔC(M2 − M1) attributable to
the pathway block** rather than to differential shrinkage of clinical terms. And it is already what
the M3 specification assumes, so B makes the three models mutually consistent.

#### Audited in v2.1 — verdict **SOUND**

M1 is **exactly** M2 with the pathway coefficients set to zero, so the models are properly nested;
had the clinical block been penalised with a model-specific λ, M1 would **not** be the limit of M2
and the nesting would be **broken** — this choice is what *makes* the comparison coherent. Clinical
coefficients are **comparable** across models because they are unpenalised in both, so differences
are genuine adjustment rather than shrinkage artefact.

> **Statement required whenever ΔC is reported.** ΔC(M2 − M1) estimates the gain from adding a
> **shrunk** pathway block to an **unshrunk** clinical model. The penalty pulls the pathway
> coefficients toward zero — toward M1 — so the comparison is **conservative for the molecular
> block**: biased *against* finding a molecular contribution. That is the right direction of error
> here, but it must be stated, not relied on silently.

#### ✅ Amendment **P6B-A1** — **APPLIED 2026-11-03**

| | |
|---|---|
| **Target** | `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` §3.1, *Penalisation* row of the M1 table |
| **From** | `ridge, λ tuned inside resampling` |
| **To** | the clinical block (four binaries + age spline basis) is **force-included and unpenalised** in M1, M2 and M3; ridge applies only to the M2 pathway scores and elastic net only to the M3 candidate genes; λ (and α for M3) tuned by inner CV inside each training fold. **M1 therefore has no penalised term and no λ to tune** |
| **Register A item?** | **No** — A-11 fixes the predictor list, A-15 the roles, A-17 the resampling |

**Limitations.** M1 becomes an unpenalised Cox model, so no λ is tuned for it — **a change from the
6B draft text**, which must be stated in any comparison with that draft. Unpenalised clinical terms
are not shrunk, so M1 carries the usual apparent optimism — which is exactly what the bootstrap
correction estimates, and the uniform shrinkage factor is reported regardless.

---

## 4. The models

### 4.1 M1 — clinical baseline (the benchmark)

| | |
|---|---|
| **Predictors** | age, sex, stage (II/III), tumour location (proximal/distal), adjuvant chemotherapy |
| **Functional form — age** | **restricted cubic spline, 3 knots (2 df)** — fully specified in §4.1.1 |
| **Functional form — all others** | binary, reference levels fixed now: sex = Female, stage = II, location = distal, chemo = N |
| **Parameters** | 6 df (age spline 2 + four binaries) against 135–136 events → **EPV ≈ 22** |
| **Missing data** | complete-case, n = 452; the single drop (GSM972262) reported, never imputed |
| **Penalisation** | ridge, λ tuned by inner CV **inside** each training fold |
| **Chemotherapy** | **prognostic covariate only — its coefficient is never interpreted as a treatment effect** (A-12; confounding by indication, confounders absent from the metadata) |
| **Stage** | included as a covariate; **stratification by stage** assessed as the pre-specified alternative (A-11) |
| **PH assessment** | scaled Schoenfeld residuals per covariate + global test; smoothed residual plots. The **diagnostic** is pre-specified (B-05); any **remedy** is reported as a protocol deviation |
| **Collinearity** | VIF on the design matrix and the smallest eigenvalue of the scaled information matrix, reported before any coefficient is interpreted |
| **Convergence** | log-likelihood change < 1e-9, no coefficient with \|β\| > 10 or SE > 10 (separation sentinel); failures reported, never silently refitted |

**M1 is built first and never re-tuned afterwards.** Every genomic model is measured against this
fixed benchmark.

#### 4.1.1 Age functional form — finalised *(new in v2.0)*

Machine-readable: `reports/phase9/age_spline_specification.json`.

| | |
|---|---|
| **Form** | restricted cubic spline, **3 knots (2 df)** |
| **Why 3** | with 6 df total against 135–136 events (EPV ≈ 22), 3 knots buys curvature without spending the budget; 4–5 knots would be justified only at a far larger event count |
| **Knot placement** | 10th, 50th and 90th percentiles of age **among the training samples of that fold** |
| **Quantile definition** | `numpy.percentile` default (linear interpolation, type 7) — named so two implementations agree |
| **Recomputed per fold** | yes; knots never come from the full cohort |
| **Outcome blindness** | knots depend on the age distribution only. The outcome is **never** inspected to place a knot or choose the form |
| **Held-out samples** | mapped through the **training fold's** knots, unchanged |

**Reference age for reporting: fixed at 69 years** (the cohort median — an outcome-blind covariate
summary), declared now and never changed. The Cox model has no intercept, so the reference does
**not affect the fit**; it affects only the reported HR(age = *a* vs *a*<sub>ref</sub>). A per-fold
reference would make reported HRs **non-comparable across folds**, which is why it is a single
constant rather than a fold quantity.

**Degenerate cases — defined, and recorded rather than silent:**

| Case | Rule |
|---|---|
| Two knots coincide (within 1e-6) | use the distinct knot values; if **fewer than 3** remain, that fold **falls back to linear age**, recorded |
| Rank-deficient basis (constant column, or collinear to 1e-10) | that fold **falls back to linear age**, recorded |

Neither is expected: the cohort has **89 distinct ages** with **14-** and **20-year gaps** between
the full-cohort 10/50/90 percentiles (p10 = 49, median = 69, p90 = 83; range 22–97). Those figures
are quoted to show duplicate knots are implausible — **not** to place the knots, which are per
training fold.

**Pre-specified alternative:** a linear age term, compared by likelihood ratio **within** the
training fold. The comparison is reported; it never selects the form using held-out data or
outcomes.

**Covariate set unchanged** — the five pre-specified clinical covariates only. T, N, MMR and
molecular covariates are **not** added (A-11).

### 4.2 M2 — primary pathway model

M1's five clinical covariates **plus the four frozen scores** — EMT, E2F_TARGETS,
WNT_BETA_CATENIN, P53_PATHWAY — scored by `src/pathway_scoring.py` with `sigma_ddof = 1` read from
the specification. Panel membership, ordering and roles are frozen and are **not** revisited.

Scaling: **gene-level** standardisation is training-fold-only and frozen (M2-A1). **Score-level**
standardisation is **D8D-1, pending approval**; if approved it is likewise training-fold-only.

Penalisation: ridge, λ by inner CV inside the training fold (B-07). Fold-safe preprocessing per
the matrix in §5. Coefficients are reported as cause-specific log-hazard ratios with CIs, and —
if D8D-1 is approved — per 1 SD of score. **No coefficient is interpreted as causal**, and no
coefficient is evidence of prognostic value absent the full validated comparison.

### 4.3 M2 sensitivity

The 5-set panel (+PI3K), compared to M2 primary under the §3.2 framework. Carries the WNT
representation comparison (§3.3, via S-7) and the batch sensitivity analyses (§3.4) **if approved**.
**Role unchanged: sensitivity-only.**

### 4.4 M3 — exploratory transcriptome-wide model

**Exploratory status is not a disclaimer; it is the model's identity.** M3 generates hypotheses
and makes **no confirmatory claim**.

Everything is **strictly nested inside the training fold**: the variance-filter cut-point (B-02,
rule frozen, value per fold) · **all** univariate screening · feature selection · scaling ·
penalty tuning (inner CV). **No full-cohort univariate screening occurs at any point** — that is
the classic selection-leakage failure, and it is prohibited here, not merely discouraged.
Held-out samples contribute **nothing** to any of it.

Selector: elastic net over the training fold's surviving genes, α and λ by inner CV. Clinical
covariates are force-included and unpenalised so the comparison with M1 is interpretable.

**Multiplicity and overfitting are reported, not minimised:** 21,815 candidate genes · the
selected set will differ across folds and seeds, and **selection stability is itself reported**
(selection frequency per gene across replicates) · optimism will be largest here, and M3 is
expected to need the largest B (MC-1…MC-4). A gene selected in one fold is **not** a finding.

---

## 5. Leakage-control matrix

Machine-readable: `reports/phase9/leakage_control_matrix.csv`.

| Operation | Permitted location | Authority |
|---|---|---|
| Frozen panel membership | OUTSIDE | A-13 |
| Gene→probe representation map | OUTSIDE | A-04/A-06 |
| Cohort membership (S3) | OUTSIDE | A-01 |
| QC exclusions | OUTSIDE (none performed) | 7C D2 |
| Gene-level μ/σ for scoring | **INSIDE** | A-14 / M2-A1 |
| Zero-variance gene exclusion | **INSIDE** (per fold, recorded) | A-14 |
| **Pathway-score standardisation** | **INSIDE** *(if D8D-1 approved)* | D8D-1 **PENDING** |
| Variance-filter cut-point (M3) | **INSIDE** (rule frozen, value per fold) | B-02 |
| Univariate screening (M3) | **INSIDE** — no full-cohort screening ever | M3 spec |
| Feature selection (M3) | **INSIDE**, strictly nested | M3 spec |
| Penalty tuning (λ, α) | **INSIDE** (inner CV) | B-07 |
| Imputation | N/A — complete-case | A-19 |
| PCA rotation, if ever a predictor | **INSIDE** (8A/8B rotations are descriptive only) | 6B §3.2.2 |
| Batch stratum definition | OUTSIDE (from counts, pre-unblinding) | D8D-4 **PENDING** |
| Choice of B | OUTSIDE (MC-driven; never lowered to suit a result) | B-06 |
| Random seeds | OUTSIDE (declared in advance) | A-17 |
| **Outcome column join** | at model fitting **only** — never during preprocessing or scoring | A-16 |

---

## 6. Validation and resampling

Machine-readable: `reports/phase9/validation_specification.json`.

| | |
|---|---|
| **Primary** | bootstrap optimism correction (Harrell); the **entire** pipeline repeated inside every replicate |
| **Cross-check** | repeated 10-fold cross-validation |
| **External** | **none** — no suitable cohort exists. **No analysis here is external validation.** |
| **Resampling unit** | the **patient**. Phase 8B established **0 declared dependent arrays** in S3, so array and patient coincide |
| **Stratification** | bootstrap and CV folds stratified by **event status** |
| **Not stratified by** | batch, partition, MMR — stratifying on a confounded covariate makes replicates unrepresentative of the estimand's target population |
| **Tuning boundary** | inner CV **inside** each training fold; never on held-out data, never on the full cohort |
| **B** | initial 500; **B is an output of the stability analysis**, increased (1000, 2000, …) until MC-3 passes, **never lowered** |
| **Seeds** | **20261102, 20261103, 20261104, 20261105, 20261106** — declared here, never changed after seeing a result |

> ⚠ **Stratifying by event status requires the outcome** and is therefore a
> **model-fitting-time** operation. It is *specified* here and *executed* only in Phase 9
> execution. No fold used in Phases 8A–8C was outcome-stratified; those folds were explicitly
> labelled illustrative.

### 6.1 Monte-Carlo error vs sampling uncertainty — never conflated

| | |
|---|---|
| **Sampling SE** | uncertainty from having 453 patients rather than infinitely many. Estimated by the bootstrap distribution. **Does not shrink** with more replicates |
| **Monte-Carlo SE** | uncertainty from having run *B* replicates rather than infinitely many. Estimated as the **SD of the corrected estimate across the five declared seeds** (MC-2). **Shrinks roughly as 1/√B** |
| **Reporting** | both reported for every corrected estimate, **labelled separately**. Never summed into one interval |
| **Criterion** | MC-SE ≤ **10 %** of sampling SE (MC-3); otherwise increase B and repeat |

### 6.2 Optimism correction, convergence, and non-estimable replicates *(new in v2.0)*

**Optimism correction (Harrell's bootstrap):**

1. Fit the **entire** pipeline on the original sample → apparent performance *P*<sub>app</sub>.
2. For *b* = 1…*B*: draw a bootstrap sample (stratified by event status), **refit the entire pipeline** on it — scoring parameters, any score standardisation, filtering, selection and tuning — then record *P*<sub>boot</sub>(*b*) on the bootstrap sample and *P*<sub>orig</sub>(*b*) on the **original** sample using the bootstrap-fitted model.
3. optimism(*b*) = *P*<sub>boot</sub>(*b*) − *P*<sub>orig</sub>(*b*).
4. corrected = *P*<sub>app</sub> − mean<sub>*b*</sub> optimism(*b*).

Apparent performance, mean optimism and corrected performance are reported **together**; apparent
performance is **never** reported alone (6B §5). A uniform shrinkage factor is estimated from the
bootstrap and reported whether or not it is applied.

**Convergence criteria:** log-likelihood change < 1e-9 between iterations · max 100 iterations ·
separation sentinel at |β| > 10 or SE > 10.

**Non-estimable replicates.** Every replicate that fails to converge, hits the iteration cap, trips
the separation sentinel, or cannot be fitted (rank-deficient design, undefined scaling) is
**recorded** with its seed, replicate index and failure reason in a per-model replicate log. Such a
replicate is **excluded from the optimism average and counted**, with the count and percentage
reported beside every corrected estimate.

> **Pre-registered threshold (fixed now):** if **more than 5 %** of a model's replicates are
> non-estimable, its corrected estimate is reported as **UNRELIABLE** and the failure pattern
> described. **Replicates are never silently re-drawn** to replace failures — that would bias the
> optimism estimate toward well-behaved samples.

**Tuning nesting:** 5-fold inner cross-validation **inside** each bootstrap replicate / outer
training fold, tuning ridge λ (M2 pathway block) and elastic-net α, λ (M3). No hyperparameter is
ever selected using held-out data or the full cohort. If D8D-1 is approved, score standardisation
is recomputed inside the **inner** loop too.

### 6.3 What is estimated — and what is not

**Estimated:** optimism-corrected Harrell C-index · calibration slope on the **log-hazard scale**
(unit-free) · calibration-in-the-large on the log-hazard scale · cause-specific hazard ratios with
CIs · ΔC vs M1 with CI and MC-SE · uniform shrinkage factor (reported whether or not applied) ·
effective degrees of freedom.

**Not estimated — and why:** t-year / fixed-horizon absolute risk · predicted survival or CIF
curves · time-dependent AUC or Brier · cloglog calibration at a horizon (Riley's estimator) — **all
blocked by U1** (A-10, Phase 6A §3.1). Any external-validation quantity — **no cohort exists**.

> ### ⛔ Explicitly prohibited while U1 is unresolved
>
> The following are **prohibited**, not merely omitted, and may not be produced, reported or
> referenced until **U1 is resolved and a valid method is separately approved**:
>
> - fixed-time absolute-risk prediction of any kind
> - **five-year** (or any named-horizon) recurrence risk
> - named-horizon calibration, including cloglog calibration at a horizon
> - time-dependent AUC, time-dependent Brier score, or any horizon-indexed discrimination measure
> - predicted survival or cumulative-incidence **curves**
> - the 60-unit endpoint (estimand S-3), except as an explicitly labelled **conditional** analysis
>
> A "60-unit" horizon is 5 years, 1.15 years or 0.16 years depending on whether the unit is months,
> weeks or days. **No source states the unit** (eight searched, Phase 6A §1).

**C-index limitation, stated plainly.** Harrell's C is **censoring-distribution dependent**: it is
computed over comparable pairs, and heavier censoring both reduces the usable pairs and biases the
estimate relative to a censoring-free population. It is reported as a **relative** discrimination
measure for comparing models **on this cohort**, not as a transportable performance figure.

> **The Phase 5 simulation is an APPROXIMATION of the Riley et al. (2022) framework, not a
> faithful implementation of it.** The ~338-event figure is an approximation and was explicitly
> **not** a cross-check of Riley's published ~365 (Phase 5 §0.1).

---

## 7. Estimand separation — three different quantities

Machine-readable: `reports/phase9/estimand_separation.csv`.

| ID | Estimand | Events | Quantity | Distinctness |
|---|---|---|---|---|
| **PRIMARY** | cause-specific recurrence hazard | 136 | cause-specific HR + discrimination | censors competing death; "rate among those still at risk" |
| **S-1** | Fine–Gray subdistribution / CIF | 46 competing | absolute recurrence probability in the presence of death | **different estimand**; coefficients **not comparable** with the primary |
| **S-2** | composite DFS | 182 | recurrence-or-death hazard | **different event definition**; a larger event count is **not** more power for the primary question |
| **S-3** | 60-unit administrative horizon | 130 | primary, truncated | **conditional on U1** — "60 units" has no defined length while the unit is unresolved |
| **S-7** | representative-IQR representation | 136 | primary | aggregation sensitivity; the vehicle for D8D-3 |

**These are not three measurements of one thing.** Their coefficients and performance measures are
reported separately and **never pooled, averaged or substituted** for one another.

---

## 8. Unresolved limitations carried forward

| Limitation | Blocks | Does **not** block |
|---|---|---|
| **U1** — `rfs.delay` unit `[PROVISIONAL: months]` | every **fixed-time absolute-risk** analysis: t-year risk, survival/CIF curves, cloglog calibration at a horizon, time-dependent AUC/Brier, the 60-unit endpoint (S-3 conditional) | **all relative-risk analyses** — HRs, C-index, log-hazard calibration slope, ΔC. The Cox partial likelihood uses only the **ordering** of times; rescaling leaves β, SE and C identical to 10 decimal places (Phase 6A §3.1) |
| **U2** — death without recurrence | treating competing death as non-informative censoring **without comment**; the censoring is **demonstrably informative** (46 competing deaths, associated with age, stage, treatment) | the primary cause-specific analysis, **provided** S-1 and S-2 are reported as distinct estimands |
| **`combat_mod` provenance** | any claim about what ComBat protected; **any imputation of missing MMR** | the primary analysis. 43 S3 samples have MMR unrecorded and are **not** imputed |
| **Eight sex-axis-discordant samples** | nothing | all remain in the cohort; **genotype-level confirmation unavailable** |
| **No confirmatory external cohort** | **any external-validation claim** | internal bootstrap optimism correction |
| **GSE14333 / GSE33113** | pooling — **counts must not be combined**; GSE33113 is a separate estimand (stage II, composite, days) | nothing in the internal analysis; both remain exploratory candidates only |

---

## 9. Readiness

Machine-readable: `reports/phase9/phase9_readiness_checklist.csv`.

**16 READY** · **5 DECISION REQUIRED** (D8D-1…D8D-4, D9-1) · **1 CONFIRM AT EXECUTION**
(complete-case event count, in [135, 136]) · **2 CONSTRAINTS** (U1, U2) · **1 NOT AVAILABLE**
(external validation) · **1 AWAITING USER** (this protocol).

**Conflict scan against frozen decisions: 0 real conflicts.** D8D-1 proposes a change to a frozen
field; it is recorded as amendment **M2-A2** and is **not applied**. The frozen records are
byte-unchanged on disk.

---

## 10. Decisions requiring approval

| # | Decision | Recommendation | Status |
|---|---|---|---|
| **D8D-1** | Score-level standardisation | **Option B** — training-fold-only standardisation; amendment **M2-A2** drafted | **PENDING** |
| **D8D-2** | 4-set vs 5-set interpretation | adopt the reporting framework; **no role change** | **PENDING** |
| **D8D-3** | Gene representation reporting | median **primary**, representative-IQR **sensitivity-only** (estimand S-7) | ✅ **APPROVED 2026-11-03** |
| **D8D-4** | Batch sensitivity | LOBO diagnostic + batch-**stratified** model; `combat_batch` only; **all nine strata, no pooling** | **PENDING** |
| **D9-1** | Penalisation scope / clinical scaling | **Option B** — clinical block force-included and **unpenalised**; amendment **P6B-A1** drafted | **PENDING** |
| **D8D-4-R1** | Estimability rule replacing the withdrawn 50 % criterion | convergence + information-matrix rank + separation sentinel; information loss **reported**, not gated | **PENDING** |

Additionally: **approval of this protocol as a whole** is required before any model is fitted.

### 10.1 What changed in v3.1

**D8D-2 and D8D-3 approved** (project owner, 2026-11-03). Both formalise **reporting protocol**;
**no specification amendment was required**, and none was made.

1. **D8D-2** — four-set panel is the **pre-specified primary** M2 analysis; five-set is **sensitivity-only**. The pre-specified **hierarchy** is now explicitly distinguished from **observed performance**, the P53–PI3K correlation (0.6314) and collinearity diagnostics are **mandatory** with any five-set result, and PI3K's sensitivity-only status is recorded as a **conditioning** matter — explicitly **not** evidence of biological irrelevance.
2. **D8D-3** — median aggregation is the **primary** representation, representative-IQR **sensitivity-only**. **WNT's weakest agreement (0.8577)** is acknowledged explicitly against E2F 0.9975, EMT 0.9920, P53 0.9720, PI3K 0.9096, with the correlations explicitly **not** ranking the representations, and **robustness kept distinct from performance**.

**No Phase 9 decision remains pending.** Panel membership, gene mappings, scoring definitions, the
frozen representations and the M2-A1/M2-A2 amendment chain are **unchanged**.

### 10.2 What changed in v3.0

Three **approved** amendments applied under user authorisation of 2026-11-03.

1. **M2-A2 APPLIED** — pathway scores standardised **training-fold only**. Both panels and the freeze record move to **specification revision 3**; implementation at `src/score_standardisation.py`. **Panel membership, roles and the gene-level `sigma_ddof = 1` are unchanged**, and `src/pathway_scoring.py` is **byte-identical**.
2. **P6B-A1 APPLIED** — Phase 6B §3.1 penalty scope amended; the **original wording is preserved verbatim** in the document and in `config/analysis/PHASE6B_AMENDMENT_RECORD.json`. **Register A untouched.** The binding ΔC interpretation is recorded.
3. **D8D-4-R1 APPLIED** — estimability on convergence, information-matrix rank and separation only. The v2.1 **squared-SE ratio is rejected** and replaced by **risk-set comparison retention**, descriptive and unthresholded.

**D8D-2 and D8D-3 remain PENDING** — not part of this authorisation.

### 10.3 What changed in v2.1

Issued by the **Phase 9 Final Protocol Audit** (`docs/PHASE9_FINAL_AUDIT.md`).

1. **Batch counts clarified — no error of fact.** The nine counts were correct for the complete-case set; both count columns are now shown side by side. Join integrity verified against the bundle *and* the SDRF: 0 missing, 0 duplicated, 0 blank, 0 ambiguous.
2. **D8D-4's 50 % rule WITHDRAWN** and replaced by **D8D-4-R1**. It was near-vacuous and statistically unfounded.
3. **D9-1 audited sound**, with amendment **P6B-A1** drafted for the Phase 6B §3.1 penalisation row, and the conservative-ΔC statement made mandatory.
4. **33/33 internal consistency checks pass.**

**No frozen record was modified by v2.1** (M2-A2 was applied later, in v3.0).

### 10.4 What changed in v2.0

**One recommendation was reversed.** v1.0's rule pooling batch strata below 20 samples is
**withdrawn** (§3.4): stratification consumes no degrees of freedom, so the threshold was a reflex
imported from the covariate case, and pooling would have concealed the very heterogeneity the
analysis probes. **One decision was added** — D9-1, because D8D-1 exposed that the clinical block
would inherit the same scale-dependence. Everything else is completion, not change.

**As of v2.1, no frozen record had been modified** — `scoring.score_standardisation` still read
`NONE` and the freeze record carried only `M2-A1`. **That changed in v3.0**, when the approved
amendment M2-A2 was applied (§"Applied amendments").

**Phase 9 protocol v3.1: all five decisions approved, three amendments applied, none pending.
Phase 9 execution has not been started — zero models fitted, zero outcome analyses performed.**
