# Portfolio presentation blueprint: GSE39582 colon-cancer recurrence study

**Date:** 8 November 2026. **Revision 3 — corrected and consolidated after independent audit.**

**Type:** presentation plan only.
- No analysis was run and no model was fitted.
- No frozen result, figure, caption or existing report changed.
- Every number below is quoted from existing project records, with its source.
- **Implementation status (8 November 2026):** the five diagrams and the presentation are built in `portfolio/` from frozen outputs only (see "Implementation" below).

**Central research question (verbatim, used on every relevant slide).**
Do four biologically selected tumour-expression pathway scores add reproducible recurrence-discrimination information beyond established clinical variables in GSE39582?

**Central conclusion (verbatim).**
- The addition of four frozen Hallmark pathway scores to the clinical Cox model did not demonstrate a reproducible improvement in recurrence discrimination in this GSE39582 analysis.
- This does not establish that the biological pathways have no role in colorectal-cancer recurrence.

**One-line summary.**
- A pre-specified, internally validated evaluation of whether four colorectal-cancer pathway-expression scores add recurrence discrimination beyond clinical variables in GSE39582.
- The analysis did not demonstrate a reproducible improvement, while documenting batch-sensitive coefficients and subtype-associated score variation that complicate biological interpretation.
- It does not establish that these pathways have no biological role.

### Evidence labels (use on every biology slide)

| Label | Meaning |
|---|---|
| **Observed** | Observed in this analysis (frozen project output) |
| **Literature** | Supported by external literature: only R1–R8, verified in `docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md` |
| **Hypothesis** | Interpretation or hypothesis; not demonstrated here |
| **Not tested** | Not assessed in this project |
| **Textbook** | Standard textbook background, not covered by R1–R8; presented as background only |

**References.**
- **R1** Guinney 2015, *Nat Med*, doi:10.1038/nm.3967
- **R2** Calon 2015, *Nat Genet*, doi:10.1038/ng.3225
- **R3** Isella 2015, *Nat Genet*, doi:10.1038/ng.3224
- **R4** Anjomshoaa 2009, *Br J Cancer*, doi:10.1038/sj.bjc.6605229
- **R5** Luo 2019, *BMC Cancer*, doi:10.1186/s12885-019-5324-y
- **R6** TCGA 2012, *Nature*, doi:10.1038/nature11252
- **R7** Liberzon 2015, *Cell Syst*, doi:10.1016/j.cels.2015.12.004
- **R8** Marisa 2013, *PLoS Med*, doi:10.1371/journal.pmed.1001453

**Figures.**
- Fig 1–9 and S1–S2 are the frozen figures in `reports/figures/final/`, with captions in `FIGURE_CAPTIONS.md`.
- Items marked **DIAGRAM D1–D5** are built in `portfolio/diagrams/` (SVG with editable text, and PNG). They use only existing numbers, are not added to `reports/figures/final/`, and carry the caveats given here.

---

## A. Biological background

**Scientific message.** The panel represents four biological axes of colorectal cancer. A bulk-expression score for each axis summarises mRNA levels of a gene set. It does not measure the axis's activity.

| Programme | Mechanism (**Textbook**) | Colorectal-cancer context (**Literature**) | What the project's score is |
|---|---|---|---|
| **Wnt/β-catenin** | Without a Wnt signal, a destruction complex (APC, AXIN, GSK3β, CK1) phosphorylates β-catenin for degradation. APC loss or activating CTNNB1 mutation stabilises β-catenin, which with TCF/LEF factors drives targets such as MYC and AXIN2 | Wnt-pathway genes are genomically altered in 93 % of tumours, with APC or CTNNB1 alterations in about 80 % (R6). Marked WNT/MYC activation characterises CMS2 (R1). **Genomic alteration is not measured pathway activation** | Mean standardised mRNA of 41 member genes |
| **E2F targets** | Cyclin D–CDK4/6 and cyclin E–CDK2 phosphorylate RB, releasing E2F factors that activate G1/S and DNA-replication genes (for example MCM family, PCNA, cyclin E) | The relation of proliferation to outcome is **mixed**: lower proliferation-signature expression in relapsed primaries (R4), but high Ki-67 associated with poorer survival in a meta-analysis (R5) | Mean standardised mRNA of 195 member genes |
| **p53 pathway** | Stress releases p53 from MDM2-mediated degradation; p53 transactivates arrest (CDKN1A/p21), feedback (MDM2) and apoptosis (for example BAX) genes | TP53 mutated in about 60 % of non-hypermutated vs about 20 % of hypermutated tumours (R6) | Mean standardised mRNA of 194 member genes, described as p53 transcriptional targets. **Not a TP53 mutation assay or a functional p53 measurement** |
| **EMT** | Epithelial cells lose polarity and E-cadherin junctions and gain mesenchymal features (vimentin, N-cadherin, fibronectin, collagens), under SNAIL, ZEB and TWIST, often TGF-β-induced | The mesenchymal CMS4 subtype had worse relapse-free survival (R1; GSE39582 contributed to the CMS derivation, so this is not independent). Mesenchymal-subtype transcripts in bulk tumours are mostly stromal (R2, R3) | Mean standardised mRNA of 198 member genes. **Not a direct measurement of EMT in malignant cells**; may also reflect fibroblast, endothelial or immune expression |
| PI3K/AKT/mTOR (five-set sensitivity only) | Growth-factor signalling through PI3K, AKT and mTOR | PIK3CA mutated in about 15 % of non-hypermutated tumours (R6) | 104 member genes |

**Textbook framing for the opening slide.**
- Colorectal cancers commonly arise through Wnt-pathway disruption, with proliferation and TP53 loss accompanying progression.
- Mesenchymal and stromal programmes characterise a poor-prognosis subtype (CMS4; R1).

**Key concept.**
- A bulk tumour sample mixes cancer epithelium with stromal, fibroblast, endothelial and immune cells.
- A bulk score may therefore combine tumour-cell programmes, contributions from these other cells, subtype-associated variation and technical batch structure.
- These are **possible contributors, not quantified components**. Cell composition was **not measured or deconvolved** in this project.

| Element | Specification |
|---|---|
| Best visual | **DIAGRAM D1** (pathway lanes; conceptual schematic) |
| Evidence | Gene counts: `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` §3.2.1 (MSigDB Hallmark v2026.2.Hs). Literature: audit §1–§2. The R7 citation for "p53 transcriptional targets" was removed by the audit, so that description is given uncited |
| Recruiter should understand | You know the biology of each axis and the difference between mutation, pathway activity and a bulk mRNA score |
| Avoid | Describing scores as activation, protein activity, flux or mutation status; describing genomic Wnt alteration as activation; ranking CMS4 beyond "poor prognosis"; citations beyond R1–R8 |

---

## B. Biological research design

**Scientific message.** The biological question was turned into a constrained, leakage-controlled transcriptomic analysis, fixed before any outcome was examined.

### 1. Programme selection, outcome-blind (`docs/PHASE6C_PATHWAY_PANEL_REVIEW.md`)
- **Candidates:** ten Hallmark sets (§4).
- **Correlation clusters:** six sets measured one mesenchymal axis (r = 0.60–0.91 with EMT; §5–§6).
- **Selection:** one representative per cluster, on biology and literature only. "No recurrence, survival or outcome variable was consulted" (§6).
- **Descriptive full-cohort score-matrix condition numbers (Phase 6C):**
  - 10-set panel: **77**;
  - four-set primary: **2.95**, with maximum \|r\| 0.318 and no pair above 0.5;
  - five-set sensitivity (+PI3K): **6.95**, with p53 ↔ PI3K r 0.632.

  These are **not** the Phase 9 model-fit values given in section C.
- **Frozen** 2 November 2026, 21:13 UTC, before any association analysis. Covered by 22 freeze tests (`tests/test_phase6c_m2.py`).

### 2. From probes to pathway scores
1. 54,675 Affymetrix probe sets; 43,163 map to one gene symbol (Fig 2a).
2. These give **21,815 unique genes** (median of probes per gene in the primary representation).
3. The four primary Hallmark sets contain **607 unique member genes** in total (`PHASE6B` §3.2.1).
4. Each pathway score is the mean of its member genes' standardised expression.
5. Gene standardisation (ddof = 1, amendment M2-A1) and every other estimated constant are re-estimated **inside each training fold**, per the frozen protocol (`PHASE6B` §3.2.2). Held-out patients never shape the features.

### 3. Why a constrained panel rather than transcriptome-wide selection
- **Overfitting:** a null simulation at 136 events gave apparent C 0.92 when the true C was 0.50 (`docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md` §6.2).
- **Parameter budget:** four scores plus about six clinical terms suits 136 events; ten correlated scores do not (`PHASE6C` §5).
- **Interpretability:** one score per axis.

### Descriptive measurement properties (`docs/PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md` §3.1)
- Implied mean inter-gene correlation (r̄): E2F 0.304, EMT 0.241, Wnt 0.041, p53 0.033, PI3K 0.025.
- These are **descriptive measures of within-set expression coherence in this tissue and processing**.
- High coherence is not evidence of pathway activation. Low coherence is not evidence that a pathway is biologically unimportant.
- The production scorer agrees with an independent re-implementation to 3×10⁻¹⁵ (`docs/PHASE8D_CONSOLIDATED_EDA_REVIEW.md`).

| Element | Specification |
|---|---|
| Best visual | **Fig 4** (score distributions, correlations, representation agreement) and **DIAGRAM D2** (panel-selection funnel) |
| Recruiter should understand | You can turn a biological hypothesis into a defensible feature set while controlling leakage, collinearity and sample size |
| Avoid | Describing the panel as outcome-optimised; describing the scores as built from the whole gene matrix rather than the 607 member genes; mixing descriptive and model-fit condition numbers |

---

## C. Pathway-by-pathway interpretation

**Scientific message.** For each pathway, connect biology, measurement properties, estimate, uncertainty and robustness. Then state what follows and what does not.

**Sources:**
- **Primary median HRs:** `reports/phase9/execution/phase9_execution_results.json` (`apparent.M2_primary`; λ 0.351).
- **S-7 HRs:** `reports/phase9/execution_s7/phase9_s7_results.json` (`apparent.M2_primary`; λ 0.152).
- **Representation agreement:** `reports/phase8c/representation_sensitivity.csv`.
- **Batch-stratified HRs:** execution JSON `batch_sensitivity.stratified_fit` (exponentiated coefficients; also SRF §11).
- **Leave-one-batch-out (LOBO):** execution JSON `batch_sensitivity.leave_one_batch_out`; FR §12.
- **MMR association:** `reports/phase8c/score_vs_annotation.csv` (`mmr.status`).

All HRs are per one training-fold SD of score.

| | **EMT** | **E2F targets** | **Wnt/β-catenin** | **p53 pathway** |
|---|---|---|---|---|
| Biological axis (Textbook) | Mesenchymal/stromal state | Proliferation | Initiating driver pathway | DNA-damage checkpoint |
| Member genes | 198 | 195 | 41 | 194 |
| Within-set coherence r̄ (Observed, descriptive) | 0.241 | 0.304 | 0.041 | 0.033 |
| **Primary HR, median representation (95 % Wald CI)** | 1.070 (0.951–1.202) | 0.906 (0.807–1.016) | 1.016 (0.906–1.138) | 0.995 (0.882–1.123) |
| **S-7 HR, representative-IQR (95 % Wald CI)** | 1.084 (0.937–1.254) | 0.872 (0.758–1.003) | 1.011 (0.879–1.162) | 0.993 (0.852–1.158) |
| Median vs representative-IQR score agreement r | 0.992 | 0.997 | 0.858 (lowest) | 0.972 |
| Batch-stratified HR (**point estimate, no interval reported**) | 1.073 | 0.912 | 1.004 | 1.022 |
| LOBO log-HR range (all-data value; sign changes) | +0.009 to +0.204 (+0.067; 0) | −0.249 to −0.014 (−0.099; 0) | −0.002 to +0.067 (+0.015; 1) | −0.113 to +0.014 (−0.005; 2) |
| MMR association, ε² (q, BH), n = 410 | 0.0394 (1.3×10⁻⁴) | 0.0251 (0.0021) | 0.0187 (0.0070) | **0.1716** (4.7×10⁻¹⁶) |
| **Observed** | Point estimate > 1 in all reported fits; every interval includes 1 | Point estimate < 1 in all reported fits; no primary-model interval excludes 1 (S-7 upper bound 1.003) | Near-zero; least representation-stable score | Near-zero; sign changes twice; largest MMR association among the non-circular tests |
| **Literature** | Compatible with reports linking stromal/mesenchymal expression to poor prognosis (R1–R3); R1 is not independent of GSE39582 | Mixed (R4 vs R5) | Alteration in about 93 % (R6); expression varies by subtype (R1) | TP53 mutation rates differ between non-hypermutated and hypermutated tumours (R6) |
| **Hypothesis / Not tested** | Whether the signal is malignant-cell EMT or non-malignant stromal expression: **not tested** | Direction is not established biology | No interpretation offered for the null estimate | Direction and mechanism of the MMR association: **not tested**; TP53 status not analysed |

**Notes on the cards.**
- **Wald intervals.** Ridge penalisation makes the Wald intervals approximate and descriptive (FR §9).
- **S-7 is a representation sensitivity analysis.** It uses the same arrays and patients with a different probe summary, so it is **not an independent replication**.
- **Batch-stratified values are point estimates only** and are not inferentially confirmed effects.
- **LOBO caveat:** omitted batches range from 14 to 103 patients, so the coefficient spread combines batch identity with unequal sample size (FR §12).
- **MMR associations** are exploratory, outcome-free Kruskal–Wallis tests: all four q < 0.05.
  - p53's ε² (0.1716) is the largest among all non-circular score–annotation tests. The circular annotation `cit.molecularsubtype` gives larger values and is excluded.
  - An ε² is an association effect size and **must not be compared with a hazard ratio**. It demonstrates no mechanism.
- **Coherence and stability (Hypothesis).** The two coherent sets (E2F, EMT) were the two with no LOBO sign changes. The two low-coherence sets (Wnt, p53) changed sign. This co-occurrence across four pathways **does not show that coherence caused stability**. Wnt and p53 coefficients were close to zero (+0.015 and −0.005), and near-zero estimates change sign easily regardless of coherence.

**Five-set sensitivity model (sensitivity only; FR §9.1).**
- E2F HR **0.862 (0.746–0.996)**, λ 0.152.
- This interval excluding 1 **coincided with a collinearity-driven shift** after PI3K was added:
  - p53 ↔ PI3K r 0.636;
  - Phase 9 model-fit condition number 2.9 → 7.1;
  - molecular SEs inflated by a factor of 1.20–1.34;
  - E2F coefficient −0.099 → −0.149.
- The panel was fixed as sensitivity-only before any fit. Its corrected C (0.6066) does not differ from the primary model's (0.6059) beyond Monte Carlo error.
- It is **not promoted**.

| Element | Specification |
|---|---|
| Best visual | **DIAGRAM D3** (pathway cards), supported by **Fig 6**, **Fig 9a–b** and **Fig 4c** |
| Recruiter should understand | You interpret molecular estimates in biological context and check measurement and robustness first |
| Avoid | Any per-pathway prognostic or significance claim; comparing an ε² with a hazard ratio; using set coherence or Wnt biology to explain the null Wnt estimate (both explanations are unsupported); promoting the five-set model because one interval excludes 1; describing S-7 as replication |

---

## D. Integrated biological interpretation

**Scientific message.** Together, the four pathways show why coefficient-level biology is hard to read from bulk transcriptomic scores in a cohort of this size. The project documents several complications; it does not measure or resolve them.

1. **Molecular subtype (Observed, descriptive).**
   - The dominant expression axis in this cohort is MMR/CIMP/CIN/BRAF/tumour location, not stage or partition (`PHASE8D` executive summary).
   - All four scores associate with MMR status, p53 most strongly (ε² 0.1716).
   - The frozen panel was not designed around this axis, and MMR status was not a model covariate.
2. **Technical batch (Observed).**
   - In exploratory tests, batch and RNA-extraction **location** associations with the scores did not survive a permutation test stratified by the ComBat model covariate `combat_mod` (0 of 10; `PHASE8C` §5.3).
   - `combat_mod` matches recorded MMR status in all 410 samples where both are available (0 discordant; `PHASE8B`). Its provenance is unresolved, and it must not be described as the MMR variable.
   - Four QC-flag **dispersion** tests survived.
   - **This does not show that technical batch effects were eliminated.** In Phase 9, LOBO refits still changed the signs of Wnt and p53, and residual batch-associated structure remains after ComBat (PC3, FR §12).
3. **Cell composition (Not tested).** Bulk expression cannot separate tumour-cell expression from stromal, fibroblast, endothelial or immune contributions. The literature attributes much bulk mesenchymal signal to stroma (R2, R3), but this project neither measured nor deconvolved cell composition. It remains an alternative explanation, especially for EMT.
4. **Score construction (Hypothesis).** Coherent and low-coherence sets behaved differently under LOBO, but causation is not shown, and near-zero coefficients are an alternative explanation (section C).
5. **Power and collinearity (Observed).**
   - With 136 events, effects of the observed size (HRs about 0.9–1.1 per SD) cannot be distinguished from zero.
   - In the five-set model, an excluding interval coincided with added collinearity.

**Conclusion to present.**
- In this cohort, bulk pathway scores varied with molecular subtype (descriptive), and their coefficients were sensitive to batch composition, with 136 events limiting precision.
- Cell composition was not assessed and remains an alternative explanation.
- Separating these influences would need deconvolution, subtype-aware models and external data.

| Element | Specification |
|---|---|
| Best visual | **DIAGRAM D4** (conceptual schematic), beside **Fig 3** (residual batch structure) and **Fig 4b** (score correlations) |
| Recruiter should understand | You treat tumour heterogeneity, stroma, subtype and batch as real obstacles to interpretation |
| Avoid | Any statement that cell composition was measured or shown to affect the results; that the EMT score measures stromal content; that batch effects were removed, resolved or attributed to MMR; that `combat_mod` is an MMR variable; reading LOBO spread as a pure batch effect (it is confounded by batch size); treating batch-stratified point estimates as confirmed effects; that the panel captured the dominant expression axis |

---

## E. Findings and descriptive observations

**Findings from this cohort:**

| # | Finding | Evidence | Figure |
|---|---|---|---|
| E1 | **No demonstrated reproducible improvement in discrimination.** 452 patients, 136 recurrences. **Optimism-corrected ΔC = +0.004316** (plan v2.1 §3.1). Five verified run-specific 95 % percentile intervals all include zero: [−0.0086, +0.0519], [−0.0090, +0.0469], [−0.0093, +0.0484], [−0.0091, +0.0518], [−0.0091, +0.0494]. The runs differ only in random seed. The intervals do **not** establish that the true effect is zero | `corrected_interval_results.json`; FR §8; verifier PASS 60/60 | **Fig 8** |
| E2 | **Most of the apparent gain was optimism.** Apparent C 0.6242 (M1) vs 0.6382 (M2), apparent ΔC +0.0140. Optimism-corrected C 0.6016 vs 0.6059 | FR §7.3 | **Fig 7** |
| E3 | **A different quantity looked positive.** The paired out-of-original ΔC interval was [+0.0002, +0.0240] (median) and [−0.0001, +0.0224] (S-7). It is not the corrected ΔC interval, and the two must not be interchanged | FR §8.1 | (table) |
| E4 | **No primary-model pathway HR interval excluded 1** | execution and S-7 JSON | **Fig 6** |
| E5 | **Molecular coefficients are batch-sensitive.** LOBO changes Wnt's sign once and p53's twice; apparent (not corrected) LOBO C 0.614–0.659; stage HR 1.94 → 1.68 with batch strata | FR §12; execution JSON | **Fig 9** |
| E6 | **Residual overfitting.** Corrected calibration slopes 0.833 (M1), 0.862 (M2) | FR §10 | **Fig 9d** |
| E7 | **Confirmatory external validation was infeasible with available public cohorts.** GSE14333 has 48 events and GSE33113 18 (exploratory only), against about 338 required at an assumed C of 0.65 (indicative; an approximation of a published framework). GSE17536 overlaps GSE14333 (129 unambiguous patient pairs) and is not independent | `PHASE4_PREACQUISITION.md` D1, D3; `PHASE5` §7–§8 | **DIAGRAM D5** |

**Descriptive observations (not findings or discoveries):**
- **O1:** Pathway sets differ in within-set coherence (section B). Its co-occurrence with sign stability is a hypothesis (section C).
- **O2:** All four scores associate with MMR status in exploratory, outcome-free tests (section C).
- **O3:** The five-set E2F interval excluding 1 coincided with collinearity (section C).
- **O4:** Informative censoring. The 46 deaths without recurrence (censored in the deposited data) were older (median 75.0 vs 68.0), less often stage III (32.6 % vs 38.7 %) and less often treated (23.9 % vs 42.1 %) (**Fig S2**). Not resolved.

**Estimand and endpoint, stated on the result slides.**
- The models estimate the **cause-specific hazard of recurrence**, not absolute recurrence probability or cumulative incidence.
- The `rfs.delay` unit is undocumented (provisional assumption: months; `PHASE1` P2), so only rank-based quantities are reported.
- Competing-risk cumulative incidence and composite-DFS analyses were specified but not performed.
- The CIT partition is internal, not external, validation.
- No external validation or clinical utility was established.

| Element | Specification |
|---|---|
| Recruiter should understand | Which results are robust (E1–E6) and which are descriptive (O1–O4) |
| Avoid | Interchanging corrected ΔC, paired out-of-original ΔC, apparent ΔC and corrected C; stating that the pathways have no role; biomarker language |

---

## F. Research decisions and audit trail

**Scope of contribution.**
- I set the scientific direction and constraints, made or approved the methodological decisions, commissioned audits, reviewed their findings and decided what could be claimed.
- The records generally attribute decisions to the project owner and do not record who first identified each problem.

**My contribution, stated with accurate verbs.**

| # | Decision | What happened | What I can say | Record |
|---|---|---|---|---|
| DC1 | Penalty scope | Clinical block left unpenalised (scale-equivariant); ridge only on pathway scores | "I **approved** amendment P6B-A1 and can explain why ridge penalties are not scale-invariant" | `PHASE6B` §3.1.1 ("approved by the project owner, 2026-11-03") |
| DC2 | CIT partition | Never described as external validation | "I **decided** the CIT partition would be treated as internal only" | `PHASE1` A6 (approved decision) |
| DC3 | Random-stream deviation | Corrected driver; full re-run; as-run results superseded | "I **adjudicated** the seed deviation and **authorised** the corrected re-run" | `PHASE9_THETA_ALIGNMENT_AUDIT.md` (owner Decisions 1–2); FR §6 |
| DC4 | Monte Carlo precision gate | Rule-driven top-up from 500 to 750 outer resamples | "I **chose** to complete the pre-specified top-up rather than report early" | `PHASE9_DECISION13_FREEZE_PROCEDURE.md` ("The owner chose to complete…") |
| DC5 | Restraint | Five-set E2F and the paired out-of-original interval kept out of the conclusions | "I **rejected** promoting a collinear sensitivity result and an interval for the wrong estimand" | FR §8.1, §9.1 |
| DC6 | Audits | θ̂ alignment, independent numerical audit, citation audit, interpretation review | "I **commissioned** these audits, **reviewed** the findings and **approved** the corrections" | `OWNER_DECISION_REGISTER.md`; `BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md` §5; `PHASE9_SCIENTIFIC_INTERPRETATION_REVIEW.md` |
| DC7 | Overall direction | Pre-specification, outcome-blind panel, restrained interpretation | "I **directed** the project to fix every choice before outcomes and not to overstate evidence" | `PHASE6B` §2, §3.2.4; `PHASE6C` header |

**Problems addressed during the project.** Collaborative; identification not attributed in the records.

| Problem | Resolution | Record |
|---|---|---|
| Undocumented time unit and death handling | Labelled as provisional; rank-based analysis only | `PHASE1` P1–P3 |
| Patient overlap between public cohorts | GSE17536 excluded; subtraction prohibited | `PHASE4_PREACQUISITION.md` D1 |
| External validation feasibility | NO-GO; an earlier overclaim withdrawn in writing | `PHASE5` §0.1, §8 |
| Library-default SD convention | Pinned (M2-A1) | `PHASE6C` §0 |
| No interval for the headline gain | Nested bootstrap planned and approved | `PHASE9_EXECUTION_REPORT.md` §15 |
| Bypassable pre-execution gate | Found by adversarial test; remediated | `PHASE9_DECISION13_FREEZE_PROCEDURE.md` |
| Sample-level sex discordance (8 arrays) | Documented, not silently removed | `PHASE8B` §7 |

**Do not claim** that I personally:
- discovered every issue or found the original RNG defect;
- designed the nested bootstrap alone;
- performed the independent audits;
- matched patient records;
- authored all scripts or figures.

Do not describe the independent numerical audit as external replication or model refitting. It recomputed reported quantities from recorded checkpoint values within the project.

| Element | Specification |
|---|---|
| Best visual | Decision table (DC1–DC7) beside the week timeline from the portfolio page |
| Recruiter should understand | Your contribution is scientific direction, methodological judgement, protocol decisions, audit commissioning, interpretation and refusal to overstate evidence |

---

## G. Future biological directions (future work, not achievements)

Each would be a **new, separately pre-specified analysis**. None has been started.

| # | Question | Rationale | Design sketch |
|---|---|---|---|
| G1 | Is the EMT-direction signal malignant-cell or stromal? | Consistent EMT direction; literature attributes bulk mesenchymal signal to stroma (R2, R3); **not tested** here | Estimate stromal and immune fractions by deconvolution; test whether any EMT association persists after adjustment |
| G2 | Does pathway information add value within molecular subtypes? | Dominant MMR/CIMP/BRAF axis; all scores associate with MMR status | Pre-specified MMR- or CMS-stratified or adjusted models, with event counts per stratum checked first |
| G3 | Does any signal transport to an independent cohort? | Original aim; blocked by event counts | A cohort with about 340 or more stage II–III recurrence events, a documented time unit and death handling, and uniform processing; frozen M1/M2 tested unchanged |

---

## Diagrams D1–D5 (originally proposed; now implemented)

Each diagram is built by `portfolio/src/make_diagrams.py` from frozen outputs, as SVG with editable text and as a 200-dpi PNG with the same content.

| # | File (`portfolio/diagrams/`) | Planned as | Status and form | Data (all existing) | Caption and limitation |
|---|---|---|---|---|---|
| D1 | `D1_biological_rationale` | Pathway lanes | **Implemented; conceptual schematic.** Four parallel lanes, no arrows between lanes | Gene counts 198/195/41/194 | "Conceptual schematic: lanes are parallel axes, not a causal sequence. Scores are bulk mRNA averages and do not measure the activity shown" |
| D2 | `D2_processing_feature_engineering` | Panel-selection funnel | **Implemented.** Processing flow plus outcome-blind selection strip | 54,675 → 43,163 → 21,815 → 607 → 4; condition numbers 77 / 2.95 / 6.95 (descriptive) | "Condition numbers are descriptive score-matrix values, not model-fit values. No outcome consulted; frozen 2 Nov 2026, 21:13 UTC" |
| D3 | `D3_pathway_interpretation` | Pathway cards | **Implemented** to the section C table: primary and S-7 intervals, stratified point estimates, LOBO ranges, MMR ε² and q, Observed / Literature / Not tested panels | All values in section C | "Stratified HRs are point estimates; LOBO mixes batch and sample size; MMR is an association, not a mechanism; coherence is not activation; not a biomarker" |
| D4 | `D4_heterogeneity_limitations` | What a bulk score may contain | **Implemented; conceptual schematic.** Four unweighted inputs with status badges, plus an observed-results panel | MMR ε² and q; `combat_mod` caveat; LOBO sign changes | "Conceptual: no contribution was estimated; lines carry no weights or causal direction" |
| D5 | `D5_model_comparison_uncertainty` | Event-count bar (extended) | **Implemented.** Model design, apparent vs corrected C, five corrected intervals, key numbers, event-count bar | C values, θ̂, five intervals, slopes; 136 / 48 / 18 / about 338 events | "Runs differ only in random seed; including zero does not establish no effect; CIT internal; required events indicative and not combined" |

None of the diagrams may imply that coherence proves activation, that MMR explains all batch effects, that batch adjustment resolved confounding, that scores identify the cellular origin of expression, or that a prognostic biomarker was demonstrated.

---

## Twenty-one-slide plan (as implemented)

Built as `portfolio/presentation/GSE39582_recurrence_portfolio.pptx` (editable) and `.pdf`: 19 main slides and 2 appendix slides (round 7 added slide 18), with speaker notes on every slide and in `portfolio/SPEAKER_NOTES.md`. This replaces the earlier twelve-slide plan.

All slide visuals are original vector drawings built by `portfolio/src/make_slide_visuals.py` (SVG with editable text, plus 200-dpi PNG in `portfolio/diagrams/slides/`). Every number is read from frozen outputs through `portfolio/src/values.py`. No external, stock or licensed images are used. Mechanism panels carry four evidence labels: **Established mechanism** (textbook background), **CRC literature** (R1–R8), **Measured in this study**, and **Not tested / hypothesis**.

| # | Slide | Content | Visual | Labels |
|---|---|---|---|---|
| 1 | **Cover** | Research question; GSE39582; Prasad Akula, MS in Bioinformatics, Northeastern University; November 2026 | V01 conceptual crypt and tumour-gland illustration | Conceptual |
| 2 | **Clinical problem and research question** | Problem → question → optimism trap → approach; 452 / 136 / 4 / 54,675; cause-specific hazard; internal validation | V02 | Observed |
| 3 | **Four biological axes** | Wnt ⊣ destruction complex ⊣ β-catenin → TCF/LEF → MYC, cyclin D; cyclin D–CDK4/6 ⊣ RB ⊣ E2F → G1/S genes; DNA damage → p53 → p21 ⊣ CDK and → MDM2 ⊣ p53; TGF-β → SNAIL/ZEB/TWIST ⊣ CDH1, → vimentin/fibronectin; CAFs also express mesenchymal genes | V03 mechanism map with compartments | Established mechanism |
| 4 | **Wnt/β-catenin** | APC intact vs APC lost; nuclear β-catenin–TCF/LEF; HR 1.016 (0.906–1.138); S-7; coherence 0.041; 1 sign change | V04 | All four |
| 5 | **E2F** | Cyclin D–CDK4/6 phosphorylates RB, releasing E2F; p21 restraint; HR 0.906 (0.807–1.016); five-set 0.862 as sensitivity only | V05 | All four |
| 6 | **p53** | ATM/ATR → p53; MDM2 feedback; p21, BAX/PUMA, GADD45; HR 0.995 (0.882–1.123); MMR ε² 0.1716 | V06 | All four |
| 7 | **EMT** | Tumour-cell EMT vs cancer-associated fibroblasts feeding one bulk score; "An elevated bulk EMT-associated score does not establish tumour-cell EMT"; HR 1.070 (0.951–1.202) | V07 | All four |
| 8 | **What a bulk sample measures** | Tumour cells not separated; stromal, immune, endothelial fractions not measured; mixed tissue → one array → one score | V08 conceptual tissue illustration | Conceptual, Observed |
| 9 | **Dataset and cohort** | 585 → 452 patients; 136 recurrences; 46 censored deaths; time unit undocumented; CIT internal | Fig 1 | Observed |
| 10 | **Preprocessing and features** | 54,675 → 43,163 → 21,815 → 607 → 4; condition numbers 77 / 6.95 / 2.95 (descriptive) | V10 | Observed |
| 11 | **Model architecture** | M1 unpenalised clinical block; M2 + ridge four-set block (λ 0.351); five-set separate, sensitivity only | V11 | Observed |
| 12 | **Nested bootstrap** | 750 outer refits; 100 inner resamples with paired optimism; corrected ΔC; 5 runs as Monte Carlo replicates; precision 4.2 % ≤ 5 %; verifier 60/60 | V12 | Observed |
| 13 | **Central result** | Apparent vs corrected C; five intervals include zero; θ̂ +0.004316; calibration slopes; conclusion sentence | V13 | Observed |
| 14 | **Pathway hazard ratios** | Primary and S-7 forest plot, log axis, line at 1; five-set E2F in a separate sensitivity-only box | V14 | Observed |
| 15 | **Robustness** | Per-batch leave-one-batch-out coefficients vs omitted batch size; MMR ε²; calibration slopes | V15 | Observed |
| 16 | **What this project established** | Four contributions: pathway-audit framework for bulk and single-cell data (13/13 tests; CI passing; Docker not built); EMT score tracks stroma in four bulk cohorts against cohort-specific random baselines and in three single-cell cohorts; p53 score tracks TP53 (2 of 3 external cohorts, robust to hypermutation); corrected ΔC +0.004316 with all intervals including zero. Compact "What remains unresolved" strip; full detail in the limitation register | V16 | Replicated (post-hoc, registered); tested software; pre-specified |
| 17 | **Decisions and audit trail** | DC1–DC7 on the project timeline; verifier 60/60; 138 quantities within 2.2×10⁻¹⁶; 1,097 tests | V17 | — |
| 18 | **From finding to solution** | Gene-disjoint EMT–fibroblast ρ in four cohorts (GSE103479 failed QC); the five pathway-audit checks; case-study results | V19 | Post-hoc, outcome-free |
| 19 | **Conclusion and the next scientific questions** | Asked → found → asked next → found → delivered; take-home; Q1–Q3 (rare or partial tumour-cell EMT states; pathway activity; independent predictive validation) with how and why, marked not performed | V18 | Not tested |
| A1 | **References and figure sources** | R1–R8 with DOIs; visual sources | Text | — |
| A2 | **Quality control and informative censoring** | Soft QC flags; batch–partition confounding; censored deaths profile | Fig 2, Fig S2 | Observed |

**Design rules.**
- One colour for the primary result, one for sensitivity, grey for apparent or diagnostic values.
- Each visual carries its caveat on the visual itself.
- Never crop a figure so that a caveat or source line disappears.

---

## Recommended final narrative (about two and a half minutes)

Exact figures in brackets are for the slide; say the rounded version aloud.

"I asked whether four predefined Hallmark pathway-expression scores add recurrence-discrimination information beyond clinical variables in stage II–III colon cancer, using 452 patients from GSE39582 with 136 recurrences.

With only 136 events, selecting genes from the whole transcriptome would overfit, so the panel was constrained and fixed before any outcome was examined: EMT, E2F targets, Wnt/β-catenin and p53, one per biological axis. Of 54,675 Affymetrix probes, 43,163 map to a single gene, giving 21,815 unique genes. The four sets contain 607 of those genes. Their scores were calculated and standardised inside each training fold only.

I compared a Cox model with age, sex, stage, tumour location and adjuvant chemotherapy against the same model plus the four ridge-penalised scores. The outcome is the cause-specific hazard of recurrence, not absolute recurrence risk.

The apparent C-index gain of about 0.014 was mostly optimism. After correction it was +0.0043. Five nested-bootstrap runs, differing only in random seed, gave 95 % intervals from about −0.009 to +0.05, and every one included zero. That does not prove the effect is zero, or that these pathways have no biological role. Calibration slopes below one showed residual overfitting in both models, and no primary pathway hazard ratio excluded one.

In exploratory, outcome-free tests on the 410 patients with recorded mismatch-repair status, all four scores were associated with it, p53 most strongly [ε² and BH q: EMT 0.0394, 1.3×10⁻⁴; E2F 0.0251, 0.0021; Wnt 0.0187, 0.0070; p53 0.1716, 4.7×10⁻¹⁶]. These are associations, not mechanisms. Batch associations with the scores did not survive stratification by the ComBat covariate `combat_mod`. That covariate matches recorded MMR status in all 410 samples, but it is not an established MMR variable. Leaving out single batches still flipped the Wnt and p53 coefficients. So batch effects were characterised, not resolved, and cell composition was never measured. Stromal content remains an untested alternative explanation for the EMT direction.

Validation was internal only. The CIT partition is not an external test, and no available external cohort had enough events.

I directed the project to fix every analytic choice before outcomes, adjudicated the random-seed deviation and authorised the corrected re-run, chose to complete the pre-specified precision top-up, kept the collinear five-set result out of the conclusions, and commissioned and reviewed the audits.

Next, I would estimate stromal and immune fractions by deconvolution and test whether the EMT signal persists, test pathway value within mismatch-repair or consensus subtypes, and validate the frozen models in an external cohort with roughly 340 or more recurrence events."

---

## Claims to avoid everywhere
- Pathways as prognostic biomarkers, or any significant predictive gain.
- "No effect" or equivalence.
- Activation, mechanism, causality or mutation inference from scores.
- That cell composition was measured, deconvolved or shown to confound.
- That batch effects were eliminated, or explained by MMR.
- That `combat_mod` is the MMR variable.
- Comparing ε² with hazard ratios.
- Clinical utility, absolute risk, external validation or transportability.
- Five independent experiments: they are five runs differing only in random seed.
- External replication or model refitting by the independent audit.
- The CIT partition as external validation.
- The five-set E2F result as a finding.
- Personal origination of work the records do not attribute.

## Integrity
- The portfolio implementation added files under `portfolio/` only; no frozen output, figure or report changed.
- No analysis was run.
- No frozen output, figure, caption, configuration, checkpoint or other document was changed.
- The public-release audit was not touched.
- Nothing was committed or pushed.
