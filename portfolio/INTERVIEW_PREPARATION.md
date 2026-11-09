# Interview preparation: GSE39582 colon-cancer recurrence study

Audience: me, preparing for recruiter screens and technical interviews. The answers are grounded in the project records and must stay consistent with `RESEARCH_REPORT.md`.

---

## Part 1. The modelling strategy at three levels

### Recruiter level (30 seconds)

I tested whether gene-expression information from tumours helps predict which colon cancers come back after surgery, beyond what doctors already know from stage and other clinical facts.

I chose four biological programmes in advance and built a careful comparison between a clinical-only model and the same model plus those programmes. I also measured how much any apparent improvement was just the model fitting its own data.

The gene information did not reliably improve prediction. The apparent improvement was mostly an artefact, and I documented exactly why the biology is hard to read from this kind of data.

### Scientific level (2 minutes)

**Question.** Do four outcome-blind Hallmark pathway scores (EMT, E2F targets, Wnt/β-catenin, p53) add recurrence discrimination beyond a clinical Cox model in 452 stage II–III colon cancers with 136 recurrences?

**Why this design:**
- **Pre-specification.** Selecting genes from about 22,000 after seeing outcomes overfits: under no true signal, a simulation gave an apparent C of 0.92 for a true 0.50. So the features were fixed in advance, one representative per biological axis, chosen from correlation structure and literature.
- **Nested comparison.** M1 is clinical only. M2 adds the scores. The clinical block is identical and unpenalised in both, so the only difference is the molecular block.
- **Optimism correction.** Apparent C rewards complexity. I estimated the gain after removing optimism, and put a valid interval on that corrected gain with a nested bootstrap.

**Findings:**
- The apparent gain of +0.014 shrank to +0.0043 after correction.
- All five 95 % intervals include zero (about −0.009 to +0.05).
- No pathway hazard ratio was distinguishable from 1.
- Coefficients changed sign when single processing batches were omitted.
- The scores tracked mismatch-repair subtype.

### Technical level (5 minutes)

**Survival model:**
- Cox proportional hazards models the instantaneous recurrence rate as a baseline hazard times exp(linear predictor).
- Because deaths without recurrence are censored, the target is the **cause-specific** hazard of recurrence. That is a valid rate, but it cannot be turned into an absolute recurrence probability without modelling the competing risk of death.
- The time unit is undocumented, so I used only rank-based measures (C-index, hazard ratios).

**Feature blocks and penalisation:**
- **Clinical block:** age as a restricted cubic spline, sex, stage, location, chemotherapy, unpenalised.
- **Molecular block:** four scores with a ridge (L2) penalty. λ is chosen by five-fold event-stratified cross-validation inside each training set.
- Unpenalised clinical terms avoid making shrinkage depend on measurement units (ridge is not scale-invariant). They also keep M1 exactly nested in M2.
- Shrinkage of the molecular block makes ΔC conservative.

**Leakage control:**
- Each pathway score is the mean of member-gene z-scores.
- The gene means and SDs are estimated on the training portion only and applied unchanged to held-out samples.
- λ tuning happens inside the training portion too.

**Optimism and the nested bootstrap:**
- Optimism is the gap between performance on the fitting data and on new data from the same population.
- **Outer loop:** draw 750 patient-level, event-stratified bootstrap resamples.
- **For each outer resample:**
  1. refit both models and compute the apparent ΔC;
  2. estimate optimism from 100 inner bootstrap resamples, refitting both models each time, paired on the same resamples;
  3. subtract to get the corrected ΔC for that replicate.
- **Interval:** the 2.5th and 97.5th percentiles across the 750 replicates.
- **Point estimate:** the full-data optimism-corrected ΔC, θ̂ = +0.004316.

**Monte Carlo precision:**
- Five runs differ only in random seed, so they measure simulation noise, not sampling uncertainty.
- A pre-specified rule required the largest endpoint Monte Carlo SE to be ≤ 5 % of interval width. 500 outer resamples failed it, so each run was topped up to 750, reaching 4.2 %.
- Every random draw has its own stream, and checkpoints are hash-linked and resumable.

**Calibration slope:**
- Regress the outcome on the model's linear predictor in new data; a slope of 1 is ideal.
- Corrected slopes of 0.833 and 0.862 mean predictions are too extreme and would need about 15 % shrinkage.

**Verification:**
- A frozen verifier (60 checks) gated reporting.
- A separate standard-library implementation recomputed 138 reported quantities from the raw checkpoints, agreeing within 2.2×10⁻¹⁶.
- That checks the arithmetic from recorded model outputs. It is not an external replication and does not refit models.

---

## Part 2. Question bank

### Study design and statistics

**1. Why not use machine learning on all 22,000 genes?**
With 136 events, high-dimensional selection overfits badly; a null simulation gave apparent C 0.92 when the truth was 0.5. A small, pre-specified biological panel keeps the number of parameters proportionate to the events and makes coefficients interpretable. The cost is that a real signal outside the four pathways would be missed; that was a deliberate trade-off.

**2. What exactly is "optimism", and why correct it?**
Optimism is how much better a model looks on the data it was fitted to than it would on new patients from the same population. Bigger models gain more from fitting noise. Here optimism was 0.023 for the clinical model and 0.032 for the pathway model, which erased most of the apparent advantage.

**3. Why a nested bootstrap rather than a single bootstrap or a train-test split?**
- A single bootstrap gives an optimism-corrected point estimate but not a valid interval for it. Each outer replicate needs its own optimism estimate, hence the inner loop.
- A train-test split wastes events. With 136 events, a 30 % test set has about 40 events, and its C would be very noisy.

**4. Your interval includes zero. Did you show the pathways are useless?**
No. An interval from −0.009 to +0.05 is compatible with a small loss and with a modest gain. Absence of evidence is not evidence of absence, and the study was not designed to show equivalence.

**5. Why five runs? Aren't five intervals that agree strong evidence?**
They differ only in random seed. They show the result is not a Monte Carlo accident, but they are replicates of one analysis of one cohort, not five independent experiments.

**6. You reported an earlier interval of [+0.0002, +0.024] that excludes zero. Why not use it?**
It is a different quantity, the paired out-of-original ΔC. It mixes in-sample and out-of-sample performance and is not an interval for the corrected ΔC. It also included zero under the alternative probe representation. Fixing the estimand before looking prevented a misleading positive claim.

**7. Why is the clinical block unpenalised?**
Ridge penalties depend on the scale of the variables. Penalising a year-scale age spline alongside binary variables would make shrinkage depend on units. Leaving clinical terms unpenalised keeps the comparison fair, and makes the clinical model exactly the pathway model with pathway coefficients at zero.

**8. Isn't penalising only the molecular block biased against it?**
Yes, slightly: shrinkage makes ΔC conservative, and I state that. But without penalisation the molecular block would be more optimistic, and the corrected comparison would still remove that. The design favours not overstating a gain.

**9. What does a calibration slope of 0.86 mean in practice?**
The model's risk predictions are too spread out. High-risk patients are predicted riskier, and low-risk patients safer, than they really are. You would shrink the linear predictor by about 14 % before using it. The clinical model shows the same problem (0.833), so even six terms overfit with 136 events.

**10. Why Harrell's C and not time-dependent AUC or Brier score?**
The recurrence-time unit is undocumented, so any horizon-based measure (5-year AUC, Brier score at a time point) would rest on an unverified assumption. Harrell's C depends only on ranking. Its limitations are that it is insensitive to small added predictors and says nothing about calibration or clinical utility.

**11. How did you handle deaths without recurrence?**
They are censored, so the primary model estimates the cause-specific hazard. That is valid for aetiology but cannot give an absolute recurrence probability. The 46 censored deaths were older and less often treated, so the censoring is informative. In the registered follow-up (X2), the Aalen–Johansen cumulative incidence of recurrence at 60 time units was 30.7 %, against 31.3 % from 1 − Kaplan–Meier, and Fine–Gray subdistribution hazard ratios were close to the cause-specific ones. So competing deaths did not change the conclusions here (see question 31).

**12. How do you know there was no data leakage?**
Every estimated quantity (gene means and SDs, score scaling, λ) is re-estimated inside each training resample and applied unchanged to held-out samples. The production scorer was checked against an independent re-implementation (agreement within 3×10⁻¹⁵) and passed fold-isolation tests in which held-out data were deleted or corrupted.

**13. What was the random-seed issue?**
An earlier version of the nested bootstrap reused a resample's random state inside cross-validation, so draws that should have been independent were linked. I adjudicated that deviation, authorised a full corrected re-run with a separate random stream per draw, and kept the earlier results only as a superseded record. The corrected intervals were close to the earlier ones, but only the corrected ones are reported.

### Biology

**14. Why these four pathways?**
They represent distinct axes of colorectal-cancer biology:

| Pathway | Axis |
|---|---|
| Wnt/β-catenin | The initiating driver |
| E2F targets | Proliferation |
| p53 pathway | DNA-damage checkpoint |
| EMT | Mesenchymal/stromal state linked to the poor-prognosis CMS4 subtype |

They were chosen as one representative per correlation cluster without outcome data. PI3K/AKT/mTOR was kept as a sensitivity addition because it was highly correlated with p53.

**15. EMT had a hazard ratio above 1. Does that mean EMT drives recurrence?**
No. The interval includes 1, so there is no demonstrated association. Even if there were, a bulk EMT score mixes tumour cells undergoing EMT with fibroblasts that express the same mesenchymal genes. Published work attributes much bulk mesenchymal signal to stroma. Cell composition was not measured, but in the registered follow-up the EMT score tracked estimated fibroblast content (ρ 0.94–0.98 with shared genes removed) in four cohorts. So in bulk tissue it behaves as a stromal index, and tumour-cell EMT cannot be inferred from it.

**16. Wnt is altered in over 90 % of colorectal cancers. Why no signal?**
The 93 % figure is genomic alteration, not measured pathway activity, and the score averages 41 target genes that are loosely co-expressed in this tissue. The estimate sat near zero and changed sign under batch omission. I did not offer an explanation for the null, because none was tested.

**17. The p53 score was strongly associated with MMR status. What does that tell you?**
The score carries molecular-subtype structure. MMR-deficient and proficient tumours differ in TP53 mutation rates and expression programmes. The association is exploratory, says nothing about mechanism, and is not an association with recurrence: p53's hazard ratio was 0.995. An ε² and a hazard ratio measure different things and should not be compared.

**18. What would you need to know whether the EMT signal is tumour-intrinsic?**
- Computational deconvolution: done post hoc (ESTIMATE, MCP-counter). It shows the score tracks stromal estimates, but cannot locate expression in cells.
- Adjusting the EMT–recurrence association for stromal fraction: tried (X7), but EMT and the stromal estimate are near-collinear (VIF > 10), so their coefficients cannot be separated.
- Single-cell RNA-seq: done, as a registered extension in three public cohorts (`docs/SINGLE_CELL_RESULTS.md`). EMT-set transcripts come mainly from stromal cells (5.1× the expression-matched share; epithelial 0.45×), and no coherent tumour-cell EMT programme was detected (see questions 39 and 40).
- Spatial transcriptomics or cell-resolved protein staining.

**19. Isn't a coherence of 0.04 for Wnt evidence that the pathway is irrelevant?**
No. Coherence describes how tightly the member genes co-vary in this tissue and processing. Low coherence means the averaged score is a weak summary. It says nothing about the pathway's biological importance, and high coherence is not evidence of activation.

### Data and batch

**20. GSE39582 was already ComBat-corrected. Why worry about batch?**
Residual batch-associated structure remained (on PC3), and leaving out single processing batches flipped the signs of two pathway coefficients. Batch is also confounded with the discovery/CIT partition. I characterised these effects but did not correct them, and I state that they are unresolved.

**21. Couldn't the leave-one-batch-out spread just reflect batch size?**
Partly, yes. Omitted batches range from 14 to 103 patients, so removing a large batch removes more data as well as a batch. The spread mixes batch identity with sample size, which is why I describe coefficients as batch-sensitive rather than attributing the variation to batch alone.

**22. Why is the CIT partition not external validation?**
It comes from the same deposit, shares processing history, and is confounded with batch. It also has only 29 events. Calling it external validation would overstate generalisability.

**23. Why not validate on another public cohort?**
I checked, twice. The registered eligibility audit (X8) screened 554 GEO hits, assessed 14 cohorts in detail and found none eligible: the largest had 124 recurrences against about 337 needed for a usefully precise estimate, and the others lacked stage, endpoint or platform compatibility. Validating on an underpowered cohort would have produced a misleading confirmation or a misleading failure.

### Ownership and process

**24. What was your personal contribution?**
Scientific direction and methodological decisions, all recorded in the project's decision records. I:
- fixed the question and constraints before outcomes;
- chose an outcome-blind panel and a clean clinical-versus-molecular comparison;
- approved keeping clinical terms unpenalised;
- restricted the CIT partition to internal validation;
- adjudicated the random-seed deviation and authorised the corrected re-run;
- chose to complete the pre-specified precision top-up;
- kept a tempting but wrong-estimand interval and a collinear sensitivity result out of the conclusions;
- commissioned and reviewed verification and citation audits.

**25. Which parts did you implement yourself, and how?**
Expect this question: the repository is large and was built in about a week. Answer accurately and specifically: which components you wrote, which you reviewed, and which tools or collaborators generated or executed code, and how you checked their output. The project's value rests on decisions and verification, so describe the workflow as it actually was, and be ready to walk through any file line by line.

**26. What would you do differently?**
- Plan a competing-risk analysis from the start.
- Obtain a documented time unit before modelling.
- Include cell-composition estimates in the protocol.
- Pre-register the batch-sensitivity analysis as part of the primary design.

**27. What is the main lesson?**
Validation design matters more than model complexity. A carefully verified null result is a real result, and separating the estimand from convenient-looking statistics is what keeps it honest.

### The follow-up work and the tool

**28. Why does the ridge penalty apply only to the pathway block?**
Two reasons. First, fairness and scale: ridge shrinkage depends on variable units, and the clinical block mixes an age spline with binary terms, so penalising it would make the result depend on coding choices. Second, nesting: with the clinical block unpenalised, M1 is exactly M2 with the pathway coefficients at zero, so any change in C is attributable to the pathways. The cost is that shrinking only the pathway block makes ΔC slightly conservative, which I state (questions 7 and 8).

**29. What exactly does optimism correction remove, and what does it not remove?**
For each bootstrap resample, the whole pipeline is refitted: spline knots, gene scaling, score scaling, λ tuning and the Cox fit. Optimism is the average of (C on the resample it was fitted to) minus (C of that model on the original data). Subtracting it removes the gain a model gets from fitting noise in this cohort. It does **not** remove sampling uncertainty (hence the nested bootstrap for an interval), and it says nothing about other populations, batches or platforms: it is internal validation only.

**30. Why is the TCGA random-gene-set baseline about 0.64, when it is near zero in the microarray cohorts?**
A score is the mean of z-scored genes. In bulk RNA-seq, a large share of genes move together along one dominant sample-wide axis, most plausibly tumour purity versus stromal and immune content, plus normalisation effects. The average of about 190 random genes therefore picks up part of that axis, which also correlates with fibroblast markers. In GSE39582 the random median was −0.14 (upper bound 0.05); in TCGA it was 0.64 (upper bound 0.76). That is why every cohort is compared with its own random baseline: TCGA's 0.98 still exceeds 0.76, but it is less informative than its size suggests. The cause of the shared axis was not formally decomposed, so I call it the most plausible explanation, not a finding.

**31. What would Fine–Gray change compared with the cause-specific model?**
The cause-specific model treats deaths without recurrence as censored and estimates the recurrence rate among patients still alive and recurrence-free: the right question for aetiology. Fine–Gray keeps patients who died in the risk set and models the subdistribution hazard, which maps directly onto the cumulative incidence of recurrence: the right question for absolute-risk prediction. With 46 competing deaths, the two can differ. Here they did not differ materially: the Aalen–Johansen incidence at 60 time units was 30.7 % against 31.3 % naive, and the subdistribution and cause-specific hazard ratios were similar. Absolute risks still are not reported, because the time unit is undocumented.

**32. Why the gene-disjoint test?**
The EMT set and the MCP-counter fibroblast markers share six genes, which would inflate their correlation by construction. Removing the shared genes and recomputing gave ρ 0.94–0.98 in all four cohorts, so the association is not an overlap artefact.

**33. You added the quality-control gate after seeing an implausible result. Isn't that cherry-picking?**
It is a declared deviation, and I treat it as one. The gate uses canonical co-expression pairs (collagens, proliferation, T-cell receptor, intestinal epithelium) that are unrelated to the hypotheses, and it would exclude a cohort whichever direction its results went. GSE103479 failed it badly (median r 0.07 against 0.78–0.94 elsewhere), probes for the same gene were uncorrelated, and adjacent probe IDs carried near-identical values. That is evidence of a defective deposit, not an inconvenient result.

**34. Why rank-biserial effect sizes for the replication?**
It is a scale-free, non-parametric effect size for two groups, so it can be compared across microarray and RNA-seq cohorts without assuming normality or a shared scale.

**35. Why were the TCGA BRAF effects judged "not established"?**
31 of the 44 BRAF-mutant tumours with mutation data were hypermutated (more than 300 mutated genes, against a cohort median of 113). After excluding hypermutated tumours, or adjusting for log mutation burden, the BRAF effects on the Wnt and p53 scores largely disappeared. The TP53 effect on the p53 score survived the same adjustment (−0.37 after exclusion).

**36. Why re-implement ESTIMATE in Python?**
So the tool runs without an R dependency. The re-implementation reads the official gene lists at run time (they are not bundled, because of their licence) and matches the package's own example output to 1.4 × 10⁻¹⁴.

**37. How do you test a tool when you do not know the true answer for real data?**
With synthetic data that has planted structure: a stromal programme, a mutation-driven programme, a programme confounded by an MSI-like factor, and a random set. The tests check that the audit labels each one correctly, that the quality gate fails when gene labels are shuffled, and that the survival screen shrinks the optimism of noise scores. The built-in demo recovers all four planted truths.

**38. What is new about `pathway-audit`?**
Not the individual methods, which are established. What is new is putting them in one tested workflow with explicit verdict rules and a data-quality gate in front, so a pathway score is checked for composition, a few dominant genes, genotype confounding and a broken dataset before anyone interprets it.

**39. What did the single-cell analysis add, and how did you keep it honest?**
It measured what the bulk work could only estimate: which cells express the EMT genes. I registered the plan and committed it before downloading anything, then committed the cell-type-to-compartment mapping after reading only the annotation files. For each patient I compared the EMT set's share of transcripts in each compartment with 1,000 random gene sets matched on expression level. Stromal cells carried 5.1 times the expected share and epithelial cells 0.45 times, in the discovery cohort (18 patients, p = 7.6 × 10⁻⁶) and in the same direction in two replication cohorts. A pseudo-bulk EMT score tracked the counted stromal fraction (ρ 0.97). I also re-read the raw matrix for two patients with separate code to check the streamed sums exactly.

**40. Does S3 show that tumour-cell EMT does not exist in colorectal cancer?**
No. It shows that, within malignant epithelial cells, the EMT genes were not more coherently co-expressed than expression-matched random genes in any of 61 patients (4 of 61 after an ambient-RNA control). That argues against a strong, uniform tumour-cell EMT programme in these samples. Rare or partial EMT states, transcripts lost to 3′ capture, or changes visible only in protein or spatial data would not be detected, so the verdict is "not detected at this resolution", which is exactly how the plan defined it in advance.

**41. How did you process a 370,000-cell atlas on a laptop?**
The matrix has about 764 million non-zero values, roughly 9 GB in memory. I streamed it in blocks of 20,000 cells straight from the HDF5 file and kept only per-patient, per-compartment sums plus a seeded subsample of up to 1,000 malignant cells per patient. The cap was declared in the plan's deviation log before any result was computed.

---

## Part 3. Numbers to know

| Quantity | Value |
|---|---|
| Cohort | 452 patients, 136 recurrences, 316 censored; 46 deaths censored |
| Probes → genes → panel | 54,675 → 43,163 → 21,815 → 607 |
| Clinical model | stage III vs II HR 1.97 (1.30–2.98) |
| C-index, apparent / corrected | M1 0.6242 / 0.6016; M2 0.6382 / 0.6059 |
| ΔC | apparent +0.0140; corrected +0.004316 |
| Intervals | five runs, about −0.009 to +0.05, all include zero |
| Calibration slopes | 0.833 (M1), 0.862 (M2) |
| Pathway HRs, primary (95 % CI) | EMT 1.070 (0.951–1.202), E2F 0.906 (0.807–1.016), Wnt 1.016 (0.906–1.138), p53 0.995 (0.882–1.123); all intervals include 1 |
| MMR ε² | EMT 0.0394, E2F 0.0251, Wnt 0.0187, p53 0.1716 (n = 410) |
| Bootstrap | 5 runs × 750 × 100; about 765,000 model fits; verifier 60/60 |
| Follow-up: EMT vs fibroblast markers, gene-disjoint ρ | 0.94 (GSE39582), 0.95 (GSE41258), 0.94 (GSE143985), 0.98 (TCGA-COAD) |
| Random-gene-set baseline, upper 97.5 % | 0.05, −0.15, 0.26, 0.76 (TCGA median 0.64) |
| TP53 → p53 score (rank-biserial) | GSE41258 −0.22; TCGA −0.42; TCGA excluding hypermutated −0.37; GSE143985 −0.06 (n.s.) |
| Competing risks (X2) | recurrence incidence at 60 units 30.7 % (Aalen–Johansen) vs 31.3 % naive |
| External audit (X8) | 554 hits, 14 assessed, 0 eligible; largest 124 recurrences vs ≈337 needed |
| Tool | pathway-audit 0.2.0; 13/13 tests; CI passes on Python 3.10 and 3.12; demo recovers 4/4 planted truths; 5 bulk cohorts (1 failed QC) + 1 single-cell cohort |
| Single cell (S1–S5) | EMT-set enrichment stromal 5.08× vs epithelial 0.45× (discovery, 18 patients); tumour-cell programme 0/61 patients; pseudo-bulk ρ 0.97 (random ≤ 0.69) |
