# Phase 10: machine-learning upgrade, methods and model card

**Status:** post-hoc and exploratory. **Phase 9 is frozen and unchanged:** the 27 hash-frozen files and 52 archive hashes are checked before and after every run (`phase10.common.phase9_integrity`, `tests/phase10`).

**Intended use:** a research and teaching demonstration of leakage-controlled ML on transcriptomic data. **It is not a clinical tool.** No model here supports a treatment decision, an absolute-risk estimate or a diagnostic claim.

## 1. Why Phase 10 exists

Phase 9 asked one pre-specified question: do four Hallmark pathway scores add recurrence discrimination beyond clinical variables? They did not (optimism-corrected ΔC +0.0043; all intervals include zero). Phase 10 adds two things.

- **Track A, survival ML.** It runs the transcriptome-wide elastic-net Cox model (M3), whose role was registered in Phase 6B (A-15) but which was never specified at execution level or run. The comparators run under one common, strictly nested cross-validation.
- **Track B, classification.** Predicting MMR/MSI status from expression is a **positive-control task** with strong, known biology. It shows whether the same machinery finds signal where signal exists. It also adds the classification metrics, calibration and independent-cohort evaluation that the recurrence endpoint cannot support.

**Declared deviations from the Phase 6B/9 protocols:**
- Repeated cross-validation replaces the A-17 nested bootstrap, for cost.
- The B-02 variance-filter rule is defined here for the first time.
- The M3 λ grid was extended once (see §4).

## 2. Cohorts and labels

**Track A**

| | |
|---|---|
| Patients | **452** (S3 complete cases; frozen loader `src/phase9_corrected/data.py`) |
| Events | **136** recurrences, 316 censored |
| Endpoint | `rfs_event_full` / `rfs_time_full` |
| Estimand | Cause-specific recurrence hazard |
| Caveat (time unit) | `rfs.delay` has **no documented time unit** |
| Caveat (competing risk) | **46 deaths without recurrence are censored** at death |

**Track B, training and cross-validation (GSE39582)**

| | |
|---|---|
| Patients | **409** = the 452 cohort with non-missing `mmr.status` |
| Classes | **59 dMMR / 350 pMMR** (14.4%) |
| Missing labels | `combat_mod` is never used to fill the 43 missing labels |

**Track B, independent evaluation (TCGA-COAD)**

| | |
|---|---|
| Samples | Primary tumours (`-01`): 286, one per patient |
| Primary population | **Stage II–III with complete age, sex and location: 182** |
| Sensitivity population | All stages, 274 (not run in this session) |
| Label | `MSI_updated_Oct62011` |
| Positive class | MSI-H |
| Negative class | MSS + MSI-L |
| Not used | The `microsatellite_instability` field |

**The two labels are different but related constructs.**

- **GSE39582 `mmr.status`** is **PCR-based MSI** (Marisa et al. 2013, *PLoS Med* 10:e1001453, Text S1).
  - Method: "a panel of five different microsatellite loci from the Bethesda reference panel".
  - Thresholds: MSI-H if ≥2 of the 5 markers are unstable; MSI-H is recoded "dMMR", and MSI-L/MSS are recoded "pMMR".
  - **No immunohistochemistry is described**, so "dMMR" here means MSI-H by PCR, not protein loss by IHC.
- **TCGA MSI** (TCGA 2012, *Nature* 487:330, Supplementary Information) is also PCR-based, but used two panels at two sites:
  - a 7-marker MDACC panel (BAT25, BAT26, BAT40, TGFBR2, D2S123, D5S346, D17S250);
  - the Promega MSI Analysis System v1.2 (5 mononucleotide markers).
  - MDACC thresholds: MSI-H at 5–7 of 7 markers (>70%), MSI-L at 1–2.
  - **Not verified:** that the Xena field `MSI_updated_Oct62011` comes from this assay. Its date and three-level coding fit the 2012 data freeze, but Xena documents no per-field provenance.
- **The panels, thresholds and laboratories differ**, so the labels are related, not identical. The TCGA analysis is therefore called **independent-cohort, cross-platform evaluation of MMR/MSI-associated expression classification**, and never validation of a dMMR classifier.

## 3. The ComBat limitation (Track B)

The deposited GSE39582 matrix was ComBat-corrected **with MMR status as a protected covariate**. Marisa et al. Text S1: "Centre and RNA extraction method as batch effects and with tumoral and MMR status as features of interest".

- The repository's X4 analysis found that `combat_mod` lines up with MMR status: modcan1 ↔ dMMR in 75/75 samples, modcan2 ↔ pMMR in 444/444.
- Protecting an unbalanced covariate can **inflate group differences** (Nygaard, Rødland & Hovig, *Biostatistics* 2016;17:29–39).
- CEL reprocessing is excluded by project decision B-03.

**Consequence:** every GSE39582 MMR cross-validation metric is **exploratory and likely biased upward**. The TCGA expression data were not corrected with MSI status, so the TCGA evaluation is free of this particular bias. It has its own limitations, listed in §6.

## 4. Validation protocol

**Common to both tracks**

- **Folds:** patient-level, stratified (event status for Track A; label for Track B). 5 outer folds × **10 repeats**; the benchmark showed 10 repeats were affordable.
- **Saved and hashed:** the fold tables are saved (`reports/phase10/*/folds.csv`) and their SHA-256 is recorded. Every model in a track reads the same table.
- **Inner tuning:** 5-fold stratified, seeded from `(seed_base, track, repeat, fold)`, so it is identical for every model.
- **Strict nesting.** Every data-dependent step is refit inside each inner-training fold during tuning, and on the outer-training set for the final fit:
  - spline knots;
  - pathway-score gene moments and score standardisation;
  - the top-2,000 variance filter;
  - gene and clinical scaling.
  
  The outer test fold never enters any fit. A spy transformer in the tests proves this.
- **Seeds:** `numpy.random.SeedSequence([2026111000, …])`, deliberately distinct from Phase 9's base.
- **Failed fits:** recorded per fold (`per_fold.csv`, `status`). Paired comparisons use only folds where both models succeeded.

**Track A models**

| Model | Definition |
|---|---|
| M1 | Frozen `fit_cox` (lifelines, Efron ties) on the 6 clinical df: age restricted cubic spline, sex, stage, location, adjuvant chemotherapy |
| M1-Breslow | Same design, scikit-survival `CoxPHSurvivalAnalysis(ties="breslow")`. A technical control for software and tie method |
| M2-strict | Same estimator as Phase 9 M2 (unpenalised clinical block, ridge on the four standardised scores, frozen 12-value λ grid, mean inner Harrell C, Phase 9 tie and skip rules). The design is rebuilt from the frozen functions inside every inner fold. `tests/phase10` proves that at a fixed λ its coefficients equal frozen `pipeline_fit` to within 1e-12 |
| M3 elastic-net Cox | Clinical 6 df **unpenalised** (`penalty_factor` 0) + top-2,000 variance genes (`penalty_factor` 1), all z-scaled on the training rows |

M3 tuning details:
- **Software:** scikit-survival `CoxnetSurvivalAnalysis`.
- **Grid:** `l1_ratio` ∈ {0.1, 0.5, 0.9} × a fixed absolute α grid.
- **Ties in tuning:** go to the larger α (sparser), then the larger `l1_ratio`.
- **Penalty rescaling:** scikit-survival rescales `penalty_factor` to sum to the number of features, so the effective gene penalty is α × (n_features / 2,000) ≈ 1.003α.
- **The one grid change:** the grid started at geomspace(1.0, 0.001, 30). The one-fold benchmark's inner-CV optimum sat at its top edge, so the pre-declared one-time contingency extended it to geomspace(10, 0.001, 40). Only training data informed this.

**Track B models**

All four take expression on the shared GSE39582/TCGA gene space.

| Model | Specification |
|---|---|
| Elastic-net logistic | `saga`; `l1_ratio` {0.1, 0.5, 0.9} × C logspace(−3, 1, 10); max_iter 5,000 |
| Random forest | 500 trees; `max_features` {sqrt, 0.1} × `min_samples_leaf` {1, 5} |
| MLH1-only | Unpenalised logistic regression |
| Clinical-only | Unpenalised logistic regression on age, sex, stage and location. **Chemotherapy is excluded**, because MMR status influences adjuvant-treatment decisions |

Shared details:
- **Preprocessing:** the gene models use VarianceFilter(top 2,000 on training rows), then StandardScaler. The forest's scaler has no effect on its splits; it is there so the TCGA adaptation applies equally to every gene model.
- **Class imbalance:** no reweighting. Folds are stratified, the primary metrics are threshold-free (AUROC, AUPRC), and predicted probabilities keep their calibration.

**Metrics**

- **Track A:** Harrell's C with lifelines 0.30.3 rules, checked against lifelines in the tests.
  - A pair is comparable when the earlier time is an event. An event tied with a censoring counts as comparable.
  - Two events at the same time are not comparable. Tied risk scores count 0.5.
  - The C is **within-fold pooled**: Σ concordant / Σ comparable, with pairs formed only inside each test fold. Every model is therefore scored on identical pairs, and ΔC is paired.
- **Track B:**
  - **AUROC:** within-fold pooled positive–negative pairs; principal.
  - **AUPRC:** principal; compare with the prevalence baseline of 0.144.
  - **Supporting:** Brier score, scaled Brier, and calibration intercept (logit offset) and slope.
  - AUPRC, Brier and calibration are computed **separately within each repeat**. Within a repeat each patient has exactly one out-of-fold prediction, so predictions are never stacked across repeats.
  - Logits are clipped at 1e-6 for calibration only. Random-forest probabilities can be exactly 0 or 1.

**Uncertainty**

- **Reported in CV:** the mean over repeats, the SD across repeats and the Monte Carlo SE (SD/√R).
  - These describe **partition variability given these data**. The repeats are independent given the data; the folds within a repeat are not.
  - **They are not population confidence intervals.**
  - No corrected *t* interval is reported.
- **Reported in TCGA:** frozen models evaluated on independent patients, so a **stratified patient bootstrap** gives legitimate 95% CIs for TCGA sampling variability. These do not cover training variability.
- **Not comparable:** Phase 9's optimism-corrected values appear only as a separate, labelled reference. They are never subtracted from Phase 10 values.

**Null diagnostic (Track A, N1)**

- **Construction:** (time, event) **pairs** are permuted together within each outer-training set. The test fold is untouched. The full M3 tuning then runs.
- **Expected:** a held-out C near 0.5.
- **Alarm, pre-registered:** mean > 0.52 *and* > 0.5 + 2·MC SE. It is a pipeline-defect diagnostic, not a significance test.

## 5. TCGA freeze and sealed labels

1. **Label-free preparation** (`python -m phase10.tcga prepare`):
   - primary tumours only;
   - clinical columns read through an explicit **allow-list** that excludes every MSI field;
   - stage mapping: Stage II*/III*;
   - location mapping: cecum, ascending colon, hepatic flexure and transverse colon are proximal; the rest are distal;
   - population files written.
2. **Freeze** (`python -m phase10.tcga freeze`):
   - final GSE39582 models are fit on all 409 patients, with inner-CV tuning;
   - `reports/phase10/tcga/freeze_manifest.json` records the models with their SHA-256, the feature lists, hyperparameters, preprocessing, the gene-space hash, cohort hashes, data and code hashes, and the environment.
3. **Sealed labels:**
   - `load_sealed_labels` reads `MSI_updated_Oct62011` only if three conditions hold: an approval file carries the manifest's SHA-256; every frozen file still matches its hash; the manifest is **committed unchanged in git**. Tests prove refusal otherwise.
   - Evaluation runs **once**; a second run is refused. No refitting, retuning, feature change or model selection happens after unblinding.
4. **Normalisation is unsupervised target-cohort adaptation.**
   - Each frozen model's genes are z-scored with the TCGA evaluation population's own mean and SD. Gene selection and coefficients stay frozen, and no labels are used.
   - **This is not a frozen single-patient prospective pipeline:**
     - it needs a batch of target samples;
     - gene means depend on the target case mix, so calibration depends on prevalence;
     - it aligns only each gene's mean and variance.
   - Clinical features keep the GSE-fitted scaler.

**Prior exposure, disclosed:** the marginal counts of `MSI_updated_Oct62011` across all TCGA samples (MSS 121, MSI-L 34, MSI-H 33) were seen during Stage 0 planning. They were never linked to expression.

## 6. Cross-platform and data limitations

- **Platforms:**
  - GSE39582 is Affymetrix HG-U133 Plus 2.0, median-of-probes per gene symbol, log₂ RMA with ComBat (as deposited).
  - TCGA is Illumina HiSeq RNASeqV2, `log2(RSEM normalized count + 1)` (UNC; Xena metadata).
  - The dynamic ranges differ: microarrays have a floor and saturation; RNA-seq has zeros. Of the shared genes, 87 are all-zero and 1,477 are more than half zero in the TCGA primary population.
- **Gene matching** is by symbol only: 16,322 shared; 5,493 GSE-only and 4,208 TCGA-only symbols are lost, including any gene renamed between annotations.
- **TCGA batch effects** (tissue source sites, plates) are untreated.
- **TCGA positives:** only 182 patients in the primary population, so CIs will be wide.
- **Track A:**
  - a single cohort and a weak signal;
  - the batch is confounded with the CIT partition;
  - the endpoint caveats in §2;
  - the ComBat-protected MMR axis is also present in the expression matrix Track A uses.
- **No external validation of the recurrence model** is possible (Phase 5).

## 7. Interpreting features (descriptive only)

**What is reported:**
- **Elastic-net models:** selection frequency over the outer fits, sign consistency and median coefficient.
- **Random forest:** how often each gene is in the top 20 by impurity importance. Impurity importance favours high-variance features and splits credit among correlated genes.

**What it doesn't mean:** outer-training sets overlap by about 60%, so frequency is **not** a replication probability. No selected gene is a biomarker, a driver or a validated finding. Coefficients do not imply mechanism.

**MLH1 context, verified:**
- Sporadic MSI colorectal cancers commonly silence MLH1 through promoter hypermethylation (Herman et al. 1998, *PNAS* 95:6870).
- EPM2AIP1 shares the bidirectional MLH1 promoter (Perera et al. 2011, *Oncol Rep*), but protein co-loss in colorectal cancer is imperfect (Challa et al. 2024, *Hum Pathol*: 64% sensitivity).
- MMR proteins act as heterodimers, MLH1/PMS2 and MSH2/MSH6 (Dedeurwaerdere et al. 2022, *Mod Pathol*). MSH2/MSH6-deficient tumours, often Lynch syndrome, would not be expected to lose MLH1. That is an inference from the heterodimer model, not a directly quoted finding.
- So the **MLH1-only model can capture only the MLH1-silenced subset**. Its misses can't be read as biology without the per-protein status, which these data lack.

## 8. Reproduce

Environment: the same pins as the Phase 9 `.venv` plus scikit-survival 0.28.0 (`requirements-phase10.txt`), in a venv outside the repo.

```bash
python3.13 -m venv ~/.venvs/coad-phase10 && ~/.venvs/coad-phase10/bin/pip install -r requirements-phase10.txt
P=~/.venvs/coad-phase10/bin/python
$P -m pytest tests/phase10 -q -o addopts=""           # safeguards
$P -m phase10.mmr --repeats 10                        # Track B cross-validation (~25 min)
$P -m phase10.recurrence --repeats 10 --null 5        # Track A + N1 diagnostic
$P -m phase10.tcga prepare && $P -m phase10.tcga freeze
# after the manifest is committed and approved:  $P -m phase10.tcga evaluate   (runs once)
```
