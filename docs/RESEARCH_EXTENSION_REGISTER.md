# Research extension register (E1–E3)

**Registered:** 8 November 2026, before any E1–E3 association was computed. The SHA-256 of this file at registration is stored in `reports/research_extensions/registration.json` and re-checked by the validation script.

**Approval:** project owner, 8 November 2026 ("proceed, approve E3 downloads").

**Relationship to earlier work:**
- The original registered study (M1/M2, corrected ΔC) remains the primary research question and is untouched.
- X1–X9 (`docs/EXPLORATORY_ANALYSIS_PLAN.md`) are not repeated. E1–E3 build on their outputs.
- **All of E1–E3 is outcome-free.** No recurrence data are used, so none of it can affect the frozen models or estimands.

**Extension research question (new, registered separately; not retrospectively part of the original study):**
> What do predefined Hallmark pathway scores measure in bulk colorectal cancer tissue? How much of their variation reflects tumour genotype, molecular subtype and estimated tissue composition, and is that structure reproducible in independent tumours?

## Common definitions

| Item | Definition |
|---|---|
| GSE39582 population | Frozen A-19 cohort, 452 patients (`load_context`) |
| Expression | `expression_median_S0.parquet` (primary). Representative-probe check: `data/processed/candidates/gene_matrix_representative_iqr.parquet` |
| Gene sets | Frozen `config/pathways/m2_panel_4set.json` `members_resolved` |
| Scores | Frozen formula: mean of gene z-scores (ddof = 1), standardised within the analysed cohort |
| Composition | X5 outputs (ESTIMATE, MCP-counter) for GSE39582. For E3, MCP-counter gene-symbol mode (platform-independent) |
| Software | Project `.venv`; seeds fixed at 20261108 |
| Outputs | `scripts/research_extensions/`, `reports/research_extensions/`, `results/research_extensions/`, `figures/research_extensions/`; downloads in `data/external/research_extensions/` with `PROVENANCE.json` |

## E1. Within-pathway gene architecture (GSE39582)

- **Question:** do the aggregate scores hide distinct gene-level programmes, and which genes drive each score?
- **Rationale:** a mean of 41–198 genes can mix biologically different programmes (for example canonical p53 targets and immune/interferon-responsive genes in the Hallmark p53 set). Bulk interpretation depends on which component varies.

**Methods (all fixed in advance):**
1. **Per-gene item–rest correlation:** Pearson correlation of each gene's z-score with the score recomputed without that gene.
2. **Leave-one-gene-out:** correlation of each leave-one-out score with the full score. Report the minimum, and the genes whose removal changes the score most.
3. **Modules:**
   - Within each set, average-linkage hierarchical clustering of genes with distance 1 − Pearson r (computed across the 452 patients), **cut at distance 0.7** (genes joined only if average r ≥ 0.3).
   - Clusters of **≥ 10 genes** are modules; the rest are "unassigned".
   - The Wnt set (41 genes) uses a minimum of 5 genes.
4. **Stability:**
   - 200 patient bootstrap resamples, reclustered with the same rule.
   - Adjusted Rand index (ARI) between the bootstrap and full-data partitions of module genes.
   - **A module structure is accepted as stable if the median ARI ≥ 0.7.**
5. **Module scores:** mean gene z-score of each accepted module.
6. **Annotation tests:**
   - Module scores vs TP53, KRAS, BRAF, CIMP, CIN and MMR: Mann–Whitney with rank-biserial effects. BH-FDR within E1.
   - Module scores vs composition (ESTIMATE stromal, immune, purity; MCP fibroblasts, cytotoxic lymphocytes): Spearman ρ, descriptive.
7. **Representation check:** recompute modules and module scores from the representative-probe matrix. Report the module-score correlation between representations and agreement in the sign of each annotation effect.

**Pre-declared hypothesis H-E1 (either outcome is reported):** the Hallmark p53 set contains at least two separable modules. The TP53-mutation association concentrates in a module containing canonical p53 transcriptional targets (any of CDKN1A, MDM2, TP53I3, SESN1, ZMAT3, FDXR, RRM2B, TNFRSF10B, BAX, DDB2, XPC). The immune-estimate association concentrates in a different module.
- **Supported** if the canonical-target module's TP53 effect (rank-biserial) is larger in magnitude than every other module's, and the module with the largest |ρ| with ImmuneScore is a different module.
- Both must hold in the representative-probe check.

**Risks:**
- Modules are data-driven; this is contained by the fixed rule, the stability threshold and the representation check.
- Multiplicity: contained by BH.
- No outcome is used.

## E2. Composition- and subtype-conditioned structure (GSE39582)

- **Question:** which score relationships survive conditioning on estimated tissue composition, and does the EMT score track estimated stroma within every molecular subtype?

**Methods:**
1. **Partial Spearman correlations** between the four scores, conditioning on ESTIMATE stromal and immune scores. Implemented as Pearson correlation of residuals of rank-transformed variables, with 2,000-resample bootstrap 95 % CIs. Purity is omitted because it is a function of stromal + immune.
2. **Genotype effects adjusted for composition:** a linear model of each score on mutation indicator + stromal + immune (+ dMMR for TP53, BRAF and CIMP). Report the genotype coefficient in score-SD units with bootstrap 95 % CI, alongside the unadjusted difference in means.
3. **Within-subtype EMT–fibroblast Spearman ρ** in each CIT subtype (n 43–122), with bootstrap 95 % CIs. The same is done for EMT–ESTIMATE stromal.

**Acceptance and decision rules (fixed):**
- "EMT tracks estimated stroma within subtypes" requires ρ(EMT, MCP fibroblasts) ≥ 0.7 in ≥ 5 of 6 subtypes.
- A genotype association is "robust to estimated composition" if its adjusted 95 % CI excludes zero with the same sign as the unadjusted effect.

**Risks:**
- Estimates and scores come from the same arrays.
- CIT subtypes were derived from these data (stated, not solved).
- Over-adjustment: composition may lie on the biological path.

## E3. Independent biological replication (outcome-free)

**Cohorts** (independent of GSE39582 and CIT, per the X8 audit):

| Cohort | Population used | Platform; deposited values | Annotations used |
|---|---|---|---|
| GSE41258 (MSKCC) | `tissue: Primary Tumor` only | GPL96 (U133A); PLIER, linear | p53 mutation status; microsatellite instability (MSI-high vs MSS) |
| GSE143985 (Japan) | all primary colon tumours | GPL570; MAS5, linear | tp53_mutation, kras_mutation, braf_mutation, msi_status |
| GSE103479 (Aberdeen, Barcelona, Florence, Dublin) | all tumours | GPL23985 (Almac Xcel); RMA + ComBat, log2 | p53, kras, braf (V600E) |

**Processing (fixed):**
- Linear-scale data are log2-transformed (log2 of max(x, 1)).
- Only probes mapping to exactly one gene symbol in the platform table are kept, then collapsed by the median per gene.
- Frozen gene sets and the frozen coverage rule apply: ≥ 80 % of source members present and ≥ 20 genes. A score failing the rule is not computed for that cohort.
- Scores and MCP-counter (gene-symbol mode) are computed within each cohort.
- Cohorts are **never pooled**.

**Pre-registered replication hypotheses** (directions taken from GSE39582 X1/X6):

| ID | Hypothesis | GSE39582 effect |
|---|---|---|
| R1 | p53 score lower in TP53-mutant tumours | rank-biserial −0.30 |
| R2 | E2F score higher in TP53-mutant tumours | +0.19 |
| R3 | Wnt score lower in BRAF-mutant tumours | −0.46 |
| R4 | p53 score higher in BRAF-mutant tumours | +0.48 |
| R5 | p53 score higher in MSI-high/dMMR tumours | ε² 0.167, dMMR higher |
| R6 | EMT score strongly correlated with MCP fibroblast estimate (ρ ≥ 0.7) | ρ 0.95 |

**Estimation and decision rules:**
- Effect: rank-biserial (R1–R5) or Spearman ρ (R6), each with a 2,000-resample bootstrap 95 % CI.
- A cohort is eligible for a hypothesis only if each group has ≥ 10 tumours (R1–R5), or n ≥ 30 (R6).
- **Replicated:** same direction with the CI excluding zero (R6: point estimate ≥ 0.7 and CI lower bound ≥ 0.5) in **≥ 2 eligible cohorts**, and in no eligible cohort a CI excluding zero in the opposite direction.
- **Partially replicated:** satisfied in exactly 1 eligible cohort, with no opposite result.
- **Not replicated:** otherwise.
- **Not testable:** fewer than 1 eligible cohort.

**Risks:**
- Platform differences: U133A lacks some genes, Almac is a different design, MAS5/PLIER are different summarisations. Coverage is reported.
- Mutation-calling methods differ between cohorts.
- Small groups.
- Literature agreement is **not** counted as replication.

## Multiplicity and reporting

- E1: BH within its annotation family.
- E2: estimation only.
- E3: six pre-registered hypotheses with a fixed decision rule, and no further testing.
- Every result is labelled "post-hoc extension, outcome-free".

## Deviations

Recorded in `docs/FINAL_RESEARCH_CONTRIBUTION.md` §Deviations.
