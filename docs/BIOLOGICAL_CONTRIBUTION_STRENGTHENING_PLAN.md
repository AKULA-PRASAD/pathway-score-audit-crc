# Biological contribution strengthening plan

**Date:** 8 November 2026.
**Scope:** how to strengthen biological interpretation and discovery through analyses, not wording.
**Execution record:** `docs/RESEARCH_EXTENSION_REGISTER.md` (plan) and `docs/FINAL_RESEARCH_CONTRIBUTION.md` (results).

## 1. Starting point (from repository inspection, not summaries)

| Area | Status before this plan | Evidence |
|---|---|---|
| Biological interpretation | Moderate to strong | X1 (TP53 anchors the p53 score), X6 (EMT ≈ stromal index), Phase 8C coherence |
| Original biological discovery | Limited | All biology was association-level within one cohort; no independent replication |
| Gene-level structure of scores | **Not examined** | Phase 8C covered set coherence, sample influence and score-level representation only |
| Composition-conditioned relationships | **Not examined** | X7 adjusted outcome models only |
| Within-subtype composition relationships | **Not examined** | — |
| Independent biological replication | **Not attempted** | X8 assessed predictive validation only |

## 2. Classification of candidate analyses

| Class | Analyses |
|---|---|
| High value, feasible with existing data | E1 within-pathway gene modules (item–rest, leave-one-gene-out, modules, stability, representation check); E2 partial correlations, composition-adjusted genotype effects, within-subtype EMT–stroma |
| High value, needing additional public data | E3 outcome-free replication of genotype–score and EMT–stroma relationships in independent cohorts with mutation annotation (owner approved downloads) |
| Potentially useful but underpowered | Subtype- or MMR-stratified recurrence models (14–35 and 10 events); gene-level outcome models |
| Interesting but not defensible here | CMS calling (CMS training included GSE39582, so circular); tumour-cell EMT inference from bulk data; claims of pathway activity |
| Redundant or low value | Re-running X1–X9; further preprocessing comparisons (representation r ≥ 0.99 already shown); single-cell integration (excluded by owner) |

## 3. Selected extensions

| Field | E1 Gene architecture | E2 Conditioned structure | E3 Independent replication |
|---|---|---|---|
| Question | Do scores hide distinct gene programmes, and which genes drive them? | Which relationships survive conditioning on estimated composition, and does EMT track stroma within subtypes? | Do genotype–score and EMT–stroma relationships replicate in independent tumours? |
| Rationale | A mean of 41–198 genes can mix biologically different programmes | Composition is the dominant source of bulk variation (X6) | Single-cohort associations may be cohort-specific |
| Inputs | `expression_median_S0.parquet`; representative-probe gene matrix; frozen panel; X5 composition; covariates | Same, plus CIT subtype from the GEO SOFT file | GSE41258, GSE143985, GSE103479 series matrices; GPL96 / GPL570 / GPL23985 tables; MCP-counter genes |
| Method | Item–rest r; leave-one-gene-out; average-linkage modules (cut 0.7); bootstrap ARI; module–annotation tests (BH); representation check | Partial Spearman given stromal + immune; OLS genotype effects ± composition ± dMMR; within-subtype ρ; bootstrap CIs | Frozen sets and coverage rule; within-cohort scores; six pre-registered hypotheses; fixed decision rule; never pooled |
| Estimand | Module membership, stability, module–annotation effects | Partial ρ; adjusted genotype difference (SD); within-subtype ρ | Rank-biserial or ρ per cohort; replication verdict |
| Sample size | 452 | 452 (subtypes 43–122) | 186 / 91 / 156 tumours; genotype groups ≥ 10 |
| Validation | Bootstrap stability; representative-probe replication | Bootstrap CIs; two composition estimators | Independent cohorts; pre-registered directions |
| Main risks | Data-driven modules; multiplicity | Shared-array circularity; CIT circularity; over-adjustment | Platform differences; annotation methods; small groups; marker overlap |
| Acceptance | Median ARI ≥ 0.7; same signs under repIQR | ρ ≥ 0.7 in ≥ 5 of 6 subtypes; adjusted CI excludes 0 with same sign | Same direction, CI excluding 0, in ≥ 2 eligible cohorts |
| Expected contribution | Shows what each score is made of | Separates genotype from composition signals | First independent evidence for this project's biology |

## 4. What this plan will not do

- Re-fit, replace or reinterpret M1/M2 or ΔC.
- Use recurrence outcomes.
- Claim pathway activity, tumour-cell-specific expression, or tumour-cell EMT.
- Pool external cohorts.
- Edit the presentation or web page.
