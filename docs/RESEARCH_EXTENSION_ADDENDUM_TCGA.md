# Registered addendum: E3 replication in TCGA-COAD

**Registered:** 8 November 2026, after the TCGA files were downloaded and their column structure inspected, and before any score or association was computed. The SHA-256 of this file is stored in `reports/research_extensions/registration_tcga.json`.

**Approval:** project owner, 8 November 2026 ("take full authority").

## Data
UCSC Xena TCGA hub (`data/external/tcga_coad_xena/`, provenance in `PROVENANCE.json`):
- `TCGA.COAD.sampleMap/HiSeqV2`: log2(RSEM + 1), gene symbols.
- `COAD_clinicalMatrix`: `sample_type`, `MSI_updated_Oct62011`.
- `mc3_gene_level/COAD_mc3_gene_level`: non-silent mutation calls (1/0) for profiled samples.

## Population
- Primary tumours (`-01` samples, `sample_type = Primary Tumor`).
- All stages, because TCGA-COAD is not stage II–III restricted. This differs from GSE39582 and is stated wherever results appear.
- Mutation status is defined only for samples present in the MC3 file. MSI-H vs MSS (MSI-L excluded).

## Hypotheses and rules
Identical to `docs/RESEARCH_EXTENSION_REGISTER.md` §E3:
- R1–R6 with the same directions, eligibility (≥ 10 per group; n ≥ 30 for R6), bootstrap CIs and decision rule.
- The frozen gene sets and coverage rule apply.
- Scores are computed within the cohort; MCP-counter fibroblast markers are used in gene-symbol mode.
- The positive-control QC (`e3b`) is applied as a gate.

TCGA is counted as an additional eligible cohort in the replication verdicts. Verdicts are reported **both** for the original three-cohort set and with TCGA added. TCGA is never pooled with any other cohort.
