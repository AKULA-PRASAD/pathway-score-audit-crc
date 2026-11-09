# pathway-audit

**Check what a pathway score actually measures in bulk tumour expression data before you interpret it.**

Gene-set ("pathway") scores, such as a Hallmark EMT score, are routinely read as biology: "EMT is high in these tumours". In bulk tissue a score can instead be driven by:
- tissue composition;
- a handful of genes;
- correlated molecular subtypes;
- a broken dataset.

`pathway-audit` runs the checks that separate these explanations and writes a plain-language report.

## What it reports

| Check | Question | Verdict examples |
|---|---|---|
| Data-quality gate | Do genes that co-express in virtually every bulk tumour (collagens, proliferation, T-cell receptor, intestinal epithelium) co-express here? | `PASS` / `FAIL` (gene-level results unsafe) |
| Coherence and gene drivers | Is the score a coordinated programme or an average of unrelated genes? Which genes drive it? | coherent / incoherent; driven by few genes |
| Gene modules | Does the set split into distinct co-expression modules, and are they stable under resampling? | stable / pattern / weak / dissolved (clusterwise Jaccard) |
| Composition | How much of the score is stromal, immune or fibroblast content (ESTIMATE, MCP-counter), and does that hold with shared genes removed? | composition-dominated / -associated / largely independent |
| Genotype effects | Is a mutation–score association robust to adjustment for composition and covariates (e.g. MMR status)? | robust / explained by adjustment covariates / not associated |
| Score relationships | Do two scores relate beyond what conditioning on composition produces for random gene sets? | exceeds random-set null: yes / no |
| Optional survival screen | Does adding the scores to a clinical Cox model improve optimism-corrected C? | corrected ΔC (screening estimate) |
| Single-cell source check (`sc-sources`) | Given single-cell counts with cell-type labels, which cell types supply each set's transcripts, compared with expression-matched random genes? | "mainly supplied by stromal (5.1x); depleted in epithelial (0.45x)" |

## Try it in 30 seconds

```bash
pip install -e .
pathway-audit demo --out demo/      # synthetic data with planted truth; writes demo/report.html and demo/REPORT.md
```

The demo plants four programmes: a stromal programme, a mutation-driven programme, a programme confounded by an MSI-like factor, and a random gene set. The audit recovers each one: composition-dominated; robust to adjustment; explained by adjustment covariates; incoherent. An example output is in `examples/demo_output/report.html`.

## Why it exists: what it found in five colorectal-cancer cohorts

Run on GSE39582, GSE41258, GSE143985, GSE103479 and TCGA-COAD (`examples/case_study_crc.py`; reports in `reports/pathway_audit_case_study/`):
- **The Hallmark "EMT" score is composition-dominated in every usable cohort.** Spearman ρ with the fibroblast estimate is 0.94–0.98, and still 0.94–0.98 with the six shared genes removed. In bulk colorectal tissue it reads stromal content, not tumour-cell EMT.
- **One public dataset failed the quality gate.** GSE103479 has a median control-pair r of 0.07 (the others are 0.78–0.94); its deposited values appear misaligned to their probe IDs. Without the gate, this would have been reported as "the EMT–stroma relationship does not replicate".
- **Genotype effects separate.** In GSE39582, the higher p53 score in BRAF-mutant tumours is *explained* by MMR status and composition, whereas the lower Wnt score is *robust* to them.
- **Hypermutation confounding is caught in TCGA-COAD.** Using log mutation burden as a covariate, both BRAF effects are classified *explained by adjustment covariates* (Wnt −0.70 → −0.11 SD). The TP53 → p53-score link stays *robust*.
- **Baselines matter.** In TCGA RNA-seq even random gene sets correlate about 0.64 with the fibroblast estimate. Compare scores with a cohort-specific random-set baseline, not with zero.
- **Survival screen.** On GSE39582 it gives corrected ΔC +0.0048 for the four scores. The full nested-bootstrap analysis gave +0.0043 with an interval spanning zero.
- **Single-cell confirmation.** On GSE132465 tumour cells (18 patients), `sc-sources` finds the EMT set's transcripts mainly supplied by stromal cells (5.1x the expression-matched share) and depleted in epithelial cells (0.45x), reproducing the registered single-cell analysis (5.08x and 0.45x; `examples/case_study_single_cell.py`).

## Install and run

```bash
pip install -e ".[survival,test]"        # Python >= 3.10
pathway-audit fetch-references --cache refs/     # downloads ESTIMATE + MCP-counter references, verifies ESTIMATE re-implementation
pathway-audit run \
  --expr expression.tsv  --gene-sets sets.gmt  --out audit/ \
  --annotations samples.csv \
  --binary "TP53 mutant=tp53:M:WT" "BRAF mutant=braf:M:WT" \
  --adjust mmr=dMMR  --references refs/ \
  --survival time,event,age,stage3          # optional
  # --composition my_cibersortx.tsv          # use your own composition estimates instead of ESTIMATE/MCP-counter

# single-cell source check on 10x output plus your cell annotations (one row per cell, in matrix order)
pathway-audit sc-sources --counts matrix.mtx.gz --genes features.tsv.gz --cells cells.csv \
  --gene-sets sets.gmt --cell-type-column cell_type --group-column patient --out sc/
```

**Inputs:**
- `expression.tsv`: genes × samples (log scale; linear MAS5/PLIER values are detected and log2-transformed).
- Gene sets: GMT or JSON.
- Annotations: one row per sample.

**Outputs:**
- `report.html` (self-contained, with charts)
- `REPORT.md`
- `audit.json`
- `scores.csv`
- `composition_estimates.csv`

Python API: `from pathway_audit import audit; res = audit(expr, gene_sets, annotations, ...)`.

## Tests

```bash
pytest                                    # synthetic data with known structure
PATHWAY_AUDIT_REFS=refs/ pytest           # also checks the ESTIMATE re-implementation against the package example
```

13 tests. Continuous integration (GitHub Actions) runs lint, the tests and the demo on Python 3.10 and 3.12 at every push, and passes. The Docker image has not yet been built.

The synthetic tests check that:
- the QC gate fails when gene labels are shuffled;
- planted modules are recovered as stable;
- a confounded genotype effect is classified "explained by adjustment covariates" while a direct one is "robust";
- gene-disjoint and random-set-null tests behave as expected;
- the survival screen shrinks the optimism of noise scores;
- the built-in demo recovers all four planted truths and the HTML report is self-contained;
- the single-cell source check finds a planted fibroblast-expressed set and leaves a random set near 1, through the API and the CLI.

## Limits (stated in every report)

- Scores are averaged bulk mRNA. They do not measure pathway activity, protein function or cell-type-specific expression.
- Composition values are computational estimates from the same data, not measured fractions.
- Adjusted effects can be over-adjusted when composition lies on the biological path.
- The survival screen is single-level and gives no interval. Use a nested bootstrap for confirmatory work.
- The single-cell source check depends on the supplied cell annotations, and tissue dissociation under-represents some cell types (often fibroblasts).

## Third-party references

Third-party reference data are not bundled:
- ESTIMATE (Yoshihara et al., *Nat Commun* 2013; GPL-2);
- MCP-counter (Becht et al., *Genome Biol* 2016; GPL-3).

`fetch-references` downloads them from their official sources and records checksums.

**Licence:** MIT (see `LICENSE`). The ESTIMATE and MCP-counter reference files are downloaded separately and keep their own licences.

Method details: [`docs/METHOD.md`](docs/METHOD.md).
