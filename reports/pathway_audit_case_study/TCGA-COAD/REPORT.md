# Pathway-score audit: TCGA-COAD

Samples: 286 · genes: 20530

## 1. Data-quality gate (positive-control co-expression)

**PASS**: median r = 0.94 over 8 control pairs (threshold 0.5).

## 2. Verdict per gene set

| Gene set | Genes present | Verdict |
|---|---|---|
| EMT | 195/198 | coherent (mean inter-gene r 0.38); composition-dominated (|rho| 0.98 with MCP_Fibroblasts; holds without shared genes) |
| E2F | 189/195 | coherent (mean inter-gene r 0.28); largely composition-independent (max |rho| 0.36) |
| Wnt | 41/41 | coherent (mean inter-gene r 0.08); composition-associated (|rho| 0.58 with MCP_Endothelial cells) |
| p53 | 190/194 | incoherent (mean inter-gene r < 0.05); composition-associated (|rho| 0.52 with ImmuneScore) |

### EMT

Top genes (item–rest r): EFEMP2 0.95, FSTL1 0.95, VIM 0.95, FBN1 0.95, COL3A1 0.95, COL6A2 0.95, PDGFRB 0.94, COL5A2 0.94
Modules: M1 (149 genes, Jaccard 0.94, stable)
Composition ρ: StromalScore +0.96, ImmuneScore +0.63, MCP_Fibroblasts +0.98, MCP_Endothelial cells +0.85, MCP_Cytotoxic lymphocytes +0.46, MCP_T cells +0.48, MCP_Monocytic lineage +0.76
Gene-disjoint fibroblast test: ρ 0.98 (0.97–0.98); 6 shared genes removed

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 162/105 | -0.02 (-0.27 to +0.24) | -0.05 (-0.13 to +0.02) | not associated |
| BRAF positive | 44/223 | +0.03 (-0.25 to +0.35) | +0.04 (-0.08 to +0.17) | not associated |

### E2F

Top genes (item–rest r): CDK1 0.83, CHEK1 0.83, ASF1B 0.82, HELLS 0.82, PLK4 0.81, RACGAP1 0.81, MAD2L1 0.81, GINS3 0.81
Modules: M1 (134 genes, Jaccard 0.83, stable); M2 (10 genes, Jaccard 0.35, dissolved)
Composition ρ: StromalScore -0.33, ImmuneScore -0.20, MCP_Fibroblasts -0.31, MCP_Endothelial cells -0.36, MCP_Cytotoxic lymphocytes -0.06, MCP_T cells -0.16, MCP_Monocytic lineage -0.19

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 162/105 | +0.13 (-0.13 to +0.38) | +0.26 (+0.02 to +0.50) | robust to adjustment |
| BRAF positive | 44/223 | +0.08 (-0.33 to +0.46) | -0.35 (-0.85 to +0.18) | not associated |

### Wnt

Top genes (item–rest r): NOTCH1 0.61, PTCH1 0.60, NOTCH4 0.55, NCOR2 0.52, LEF1 0.51, NKD1 0.51, WNT6 0.49, WNT5B 0.48
Modules: M1 (8 genes, Jaccard 0.69, pattern); M2 (7 genes, Jaccard 0.58, weak)
Composition ρ: StromalScore +0.47, ImmuneScore +0.14, MCP_Fibroblasts +0.48, MCP_Endothelial cells +0.58, MCP_Cytotoxic lymphocytes +0.02, MCP_T cells +0.22, MCP_Monocytic lineage +0.27

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 162/105 | +0.28 (+0.03 to +0.51) | +0.03 (-0.17 to +0.23) | explained by adjustment covariates |
| BRAF positive | 44/223 | -0.70 (-1.00 to -0.38) | -0.11 (-0.43 to +0.24) | explained by adjustment covariates |

### p53

Top genes (item–rest r): PLK3 0.59, CDKN1A 0.57, PHLDA3 0.54, RNF19B 0.52, TNFSF9 0.51, FAS 0.50, HMOX1 0.46, CTSD 0.46
Modules: M1 (25 genes, Jaccard 0.76, stable); M2 (10 genes, Jaccard 0.51, weak)
Composition ρ: StromalScore +0.36, ImmuneScore +0.52, MCP_Fibroblasts +0.36, MCP_Endothelial cells +0.32, MCP_Cytotoxic lymphocytes +0.45, MCP_T cells +0.50, MCP_Monocytic lineage +0.40

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 162/105 | -0.71 (-0.93 to -0.47) | -0.58 (-0.81 to -0.33) | robust to adjustment |
| BRAF positive | 44/223 | +0.40 (+0.10 to +0.70) | +0.01 (-0.38 to +0.40) | explained by adjustment covariates |

## 3. Relationships between scores (partial Spearman given stromal + immune estimates)

| Pair | Spearman | Partial | Random-set null 95 % | Exceeds null |
|---|---|---|---|---|
| EMT | E2F | -0.36 | -0.15 | +0.03 to +0.31 | yes |
| EMT | Wnt | +0.51 | +0.12 | +0.03 to +0.31 | no |
| EMT | p53 | +0.42 | +0.50 | +0.03 to +0.31 | yes |
| E2F | Wnt | -0.38 | -0.26 | -0.20 to +0.15 | yes |
| E2F | p53 | -0.31 | -0.26 | -0.20 to +0.15 | yes |
| Wnt | p53 | +0.10 | +0.07 | +0.24 to +0.55 | yes |

## Interpretation limits

- Scores are averaged bulk mRNA; they do not measure pathway activity, protein function or cell-type-specific expression.
- Composition values are computational estimates (ESTIMATE, MCP-counter) from the same expression data, not measured fractions.
- Adjusted genotype effects may be over-adjusted if composition lies on the biological path.
