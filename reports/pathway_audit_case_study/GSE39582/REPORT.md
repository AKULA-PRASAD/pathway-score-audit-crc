# Pathway-score audit: GSE39582

Samples: 452 · genes: 21815

## 1. Data-quality gate (positive-control co-expression)

**PASS**: median r = 0.84 over 8 control pairs (threshold 0.5).

## 2. Verdict per gene set

| Gene set | Genes present | Verdict |
|---|---|---|
| EMT | 198/198 | coherent (mean inter-gene r 0.24); composition-dominated (|rho| 0.95 with MCP_Fibroblasts; holds without shared genes) |
| E2F | 195/195 | coherent (mean inter-gene r 0.30); largely composition-independent (max |rho| 0.35) |
| Wnt | 41/41 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.19) |
| p53 | 194/194 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.39) |

### EMT

Top genes (item–rest r): SPARC 0.93, EFEMP2 0.92, COL6A3 0.91, COL5A2 0.91, LGALS1 0.90, THY1 0.89, PCOLCE 0.89, COL3A1 0.89
Modules: M1 (117 genes, Jaccard 0.92, stable)
Composition ρ: StromalScore +0.91, ImmuneScore +0.57, MCP_Fibroblasts +0.95, MCP_Endothelial cells +0.69, MCP_Cytotoxic lymphocytes +0.34, MCP_T cells +0.23, MCP_Monocytic lineage +0.64
Gene-disjoint fibroblast test: ρ 0.94 (0.93–0.96); 6 shared genes removed

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 mutant | 154/128 | -0.15 (-0.38 to +0.09) | -0.06 (-0.14 to +0.03) | not associated |
| BRAF mutant | 35/320 | +0.32 (+0.03 to +0.61) | +0.01 (-0.15 to +0.17) | explained by adjustment covariates |
| KRAS mutant | 149/239 | +0.06 (-0.17 to +0.28) | +0.14 (+0.06 to +0.21) | robust to adjustment |
| dMMR | 59/350 | +0.47 (+0.26 to +0.71) | +0.12 (+0.06 to +0.17) | robust to adjustment |

### E2F

Top genes (item–rest r): MELK 0.84, SPC25 0.84, BUB1B 0.84, KIF4A 0.84, SPAG5 0.83, KIF2C 0.83, CHEK1 0.82, MCM2 0.82
Modules: M1 (142 genes, Jaccard 0.86, stable); M2 (16 genes, Jaccard 0.48, dissolved)
Composition ρ: StromalScore -0.32, ImmuneScore -0.25, MCP_Fibroblasts -0.25, MCP_Endothelial cells -0.35, MCP_Cytotoxic lymphocytes -0.04, MCP_T cells -0.19, MCP_Monocytic lineage -0.13

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 mutant | 154/128 | +0.28 (+0.04 to +0.51) | +0.31 (+0.10 to +0.54) | robust to adjustment |
| BRAF mutant | 35/320 | +0.32 (+0.01 to +0.61) | +0.03 (-0.30 to +0.34) | explained by adjustment covariates |
| KRAS mutant | 149/239 | -0.17 (-0.37 to +0.04) | -0.15 (-0.34 to +0.04) | not associated |
| dMMR | 59/350 | +0.51 (+0.30 to +0.75) | +0.37 (+0.24 to +0.50) | robust to adjustment |

### Wnt

Top genes (item–rest r): NOTCH1 0.56, NKD1 0.49, AXIN2 0.40, CSNK1E 0.36, TCF7 0.35, PTCH1 0.35, WNT6 0.34, AXIN1 0.34
Modules: M1 (5 genes, Jaccard 0.84, stable)
Composition ρ: StromalScore -0.04, ImmuneScore -0.18, MCP_Fibroblasts +0.04, MCP_Endothelial cells +0.08, MCP_Cytotoxic lymphocytes -0.10, MCP_T cells -0.05, MCP_Monocytic lineage -0.19

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 mutant | 154/128 | +0.21 (-0.03 to +0.45) | +0.09 (-0.14 to +0.33) | not associated |
| BRAF mutant | 35/320 | -0.83 (-1.12 to -0.54) | -0.65 (-1.02 to -0.29) | robust to adjustment |
| KRAS mutant | 149/239 | +0.14 (-0.06 to +0.34) | +0.10 (-0.10 to +0.29) | not associated |
| dMMR | 59/350 | -0.41 (-0.64 to -0.16) | -0.11 (-0.23 to +0.01) | explained by adjustment covariates |

### p53

Top genes (item–rest r): RNF19B 0.62, UPP1 0.60, TNFSF9 0.58, EPHA2 0.57, CDKN1A 0.55, FDXR 0.55, DDB2 0.55, PLK2 0.52
Modules: M1 (14 genes, Jaccard 0.75, pattern); M2 (14 genes, Jaccard 0.60, weak); M3 (10 genes, Jaccard 0.60, weak)
Composition ρ: StromalScore +0.16, ImmuneScore +0.32, MCP_Fibroblasts +0.18, MCP_Endothelial cells +0.12, MCP_Cytotoxic lymphocytes +0.39, MCP_T cells +0.24, MCP_Monocytic lineage +0.32

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 mutant | 154/128 | -0.43 (-0.64 to -0.22) | -0.29 (-0.51 to -0.08) | robust to adjustment |
| BRAF mutant | 35/320 | +0.89 (+0.52 to +1.28) | +0.03 (-0.41 to +0.48) | explained by adjustment covariates |
| KRAS mutant | 149/239 | +0.16 (-0.05 to +0.37) | +0.27 (+0.10 to +0.44) | robust to adjustment |
| dMMR | 59/350 | +1.27 (+0.98 to +1.54) | +0.52 (+0.37 to +0.65) | robust to adjustment |

## 3. Relationships between scores (partial Spearman given stromal + immune estimates)

| Pair | Spearman | Partial | Random-set null 95 % | Exceeds null |
|---|---|---|---|---|
| EMT | E2F | -0.22 | +0.17 | +0.06 to +0.29 | no |
| EMT | Wnt | +0.08 | +0.23 | +0.06 to +0.29 | no |
| EMT | p53 | +0.30 | +0.55 | +0.06 to +0.29 | yes |
| E2F | Wnt | +0.06 | +0.04 | -0.17 to +0.25 | no |
| E2F | p53 | +0.17 | +0.25 | -0.17 to +0.25 | no |
| Wnt | p53 | +0.16 | +0.26 | +0.21 to +0.44 | no |

## 4. Optional survival screen

Apparent ΔC +0.0174; optimism +0.0126; **corrected ΔC +0.0048** (n 452, events 136, B 200; screening estimate: single-level bootstrap, no interval).

## Interpretation limits

- Scores are averaged bulk mRNA; they do not measure pathway activity, protein function or cell-type-specific expression.
- Composition values are computational estimates (ESTIMATE, MCP-counter) from the same expression data, not measured fractions.
- Adjusted genotype effects may be over-adjusted if composition lies on the biological path.
