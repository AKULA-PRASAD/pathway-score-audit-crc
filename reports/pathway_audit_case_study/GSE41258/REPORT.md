# Pathway-score audit: GSE41258

Samples: 183 · genes: 12548

## 1. Data-quality gate (positive-control co-expression)

**PASS**: median r = 0.79 over 8 control pairs (threshold 0.5).

## 2. Verdict per gene set

| Gene set | Genes present | Verdict |
|---|---|---|
| EMT | 192/198 | coherent (mean inter-gene r 0.24); composition-dominated (|rho| 0.95 with MCP_Fibroblasts; holds without shared genes) |
| E2F | 187/195 | coherent (mean inter-gene r 0.29); composition-associated (|rho| 0.49 with StromalScore) |
| Wnt | 37/41 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.32) |
| p53 | 194/194 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.21) |

### EMT

Top genes (item–rest r): SPARC 0.93, EFEMP2 0.92, COL6A2 0.92, COL5A2 0.92, LGALS1 0.92, FSTL1 0.91, COL6A3 0.91, FAP 0.90
Modules: M1 (117 genes, Jaccard 0.90, stable)
Composition ρ: StromalScore +0.91, ImmuneScore +0.52, MCP_Fibroblasts +0.95, MCP_Endothelial cells +0.66, MCP_Cytotoxic lymphocytes +0.21, MCP_T cells -0.21, MCP_Monocytic lineage +0.61
Gene-disjoint fibroblast test: ρ 0.95 (0.92–0.96); 6 shared genes removed

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 90/53 | -0.11 (-0.48 to +0.24) | -0.01 (-0.12 to +0.11) | not associated |
| MSI positive | 35/116 | +0.01 (-0.32 to +0.34) | -0.01 (-0.15 to +0.12) | not associated |

### E2F

Top genes (item–rest r): KIF2C 0.87, CHEK1 0.86, MAD2L1 0.85, E2F8 0.85, MCM2 0.85, PLK4 0.85, ASF1B 0.84, MCM6 0.84
Modules: M1 (130 genes, Jaccard 0.84, stable)
Composition ρ: StromalScore -0.49, ImmuneScore -0.25, MCP_Fibroblasts -0.48, MCP_Endothelial cells -0.43, MCP_Cytotoxic lymphocytes +0.13, MCP_T cells -0.07, MCP_Monocytic lineage -0.22

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 90/53 | -0.21 (-0.53 to +0.15) | -0.25 (-0.52 to +0.05) | not associated |
| MSI positive | 35/116 | +0.67 (+0.30 to +1.05) | +0.66 (+0.32 to +0.96) | robust to adjustment |

### Wnt

Top genes (item–rest r): DVL2 0.32, AXIN1 0.30, JAG2 0.27, HDAC5 0.25, HDAC11 0.23, PTCH1 0.23, MAML1 0.18, NCOR2 0.18
Modules: M1 (5 genes, Jaccard 0.72, pattern)
Composition ρ: StromalScore -0.07, ImmuneScore -0.32, MCP_Fibroblasts +0.02, MCP_Endothelial cells +0.07, MCP_Cytotoxic lymphocytes -0.29, MCP_T cells -0.07, MCP_Monocytic lineage -0.29

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 90/53 | +0.25 (-0.10 to +0.60) | +0.11 (-0.29 to +0.46) | not associated |
| MSI positive | 35/116 | -0.44 (-0.77 to -0.08) | -0.29 (-0.58 to +0.07) | inconclusive |

### p53

Top genes (item–rest r): SFN 0.56, FDXR 0.54, CDKN1A 0.54, LIF 0.53, EPHA2 0.48, DGKA 0.46, UPP1 0.45, RNF19B 0.44
Modules: M1 (12 genes, Jaccard 0.44, dissolved)
Composition ρ: StromalScore -0.04, ImmuneScore +0.13, MCP_Fibroblasts -0.05, MCP_Endothelial cells -0.21, MCP_Cytotoxic lymphocytes +0.10, MCP_T cells -0.00, MCP_Monocytic lineage +0.06

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 90/53 | -0.39 (-0.78 to -0.02) | -0.27 (-0.60 to +0.10) | inconclusive |
| MSI positive | 35/116 | +0.39 (-0.04 to +0.82) | +0.29 (-0.14 to +0.71) | not associated |

## 3. Relationships between scores (partial Spearman given stromal + immune estimates)

| Pair | Spearman | Partial | Random-set null 95 % | Exceeds null |
|---|---|---|---|---|
| EMT | E2F | -0.50 | -0.06 | -0.26 to +0.12 | no |
| EMT | Wnt | +0.03 | +0.05 | -0.26 to +0.12 | no |
| EMT | p53 | +0.04 | +0.32 | -0.26 to +0.12 | yes |
| E2F | Wnt | -0.23 | -0.25 | -0.26 to +0.38 | no |
| E2F | p53 | -0.16 | -0.26 | -0.26 to +0.38 | yes |
| Wnt | p53 | -0.05 | +0.05 | -0.06 to +0.31 | no |

## Interpretation limits

- Scores are averaged bulk mRNA; they do not measure pathway activity, protein function or cell-type-specific expression.
- Composition values are computational estimates (ESTIMATE, MCP-counter) from the same expression data, not measured fractions.
- Adjusted genotype effects may be over-adjusted if composition lies on the biological path.
