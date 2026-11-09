# Pathway-score audit: GSE143985

Samples: 91 · genes: 21655

## 1. Data-quality gate (positive-control co-expression)

**PASS**: median r = 0.78 over 8 control pairs (threshold 0.5).

## 2. Verdict per gene set

| Gene set | Genes present | Verdict |
|---|---|---|
| EMT | 198/198 | coherent (mean inter-gene r 0.30); composition-dominated (|rho| 0.94 with MCP_Fibroblasts; holds without shared genes) |
| E2F | 195/195 | coherent (mean inter-gene r 0.24); composition-associated (|rho| 0.46 with StromalScore) |
| Wnt | 41/41 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.11) |
| p53 | 194/194 | incoherent (mean inter-gene r < 0.05); composition-associated (|rho| 0.44 with MCP_Endothelial cells) |

### EMT

Top genes (item–rest r): EFEMP2 0.95, COL5A2 0.94, BGN 0.94, THY1 0.94, LGALS1 0.93, INHBA 0.93, COL5A1 0.93, COL6A2 0.93
Modules: M1 (136 genes, Jaccard 0.90, stable)
Composition ρ: StromalScore +0.90, ImmuneScore +0.71, MCP_Fibroblasts +0.94, MCP_Endothelial cells +0.69, MCP_Cytotoxic lymphocytes +0.44, MCP_T cells +0.42, MCP_Monocytic lineage +0.74
Gene-disjoint fibroblast test: ρ 0.94 (0.90–0.96); 6 shared genes removed

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 53/38 | +0.16 (-0.28 to +0.59) | +0.05 (-0.10 to +0.18) | not associated |
| BRAF positive | 2/89 | — | — | not testable (groups 2/89 < 10) |
| MSI positive | 5/85 | — | — | not testable (groups 5/85 < 10) |

### E2F

Top genes (item–rest r): KIF18B 0.86, TOP2A 0.85, CDC20 0.84, BUB1B 0.82, KPNA2 0.82, KIF2C 0.82, MELK 0.82, CDK1 0.82
Modules: M1 (110 genes, Jaccard 0.80, stable); M2 (11 genes, Jaccard 0.32, dissolved)
Composition ρ: StromalScore -0.46, ImmuneScore -0.41, MCP_Fibroblasts -0.43, MCP_Endothelial cells -0.35, MCP_Cytotoxic lymphocytes -0.35, MCP_T cells -0.29, MCP_Monocytic lineage -0.35

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 53/38 | +0.07 (-0.32 to +0.49) | +0.08 (-0.33 to +0.47) | not associated |
| BRAF positive | 2/89 | — | — | not testable (groups 2/89 < 10) |
| MSI positive | 5/85 | — | — | not testable (groups 5/85 < 10) |

### Wnt

Top genes (item–rest r): NKD1 0.46, KAT2A 0.45, PTCH1 0.43, AXIN2 0.36, TCF7 0.36, NUMB 0.26, JAG2 0.23, CTNNB1 0.22
Composition ρ: StromalScore -0.06, ImmuneScore -0.11, MCP_Fibroblasts -0.02, MCP_Endothelial cells +0.07, MCP_Cytotoxic lymphocytes -0.00, MCP_T cells +0.09, MCP_Monocytic lineage -0.08

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 53/38 | +0.26 (-0.18 to +0.69) | +0.24 (-0.18 to +0.68) | not associated |
| BRAF positive | 2/89 | — | — | not testable (groups 2/89 < 10) |
| MSI positive | 5/85 | — | — | not testable (groups 5/85 < 10) |

### p53

Top genes (item–rest r): PHLDA3 0.68, CDKN1A 0.60, LIF 0.56, TNFSF9 0.56, PLK3 0.47, PPP1R15A 0.46, RGS16 0.44, PLK2 0.43
Modules: M1 (15 genes, Jaccard 0.45, dissolved); M2 (15 genes, Jaccard 0.41, dissolved)
Composition ρ: StromalScore +0.30, ImmuneScore +0.31, MCP_Fibroblasts +0.35, MCP_Endothelial cells +0.44, MCP_Cytotoxic lymphocytes +0.33, MCP_T cells +0.35, MCP_Monocytic lineage +0.29

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 53/38 | -0.05 (-0.51 to +0.37) | -0.05 (-0.44 to +0.29) | not associated |
| BRAF positive | 2/89 | — | — | not testable (groups 2/89 < 10) |
| MSI positive | 5/85 | — | — | not testable (groups 5/85 < 10) |

## 3. Relationships between scores (partial Spearman given stromal + immune estimates)

| Pair | Spearman | Partial | Random-set null 95 % | Exceeds null |
|---|---|---|---|---|
| EMT | E2F | -0.46 | -0.16 | +0.19 to +0.44 | yes |
| EMT | Wnt | +0.10 | +0.33 | +0.19 to +0.44 | no |
| EMT | p53 | +0.41 | +0.40 | +0.19 to +0.44 | no |
| E2F | Wnt | -0.01 | -0.04 | -0.24 to +0.14 | no |
| E2F | p53 | -0.34 | -0.24 | -0.24 to +0.14 | yes |
| Wnt | p53 | +0.30 | +0.34 | +0.28 to +0.52 | no |

## Interpretation limits

- Scores are averaged bulk mRNA; they do not measure pathway activity, protein function or cell-type-specific expression.
- Composition values are computational estimates (ESTIMATE, MCP-counter) from the same expression data, not measured fractions.
- Adjusted genotype effects may be over-adjusted if composition lies on the biological path.
