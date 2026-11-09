# Pathway-score audit: GSE103479

Samples: 156 · genes: 24457

## 1. Data-quality gate (positive-control co-expression)

**FAIL**: median r = 0.07 over 8 control pairs (threshold 0.5).

> Gene-level results below are **not interpretable**: control genes that co-express in virtually every bulk tumour dataset do not co-express here, which indicates a processing or probe-labelling defect.

## 2. Verdict per gene set

| Gene set | Genes present | Verdict |
|---|---|---|
| EMT | 197/198 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.14) |
| E2F | 194/195 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.31) |
| Wnt | 39/41 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.19) |
| p53 | 192/194 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.30) |

### EMT

Top genes (item–rest r): SAT1 0.62, FMOD 0.53, ID2 0.52, CAPG 0.52, IL32 0.52, FBLN1 0.52, SNTB1 0.51, MYLK 0.49
Modules: M1 (15 genes, Jaccard 0.46, dissolved); M2 (12 genes, Jaccard 0.64, pattern); M3 (11 genes, Jaccard 0.60, weak)
Composition ρ: StromalScore -0.12, ImmuneScore +0.14, MCP_Fibroblasts +0.09, MCP_Endothelial cells -0.02, MCP_Cytotoxic lymphocytes -0.14, MCP_T cells +0.14, MCP_Monocytic lineage +0.12

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 78/59 | +0.36 (+0.00 to +0.68) | +0.33 (+0.01 to +0.66) | robust to adjustment |
| BRAF positive | 15/122 | -0.69 (-1.33 to -0.02) | -0.67 (-1.34 to +0.03) | inconclusive |

### E2F

Top genes (item–rest r): SMC6 0.38, RAN 0.36, KIF18B 0.36, RFC1 0.32, ASF1B 0.31, BRMS1L 0.30, BUB1B 0.30, CDCA8 0.29
Modules: M1 (17 genes, Jaccard 0.61, pattern); M2 (15 genes, Jaccard 0.52, weak); M3 (14 genes, Jaccard 0.40, dissolved); M4 (12 genes, Jaccard 0.43, dissolved)
Composition ρ: StromalScore -0.07, ImmuneScore +0.12, MCP_Fibroblasts -0.08, MCP_Endothelial cells +0.08, MCP_Cytotoxic lymphocytes -0.19, MCP_T cells +0.31, MCP_Monocytic lineage +0.16

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 78/59 | +0.01 (-0.32 to +0.36) | +0.02 (-0.35 to +0.38) | not associated |
| BRAF positive | 15/122 | +0.16 (-0.42 to +0.72) | +0.13 (-0.41 to +0.72) | not associated |

### Wnt

Top genes (item–rest r): DLL1 0.43, ADAM17 0.42, HDAC11 0.40, TCF7 0.36, PSEN2 0.36, WNT5B 0.34, PPARD 0.32, SKP2 0.31
Modules: M1 (6 genes, Jaccard 0.65, pattern); M2 (5 genes, Jaccard 0.65, pattern)
Composition ρ: StromalScore +0.10, ImmuneScore +0.01, MCP_Fibroblasts +0.02, MCP_Endothelial cells +0.06, MCP_Cytotoxic lymphocytes +0.19, MCP_T cells +0.07, MCP_Monocytic lineage -0.16

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 78/59 | -0.16 (-0.50 to +0.16) | -0.13 (-0.44 to +0.19) | not associated |
| BRAF positive | 15/122 | +0.39 (-0.08 to +0.91) | +0.37 (-0.16 to +0.88) | not associated |

### p53

Top genes (item–rest r): DDIT4 0.37, DDB2 0.34, STOM 0.34, ABHD4 0.32, EPS8L2 0.30, PRMT2 0.30, RXRA 0.30, AEN 0.29
Modules: M1 (21 genes, Jaccard 0.75, stable); M2 (14 genes, Jaccard 0.48, dissolved); M3 (11 genes, Jaccard 0.46, dissolved)
Composition ρ: StromalScore -0.08, ImmuneScore +0.26, MCP_Fibroblasts +0.16, MCP_Endothelial cells -0.06, MCP_Cytotoxic lymphocytes +0.06, MCP_T cells +0.30, MCP_Monocytic lineage +0.06

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| TP53 positive | 78/59 | -0.11 (-0.45 to +0.21) | -0.09 (-0.42 to +0.20) | not associated |
| BRAF positive | 15/122 | -0.37 (-1.20 to +0.49) | -0.45 (-1.28 to +0.36) | not associated |

## 3. Relationships between scores (partial Spearman given stromal + immune estimates)

| Pair | Spearman | Partial | Random-set null 95 % | Exceeds null |
|---|---|---|---|---|
| EMT | E2F | +0.43 | +0.42 | -0.54 to -0.09 | yes |
| EMT | Wnt | -0.44 | -0.44 | -0.54 to -0.09 | no |
| EMT | p53 | +0.35 | +0.33 | -0.54 to -0.09 | yes |
| E2F | Wnt | -0.12 | -0.12 | -0.27 to +0.07 | no |
| E2F | p53 | +0.42 | +0.40 | -0.27 to +0.07 | yes |
| Wnt | p53 | +0.02 | +0.03 | +0.27 to +0.63 | yes |

## Interpretation limits

- Scores are averaged bulk mRNA; they do not measure pathway activity, protein function or cell-type-specific expression.
- Composition values are computational estimates (ESTIMATE, MCP-counter) from the same expression data, not measured fractions.
- Adjusted genotype effects may be over-adjusted if composition lies on the biological path.
