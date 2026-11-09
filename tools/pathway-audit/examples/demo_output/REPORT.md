# pathway-audit demo (synthetic data with planted truth)

Samples: 300 · genes: 583

## 1. Data-quality gate (positive-control co-expression)

**PASS**: median r = 0.86 over 8 control pairs (threshold 0.5).

## 2. Verdict per gene set

| Gene set | Genes present | Verdict |
|---|---|---|
| STROMAL_PROGRAMME | 42/42 | coherent (mean inter-gene r 0.84); composition-dominated (|rho| 0.95 with StromalScore; holds without shared genes) |
| TARGET_PROGRAMME | 40/40 | coherent (mean inter-gene r 0.33); largely composition-independent (max |rho| 0.01) |
| CONFOUNDED_PROGRAMME | 40/40 | coherent (mean inter-gene r 0.34); largely composition-independent (max |rho| 0.02) |
| RANDOM_SET | 40/40 | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.05) |

### STROMAL_PROGRAMME

Top genes (item–rest r): STR32 0.93, STR14 0.93, STR0 0.93, STR36 0.92, STR16 0.92, STR19 0.92, STR13 0.92, STR3 0.92
Modules: M1 (42 genes, Jaccard 1.00, stable)
Composition ρ: StromalScore +0.95, ImmuneScore -0.04
Gene-disjoint fibroblast test: ρ 0.97 (0.96–0.98); 2 shared genes removed

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| Target mutation | 136/164 | -0.06 (-0.26 to +0.17) | -0.05 (-0.11 to +0.01) | not associated |
| Driver mutation | 82/218 | +0.03 (-0.23 to +0.26) | +0.01 (-0.07 to +0.08) | not associated |

### TARGET_PROGRAMME

Top genes (item–rest r): TGT20 0.64, TGT9 0.62, TGT22 0.61, TGT2 0.60, TGT26 0.60, TGT16 0.59, TGT3 0.59, TGT14 0.58
Modules: M1 (40 genes, Jaccard 0.92, stable)
Composition ρ: StromalScore +0.01, ImmuneScore +0.01

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| Target mutation | 136/164 | -1.96 (-2.01 to -1.91) | -1.96 (-2.00 to -1.91) | robust to adjustment |
| Driver mutation | 82/218 | -0.22 (-0.47 to +0.04) | -0.36 (-0.63 to -0.04) | robust to adjustment |

### CONFOUNDED_PROGRAMME

Top genes (item–rest r): CNF17 0.65, CNF32 0.64, CNF7 0.63, CNF20 0.61, CNF6 0.61, CNF24 0.61, CNF2 0.60, CNF15 0.60
Modules: M1 (38 genes, Jaccard 0.94, stable)
Composition ρ: StromalScore -0.02, ImmuneScore +0.02

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| Target mutation | 136/164 | -0.06 (-0.29 to +0.18) | +0.00 (-0.04 to +0.05) | not associated |
| Driver mutation | 82/218 | +1.12 (+0.87 to +1.39) | +0.02 (-0.03 to +0.08) | explained by adjustment covariates |

### RANDOM_SET

Top genes (item–rest r): G17 0.12, G1 0.09, G0 0.07, G34 0.05, G2 0.05, G28 0.05, G27 0.03, G8 0.03
Composition ρ: StromalScore -0.04, ImmuneScore +0.05

| Annotation | n (+/−) | Unadjusted SD diff (95 % CI) | Adjusted (95 % CI) | Verdict |
|---|---|---|---|---|
| Target mutation | 136/164 | -0.11 (-0.33 to +0.07) | -0.11 (-0.34 to +0.11) | not associated |
| Driver mutation | 82/218 | +0.18 (-0.05 to +0.47) | +0.15 (-0.18 to +0.46) | not associated |

## 3. Relationships between scores (partial Spearman given stromal + immune estimates)

| Pair | Spearman | Partial | Random-set null 95 % | Exceeds null |
|---|---|---|---|---|
| STROMAL_PROGRAMME | TARGET_PROGRAMME | +0.02 | +0.05 | +0.02 to +0.37 | no |
| STROMAL_PROGRAMME | CONFOUNDED_PROGRAMME | -0.01 | +0.01 | +0.02 to +0.37 | yes |
| STROMAL_PROGRAMME | RANDOM_SET | -0.03 | +0.03 | +0.02 to +0.37 | no |
| TARGET_PROGRAMME | CONFOUNDED_PROGRAMME | +0.02 | +0.02 | -0.01 to +0.50 | no |
| TARGET_PROGRAMME | RANDOM_SET | +0.05 | +0.05 | -0.01 to +0.50 | no |
| CONFOUNDED_PROGRAMME | RANDOM_SET | -0.03 | -0.04 | -0.01 to +0.44 | yes |

## Interpretation limits

- Scores are averaged bulk mRNA; they do not measure pathway activity, protein function or cell-type-specific expression.
- Composition values are computational estimates (ESTIMATE, MCP-counter) from the same expression data, not measured fractions.
- Adjusted genotype effects may be over-adjusted if composition lies on the biological path.
