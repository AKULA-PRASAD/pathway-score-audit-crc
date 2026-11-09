# pathway-audit case study: five colorectal cancer cohorts

Gene sets: the four frozen Hallmark sets of the GSE39582 study. Composition: ESTIMATE + MCP-counter (computational estimates).

| Cohort | Samples | QC (median r) | EMT | p53 |
|---|---|---|---|---|
| GSE39582 | 452 | PASS (0.84) | coherent (mean inter-gene r 0.24); composition-dominated (|rho| 0.95 with MCP_Fibroblasts; holds without shared genes) | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.39) |
| GSE41258 | 183 | PASS (0.79) | coherent (mean inter-gene r 0.24); composition-dominated (|rho| 0.95 with MCP_Fibroblasts; holds without shared genes) | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.21) |
| GSE143985 | 91 | PASS (0.78) | coherent (mean inter-gene r 0.30); composition-dominated (|rho| 0.94 with MCP_Fibroblasts; holds without shared genes) | incoherent (mean inter-gene r < 0.05); composition-associated (|rho| 0.44 with MCP_Endothelial cells) |
| GSE103479 | 156 | FAIL (0.07) | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.14) | incoherent (mean inter-gene r < 0.05); largely composition-independent (max |rho| 0.30) |
| TCGA-COAD | 286 | PASS (0.94) | coherent (mean inter-gene r 0.38); composition-dominated (|rho| 0.98 with MCP_Fibroblasts; holds without shared genes) | incoherent (mean inter-gene r < 0.05); composition-associated (|rho| 0.52 with ImmuneScore) |

Per-cohort reports: `<cohort>/REPORT.md`.