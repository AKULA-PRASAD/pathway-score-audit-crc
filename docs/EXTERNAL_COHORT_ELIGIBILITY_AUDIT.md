# External-cohort eligibility audit (registered analysis X8)

**Date:** 8 November 2026.
**Code:** `scripts/exploratory/x8_external_audit.py`.
**Outputs:** `reports/exploratory/x8_external_audit/eligibility_table.csv` and `x8_summary.json`, plus figure `reports/exploratory/figures/XF6_external_audit.png`.
**Inputs:** the recorded search and per-sample characteristics text in `data/external/exploratory/geo_search/`. No expression data were downloaded.

**Result: no independent public cohort is compatible with frozen-model validation, so none was attempted.** The CIT partition of GSE39582 remains internal validation and is not used here.

## 1. Search

| Step | Detail |
|---|---|
| GEO DataSets (E-utilities, 2026-11-08) | Query in `geo_search/query.txt`: colorectal/colon/rectal × cancer/carcinoma × recurrence/relapse/disease-free/relapse-free/recurrence-free/prognosis × *Homo sapiens* × series × array or sequencing expression profiling. **554 hits** |
| Screen | ≥ 80 samples, human tumour expression, outcome keywords: 85 series. Non-colorectal series, cell lines, non-mRNA assays and metastasis-only series were then set aside on reading |
| Assessed | 14 cohorts whose per-sample characteristics include a recurrence or disease-free outcome (table below). Their characteristics text was downloaded (`geo_search/chars/`) |
| ArrayExpress / BioStudies | Three queries. Non-GEO hits were stage-II-only Almac FFPE series (E-MTAB-863/864), exon-array or methodological series, and cell lines. None qualified. A "PETACC" query returned an error page, so PETACC-3 was **not assessed** from the repository |
| Prior review | Phase 4/5 cohorts (GSE14333, GSE33113, GSE17536 overlap) re-derived. GSE14333 (48 events) and GSE33113 (18 events) reproduce Phase 5 exactly |

## 2. Basis for the precision requirement (no unexplained threshold)

- External validation of a prognostic model is limited mainly by the **number of events**.
- The calibration slope is the measure that needs most events (Riley et al., *Stat Med* 2022;41:1280–1295). Their relaxed target is a slope 95 % CI width of **0.4** (SE ≈ 0.10).
- Phase 5 simulated this for the present design: at C = 0.65, 136 events give a slope CI width of 0.63 and a C-index CI width of 0.0945.
- Width scales with 1/√events (the scaling Phase 5 showed), so for each cohort here:
  - slope width ≈ 0.63 × √(136 / events);
  - C width ≈ 0.0945 × √(136 / events).
- The relaxed target therefore needs **≈ 337 events**.
- **This number is a precision consequence, not an exclusion rule.** Every cohort's expected precision is reported, and eligibility rests on the compatibility criteria in §3.

**What smaller cohorts could and could not show.** A ΔC of +0.004 cannot be distinguished from zero with any cohort here: the C-index CI half-width alone is ≥ 0.05 in every one.

## 3. Compatibility criteria (fixed in the registered plan)

1. Independent of GSE39582 and the CIT programme.
2. A recurrence endpoint compatible with "time from surgery to first recurrence, deaths without recurrence censored".
3. All five M1 predictors deposited: age, sex, stage, tumour location, adjuvant chemotherapy.
4. Exact reproduction of the frozen scores from the deposit: HG-U133 Plus 2.0 with RMA-type processing. A different platform or summarisation changes the measurement itself.

Frozen-model validation is run only if criteria 1–4 all hold. Precision (§2) is reported in every case.

## 4. Assessment

| Cohort | Institutions | Platform (deposited processing) | Endpoint | Stage II–III n / events | Approx. C CI width | Approx. slope CI width | 1 Indep. | 2 Endpoint | 3 All M1 predictors | 4 Exact scores |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE24551 | Oslo | HuEx exon (GPL5175/11028) | DFS (event/time fields swapped in deposit) | 320 / 124 | 0.099 | 0.66 | yes | unverified | no (stage, MSI only) | no |
| GSE58394 | Japan | Agilent (GPL18709) | recurrence status, no time | 300 / 62 (stage II only) | 0.140 | 0.93 | yes | no (no time) | no | no |
| GSE106584 | US | HTA 2.0 (GPL17586) | recurrence, DFS days | 131 / 52 | 0.153 | 1.02 | yes | partly | no (stage, adjuvant only) | no |
| GSE14333 | Melbourne + Moffitt | U133 Plus 2.0 (MAS5) | DFS (DFS_Cens 0 = event) | 185 / 48 (colon only 163 / 42) | 0.159 | 1.06 | yes (overlaps GSE17536, not GSE39582) | unverified (death handling not stated) | **yes** | no (MAS5; includes rectum) |
| GSE103479 | Aberdeen, Barcelona, Florence, Dublin | Almac Xcel (RMA + ComBat) | recurrence + PFS time | 156 / 42 | 0.170 | 1.13 | yes | unverified | **yes** | no (different array) |
| GSE161158 | Moffitt + Vanderbilt | U133 Plus 2.0 (RMA + ComBat) | DFS | 158 / 41 | 0.172 | 1.15 | yes | unverified | no (age, stage only) | yes |
| GSE37892 | six French centres | U133 Plus 2.0 | distant metastasis (dates) | 130 / 37 | 0.181 | 1.21 | cannot be excluded | no (distant metastasis only) | no (no chemotherapy) | no |
| GSE17536 | Moffitt | U133 Plus 2.0 (RMA) | DFS (recurrence) | 114 / 31 | 0.198 | 1.32 | yes | partly | no | yes |
| GSE41258 | MSKCC | U133A (PLIER) | recurrence with interval | 99 / 21 | 0.240 | 1.60 | yes | partly | no (no chemotherapy) | no |
| GSE33113 | AMC Amsterdam | U133 Plus 2.0 (MAS5) | metastasis or recurrence within 3 years | 90 / 18 (stage II) | 0.260 | 1.73 | yes | no (composite, truncated) | no | no |
| GSE143985 | Japan | U133 Plus 2.0 (MAS5) | DFS | 91 / 15 | 0.285 | 1.90 | yes | unverified | no | no |
| GSE38832 | Vanderbilt and others | U133 Plus 2.0 (MAS5) | DFS | 74 / 9 | 0.367 | 2.45 | yes | unverified | no | no |
| GSE17537 | Vanderbilt | U133 Plus 2.0 (RMA) | DFS (recurrence) | 34 / 3 | 0.636 | 4.24 | yes | partly | no | yes |
| GSE81653 | not stated | HuGene 2.0 | recurrence status only (234 events, all stages) | — | — | — | not stated | no (no time, no stage) | no | no |

Event counts come from the deposited characteristics. GSE161158 is a re-deposit spanning the GSE17536/GSE17537 institutions, so those three are not independent of one another and must never be pooled.

## 5. The two closest candidates

- **GSE14333:**
  - **Meets:** all five M1 predictors are deposited; it is independent of GSE39582; it uses the same array.
  - **Fails or unresolved:**
    - the deposited values are MAS5, not RMA + ComBat, so the frozen score scale is not reproduced;
    - it includes rectal tumours;
    - "DFS" may count deaths, and the deposit does not say (Phase 4).
  - **Precision:** 48 events (42 colon), giving a slope CI width of ≈ 1.06.
  - **What it would take:** CEL-level RMA reprocessing (raw CEL files are available on GEO, about 185 arrays) and written confirmation of the DFS event definition. Even then, the precision would support only a qualitative direction check.
- **GSE103479:**
  - **Meets:** all five predictors are deposited.
  - **Fails:** it uses a different array (Almac Xcel), so the frozen scores cannot be reproduced exactly.
  - **Precision:** 42 events.

## 6. Decision

- No cohort meets criteria 1–4, so **no frozen-model validation was performed**.
- No exploratory "direction check" was run either. The best candidate (GSE14333) would need CEL reprocessing to reproduce the measurement. Its 42–48 events give a C-index CI half-width of about 0.08, about twenty times the effect in question, so the check could not inform the scientific question.
- The owner may still choose to commission it as a clearly labelled exploratory transport check.
- Events must never be pooled across cohorts with different endpoints, platforms or overlapping patients.

**Status:** external validity **NOT ESTABLISHED**. A defensible validation is **NOT FEASIBLE** with public data identified on 8 November 2026.

**Evidence that would resolve it:** an independent stage II–III colon cancer cohort with:
- about 340 or more recurrence events (precision-driven);
- recurrence-free survival with documented time origin, unit and death handling;
- age, sex, stage, tumour location and adjuvant chemotherapy recorded;
- fresh-frozen tumour RNA profiled on HG-U133 Plus 2.0 (or raw data reprocessable to the same RMA scale).

Candidate sources are trial biobanks (for example PETACC-3, NSABP C-07/C-08, QUASAR) under data-access agreements, or prospective consortium data.
