# G1–G13 implementation record

**Approved by the owner:** 8 November 2026, covering the gaps listed in `docs/PHASE9_SCIENTIFIC_INTERPRETATION_REVIEW.md` §10.
**Scope:** documentation only. No number, model, endpoint, cohort, configuration, checkpoint, frozen code or historical analysis was changed. No analysis was run. The five corrected-ΔC intervals are unchanged.

## Completion table

FR = `docs/PHASE9_CORRECTED_FINAL_REPORT.md`; ER = `docs/PHASE9_EXECUTION_REPORT.md`; S7 = `docs/PHASE9_S7_SENSITIVITY_REPORT.md`; P6A = `docs/PHASE6A_ENDPOINT_RESOLUTION.md`; EX = `reports/phase9/execution/phase9_execution_results.json`.

| Gap | Change | Location | Source verified |
|---|---|---|---|
| G1 | Added a leave-one-batch-out coefficient table (EMT, E2F, Wnt, p53, stage: all-data value, range, sign changes), the batch size and event range, and the sample-size caveat. Added an "Implication" paragraph stating that batch effects are characterised, not resolved. Added a summary row | FR Summary ("Batch"); FR §12; README Interpretation and Limitations | EX `batch_sensitivity.leave_one_batch_out` (all 5 ranges reproduced to 3 dp); ER §10 table (14–103 patients, 4–52 events) |
| G2 | Batch-stratified fit: molecular Δ ≤ 0.027; stage 0.662 → 0.520 (HR 1.94 → 1.68); chemotherapy 0.137 → 0.029 | FR §12 | EX `batch_sensitivity.stratified_fit` (max molecular Δ 0.0266); ER §10 |
| G3 | Restricted "every pathway CI includes 1" to the primary model. Added §9.1 with the five-set HR table, including E2F 0.862 (0.746–0.996) labelled a sensitivity-model result, the collinearity evidence and the reason it is not promoted | FR Summary; §9; new §9.1; §14; README findings table and Interpretation; `docs/FINAL_PROJECT_HANDOVER.md` | EX `apparent.M2_sensitivity` (all 5 HRs and λ 0.152); ER §5 (r 0.636, condition 2.909 → 7.084, SE ×1.20–1.34, −0.099 → −0.149) |
| G4 | New §8.1 on the paired out-of-original ΔC (median +0.0121 [+0.0002, +0.0240]; S-7 +0.0105 [−0.0001, +0.0224]): why it is a different quantity and must not replace the corrected intervals. Added a README bullet | FR §8.1; README Interpretation | ER §7; S7 §4 |
| G5 | Added the penalisation-asymmetry statement (conservative by design; not an unbiased molecular estimate) | FR §5; README Limitations ("Modelling") | ER §15 |
| G6 | S-7 table rows for the five-set corrected C (0.6066 → 0.6041, exceeds MC-SE 0.00041) and p53↔PI3K correlation (0.636 → 0.688). Bullets on collinearity (condition 7.08 → 9.55; E2F↔PI3K r 0.545; PI3K +0.0417 → +0.0103) and on S-7 not being replication | FR §11 | S7 §4, §8, §10; EX `monte_carlo.M2_sensitivity.monte_carlo_SE` |
| G7 | Informative-censoring profile (median age 75.0 vs 68.0; stage III 32.6 % vs 38.7 %; chemotherapy 23.9 % vs 42.1 %). New "Estimand" row distinguishing the cause-specific hazard, cumulative incidence (Fine–Gray / CIF) and composite DFS (182 events). New limitation 8 on competing risks | FR §2; FR §13 item 8; README Limitations ("Endpoint") | P6A §2.3 table; §3.2.2–3.2.3 |
| G8 | New limitation 9 on what C cannot show (ranking only; no decision-curve, reclassification or time-dependent analysis). New limitation 10 on calibration (slopes only) | FR §13 items 9–10; README Limitations ("Discrimination and calibration measures") | FR §10; P6A U1 |
| G9 | Replaced the unreferenced "consistent with some published reports" with new §9.2: a plausibility-vs-evidence table for all five pathways and an explicit statement that scores are mRNA averages, not activation, protein activity, mutation status or mechanism. CMS citation: Guinney et al., *Nat Med* 2015. Removed the unsupported §14 phrase "largely captured by stage" | FR §9 (Direction bullet), new §9.2, §14 "Biological hypothesis"; README Interpretation | Observed values from EX and S7 §7; the background is labelled as not tested |
| G10 | Removed the comparison of the leave-one-batch-out apparent-C range with the ΔC interval distance | FR §12 | — |
| G11 | Added the status of the three ER §15 open items: S-7 executed; corrected-ΔC interval obtained; batch instability characterised but not resolved | FR §14 | ER §15; S7; FR §8 |
| G12 | Events per coefficient (136/10 ≈ 13.6; protocol EPV 15.0–15.1) | FR §2 ("Events per coefficient"); FR §13 item 11; README Limitations | `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` §2.1 |
| G13 | M2 stage HR 1.94 (1.28–2.94) attenuating to 1.68 (1.07–2.65) under batch stratification | FR §9 "Clinical covariates" | EX `apparent.M2_primary`, `stratified_fit` (recomputed HR and CI from the stored coefficient and SE) |

**Also corrected (item 8 of the approval):** the θ̂ RNG disclosure in FR §6 now states:
- the strict analytical worst-case bound of 1.6×10⁻³, with its derivation;
- the typical-size estimate of about 1.9×10⁻⁵, conditional on the linked replicate behaving typically;
- that neither figure is an empirical measurement;
- that the primary percentile interval does not depend on θ̂.

## Quality control

| Check | Result |
|---|---|
| Inserted numbers vs sources | All EX-derived values (five-set HRs and λ; M2 stage HR; stratified HR, CI and coefficient change; leave-one-batch-out ranges and sign counts; S-7 five-set C and MC-SE) were recomputed from the JSON and matched as printed. All document-derived values (P6A, ER, S7, protocol) were found verbatim in the cited source |
| Corrected-ΔC intervals | Unchanged. Report, README and handover strings equal the values in `corrected_interval_results.json`. The per-run 5-dp endpoints in FR §8 equal the independent-audit values exactly |
| Primary vs sensitivity vs S-7 vs historical | Separate labelled subsections (FR §8.1, §9.1, §11); every historical single-layer result labelled |
| Conclusion | Unchanged and restrained: "not demonstrated", "inconclusive for a modest gain", "not evidence of equivalence"; no utility, absolute-risk, external-validity, transportability or causal claim |
| Cross-document consistency | Final report, README, handover, inventory and independent audit agree on θ̂, corrected C, intervals, between-run SD, verification and audit classification. Handover and release-readiness wording updated for G3 |
| Test suite | see the final response |

## Files changed in this step
- `docs/PHASE9_CORRECTED_FINAL_REPORT.md`: G1–G13 and the θ̂ disclosure.
- `README.md`: findings table, Interpretation and Limitations (G1, G3, G4, G5, G7, G8, G9, G12).
- `docs/FINAL_PROJECT_HANDOVER.md`: G3 wording in the conclusion row.
- `docs/PHASE9_SCIENTIFIC_INTERPRETATION_REVIEW.md`: approval and applied note at the end of §10.
- `docs/release/RELEASE_READINESS.md`: status of discrepancy 3, the checklist and the hash note.
- `docs/release/G1_G13_IMPLEMENTATION.md`: this file.
