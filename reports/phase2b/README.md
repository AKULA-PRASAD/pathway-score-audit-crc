# Phase 2B evidence — file status

## Superseded files (preserved, do not use)

**`gene_reliability_alpha.csv`** — produced by revision 1 of
`scripts/21_phase2b_probe_quality.py`. It contains a `reliability_gain`
column whose interpretation has been **WITHDRAWN**.

That column was computed as `alpha - rbar`, where `rbar` was back-solved from
`alpha` and `k` via the Spearman-Brown identity. It is therefore a
deterministic algebraic function of `alpha` and `k` (verified: recomputation
from alpha and k alone reproduces it to within 8.9e-16). It carries no
information beyond alpha and k and measures no external accuracy.

The file is retained unaltered for audit. **Do not use `reliability_gain`.**

## Current file

**`gene_probe_agreement_alpha.csv`** — current output. Reports `alpha`,
`k_probes`, `implied_mean_inter_probe_r` and `within_gene_rho`.

`implied_mean_inter_probe_r` is included ONLY to expose alpha's mechanical
dependence on probe count; it is an algebraic restatement, not an independent
measurement.

Cronbach's alpha here is **descriptive evidence of probe co-variation only**.
It is not a reliability coefficient in the psychometric sense (the
uncorrelated-error assumption is violated by shared array, sample prep,
hybridisation and the depositors' joint RMA + ComBat), it is not an accuracy
measure, and it is not an aggregation decision rule. See
`docs/PHASE2B_SELECTION_SENSITIVITY.md` section 4.

## Approved representation

- **Primary:** `median` aggregation (provisional, uniform rule; not claimed biologically optimal)
- **Sensitivity:** `representative_iqr` selected on `discovery_443`
- **Historical only:** `mean`
- **Not approved:** alpha-stratified hybrid (no alpha threshold exists in any code path)

Frozen artifacts: `data/processed/frozen/`
