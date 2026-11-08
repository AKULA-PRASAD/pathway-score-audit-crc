# Optimism-corrected value of transcriptomic pathway scores for colon-cancer recurrence

A pre-specified, audited analysis of whether four Hallmark pathway scores improve recurrence-risk discrimination beyond clinical covariates in stage II–III colon cancer (GEO GSE39582). The headline result is a nested-bootstrap confidence interval for the optimism-corrected gain in Harrell's C.

**Status (8 November 2026):**
- **Scientifically complete:** the frozen verification passed (60 of 60 checks), and an independent checkpoint-level audit reproduced every reported number.
- **Repository disposition:** owner decisions approved, and the project committed locally. Push, the checkpoint-archive release, the code licence and the public-release review are pending (`docs/PROJECT_COMPLETION_INVENTORY.md`).

## Research question

Does adding reproducible transcriptomic pathway information improve recurrence-risk prediction beyond clinical characteristics? Are any improvements stable enough to be statistically and biologically credible?

Recurrence risk after surgery for stage II–III colon cancer varies widely, and better risk discrimination could in principle inform decisions about adjuvant chemotherapy; this study does not show such utility. Expression signatures are often reported to add prognostic value in colorectal cancer, but part of an apparent gain is optimism: a larger model fits its own data better than new patients. This project therefore:
- fixes a small, biologically motivated panel before looking at outcomes;
- corrects for optimism;
- puts an interval on the corrected gain, so that both its size and its stability are measured.

## Main findings

> Conditional on the observed 136 recurrence events and 316 censored observations in the 452-patient A-19 population.

| | |
|---|---|
| Optimism-corrected C | M1 (clinical) **0.6016**; M2 (clinical + 4 pathways) **0.6059** |
| ΔC_corrected (pre-specified point estimate) | **+0.0043** |
| 95 % percentile intervals, five independent runs (B_out = 750, B_in = 100) | [−0.0086, +0.0519] · [−0.0090, +0.0469] · [−0.0093, +0.0484] · [−0.0091, +0.0518] · [−0.0091, +0.0494] |
| Run-to-run SD of the interval endpoints | 0.00025 (lower), 0.00217 (upper) |
| Pathway hazard ratios, primary four-set model (per SD of score) | E2F 0.91, EMT 1.07, Wnt 1.02, p53 1.00; every 95 % CI includes 1. In the five-set sensitivity model only, E2F is 0.86 (0.75–0.996), with correlated scores; not a primary finding |
| Corrected calibration slope | M1 0.83, M2 0.86 (residual overfitting) |
| Reliability | no hard failure. Warnings: the bootstrap distribution sits above the point estimate (z₀ ≈ −0.5), and the percentile and basic intervals diverge |

**Interpretation:**
- Adding the four pathway scores was **not shown** to improve discrimination. The data are compatible with no gain and with a gain of up to about 0.05 in C.
- All percentile, basic and bias-corrected intervals include zero. That does not establish equivalence.
- No pathway shows an independent association with recurrence in the primary model.
- E2F (< 1) and EMT (> 1) keep consistent directions, which is biologically plausible (proliferation; mesenchymal/stromal programmes) but not demonstrated here. Pathway scores are averaged mRNA levels, not measures of pathway activation, protein activity or mutation status.
- The corrected ΔC is stable across Monte Carlo runs and expression representations. Molecular coefficients are not stable across processing batches: leaving out one batch at a time changes the sign of the Wnt and p53 coefficients. Batch effects are characterised, not resolved.
- A historical "paired out-of-original" ΔC interval that narrowly excluded zero is a different quantity, not an interval for the corrected ΔC, and it included zero under the alternative expression representation.
- There is no claim of clinical utility, absolute-risk prediction, external validity, transportability or causal mechanism.

**Full report:** `docs/PHASE9_CORRECTED_FINAL_REPORT.md`.

## Data, cohort and endpoint

| | |
|---|---|
| Dataset | GEO GSE39582 (Marisa et al., *PLoS Med* 2013): 585 Affymetrix HG-U133 Plus 2.0 arrays, deposited as RMA + ComBat |
| Cohort | 453 stage II–III colon tumours with recurrence follow-up. QC retained all of them |
| Analysis population | 452 patients with complete clinical covariates: 136 recurrences, 316 censored |
| Endpoint | cause-specific recurrence as deposited. Deaths without recurrence are censored, and the recurrence-time unit is undocumented, so only rank-based quantities (C, hazard ratios) are reported |
| Expression representation | gene level, median of uniquely mapped probes (primary); representative-IQR probe (sensitivity). Ambiguous and unmapped probes excluded |

## Methods

- **M1:** Cox model with age (restricted cubic spline), sex, TNM stage, tumour location and adjuvant chemotherapy.
- **M2:** M1 plus four frozen MSigDB Hallmark scores (E2F targets, EMT, p53, Wnt/β-catenin; 607 genes). Each score is the mean of training-fold-standardised member genes, then standardised again. The molecular block carries a ridge penalty, with λ tuned by five-fold event-stratified cross-validation.
- **No leakage:** every estimated quantity (spline knots, gene scaling, score scaling, λ) is re-estimated inside every fit from that fit's training data.
- **Estimand:** ΔC_corrected = C_corrected(M2) − C_corrected(M1), using Harrell's bootstrap optimism correction.
- **Interval (nested bootstrap):**
  ```
  for k in 1..B_out:                       outer: event-stratified bootstrap D*_k
      fit M1, M2 on D*_k                   -> apparent ΔC(k)
      for b in 1..B_in:                    inner: bootstrap of D*_k, evaluated on D*_k
          fit M1, M2 (paired)              -> optimism(M, k, b)
      ΔC_corrected(k) = ΔC_app(k) − mean_b[opt(M2) − opt(M1)]
  interval = 2.5th and 97.5th percentiles, in each of five independent runs
  ```
- **Monte Carlo control:** kernel quantile SE at three bandwidths. B_out grows along a fixed ladder until the endpoint SE is ≤ 5 % of the width (EP-1). Here 500 was topped up to 750.

## Validation and audit

- **Internal validation only.** No external validation was performed. The CIT partition of GSE39582 is confounded with batch and has 29 events, so it is not an external test.
- **Corrected pipeline** (`src/phase9_corrected/`). It fixes a random-stream deviation in an earlier run, where cross-validation reused the bootstrap resample's random state, and calls the original fitting code unchanged.
- **Code identity:** 27 source files are frozen by SHA-256, with a launcher that checks them before Python starts.
- **Execution integrity:** a pre-execution gate, append-only hash-linked checkpoints, and crash-consistent invocation records.
- **Frozen post-run verifier** (`scripts/43`): PASS, 60/60.
  - It regenerates every resample (3,750 outer, 375,000 inner) and checks stream assignment.
  - It recomputes the sizing decisions.
  - It recomputes 20 replicates bit for bit.
- **Independent audit** (`scripts/48`, standard library only): it rebuilt every interval quantity from raw checkpoint bytes and matched the official results within 2.2×10⁻¹⁶ (`reports/phase9/independent_audit/`).
- **Point-estimate alignment:** shown in `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md`.

## Limitations

- **Single cohort, single platform:** internal validation only.
- **Few events:** 136 events limit precision.
- **Batch:** processing batch is confounded with the data partition, and molecular coefficients are unstable across batches; batch was examined but not perturbed.
- **Endpoint (U1, U2):** the recurrence-time unit is undocumented, and deaths without recurrence (10 %) are censored, which is informative. The models estimate the cause-specific recurrence hazard, not the cumulative incidence of recurrence (competing-risk analysis) or composite disease-free survival. No absolute-risk or calibration-in-the-large analysis is possible while the time unit is undocumented, and no competing-risk analysis was performed.
- **Discrimination and calibration measures:** the C-index reflects ranking only and says nothing about clinical usefulness; no decision-curve or reclassification analysis was done. Calibration is assessed only by log-hazard slopes.
- **Modelling:** only the molecular block is penalised, so the comparison is conservative by design. Overfitting is visible even in the clinical model (about 13.6 events per M2 coefficient).
- **Bootstrap bias:** diagnostics are triggered in every run, in the opposite direction to the protocol's anticipated mechanism, and are unexplained.
- **Partly historical:** calibration, the representation sensitivity (S-7) and the five-set sensitivity panel come from the earlier single-layer execution and were not re-run under the corrected streams.

## Repository layout

| Path | Contents |
|---|---|
| `src/phase9_execution.py`, `src/pathway_scoring.py`, `src/score_standardisation.py` | model pipeline and pathway scoring (frozen) |
| `src/phase9_corrected/` | nested-bootstrap engine, random streams, checkpointing, orchestration, statistics, gate, verification |
| `scripts/01`–`14` | data download, metadata and probe audits, probe cache |
| `scripts/15`–`39` | Phases 1–9: cohort, gene mapping, bundle, QC, EDA, protocol, Phase 9 primary execution (`38`) and S-7 (`39`) |
| `scripts/40`–`48` | as-run nested driver (historical), corrected driver (`41`), launcher (`42`), verifier (`43`), freeze tooling (`44`), sizing writer (`45`), gated report (`46`), θ̂ alignment check (`47`), independent audit (`48`) |
| `config/analysis/`, `config/pathways/` | run configuration with owner decisions, freeze record, QC and amendment records, frozen pathway panels |
| `docs/` | phase reports, protocols (`PHASE9_MODELING_PROTOCOL_DRAFT.md` v3.1, `PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1), audits, final report |
| `reports/phase9/corrected_rerun/` | corrected campaign: checkpoints, invocation records, gate evidence, sizing, verification reports |
| `reports/phase9/corrected_rerun_report/` | final interval results (JSON) and tables |
| `reports/phase9/independent_audit/` | independent audit results and report |
| `reports/phase9/execution*/`, `dc_interval/`, `final_project/` | historical Phase 9 outputs, kept unchanged; the as-run intervals are superseded |
| `tests/` | 1,085 tests |

## Environment

- **Python:** 3.13.13 with NumPy 2.5.4, pandas, lifelines, SciPy, pyarrow and pytest. The exact versions are in `requirements.txt` (original capture) and `requirements-analysis.txt` (additions).
- **Platform:** macOS arm64. The freeze record pins the environment fingerprint, including the BLAS/LAPACK build, and the launcher needs system `jq` and `shasum`.

```bash
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt -r requirements-analysis.txt
.venv/bin/python -m pytest -q tests/
```

Some tests read the locally built study inputs and historical outputs. They pass only after the data pipeline below has been run.

## Reproduction

1. **Data** (public, about 356 MB): `scripts/01_download_geo.py` and `scripts/09_download_gpl570_annotation.py`. Then build the caches and the analysis bundle with scripts 13, 16–23 in order. `data/raw/CHECKSUMS.sha256` and the manifests pin every input.
2. **Historical Phase 9 primary estimate:** `scripts/38_phase9_execute.py`.
3. **Corrected campaign:** see `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md`.
   ```bash
   unset OMP_NUM_THREADS MKL_NUM_THREADS OPENBLAS_NUM_THREADS VECLIB_MAXIMUM_THREADS NUMEXPR_NUM_THREADS
   CFG=config/analysis/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json
   for r in 0 1 2 3 4; do scripts/42_phase9_corrected_launch.sh run $CFG $r 500; done
   # a top-up beyond 500 requires the sizing decision (scripts/45) and an explicit authorisation reference
   scripts/42_phase9_corrected_launch.sh verify $CFG
   P=$(mktemp -d); .venv/bin/python -I -B -X pycache_prefix=$P scripts/46_phase9_corrected_report.py $CFG <out_dir>
   python3 -I -B scripts/48_phase9_independent_audit.py <audit_dir>
   ```
   One run takes about 3–5 h at 500 replicates on 8 cores.

**Without re-running, from the archived checkpoints:**
- restore the archive (below);
- run the verifier, the report tool and the independent audit. These recompute every reported number in under an hour.

## Checkpoint archive

The five corrected checkpoints (113 MB each) exceed GitHub's 100 MB file limit.
- `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt` pins all 52 campaign files.
- `docs/LARGE_ARTIFACT_ARCHIVE.md` gives the archive, restore and integrity-check commands:
  ```bash
  shasum -a 256 -c reports/phase9/ARCHIVE_MANIFEST_SHA256.txt
  ```

The archive destination is pending the owner's decision.

## Data access and licensing

- **GSE39582:** public GEO data. Cite Marisa L. et al., *PLoS Med* 2013;10:e1001453. Raw and processed data are not redistributed here: download them from GEO.
- **Hallmark gene sets:** from MSigDB (Liberzon et al., *Cell Syst* 2015). Use is subject to MSigDB's terms. The panel files record set identifiers and resolved members.
- **Code licence:** no licence file has been chosen yet.
