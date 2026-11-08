# Final project handover — GSE39582 colon-cancer recurrence-risk modelling

**Date:** 8 November 2026. **Branch:** `phase9-corrected-pipeline`.

**Status: PROJECT SCIENTIFICALLY AND TECHNICALLY COMPLETE — PENDING OWNER SIGN-OFF.**

## 1. Scientific completion

| | |
|---|---|
| Research question | Does adding reproducible transcriptomic pathway information improve recurrence-risk prediction beyond clinical characteristics, and is any improvement stable enough to be credible? |
| Cohort | GSE39582, stage II–III colon cancer; A-19 analysis population of 452 patients, 136 recurrences, 316 censored |
| Endpoint | cause-specific recurrence as deposited (rank-based analysis only: U1 time unit undocumented; U2 deaths censored) |
| Models | M1 clinical Cox (6 coefficients) vs M2 = M1 + four frozen Hallmark scores (E2F, EMT, p53, Wnt; ridge on the molecular block; 10 coefficients). Five-set panel is sensitivity-only |
| Corrected C | M1 0.6016, M2 0.6059 |
| **ΔC_corrected** | **+0.0043**. 95 % percentile intervals in five independent runs: [−0.0086, +0.0519], [−0.0090, +0.0469], [−0.0093, +0.0484], [−0.0091, +0.0518], [−0.0091, +0.0494] |
| Frozen verification | PASS, 60/60; reporting permitted |
| Independent audit | **INDEPENDENTLY REPRODUCED**: 3,750 outer and 375,000 inner pairs; 138 quantities within 2.2×10⁻¹⁶ |
| Conclusion | No demonstrated improvement. The data are compatible with zero and with a gain of up to about 0.05 in C. Primary-model pathway HRs are all compatible with 1; the only excluding interval (five-set E2F, 0.746–0.996) is in a collinear sensitivity model. Molecular coefficients are unstable across batches. Corrected calibration slopes of 0.83/0.86 indicate residual overfitting. The result is stable across Monte Carlo runs and representations, but apparent performance varies with batch. No equivalence, clinical-utility, absolute-risk, external-validity, transportability or causal claim |
| Main limitations | internal validation only; 136 events; batch confounding; U1/U2; unexplained bootstrap-bias warnings (A7/A8) in every run; calibration, S-7 and the five-set results are historical |

Full report: `docs/PHASE9_CORRECTED_FINAL_REPORT.md`.

## 2. Technical completion

**Pipeline architecture:**
1. **Data:** GEO download, then probe cache, then outcome-blind probe-to-gene mapping (43,163 unique probes → 21,815 genes; median aggregation), then the hash-recorded analysis bundle.
2. **Model pipeline** (`src/phase9_execution.py`, `src/pathway_scoring.py`, `src/score_standardisation.py`):
   - spline knots, gene scaling, score scaling and λ (five-fold event-stratified CV) are all re-estimated inside each fit;
   - Cox fitting uses lifelines.
3. **Corrected nested bootstrap** (`src/phase9_corrected/`):
   - per-replicate `SeedSequence` streams;
   - paired M1/M2 inner resamples;
   - append-only hash-linked checkpoints;
   - runtime and non-estimability stop rules;
   - adaptive B_out sizing (EP-1).
4. **Controls:**
   - code freeze (27 files) checked by the launcher before Python starts;
   - a pre-execution gate;
   - a post-run verifier;
   - a gated report tool;
   - an independent audit.

**Key scripts:**

| Script | Role |
|---|---|
| `38` | Phase 9 primary estimate (θ̂) |
| `41` | corrected driver |
| `42` | launcher (`run` / `verify`) |
| `43` | verifier |
| `45` | sizing-decision writer |
| `46` | gated interval report |
| `47` | θ̂ alignment check |
| `48` | independent audit |

**Software environment:**
- Python 3.13.13; NumPy 2.5.4; pandas 2.3.3; SciPy 1.18.1; lifelines 0.30.3; pyarrow; pytest. Versions are in `requirements*.txt`.
- macOS arm64. The environment fingerprint is in the freeze record.

**Tests:** full suite on 8 November 2026 (`.venv/bin/python -B -m pytest -q -p no:cacheprovider -rs tests/`): **1,085 passed, 0 failed, 0 skipped**, 42 warnings, in 3 min 6 s. The warnings are NumPy runtime warnings inside deliberate non-finite tests and a fork deprecation notice. Most corrected-pipeline tests use synthetic data, and passing tests establish software conformance, not scientific validity.

**Reproduction:**
- README "Reproduction" section, and `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md`.
- **From archived checkpoints:** verifier + `scripts/46` + `scripts/48` recompute every reported number in under an hour.
- **From raw GEO data:** the full campaign is about 5 runs × 5.5–7.5 h on 8 cores.

## 3. Evidence and reproducibility

| Item | Value |
|---|---|
| Frozen verification at 750 | `reports/phase9/corrected_rerun/verification/L00750.report.json` (`9558aaa5…ff25`): PASS 60, FAIL 0, NOT_APPLICABLE 0, NOT_COMPARABLE 0; H5 false |
| Freeze record | `config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json` (`618f20eb…a790`); 27/27 file hashes re-verified on 8 Nov 2026 |
| Configuration | `config/analysis/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json`, config hash `bf939a4e…6a56` |
| Input manifest | `reports/phase9/corrected_rerun/manifest/input_manifest.json` (`1ff24387…3f27`); 29/29 entries re-verified |
| Protocols | `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` v3.1; `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1 |
| Checkpoint SHA-256 | run0 `0e27494f…5dab15`, run1 `9dfe640b…af7f43`, run2 `0204da98…0e91a3`, run3 `05274f0e…fc1e70`, run4 `af753007…58e326` (full values in `docs/LARGE_ARTIFACT_ARCHIVE.md`) |
| Archive manifest | `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt` (`c45bb806…6b6f`), 52 files, all OK |
| Interval results | `reports/phase9/corrected_rerun_report/corrected_interval_results.json` (`1c504b13…eb60`) and tables |
| Independent audit | `reports/phase9/independent_audit/INDEPENDENT_AUDIT_REPORT.md`, results JSON `b67560a8…e58d` |
| θ̂ alignment | `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md` |

**What can be reproduced from the repository alone:**
- the code;
- the tests that use synthetic data;
- the documentation and decision trail.

**What needs the archived checkpoints:**
- the verification;
- the interval report;
- the independent audit.

**What needs the public GEO data:**
- the analysis bundle;
- any refit, including θ̂ and a full rerun.

## 4. Repository status

| | |
|---|---|
| Branch | `phase9-corrected-pipeline`, tracking `origin/phase9-corrected-pipeline` (`6bb3f26`). **10 local commits not pushed**: `c432eee` (kept, decision D1), the eight approved groups `7e0df40`…`97ed2a0`, and a documentation-only follow-up commit |
| Working tree | clean: no modified or untracked files after the commits |
| Committed content | all project source, documentation, configuration records, the freeze record and `reports/` outputs (487 tracked files after G8), excluding the five checkpoints |
| Large files | five checkpoints of 113 MB, excluded by `.gitignore` and pinned by the archive manifest. The archive is not yet built or uploaded |
| Local paths | `/Users/kalki/…` appears in frozen or historical JSON and in three historical script docstrings. Not sensitive. The frozen and historical files are left unchanged |
| Secrets | none found |
| Owner decisions | D1–D5 approved on 8 Nov 2026 and recorded in `docs/OWNER_DECISION_REGISTER.md`. **Still pending (separate authorisation):** push to the remote, building and publishing the checkpoint archive as a GitHub Release, choice of a code licence, and the public-release review (including per-sample tables). The repository is not publicly released |

## 5. Portfolio readiness

**What the project demonstrates:**
- **Biological data preprocessing:** GEO SOFT parsing; probe-to-gene mapping with explicit handling of ambiguous probes; representation sensitivity; outcome-blind QC of 453 arrays; endpoint forensics (U1/U2).
- **Transcriptomic pathway analysis:** an outcome-blind MSigDB Hallmark panel; Entrez-primary mapping; training-fold standardised pathway scores with coverage rules.
- **Survival analysis:** Cox proportional hazards with spline terms, cause-specific endpoint handling, and Harrell's C.
- **Machine-learning model development:** penalised (ridge) Cox, leakage-free pipelines that refit every preprocessing step per fold, and λ tuning by event-stratified CV.
- **Nested resampling:** CV inside bootstrap inside bootstrap, about 1.5 × 10⁵ model fits per run.
- **Bootstrap optimism correction**, with a nested-bootstrap interval for the corrected ΔC.
- **Calibration:** optimism-corrected calibration slopes, interpreted as overfitting.
- **Statistical uncertainty:**
  - kernel quantile Monte Carlo error and adaptive sizing;
  - between-run reproducibility;
  - bias diagnostics (z₀, percentile vs basic);
  - explicit separation of hard failures from warnings.
- **Reproducibility and audit:**
  - hash-frozen code and inputs;
  - per-stream RNG design with collision checks;
  - crash-consistent checkpoints;
  - a post-run verifier;
  - an independent re-implementation audit;
  - a documented protocol deviation and its correction;
  - 1,085 tests.

**Project description (GitHub / interviews):**

> Pre-specified survival-modelling study on GEO GSE39582 (452 stage II–III colon cancers, 136 recurrences). It asks whether four frozen Hallmark pathway scores improve recurrence discrimination beyond clinical covariates. I built a leakage-free penalised-Cox pipeline (per-fold gene and score standardisation, CV-tuned ridge) and quantified the optimism-corrected gain in Harrell's C with a nested bootstrap: five independent runs of 750 outer × 100 inner replicates, about 765,000 model fits. The pipeline is hash-frozen, uses per-replicate random streams, checkpoints crash-consistently, and was checked by a post-run verifier and an independent re-implementation audit. Result: ΔC = +0.004, with 95 % intervals of about −0.009 to +0.05 in every run. The pathways add no demonstrable discrimination, and residual overfitting and batch sensitivity are documented. The project shows how to evaluate an omics signature honestly rather than report its apparent performance.
