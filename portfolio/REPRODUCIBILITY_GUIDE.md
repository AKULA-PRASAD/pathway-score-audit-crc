# Reproducibility guide

Audience: a reviewer who wants to check the results without contacting me.

Each command below is marked:
- **TESTED** if it was run on 8 November 2026 with the result shown;
- **DOCUMENTED, NOT RE-RUN** if it comes from the project instructions but was not executed in that check.

## 1. Requirements
- **Platform:** macOS arm64 (the frozen environment fingerprint); Python 3.13.13.
- **Analysis packages:** numpy 2.5.4, pandas 2.3.3, scipy 1.18.1, scikit-learn 1.9.1, lifelines 0.30.3, pyarrow 26.0.0, pytest 9.2.0. These are recorded in `reports/phase9/execution/phase9_execution_results.json` (`environment`) and pinned in `requirements.txt` and `requirements-analysis.txt`.
- **Portfolio assets:** `portfolio/requirements-portfolio.txt` (matplotlib, numpy, python-pptx, pypdf).
- **System tools:** `jq` and `shasum`, required by the frozen launcher.

## 2. Data acquisition
- **Expression data.** GSE39582 is public on NCBI GEO and needs no access request. Raw and processed data are not redistributed in this repository.
- **Download scripts. DOCUMENTED, NOT RE-RUN.**
  ```bash
  .venv/bin/python scripts/01_download_geo.py               # GSE39582 family SOFT via GEOparse (network)
  .venv/bin/python scripts/09_download_gpl570_annotation.py # GPL570 annotation from NCBI GEO (network)
  ```
- **Expected inputs:**
  - `data/raw/gse39582/GSE39582_family.soft.gz`
  - `data/raw/gpl570/GPL570_platform.soft`
- **Validation. TESTED:** both files match the recorded checksums.
  ```bash
  cd data/raw && shasum -a 256 -c CHECKSUMS.sha256
  # GSE39582_family.soft.gz: OK
  # GPL570_platform.soft: OK
  ```
- **Gene sets.** MSigDB Hallmark gene sets are used under MSigDB's terms. The resolved memberships are recorded in `config/pathways/m2_panel_4set.json` and `m2_panel_5set.json`.

## 3. Building the analysis inputs (DOCUMENTED, NOT RE-RUN)
Build the caches and the analysis bundle with scripts 13 and 16–23, in order (see the repository `README.md`, "Reproduction"). The bundle manifests pin every input hash.

## 4. Test suite (TESTED)
```bash
.venv/bin/python -m pytest -q -p no:cacheprovider
```
- **Result:** 1,097 passed (about 3.3 minutes).
- Some tests read the locally built inputs and frozen outputs, so they need steps 2–3 or the existing local files.

## 5. Re-checking the headline result from the archived checkpoints (TESTED)
- **Restoring checkpoints.** The five checkpoint files are in the GitHub Release asset `phase9_campaign_artifacts.tar` (579,716,096 bytes). The restore commands are in `docs/LARGE_ARTIFACT_ARCHIVE.md`.
- **Checking the archive manifest:**
  ```bash
  shasum -a 256 -c reports/phase9/ARCHIVE_MANIFEST_SHA256.txt   # 52 files: OK
  ```
- **Running the independent audit**, which uses the standard library only and is portable:
  ```bash
  python3 -I -B scripts/48_phase9_independent_audit.py /tmp/audit_check
  ```
  - **Expected output:** `independent_audit_results.json` with `"classification": "INDEPENDENTLY REPRODUCED"` and `"max_abs_diff": 2.220446049250313e-16`.
  - **Runtime:** about 7 seconds.
  - This recomputes every reported interval quantity from the raw checkpoint bytes.

## 6. Re-running the full campaign (DOCUMENTED, NOT RE-RUN)
- **Instructions:** `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md` and the repository `README.md`.
- **Runtime:** each run takes 3–5 hours at 500 outer resamples on 8 cores, plus about 2 hours for the top-up to 750.
- **Portability:** the frozen launcher and verifier check the interpreter path and virtual-environment prefix, so they reproduce only in the original environment. The independent audit (step 5) is the portable check.

## 7. Portfolio materials (TESTED)
```bash
VENV=../portfolio-venv   # any location outside the repository
python3 -m venv $VENV && $VENV/bin/pip install -r portfolio/requirements-portfolio.txt
$VENV/bin/python portfolio/src/make_diagrams.py   # 10 files in portfolio/diagrams (byte-identical on rebuild)
$VENV/bin/python portfolio/src/make_slide_visuals.py  # 34 files in portfolio/diagrams/slides (byte-identical on rebuild)
$VENV/bin/python portfolio/src/make_deck.py       # PPTX, PDF, SPEAKER_NOTES.md
$VENV/bin/python portfolio/src/qc_check.py        # writes PORTFOLIO_QC_REPORT.md; exits non-zero on any FAIL
$VENV/bin/python portfolio/src/make_inventory.py  # writes ASSET_INVENTORY.md
```
Tested from a fresh environment built from `portfolio/requirements-portfolio.txt`: the diagrams rebuilt byte-identical and QC passed. The slide visuals (`make_slide_visuals.py`, added with the 20-slide presentation) were rebuilt twice in the existing portfolio environment and were byte-identical; the PDF also rebuilds byte-identical, while the PPTX differs only in its package metadata. They have not yet been rebuilt from a fresh environment. None of these fits a model; they read frozen outputs only.

## 8. Expected key outputs

| Output | Path | Key value |
|---|---|---|
| Corrected intervals | `reports/phase9/corrected_rerun_report/corrected_interval_results.json` | θ̂ = 0.004316342804030704; five percentile intervals |
| Verification | `reports/phase9/corrected_rerun/verification/L00750.report.json` | PASS, 60/60 |
| Primary execution | `reports/phase9/execution/phase9_execution_results.json` | corrected C 0.6016 / 0.6059 |
| Figures | `reports/figures/final/` | 11 figures, captions, values, sources |
