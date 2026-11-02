# Phase 0 — Infrastructure

**Date:** 2026-11-02
**Scope:** infrastructure only. No scientific decision was made or recorded.

## What Phase 0 establishes

A fast, verified, reproducible foundation:

- one shared SOFT parser, replacing three near-duplicate inline copies
- a Parquet cache of the probe matrix (~60 s re-parse → ~2 s load)
- CWD-independent path handling
- SHA-256 integrity tracking for originals and raw data
- a 21-test verification suite

## What Phase 0 explicitly does NOT establish

Phase 0 verifies **computational fidelity**: the cache is bit-identical to the
deposited values, and the new code reproduces the original scripts' results
exactly.

It says nothing about **biological or normalization validity**. GSE39582 was
deposited as:

```
RMA normalisation (Bioconductor R affy package) + ComBat (R sva package)
```

applied jointly across all 585 arrays. Whether those values support
cross-cohort transport is a Phase 4 question. **A pipeline can pass every
Phase 0 test and still be scientifically invalid.** The two properties are
tracked separately throughout.

No endpoint definition, cohort eligibility rule, expression filter, probe
aggregation rule or normalization strategy exists anywhere in this phase.
`config/analysis_config.yaml` documents their intentional absence.

## Files created (18)

### Library — `src/`
| File | Purpose |
|---|---|
| `__init__.py` | Package marker |
| `paths.py` | Canonical paths anchored to `__file__`, not CWD |
| `soft_parser.py` | Single streaming SOFT reader |
| `annotation.py` | GPL570 reader; resolves the script 11/12 header inconsistency |
| `cache.py` | Build/load/validate the probe matrix cache |
| `reproducibility.py` | Seeding and environment capture |

### Configuration
| File | Purpose |
|---|---|
| `config/analysis_config.yaml` | Infrastructure parameters only |
| `pytest.ini` | Test discovery and import path |
| `.gitignore` | Excludes venv, raw data, caches, credentials |
| `requirements-analysis.txt` | Added packages; `requirements.txt` untouched |

### Entry points — `scripts/`
| File | Purpose |
|---|---|
| `13_build_probe_cache.py` | Build cache + manifest |
| `14_verify_cache.py` | Independent verification across all 585 samples |

### Tests
| File | Purpose |
|---|---|
| `tests/__init__.py` | Package marker |
| `tests/test_phase0.py` | 21 tests |
| `tests/baseline_checksums.json` | SHA-256 of all 19 pre-Phase-0 files |

### Integrity / docs
| File | Purpose |
|---|---|
| `data/raw/CHECKSUMS.sha256` | SHA-256 of both raw downloads |
| `docs/PHASE0_CHANGELOG.md` | This file |

### Generated (gitignored)
| File | Size |
|---|---|
| `data/processed/probe_matrix.parquet` | 54,675 × 585, float64 |
| `data/processed/probe_matrix_manifest.json` | provenance record |

## Files modified

**None.** Scripts 01–12, `requirements.txt`, `metadata/`, `reports/`,
`results/` and `data/raw/` are byte-identical to their pre-Phase-0 state,
with original permissions preserved. Verified by
`test_originals_unmodified`, `test_original_permissions_unchanged` and
`test_raw_data_unmodified`.

## Findings carried forward

1. **Probe-gene mapping is sound.** Script 12's median ρ = 0.33 is driven by
   unexpressed, low-variance probes; restricted to expressed and variable
   probes, median ρ = 0.82. The corrective action is an expression filter
   (Phase 2), not distrust of GPL570 annotation. No filter applied in Phase 0.

2. **Discovery/validation are not independent.** Joint RMA + ComBat, plus
   markedly different RNA extraction method distributions between subsets.
   Phase 4 decides the normalization strategy.

3. **Endpoint ambiguities remain open.** `rfs.delay` carries no unit label in
   GEO (unlike `os.delay (months)`); 8 stage II–III samples have
   `rfs.delay = 0`, some with `rfs.event = 1`. Phase 1 resolves these.
