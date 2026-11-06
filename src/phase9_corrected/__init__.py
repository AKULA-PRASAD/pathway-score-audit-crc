"""Phase 9 corrected nested-bootstrap implementation (independently versioned).

A NEW namespace. Nothing here modifies, or writes beside, the historical analysis:
  * historical driver  scripts/40_phase9_dc_interval.py   (frozen, never imported for execution)
  * historical outputs reports/phase9/dc_interval/         (read-only; verification only)

Rule classes used in docstrings and comments:
  [PLAN]  requirement of docs/PHASE9_DELTAC_INTERVAL_PLAN.md v2.1 (line numbers cited)
  [IMPL]  existing (historical) implementation behaviour, reproduced deliberately
  [PROP]  proposed clarification from the DRAFT adjudication record -- NOT approved
  [OWNER] unresolved owner decision -- must be supplied explicitly in configuration

The package refuses to execute unless every [OWNER] field is supplied and marked
approved in the configuration (config.validate_for_execution).
"""

SCHEMA_VERSION = "p9c-1"
