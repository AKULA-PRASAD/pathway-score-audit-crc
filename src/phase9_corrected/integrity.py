"""H4 frozen-input integrity and historical-artifact protection; input manifest.

H4 (plan line 457) is a hard check of its own, separate from V1-V5.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

from .config import canonical_json

# [IMPL] the frozen driver's FROZEN set (scripts/40_phase9_dc_interval.py lines 38-51) ...
FROZEN_INPUTS = {
    "docs/PHASE9_MODELING_PROTOCOL_DRAFT.md": "ab56708919f3746f",
    "config/pathways/m2_panel_4set.json": "43ce15b2f0e06f75",
    "config/pathways/M2_FREEZE_RECORD.json": "64a86fa8ed0d6a8a",
    "config/analysis/A19_FREEZE_RECORD.json": "a77de7b37fe09eeb",
    "src/pathway_scoring.py": "38cac035c0b2d450",
    "src/score_standardisation.py": "41b5c091125e1b11",
    "src/phase9_execution.py": "571d0c233f367be6",
    "data/processed/analysis/expression_median_S0.parquet": "a1ca41807c534584",
    "data/processed/analysis/covariates.csv": "47b2d9829a4fccec",
    "data/processed/analysis/endpoint.csv": "802a68575a34a3b6",
    "reports/phase9/execution/phase9_execution_results.json": "05d13b8b339fa7d4",
    "reports/phase9/execution_s7/phase9_s7_results.json": "b9202305bb6416e1",
    # ... plus cohort_membership.csv, which the driver READS (line 73) but never hashed.
    "data/processed/analysis/cohort_membership.csv": "bd3671f8df64c594",
}

# Historical artifacts that must remain byte-identical (preservation, not inputs).
PROTECTED_HISTORICAL = {
    "scripts/40_phase9_dc_interval.py": "4b1738a4d0ee0e39",
    "docs/PHASE9_DELTAC_INTERVAL_PLAN.md": "d27c52a29930a338",
    "scripts/38_phase9_execute.py": "82afbcf05f3b7514",
    "reports/phase9/dc_interval/run0_outer.jsonl": "885bbd1bd891b570",
    "reports/phase9/dc_interval/run1_outer.jsonl": "1bb27c56612a5f9d",
    "reports/phase9/dc_interval/run2_outer.jsonl": "0f12dce12bfce95e",
    "reports/phase9/dc_interval/run3_outer.jsonl": "4115bd4cd014d345",
    "reports/phase9/dc_interval/run4_outer.jsonl": "2202d559b0418392",
    "reports/phase9/dc_interval/HALT_RECORD.json": "c1b404309bd63744",
    "reports/phase9/dc_interval/RESUME_AUTHORIZATION.json": "0571c10e02e9059d",
    "reports/phase9/h5_spotcheck/h5_comparison.json": "35e2ed1e8c271da5",
    "reports/phase9/final_project/results/s1_interval_results.json": "99ad5e20c1485425",
    "reports/phase9/final_project/documentation/CV_SEED_DEVIATION_DRAFT_v2.md": "1e4200c100358161",
    "reports/phase9/final_project/documentation/METHODOLOGICAL_DECISION_REVIEW.md": "ade6627a786c113f",
    "reports/phase9/final_project/documentation/PRE_AUTHORIZATION_METHODOLOGY_AUDIT.md": "0d15990ad82e2006",
    "reports/phase9/final_project/documentation/PHASE9_METHODOLOGICAL_CLARIFICATION_AND_OWNER_ADJUDICATION_DRAFT.md": "fe81d8deea0400cd",
    "reports/phase9/final_project/documentation/STATISTICAL_CONVENTIONS.md": "fe0845dbce0c563c",
    "reports/phase9/final_project/audit/EXECUTION_LEDGER.md": "03aad17587478eb5",
}


class HardFailureH4(SystemExit):
    pass


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_hashes(registry: dict, root: Path) -> dict:
    rows = {}
    for rel, exp in registry.items():
        p = Path(root) / rel
        full = sha256_file(p) if p.exists() else None
        rows[rel] = {"expected_prefix": exp, "sha256": full,
                     "ok": full is not None and full.startswith(exp)}
    return {"files": rows, "all_ok": all(r["ok"] for r in rows.values())}


def h4_check(root: Path, stage: str, raise_on_fail=True) -> dict:
    res = check_hashes(FROZEN_INPUTS, root)
    res["check"], res["stage"] = "H4", stage
    res["status"] = "PASS" if res["all_ok"] else "FAIL"
    if raise_on_fail and not res["all_ok"]:
        bad = [k for k, v in res["files"].items() if not v["ok"]]
        raise HardFailureH4(f"HARD FAILURE H4 ({stage}) -- frozen input changed: {bad}")
    return res


def historical_preservation_check(root: Path) -> dict:
    res = check_hashes(PROTECTED_HISTORICAL, root)
    res["check"] = "historical_preservation"
    res["status"] = "PASS" if res["all_ok"] else "FAIL"
    return res


def code_hashes(root: Path) -> dict:
    root = Path(root)
    files = sorted((root / "src" / "phase9_corrected").glob("*.py"))
    files.append(root / "scripts" / "41_phase9_corrected_driver.py")
    return {str(f.relative_to(root)): sha256_file(f) for f in files if f.exists()}


# Decision 13 (final remediation): the manifest's BLOCKING platform identity is the OS family
# and architecture; the full platform.platform() string (which embeds the macOS point
# version) is RECORD-ONLY metadata under MANIFEST_RECORD_ONLY_KEY. Every other field is
# unchanged and remains blocking.
MANIFEST_RECORD_ONLY_KEY = "record_only"


def manifest_platform_identity() -> dict:
    return {"system": platform.system(), "machine": platform.machine(),
            "architecture": platform.architecture()[0]}


def build_manifest(config, root: Path) -> dict:
    import numpy, pandas, scipy, lifelines
    return {"schema": "p9c-manifest-2",
            "config_sha256": config.sha256(),
            "inputs": {k: v["sha256"] for k, v in check_hashes(FROZEN_INPUTS, root)["files"].items()},
            "code": code_hashes(root),
            "environment": {"python": sys.version.split()[0], "platform": manifest_platform_identity(),
                            "numpy": numpy.__version__, "pandas": pandas.__version__,
                            "scipy": scipy.__version__, "lifelines": lifelines.__version__},
            MANIFEST_RECORD_ONLY_KEY: {"platform_string": platform.platform()}}


def _manifest_blocking(m: dict) -> dict:
    return {k: v for k, v in m.items() if k != MANIFEST_RECORD_ONLY_KEY}


def write_manifest_once(path: Path, manifest: dict) -> str:
    """Write the manifest if absent; if present, every BLOCKING field must be identical
    (only the `record_only` section may differ) and the stored file must be the canonical
    serialisation of what it contains. Returns the sha256 of the STORED manifest bytes."""
    path = Path(path)
    data = canonical_json(manifest)
    if path.exists():
        stored = path.read_bytes()
        try:
            old = json.loads(stored)
        except ValueError:
            old = None
        if not isinstance(old, dict) or canonical_json(old) != stored or \
                canonical_json(_manifest_blocking(old)) != canonical_json(_manifest_blocking(manifest)):
            raise HardFailureH4(f"{path} exists with different content: inputs, code, config or "
                                "environment changed since the run started")
        return hashlib.sha256(stored).hexdigest()
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_bytes(data)
        tmp.replace(path)
    return hashlib.sha256(data).hexdigest()


def load_manifest(path: Path) -> dict:
    return json.loads(Path(path).read_text())
