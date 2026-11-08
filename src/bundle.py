"""Analysis-bundle assembly.

Storage design: S3 subset S1 subset S0 and S3 subset S2 subset S0, so the
bundle stores S0-scoped expression matrices plus a membership table. Each
cohort materialises as an exact sample-id subset, verified against the Phase 1
manifest. One superset cannot drift the way four near-duplicate copies can.

Alignment is ALWAYS by sample identifier. No positional indexing anywhere.

Phase 3 fits NOTHING: no scaler, no filter, no correction, no batch
adjustment. There is therefore no fitted object that could carry information
between partitions.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from src import paths
from src.provenance import (PARTITION_INDEPENDENCE, PROCESSING_PROVENANCE)
from src.reproducibility import environment_report

ANALYSIS_DIR = paths.PROCESSED_DIR / "analysis"

CARRIED_COVARIATES = [
    "age.at.diagnosis_(year)", "sex", "tnm.stage", "tnm.t", "tnm.n", "tnm.m",
    "tumor.location", "chemotherapy.adjuvant", "mmr.status", "cimp.status",
    "cin.status", "kras.mutation", "braf.mutation", "tp53.mutation",
]
EXCLUDED_COVARIATES = {
    "cit.molecularsubtype":
        "Derived FROM the expression data. Including it in a clinical "
        "baseline would absorb the signal being measured against.",
    "rna.extraction.method":
        "Technical variable. A batch diagnostic, never a predictor.",
    "chemotherapy.adjuvant.type":
        "69% missing, 5 levels on 51 events.",
    "dependancy_sample":
        "Tumour/normal pairing key; not a clinical covariate.",
}


def sha256_file(path: Path) -> str:
    d = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            d.update(c)
    return d.hexdigest()


def align_by_id(matrix: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """Label-based column alignment with hard assertions. Never positional."""
    missing = [s for s in sample_ids if s not in matrix.columns]
    if missing:
        raise KeyError(f"{len(missing)} sample ids absent from the matrix: "
                       f"{missing[:5]}")
    out = matrix.reindex(columns=sample_ids)
    if list(out.columns) != list(sample_ids):
        raise AssertionError("column alignment did not preserve the requested order")
    if out.isna().to_numpy().any():
        raise AssertionError("reindex introduced NaN")
    return out


def build_covariates(meta: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    idx = meta.set_index("sample_id").loc[sample_ids]
    out = idx[CARRIED_COVARIATES].copy()
    out.index.name = "sample_id"
    # CIT446 stage/M conflict: flagged, never recoded.
    out["flag_stage_m_conflict"] = (
        (idx["tnm.m"] == "M1") & (idx["tnm.stage"] != 4)).to_numpy()
    return out.reset_index()


def build_partition_and_batch(meta: pd.DataFrame, batch: pd.DataFrame | None,
                              sample_ids: list[str]) -> pd.DataFrame:
    idx = meta.set_index("sample_id").loc[sample_ids]
    out = pd.DataFrame({"sample_id": sample_ids,
                        "dataset": idx["dataset"].to_numpy()})
    if batch is not None:
        b = batch.set_index("sample_id")
        out["combat_batch"] = b["combat_batch"].reindex(sample_ids).to_numpy()
        out["combat_mod"] = b["combat_mod"].reindex(sample_ids).to_numpy()
    else:
        out["combat_batch"] = None
        out["combat_mod"] = None
    return out


def batch_composition(part: pd.DataFrame,
                      membership: pd.DataFrame, cohort: str) -> dict:
    ids = set(membership.loc[membership[f"in_{cohort}"].astype(bool),
                             "sample_id"])
    p = part[part["sample_id"].isin(ids)]
    ct = pd.crosstab(p["combat_batch"], p["dataset"])
    for col in ("discovery", "validation"):
        if col not in ct.columns:
            ct[col] = 0
    return {str(k): {"discovery": int(v["discovery"]),
                     "validation": int(v["validation"])}
            for k, v in ct.iterrows()}


def write_manifest(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n")


def base_manifest(inputs: dict[str, Path]) -> dict:
    return {
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "phase": 3,
        "input_hashes": {k: sha256_file(v) for k, v in inputs.items()},
        "input_sizes_bytes": {k: v.stat().st_size for k, v in inputs.items()},
        "processing_provenance": PROCESSING_PROVENANCE,
        "partition_independence": PARTITION_INDEPENDENCE,
        "normalization_applied": "none",
        "batch_correction_applied": "none",
        "batch_covariate_auto_included": False,
        "transformation_fitted_in_phase3": "none",
        "outcome_based_gene_selection": False,
        "differential_expression_filtering": False,
        "environment": environment_report(),
    }
