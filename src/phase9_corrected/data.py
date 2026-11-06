"""Read-only loading of the frozen analysis inputs.

[IMPL] identical selection logic to the frozen driver's `init_worker`
(scripts/40_phase9_dc_interval.py lines 72-84); equivalence is tested.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

EV, TM = "rfs_event_full", "rfs_time_full"
AGE_COL = "age.at.diagnosis_(year)"            # src/phase9_execution.py:36
M1_VARS = [AGE_COL, "sex", "tnm.stage", "tumor.location", "chemotherapy.adjuvant"]
LAM_GRID = np.geomspace(1e-3, 10.0, 12)        # driver line 33


@dataclass
class AnalysisContext:
    ids: list
    cov: pd.DataFrame
    t: np.ndarray
    e: np.ndarray
    X: object
    panel: object

    @property
    def models(self) -> dict:
        """M1 = clinical only (no panel); M2 = clinical + frozen four-set panel."""
        return {"M1": None, "M2": self.panel}


def load_context(root: Path) -> AnalysisContext:
    analysis = Path(root) / "data" / "processed" / "analysis"
    mem = pd.read_csv(analysis / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    cov = pd.read_csv(analysis / "covariates.csv", keep_default_na=False,
                      dtype=str).set_index("sample_id")
    ids = [s for s in s3 if all(cov.loc[s, v] != "" for v in M1_VARS)]
    ep = pd.read_csv(analysis / "endpoint.csv").set_index("sample_id").reindex(ids)
    X = pd.read_parquet(analysis / "expression_median_S0.parquet")[ids].astype(np.float64)
    panel = json.loads((Path(root) / "config" / "pathways" / "m2_panel_4set.json").read_text())
    return AnalysisContext(ids, cov, ep[TM].to_numpy(float), ep[EV].to_numpy(float), X, panel)


def load_endpoint_only(root: Path):
    """ids and event vector only (enough to regenerate resamples; no expression data)."""
    analysis = Path(root) / "data" / "processed" / "analysis"
    mem = pd.read_csv(analysis / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    cov = pd.read_csv(analysis / "covariates.csv", keep_default_na=False,
                      dtype=str).set_index("sample_id")
    ids = [s for s in s3 if all(cov.loc[s, v] != "" for v in M1_VARS)]
    ep = pd.read_csv(analysis / "endpoint.csv").set_index("sample_id").reindex(ids)
    return ids, ep[EV].to_numpy(float)
