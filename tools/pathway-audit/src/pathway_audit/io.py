"""Input readers. Expression is genes (rows) × samples (columns)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def read_expression(path, log2: str = "auto") -> pd.DataFrame:
    """Read a genes × samples matrix (CSV, TSV, optionally gzipped, or parquet).

    log2: "auto" applies log2(max(x, 1)) when the 99th percentile exceeds 100 (linear-scale arrays such as MAS5/PLIER);
    "yes" always; "no" never. Duplicate gene rows are collapsed by the median.
    """
    p = Path(path)
    if p.suffix == ".parquet":
        df = pd.read_parquet(p)
    else:
        sep = "," if ".csv" in p.name else "\t"
        df = pd.read_csv(p, sep=sep, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all")
    df.index = df.index.astype(str)
    if log2 == "yes" or (log2 == "auto" and np.nanpercentile(df.to_numpy(), 99) > 100):
        df = np.log2(df.clip(lower=1))
    if df.index.duplicated().any():
        df = df.groupby(level=0).median()
    return df


def read_gene_sets(path) -> dict[str, list[str]]:
    """GMT (name, description, genes...) or JSON {name: [genes]}."""
    p = Path(path)
    if p.suffix == ".json":
        obj = json.loads(p.read_text())
        return {k: list(v) for k, v in obj.items()}
    sets = {}
    for line in p.read_text().splitlines():
        parts = [x for x in line.rstrip("\n").split("\t")]
        if len(parts) >= 3:
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def read_annotations(path) -> pd.DataFrame:
    p = Path(path)
    sep = "," if p.suffix == ".csv" else "\t"
    return pd.read_csv(p, sep=sep, index_col=0, dtype=str, keep_default_na=False)
