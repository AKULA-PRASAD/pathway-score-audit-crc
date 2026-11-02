"""Streaming reader for GEO SOFT series files.

Replaces three near-duplicate inline parsers (scripts 07, 08 and 12). Those
scripts are left untouched; this module is the single implementation for all
new code.

Scope note (Phase 0): this module transports bytes from the SOFT file into
memory without interpretation. It applies no sample filter, no probe filter,
no value transformation and no aggregation. Parsing fidelity is the only
guarantee it makes -- it says nothing about whether the deposited values are
biologically or normalizationally fit for any particular analysis.
"""

from __future__ import annotations

import gzip
import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

_SAMPLE_MARKER = "^SAMPLE = "
_TABLE_HEADER = "ID_REF"
_TABLE_END = "!sample_table_end"


def _open_soft(path: Path):
    """Open a SOFT file, transparently handling gzip. Read-only, text mode."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


@dataclass
class SampleTable:
    """One sample's expression table, exactly as deposited."""

    sample_id: str
    probe_ids: list[str] = field(default_factory=list)
    values: list[float] = field(default_factory=list)
    non_numeric: int = 0

    def __len__(self) -> int:
        return len(self.probe_ids)


def iter_sample_tables(
    soft_path: Path,
    sample_ids: set[str] | None = None,
) -> Iterator[SampleTable]:
    """Yield one SampleTable per sample in file order.

    Parameters
    ----------
    soft_path
        Path to the SOFT file (``.gz`` supported).
    sample_ids
        Optional restriction. ``None`` (default) yields every sample in the
        file. Passing a subset is a caller-side convenience and must not be
        used to encode cohort eligibility.

    Notes
    -----
    Values that cannot be parsed as floats are counted in
    ``SampleTable.non_numeric`` rather than silently dropped, so callers can
    assert on them. The probe id is retained for every row that has at least
    two tab-separated fields, keeping probe ordering faithful to the file.
    """
    current: SampleTable | None = None
    inside_table = False

    with _open_soft(soft_path) as handle:
        for line in handle:
            if line.startswith(_SAMPLE_MARKER):
                if current is not None and len(current):
                    yield current
                sid = line.split("=", 1)[1].strip()
                wanted = sample_ids is None or sid in sample_ids
                current = SampleTable(sample_id=sid) if wanted else None
                inside_table = False
                continue

            if current is None:
                continue

            if line.startswith(_TABLE_HEADER):
                inside_table = True
                continue

            if not inside_table:
                continue

            if line.startswith(_TABLE_END):
                inside_table = False
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            current.probe_ids.append(parts[0])
            try:
                current.values.append(float(parts[1]))
            except ValueError:
                # Keep probe/value arrays aligned; record the anomaly.
                current.values.append(float("nan"))
                current.non_numeric += 1

    if current is not None and len(current):
        yield current


def list_sample_ids(soft_path: Path) -> list[str]:
    """Return every sample accession in the SOFT file, in file order."""
    ids: list[str] = []
    with _open_soft(soft_path) as handle:
        for line in handle:
            if line.startswith(_SAMPLE_MARKER):
                ids.append(line.split("=", 1)[1].strip())
    return ids


def series_field(soft_path: Path, key: str) -> list[str]:
    """Return every value of a ``!`` metadata key (e.g. ``Sample_data_processing``).

    Used to surface provenance such as the deposited normalization pipeline.
    """
    prefix = f"!{key}"
    out: list[str] = []
    with _open_soft(soft_path) as handle:
        for line in handle:
            if line.startswith(prefix):
                out.append(line.split("=", 1)[1].strip())
    return out


def sha256(path: Path, chunk_size: int = 1 << 20) -> str:
    """SHA-256 of a file, streamed. Opens read-only in binary mode."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()
