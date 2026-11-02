"""Reader for the GPL570 (Affymetrix HG-U133 Plus 2.0) platform annotation.

Resolves an inconsistency found during the Phase 0 audit: script 11 locates
the annotation table by matching ``ID\\t`` while script 12 matches
``ID\\tGB_ACC``. Both happen to find the same row in the current file, but the
two strategies are not equivalent in general. This module uses the stricter
match, verifies it against the looser one, and raises if they ever diverge.

Scope note (Phase 0): this module parses annotation only. It performs NO
probe filtering, NO ambiguity resolution policy and NO probe-to-gene
collapsing -- those are Phase 2 scientific decisions. ``split_symbols`` is
offered as a pure string utility; deciding what to do with multi-symbol
probes is explicitly not made here.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

_STRICT_PREFIX = "ID\tGB_ACC"
_LOOSE_PREFIX = "ID\t"

# Placeholders Affymetrix uses for "no annotation".
NULL_SYMBOLS = {"", "---", "NA", "N/A"}


def find_table_header(path: Path) -> int:
    """Return the 0-based line index of the annotation table header.

    Raises
    ------
    ValueError
        If the strict and loose header-detection strategies disagree, or if
        no header can be found at all.
    """
    strict_idx = loose_idx = None
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for i, line in enumerate(handle):
            if loose_idx is None and line.startswith(_LOOSE_PREFIX):
                loose_idx = i
            if strict_idx is None and line.startswith(_STRICT_PREFIX):
                strict_idx = i
            if strict_idx is not None and loose_idx is not None:
                break

    if strict_idx is None:
        raise ValueError(f"GPL570 annotation table header not found in {path}")
    if strict_idx != loose_idx:
        raise ValueError(
            "Ambiguous annotation header: strict match 'ID\\tGB_ACC' at line "
            f"{strict_idx} but loose match 'ID\\t' at line {loose_idx}. "
            "Resolve before trusting the annotation table."
        )
    return strict_idx


def load_annotation(path: Path | None = None) -> pd.DataFrame:
    """Load the GPL570 annotation table verbatim as strings.

    Returns every annotated probe with no filtering applied.
    """
    from src import paths

    path = Path(path) if path is not None else paths.GPL570_ANNOTATION
    header_index = find_table_header(path)

    table = pd.read_csv(
        path,
        sep="\t",
        skiprows=header_index,
        dtype=str,
        low_memory=False,
        on_bad_lines="warn",
    )

    if "ID" not in table.columns:
        raise ValueError("GPL570 annotation table has no 'ID' column.")

    # Drop the trailing SOFT control records (e.g. '!platform_table_end').
    table = table[
        table["ID"].notna() & ~table["ID"].str.startswith("!")
    ].copy()
    table["ID"] = table["ID"].str.strip()

    return table.reset_index(drop=True)


def split_symbols(value: object) -> list[str]:
    """Split an Affymetrix ``///``-delimited Gene Symbol field.

    Pure string utility: returns the de-duplicated symbol list with
    placeholders removed. Makes no decision about multi-symbol probes.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    parts = [p.strip() for p in str(value).split("///")]
    parts = [p for p in parts if p not in NULL_SYMBOLS]
    return list(dict.fromkeys(parts))
