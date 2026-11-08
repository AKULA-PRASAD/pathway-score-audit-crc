"""Probe-to-gene mapping classification for GPL570.

Classifies every probe on the array into exactly one of three disjoint
classes, so that nothing is silently discarded:

    unmapped    -- the annotation gives no gene symbol
    unique      -- exactly one gene symbol
    ambiguous   -- two or more symbols (Affymetrix '///' multi-mapping)

Gene-level aggregation uses UNIQUE probes only. Ambiguous and unmapped probes
cannot be assigned to a single gene without an arbitrary choice, so they are
excluded from gene-level matrices BY CONSTRUCTION -- but they are counted,
listed, and left untouched in the probe-level cache, which remains the
canonical representation.

This module makes no use of outcome data of any kind.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from src.annotation import load_annotation, split_symbols

UNMAPPED = "unmapped"
UNIQUE = "unique"
AMBIGUOUS = "ambiguous"


@dataclass
class MappingTable:
    """Full per-probe classification. `table` covers every probe on the array."""

    table: pd.DataFrame          # ID, n_symbols, mapping_class, gene_symbol
    n_probes: int

    def probes(self, cls: str) -> pd.Index:
        return pd.Index(self.table.loc[self.table["mapping_class"] == cls, "ID"])

    @property
    def unique_map(self) -> pd.Series:
        """probe_id -> gene_symbol, for uniquely mapped probes only."""
        u = self.table[self.table["mapping_class"] == UNIQUE]
        return pd.Series(u["gene_symbol"].values, index=u["ID"].values,
                         name="gene_symbol")

    def summary(self) -> dict:
        counts = self.table["mapping_class"].value_counts().to_dict()
        um = self.unique_map
        per_gene = um.groupby(um.values).size()
        return {
            "n_probes_total": self.n_probes,
            "n_unmapped": int(counts.get(UNMAPPED, 0)),
            "n_unique": int(counts.get(UNIQUE, 0)),
            "n_ambiguous": int(counts.get(AMBIGUOUS, 0)),
            "n_genes_unique": int(per_gene.size),
            "n_genes_single_probe": int((per_gene == 1).sum()),
            "n_genes_multi_probe": int((per_gene > 1).sum()),
            "n_probes_in_multi_probe_genes": int(per_gene[per_gene > 1].sum()),
            "max_probes_per_gene": int(per_gene.max()),
        }


def build_mapping(annotation: pd.DataFrame | None = None) -> MappingTable:
    """Classify every probe. Pure annotation work; no expression, no outcome."""
    ann = load_annotation() if annotation is None else annotation
    syms = ann["Gene Symbol"].apply(split_symbols)
    n = syms.apply(len)

    cls = pd.Series(AMBIGUOUS, index=ann.index, dtype=object)
    cls[n == 0] = UNMAPPED
    cls[n == 1] = UNIQUE

    table = pd.DataFrame({
        "ID": ann["ID"].values,
        "n_symbols": n.values,
        "mapping_class": cls.values,
        "gene_symbol": [s[0] if len(s) == 1 else None for s in syms],
        "all_symbols": ["///".join(s) if s else None for s in syms],
    })
    return MappingTable(table=table, n_probes=len(table))


def reconcile(mapping: MappingTable, expression_probes: pd.Index) -> dict:
    """Check the annotation's probe set against the expression matrix's."""
    ann_ids = pd.Index(mapping.table["ID"])
    return {
        "annotation_probes": len(ann_ids),
        "expression_probes": len(expression_probes),
        "in_both": len(ann_ids.intersection(expression_probes)),
        "annotation_only": sorted(ann_ids.difference(expression_probes))[:20],
        "expression_only": sorted(expression_probes.difference(ann_ids))[:20],
        "identical_sets": bool(
            set(ann_ids) == set(expression_probes)),
    }
