"""Case study: `pathway-audit` single-cell source check on GSE132465 (tumour cells of the 18 QC-passing patients).

Run from the repository root with the single-cell environment:
    .venv-sc/bin/python tools/pathway-audit/examples/case_study_single_cell.py
Compares the tool's result with the registered single-cell analysis (results/single_cell/). Outputs:
reports/pathway_audit_case_study/single_cell_GSE132465/{sc_sources.csv, SC_SOURCES.md, comparison.json}.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "tools" / "pathway-audit" / "src"))
sys.path.insert(0, str(REPO / "scripts" / "single_cell"))

import sccommon as C  # noqa: E402

from pathway_audit.single_cell import cell_type_sources, verdicts  # noqa: E402

OUT = REPO / "reports" / "pathway_audit_case_study" / "single_cell_GSE132465"


def main():
    acc = "GSE132465"
    gate = json.loads((C.RES / "qc_gate.json").read_text())[acc]
    cells = C.cell_table(acc)
    use = cells[(cells["cls"] == "tumour") & cells["patient"].isin(gate["retained"]) & cells["compartment"].isin(C.COMPARTMENTS)]
    f = next((C.EXT / acc).glob("*raw_UMI*"))
    keep = set(use["cell"])
    parts, genes, cols = [], [], None
    for chunk in pd.read_csv(f, sep="\t", index_col=0, chunksize=2000, dtype=defaultdict(lambda: np.float32, Index=str)):
        if cols is None:
            cols = [c for c in chunk.columns if c in keep]
        genes += list(chunk.index)
        parts.append(sp.csr_matrix(chunk[cols].to_numpy()))
    X = sp.vstack(parts).T.tocsr()
    lab = use.set_index("cell").loc[cols]
    sets, _ = C.gene_sets()
    table = cell_type_sources(X, genes, lab["compartment"].to_numpy(), sets, groups=lab["patient"].to_numpy(), n_random=1000, seed=C.SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT / "sc_sources.csv", index=False)
    v = verdicts(table)
    body = "\n".join(f"- **{k}:** {t}" for k, t in v.items())
    (OUT / "SC_SOURCES.md").write_text("# Single-cell source check: GSE132465 tumour cells (18 patients)\n\n" + body + "\n")
    study = pd.read_csv(C.RES / "s1_s2_s5_enrichment.csv", dtype={"patient": str})
    study = study[study["accession"] == acc].groupby(["set", "compartment"])["enrichment"].median()
    tool = table.set_index(["set", "cell_type"])["enrichment"]
    comp = [{"set": s, "compartment": c, "study_median_enrichment": float(study[(s, c)]), "tool_median_enrichment": float(tool[(s, c)])}
            for s, c in study.index if (s, c) in tool.index]
    C.write_json(OUT / "comparison.json", {"note": "Tool (pathway_audit.single_cell) vs registered study analysis; random draws differ, so values agree "
                                                  "closely but not exactly.", "rows": comp})
    for r in comp:
        if r["set"] in ("EMT", "EMT-ex6"):
            print(r)
    for k, t in v.items():
        print(f"  {k}: {t}")


if __name__ == "__main__":
    main()
