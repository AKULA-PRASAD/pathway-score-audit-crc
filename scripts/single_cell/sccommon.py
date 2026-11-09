"""Shared helpers for the registered single-cell extension S1-S5 (docs/SINGLE_CELL_EXTENSION_PLAN.md).

Outcome-free. Reads the public processed GEO files under data/external/single_cell/ and the frozen
Hallmark gene sets; writes compact per-patient aggregates to data/processed/single_cell/ and results
to results/single_cell/.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "data" / "external" / "single_cell"
PROC = ROOT / "data" / "processed" / "single_cell"
RES = ROOT / "results" / "single_cell"
FIG = ROOT / "figures" / "single_cell"
SEED = 20261108
COMPARTMENTS = ["epithelial", "stromal", "immune", "endothelial"]
DATASETS = ["GSE132465", "GSE144735", "GSE178341"]
ROLE = {"GSE132465": "discovery", "GSE144735": "replication 1", "GSE178341": "replication 2"}
SET_LABELS = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT", "HALLMARK_E2F_TARGETS": "E2F",
              "HALLMARK_WNT_BETA_CATENIN_SIGNALING": "Wnt", "HALLMARK_P53_PATHWAY": "p53"}
# QC positive controls (plan section 3): each must be most expressed in its own compartment
QC_MARKERS = {"epithelial": ["EPCAM", "KRT8", "KRT18"], "stromal": ["COL1A1", "DCN", "LUM"],
              "immune": ["PTPRC"], "endothelial": ["PECAM1", "VWF"]}
MIN_EPI, MIN_STROMAL = 50, 20
N_RANDOM = 1000
MAX_CELLS_S3 = 1000  # per-patient cap on malignant epithelial cells for S3 (seeded subsample; computational bound)

LEE_STROMAL = ["Stromal 1", "Stromal 2", "Stromal 3", "Myofibroblasts", "Pericytes", "Smooth muscle cells"]
LEE_ENDO = ["Tip-like ECs", "Stalk-like ECs", "Lymphatic ECs", "Proliferative ECs"]
LEE_IMMUNE = ["B cells", "T cells", "Myeloids", "Mast cells"]


def gene_sets():
    """The four frozen Hallmark sets (resolved members) and the gene-disjoint EMT-ex6 variant."""
    panel = json.loads((ROOT / "config/pathways/m2_panel_4set.json").read_text())
    sets = {SET_LABELS[k]: list(v["members_resolved"]) for k, v in panel["gene_sets"].items()}
    mcp = pd.read_csv(ROOT / "data/external/exploratory/MCPcounter_genes.txt", sep="\t")
    fib = set(mcp.loc[mcp["Cell population"] == "Fibroblasts", "HUGO symbols"])
    sets["EMT-ex6"] = [g for g in sets["EMT"] if g not in fib]
    return sets, sorted(fib)


def lee_compartment(cell_type, subtype):
    if cell_type == "Epithelial cells":
        return "epithelial"
    if cell_type in LEE_IMMUNE:
        return "immune"
    if cell_type == "Stromal cells":
        if subtype in LEE_ENDO:
            return "endothelial"
        if subtype in LEE_STROMAL:
            return "stromal"
    return "excluded"


def pelka_compartment(top, mid):
    if top == "Epi":
        return "epithelial"
    if top in ("B", "Plasma", "Myeloid", "Mast", "TNKILC"):
        return "immune"
    if top == "Strom":
        if mid == "Endo":
            return "endothelial"
        if mid in ("Fibro", "Peri", "SmoothMuscle"):
            return "stromal"
    return "excluded"


def cell_table(acc):
    """One row per cell: cell id, patient, sample class (tumour/normal/other), compartment, malignant-epithelial flag."""
    if acc in ("GSE132465", "GSE144735"):
        f = next((EXT / acc).glob("*annot*"))
        a = pd.read_csv(f, sep="\t")
        cls = a["Class"].map({"Tumor": "tumour", "Normal": "normal"}).fillna("other")
        comp = [lee_compartment(t, s) for t, s in zip(a["Cell_type"], a["Cell_subtype"])]
        mal = (cls == "tumour") & (a["Cell_type"] == "Epithelial cells") & a["Cell_subtype"].astype(str).str.startswith("CMS")
        return pd.DataFrame({"cell": a["Index"], "patient": a["Patient"], "cls": cls, "compartment": comp, "malignant_epi": mal})
    cl = pd.read_csv(next((EXT / acc).glob("*cluster*")))
    meta = pd.read_csv(next((EXT / acc).glob("*metatables*")), usecols=["cellID", "SPECIMEN_TYPE", "PID"])
    assert (cl["sampleID"].values == meta["cellID"].values).all(), "cluster and metadata tables are not aligned"
    cls = meta["SPECIMEN_TYPE"].map({"T": "tumour", "N": "normal"}).fillna("other")
    comp = [pelka_compartment(t, m) for t, m in zip(cl["clTopLevel"], cl["clMidwayPr"])]
    mal = (cls == "tumour") & (cl["clTopLevel"] == "Epi")
    return pd.DataFrame({"cell": cl["sampleID"], "patient": meta["PID"], "cls": cls, "compartment": comp, "malignant_epi": mal})


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
