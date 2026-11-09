"""Single-cell step 0: stream each dataset once into compact per-patient aggregates (outcome-free).

For every (sample class, patient, compartment) group it stores summed UMI counts per gene and the
number of cells; it also stores per-gene detection counts in tumour cells, and a seeded subsample of
malignant epithelial cells (raw counts and library sizes) for the S3 coherence analysis.
Matrices are processed in blocks of cells so the 370,000-cell atlas never has to fit in memory.

Usage: .venv-sc/bin/python scripts/single_cell/s0_build.py [ACCESSION ...]
"""
import sys
from collections import defaultdict

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import sccommon as C

BLOCK = 20000


def lee_matrix(acc):
    """Dense text (genes x cells) -> cells x genes CSR, read in chunks of genes."""
    f = next((C.EXT / acc).glob("*raw_UMI*"))
    parts, genes, cells = [], [], None
    for chunk in pd.read_csv(f, sep="\t", index_col=0, chunksize=1000, dtype=defaultdict(lambda: np.float32, Index=str), engine="c"):
        if cells is None:
            cells = list(chunk.columns)
        genes += list(chunk.index)
        parts.append(sp.csr_matrix(chunk.to_numpy()))
        print(f"  {acc}: {len(genes)} genes read", flush=True)
    M = sp.vstack(parts).T.tocsr()
    return M, np.array(genes), np.array(cells)


def blocks_from_csr(M):
    for a in range(0, M.shape[0], BLOCK):
        yield a, M[a:a + BLOCK]


def pelka_blocks(f, n_genes):
    indptr = f["matrix/indptr"][:]
    n = len(indptr) - 1
    for a in range(0, n, BLOCK):
        b = min(a + BLOCK, n)
        lo, hi = indptr[a], indptr[b]
        data = f["matrix/data"][lo:hi].astype(np.float32)
        idx = f["matrix/indices"][lo:hi]
        yield a, sp.csr_matrix((data, idx, indptr[a:b + 1] - lo), shape=(b - a, n_genes))
        print(f"  cells {b}/{n}", flush=True)


def build(acc):
    print(acc, flush=True)
    cells = C.cell_table(acc)
    if acc == "GSE178341":
        f = h5py.File(next((C.EXT / acc).glob("*.h5")), "r")
        genes = np.array([g.decode() for g in f["matrix/features/name"][:]])
        bc = np.array([b.decode() for b in f["matrix/barcodes"][:]])
        n_genes = len(genes)
        source = pelka_blocks(f, n_genes)
    else:
        M, genes, bc = lee_matrix(acc)
        n_genes = len(genes)
        source = blocks_from_csr(M)
    # duplicate gene symbols are kept here and summed when the aggregates are loaded (s1_tests.py)
    tab = cells.set_index("cell").reindex(bc)
    known = tab["patient"].notna().to_numpy()
    print(f"  matrix cells {len(bc)}, annotated {known.sum()}, annotation-only {len(cells) - known.sum()}", flush=True)
    tab = tab.reset_index().rename(columns={"index": "cell"})
    tab.loc[~known, ["cls", "compartment"]] = ["unannotated", "excluded"]
    tab["malignant_epi"] = tab["malignant_epi"].fillna(False).astype(bool)
    keep_cls = tab["cls"].isin(["tumour", "normal"]).to_numpy()
    grp = tab["cls"].astype(str) + "|" + tab["patient"].astype(str) + "|" + tab["compartment"].astype(str)
    groups = sorted(grp[keep_cls].unique())
    gid = pd.Series(range(len(groups)), index=groups)
    cell_gid = grp.map(gid).fillna(-1).astype(int).to_numpy()
    # seeded S3 subsample of malignant epithelial cells, per patient
    rng = np.random.default_rng(C.SEED)
    sel = []
    for _, sub in tab[tab["malignant_epi"]].groupby("patient", sort=True):
        idx = sub.index.to_numpy()
        sel.append(np.sort(rng.choice(idx, C.MAX_CELLS_S3, replace=False)) if len(idx) > C.MAX_CELLS_S3 else idx)
    sel = np.sort(np.concatenate(sel)) if sel else np.array([], int)
    sel_set = np.zeros(len(tab), bool)
    sel_set[sel] = True
    sums = np.zeros((len(groups), n_genes))
    det = np.zeros(n_genes)
    epi_rows, epi_lib = [], []
    tumour = (tab["cls"] == "tumour").to_numpy()
    for a, B in source:
        b = a + B.shape[0]
        g = cell_gid[a:b]
        ok = g >= 0
        if ok.any():
            ind = sp.csr_matrix((np.ones(ok.sum()), (g[ok], np.nonzero(ok)[0])), shape=(len(groups), B.shape[0]))
            sums += (ind @ B).toarray()
        t = tumour[a:b]
        if t.any():
            det += np.asarray((B[t] > 0).sum(axis=0)).ravel()
        s = sel_set[a:b]
        if s.any():
            epi_rows.append(B[s])
            epi_lib.append(np.asarray(B[s].sum(axis=1)).ravel())
    out = C.PROC / acc
    out.mkdir(parents=True, exist_ok=True)
    gtab = pd.DataFrame([g.split("|") for g in groups], columns=["cls", "patient", "compartment"])
    gtab["n_cells"] = pd.Series(cell_gid[cell_gid >= 0]).value_counts().reindex(range(len(groups))).fillna(0).astype(int).to_numpy()
    gtab.to_csv(out / "groups.csv", index=False)
    np.save(out / "group_sums.npy", sums)
    np.save(out / "tumour_detected_cells.npy", det)
    pd.Series(genes).to_csv(out / "genes.csv", index=False, header=["gene"])
    E = sp.vstack(epi_rows).tocsr() if epi_rows else sp.csr_matrix((0, n_genes))
    sp.save_npz(out / "s3_malignant_epi_counts.npz", E)
    pd.DataFrame({"cell": tab["cell"].to_numpy()[sel], "patient": tab["patient"].to_numpy()[sel],
                  "library_size": np.concatenate(epi_lib) if epi_lib else []}).to_csv(out / "s3_malignant_epi_cells.csv", index=False)
    C.write_json(out / "build_summary.json", {"accession": acc, "matrix_cells": int(len(bc)), "annotated_cells": int(known.sum()),
                                               "n_tumour_cells": int(tumour.sum()), "n_genes": int(n_genes), "n_groups": len(groups),
                                               "s3_cells": int(len(sel)), "max_cells_s3": C.MAX_CELLS_S3})
    print(f"  done: {len(groups)} groups, {len(sel)} S3 cells", flush=True)


if __name__ == "__main__":
    for acc in sys.argv[1:] or C.DATASETS:
        build(acc)
