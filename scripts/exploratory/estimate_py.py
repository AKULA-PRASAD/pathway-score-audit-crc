"""Python re-implementation of ESTIMATE v1.0.13 (estimateScore.R + filterCommonGenes.R), read from the package source.

Validated against the package's bundled example (see x5_composition.py, validation gate).
"""
import numpy as np
import pandas as pd

import xcommon as X

PKG = X.EXT / "estimate_pkg" / "estimate" / "inst" / "extdata"


def gene_sets():
    gs = {}
    for line in (PKG / "SI_geneset.gmt").read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        gs[parts[0]] = [g for g in parts[2:] if g]
    return gs


def common_genes():
    df = pd.read_csv(PKG / "common_genes.txt", sep="\t")
    return df


def filter_common(expr: pd.DataFrame) -> pd.DataFrame:
    """filterCommonGenes(id='GeneSymbol'): inner join of expression rows with the common-gene list."""
    cg = common_genes()
    keep = [g for g in cg["GeneSymbol"] if g in expr.index]
    return expr.loc[sorted(set(keep))]


def estimate_score(expr: pd.DataFrame, affymetrix=True) -> pd.DataFrame:
    """estimateScore: per-sample rank normalisation, ssGSEA-style ES (weight 0.25, sum of RES), Affymetrix purity."""
    m = expr.to_numpy(float)
    genes = np.array(expr.index)
    Ng = m.shape[0]
    ranks = np.apply_along_axis(lambda c: pd.Series(c).rank(method="average").to_numpy(), 0, m)
    m = 10000 * ranks / Ng
    out = {}
    for name, gs in gene_sets().items():
        overlap = [g for g in gs if g in set(genes)]
        tagset = np.isin(genes, overlap)
        es = np.empty(m.shape[1])
        for j in range(m.shape[1]):
            order = np.argsort(-m[:, j], kind="stable")     # R order(decreasing=TRUE): ties keep original order
            tag = tagset[order]
            corr = np.abs(m[order, j]) ** 0.25
            nh = tag.sum()
            nm = Ng - nh
            Fn = np.cumsum(np.where(tag, corr, 0.0)) / corr[tag].sum()
            F0 = np.cumsum(~tag) / nm
            es[j] = np.sum(Fn - F0)
        out[name] = es
        out[f"_overlap_{name}"] = len(overlap)
    res = pd.DataFrame({"StromalScore": out["StromalSignature"], "ImmuneScore": out["ImmuneSignature"]}, index=expr.columns)
    res["ESTIMATEScore"] = res.StromalScore + res.ImmuneScore
    if affymetrix:
        pur = np.cos(0.6049872018 + 0.0001467884 * res.ESTIMATEScore)
        res["TumorPurity"] = np.where(pur >= 0, pur, np.nan)
    res.attrs["overlap"] = {k[9:]: v for k, v in out.items() if k.startswith("_overlap_")}
    return res
