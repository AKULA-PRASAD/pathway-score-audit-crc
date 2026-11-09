"""Computational tissue-composition estimates (not measurements).

References are NOT bundled (third-party licences: ESTIMATE GPL-2, MCP-counter GPL-3). `fetch_references` downloads them
from their official sources into a local cache and records SHA-256 checksums.

* ESTIMATE (Yoshihara et al., Nat Commun 2013;4:2612): stromal and immune enrichment scores, re-implemented from the
  v1.0.13 package source (rank normalisation, weight 0.25, sum of the running-sum statistic). `validate_estimate`
  reproduces the package's bundled example output.
* MCP-counter (Becht et al., Genome Biol 2016;17:218): mean expression of population marker genes (gene-symbol mode).
"""
from __future__ import annotations

import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

SOURCES = {
    "estimate_1.0.13.tar.gz": "https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz",
    "MCPcounter_genes.txt": "https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt",
}


def fetch_references(cache) -> Path:
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    rec = {}
    for name, url in SOURCES.items():
        dest = cache / name
        if not dest.exists():
            urllib.request.urlretrieve(url, dest)
        rec[name] = {"url": url, "sha256": hashlib.sha256(dest.read_bytes()).hexdigest()}
    with tarfile.open(cache / "estimate_1.0.13.tar.gz") as t:
        t.extractall(cache / "estimate_pkg", filter="data")
    (cache / "references.json").write_text(json.dumps(rec, indent=2))
    return cache


def _estimate_dir(cache):
    return Path(cache) / "estimate_pkg" / "estimate" / "inst" / "extdata"


def estimate_gene_sets(cache) -> dict:
    out = {}
    for line in (_estimate_dir(cache) / "SI_geneset.gmt").read_text().splitlines():
        p = line.split("\t")
        out[p[0]] = [g for g in p[2:] if g]
    return out


def estimate_scores(expr: pd.DataFrame, cache) -> pd.DataFrame:
    common = pd.read_csv(_estimate_dir(cache) / "common_genes.txt", sep="\t")["GeneSymbol"]
    e = expr.loc[sorted(set(common) & set(expr.index))]
    m = e.to_numpy(float)
    genes = np.array(e.index)
    ng = m.shape[0]
    ranks = np.apply_along_axis(lambda c: pd.Series(c).rank(method="average").to_numpy(), 0, m)
    m = 10000 * ranks / ng
    res = {}
    for name, gs in estimate_gene_sets(cache).items():
        tagset = np.isin(genes, gs)
        es = np.empty(m.shape[1])
        for j in range(m.shape[1]):
            order = np.argsort(-m[:, j], kind="stable")
            tag = tagset[order]
            w = np.abs(m[order, j]) ** 0.25
            fn = np.cumsum(np.where(tag, w, 0.0)) / w[tag].sum()
            f0 = np.cumsum(~tag) / (ng - tag.sum())
            es[j] = np.sum(fn - f0)
        res[name] = es
    out = pd.DataFrame({"StromalScore": res["StromalSignature"], "ImmuneScore": res["ImmuneSignature"]}, index=expr.columns)
    out["ESTIMATEScore"] = out.StromalScore + out.ImmuneScore
    out.attrs["n_common_genes"] = int(e.shape[0])
    return out


def validate_estimate(cache) -> float:
    """Max relative error against the package's bundled example (sample_input.gct → sample_estimate.gct)."""
    d = _estimate_dir(cache)
    inp = pd.read_csv(d / "sample_input.gct", sep="\t", skiprows=2, index_col=0).drop(columns="Description")
    ref = pd.read_csv(d / "sample_estimate.gct", sep="\t", skiprows=2, index_col=0).drop(columns="Description").T
    r = estimate_scores(inp, cache)
    return float(np.abs((r[["StromalScore", "ImmuneScore", "ESTIMATEScore"]].to_numpy() - ref.to_numpy()) / ref.to_numpy()).max())


def mcp_markers(cache) -> dict:
    g = pd.read_csv(Path(cache) / "MCPcounter_genes.txt", sep="\t")
    return g.groupby("Cell population")["HUGO symbols"].apply(list).to_dict()


def mcp_scores(expr: pd.DataFrame, cache, populations=("Fibroblasts", "Endothelial cells", "Cytotoxic lymphocytes", "T cells", "Monocytic lineage")):
    mk = mcp_markers(cache)
    out = {}
    for pop in populations:
        genes = [g for g in mk.get(pop, []) if g in expr.index]
        if genes:
            out[f"MCP_{pop}"] = expr.loc[genes].mean(axis=0)
    return pd.DataFrame(out)
