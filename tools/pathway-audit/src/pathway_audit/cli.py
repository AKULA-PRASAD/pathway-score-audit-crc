"""Command-line interface.

  pathway-audit demo --out demo/            (synthetic data with planted truth; no downloads)
  pathway-audit fetch-references --cache refs/
  pathway-audit run --expr expr.tsv --gene-sets sets.gmt --out report/ [--annotations ann.csv]
                    [--binary TP53=tp53:M:WT ...] [--adjust mmr=dMMR ...] [--references refs/]
                    [--survival time,event,age,stage3]
"""
from __future__ import annotations

import argparse
import sys

from . import composition, io, report
from .audit import audit
from .report_html import write_html


def main(argv=None):
    p = argparse.ArgumentParser(prog="pathway-audit", description="Audit what gene-set scores measure in bulk expression data.")
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("demo", help="run on built-in synthetic data with planted truth (no downloads)")
    d.add_argument("--out", required=True)
    f = sub.add_parser("fetch-references", help="download ESTIMATE and MCP-counter reference files")
    f.add_argument("--cache", required=True)
    r = sub.add_parser("run", help="run an audit")
    r.add_argument("--expr", required=True, help="genes × samples matrix (csv/tsv/parquet)")
    r.add_argument("--gene-sets", required=True, help="GMT or JSON {name: [genes]}")
    r.add_argument("--out", required=True)
    r.add_argument("--annotations", help="sample-indexed CSV/TSV")
    r.add_argument("--binary", nargs="*", default=[], help="LABEL=column:positive:negative")
    r.add_argument("--adjust", nargs="*", default=[], help="column=positive_value (binary) or column (numeric)")
    r.add_argument("--references", help="reference cache from fetch-references (enables composition analysis)")
    r.add_argument("--composition", help="your own composition estimates (sample-indexed CSV/TSV, numeric columns), e.g. CIBERSORTx output")
    r.add_argument("--fibroblast-markers", help="comma-separated marker genes for the gene-disjoint test (default: MCP-counter fibroblasts)")
    r.add_argument("--survival", help="time_col,event_col,clinical_col1,... (optional ΔC screen; needs lifelines)")
    r.add_argument("--log2", default="auto", choices=["auto", "yes", "no"])
    r.add_argument("--title", default="Pathway-score audit")
    r.add_argument("--seed", type=int, default=0)
    a = p.parse_args(argv)
    if a.cmd == "demo":
        from .demo import run_demo
        res, path = run_demo(a.out)
        print(f"demo report: {path}")
        for name, s in res["gene_sets"].items():
            print(f"  {name}: {'; '.join(s['verdict'])}")
        return 0
    if a.cmd == "fetch-references":
        c = composition.fetch_references(a.cache)
        err = composition.validate_estimate(c)
        print(f"references in {c}; ESTIMATE re-implementation max relative error vs package example: {err:.2e}")
        return 0 if err < 1e-6 else 1
    expr = io.read_expression(a.expr, a.log2)
    sets = io.read_gene_sets(a.gene_sets)
    ann = io.read_annotations(a.annotations) if a.annotations else None
    binary = {}
    for spec in a.binary:
        lab, _, rest = spec.partition("=")
        col, pos, neg = rest.split(":")
        binary[lab] = (col, pos, neg)
    surv = None
    if a.survival:
        cols = a.survival.split(",")
        surv = {"time": cols[0], "event": cols[1], "clinical": cols[2:]}
    comp = None
    if a.composition:
        comp = io.read_annotations(a.composition).apply(lambda c: c.astype(float))
    fibm = a.fibroblast_markers.split(",") if a.fibroblast_markers else None
    res = audit(expr, sets, ann, binary or None, a.adjust or None, a.references, surv, composition=comp, fibroblast_markers=fibm, seed=a.seed)
    path = report.write(res, a.out, a.title)
    hpath = write_html(res, a.out, a.title)
    print(f"report: {path}\nhtml:   {hpath}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
