"""Portfolio quality control: check every portfolio number and claim against the research record.

Usage: python portfolio/src/qc_check.py
Re-reads the frozen sources independently of values.py where possible, checks the slide text,
diagram sources, README, speaker notes and blueprint, validates the PPTX/PDF/SVG/PNG outputs,
and writes portfolio/PORTFOLIO_QC_REPORT.md. Exits non-zero on any FAIL. Fits no model.
"""
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_deck  # noqa: E402
import values  # noqa: E402

P = ROOT / "portfolio"
rows = []


def check(area, item, reported, source, verified, ok):
    rows.append((area, item, str(reported), source, str(verified), "PASS" if ok else "FAIL"))


def jl(rel):
    return json.loads((ROOT / rel).read_text())


def text(rel):
    return (ROOT / rel).read_text(encoding="utf-8")


V = values.load()
slides = make_deck.build_slides(V)
deck_text = "\n".join(" ".join([s.get("title", ""), s.get("subtitle", ""), s.get("left_title", ""), s.get("right_title", ""), *s.get("left", []), *s.get("right", []), *s.get("captions", []), s["source"], s["notes"]]) for s in slides)
diagram_src = text("portfolio/src/make_diagrams.py") + "\n" + text("portfolio/src/make_slide_visuals.py")
N_SLIDES = len(slides)
N_PICS = sum(len(make_deck.images(sp)) for sp in slides)
d3svg = text("portfolio/diagrams/D3_pathway_interpretation.svg")
readme = text("portfolio/README.md")
notes_md = text("portfolio/SPEAKER_NOTES.md")
blueprint = text("docs/PORTFOLIO_PRESENTATION_BLUEPRINT.md")
report = text("portfolio/RESEARCH_REPORT.md")
interview = text("portfolio/INTERVIEW_PREPARATION.md")
skills = text("portfolio/SKILLS_EVIDENCE_AND_ROLE_FIT.md")
repro = text("portfolio/REPRODUCIBILITY_GUIDE.md")
import html as _html
web = _html.unescape(re.sub(r"<[^>]+>", " ", text("portfolio/web/index.html")))
web_src_compact = re.sub(r"\s+", "", text("portfolio/web/index.html"))
facing = {k: v.replace("**", "") for k, v in {"slides": deck_text, "README": readme, "speaker notes": notes_md, "blueprint": blueprint,
          "research report": report, "interview prep": interview, "skills doc": skills, "reproducibility guide": repro, "web page": web}.items()}

# ---------------------------------------------------------------- numbers vs authoritative sources
fv = jl("reports/figures/final/figure_values.json")
ex = jl("reports/phase9/execution/phase9_execution_results.json")
s7 = jl("reports/phase9/execution_s7/phase9_s7_results.json")
iv = jl("reports/phase9/corrected_rerun_report/corrected_interval_results.json")
p4 = jl("config/pathways/m2_panel_4set.json")
p5 = jl("config/pathways/m2_panel_5set.json")

f1 = fv["fig1"]
for item, rep, src_v in [("A-19 patients", 452, f1["A19"]["n"]), ("A-19 recurrences", 136, f1["A19"]["events"]),
                         ("A-19 censored", 316, f1["A19"]["censored"]), ("Arrays deposited", 585, f1["arrays_total"]),
                         ("Deaths without recurrence (censored)", 46, f1["deaths_censored"]),
                         ("CIT partition patients", 99, f1["A19_by_partition"]["validation"]["n"]),
                         ("CIT partition events", 29, f1["A19_by_partition"]["validation"]["events"]),
                         ("Execution population n", 452, ex["population"]["n"]), ("Execution population events", 136, ex["population"]["events"])]:
    check("Cohort", item, rep, "figure_values.json fig1 / execution JSON", src_v, rep == src_v)
pr = fv["fig2"]["probes"]
for item, rep, src_v in [("Probe sets", 54675, pr["n_probes_total"]), ("Single-gene probes", 43163, pr["n_unique"]),
                         ("Unique genes", 21815, pr["n_genes_unique"]), ("Unmapped probes", 8893, pr["n_unmapped"]),
                         ("Ambiguous probes", 2619, pr["n_ambiguous"]), ("Four-set member genes", 607, p4["union_unique_genes"]),
                         ("Five-set member genes", 702, p5["union_unique_genes"])]:
    check("Probes and genes", item, rep, "figure_values.json fig2.probes / m2_panel JSON", src_v, rep == src_v)
sizes = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": 198, "HALLMARK_E2F_TARGETS": 195, "HALLMARK_WNT_BETA_CATENIN_SIGNALING": 41, "HALLMARK_P53_PATHWAY": 194}
for pid, n in sizes.items():
    check("Pathways", f"{pid} member genes", n, "m2_panel_4set.json", p4["gene_sets"][pid]["n_resolved_unique"], p4["gene_sets"][pid]["n_resolved_unique"] == n)
check("Pathways", "Four primary set names", ", ".join(sizes), "m2_panel_4set.json gene_sets", ", ".join(p4["gene_sets"]), list(p4["gene_sets"]) == list(sizes))
terms = list(ex["apparent"]["M1"]["hazard_ratios"])
check("Model", "M1 covariates", "age (2 spline terms), sex, stage III vs II, tumour location, adjuvant chemotherapy", "execution JSON apparent.M1",
      ", ".join(terms), terms == ["age_rcs1", "age_rcs2", "sex=Male", "tnm.stage=3.0", "tumor.location=proximal", "chemotherapy.adjuvant=Y"])
st = ex["apparent"]["M1"]
check("Model", "Stage III vs II HR (M1)", "1.97 (1.30–2.98)", "execution JSON", f"{st['hazard_ratios']['tnm.stage=3.0']:.2f} ({st['hr_ci_lower']['tnm.stage=3.0']:.2f}–{st['hr_ci_upper']['tnm.stage=3.0']:.2f})",
      f"{st['hazard_ratios']['tnm.stage=3.0']:.2f} ({st['hr_ci_lower']['tnm.stage=3.0']:.2f}–{st['hr_ci_upper']['tnm.stage=3.0']:.2f})" == "1.97 (1.30–2.98)")
exp_hr = {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": ("1.070 (0.951–1.202)", "1.084 (0.937–1.254)"),
          "HALLMARK_E2F_TARGETS": ("0.906 (0.807–1.016)", "0.872 (0.758–1.003)"),
          "HALLMARK_WNT_BETA_CATENIN_SIGNALING": ("1.016 (0.906–1.138)", "1.011 (0.879–1.162)"),
          "HALLMARK_P53_PATHWAY": ("0.995 (0.882–1.123)", "0.993 (0.852–1.158)")}
for pid, (prim, sens) in exp_hr.items():
    for lab, src, want in (("primary median", ex, prim), ("S-7 representative-IQR", s7, sens)):
        a = src["apparent"]["M2_primary"]
        got = f"{a['hazard_ratios'][pid]:.3f} ({a['hr_ci_lower'][pid]:.3f}–{a['hr_ci_upper'][pid]:.3f})"
        check("Hazard ratios", f"{pid} {lab} (also in D3 SVG and blueprint)", want, "execution / S-7 JSON apparent.M2_primary", got, got == want and want in d3svg and want in blueprint)
a5 = ex["apparent"]["M2_sensitivity"]
got = f"{a5['hazard_ratios']['HALLMARK_E2F_TARGETS']:.3f} ({a5['hr_ci_lower']['HALLMARK_E2F_TARGETS']:.3f}–{a5['hr_ci_upper']['HALLMARK_E2F_TARGETS']:.3f})"
check("Hazard ratios", "Five-set E2F (sensitivity only)", "0.862 (0.746–0.996)", "execution JSON apparent.M2_sensitivity", got, got == "0.862 (0.746–0.996)")
for m, app, cor in (("M1", "0.6242", "0.6016"), ("M2_primary", "0.6382", "0.6059"), ("M2_sensitivity", "0.6413", "0.6066")):
    ga, gc = f"{ex['apparent'][m]['C']:.4f}", f"{ex['monte_carlo'][m]['corrected_C_mean']:.4f}"
    check("Discrimination", f"{m} apparent / corrected C", f"{app} / {cor}", "execution JSON apparent.C / monte_carlo.corrected_C_mean", f"{ga} / {gc}", (ga, gc) == (app, cor))
check("Discrimination", "Apparent ΔC", "+0.0140", "execution JSON", f"{ex['apparent']['M2_primary']['C'] - ex['apparent']['M1']['C']:+.4f}",
      f"{ex['apparent']['M2_primary']['C'] - ex['apparent']['M1']['C']:+.4f}" == "+0.0140")
check("Discrimination", "Corrected ΔC θ̂", "+0.004316", "corrected_interval_results.json theta_hat", f"{iv['theta_hat']:+.6f}", f"{iv['theta_hat']:+.6f}" == "+0.004316")
want_runs = ["[−0.0086, +0.0519]", "[−0.0090, +0.0469]", "[−0.0093, +0.0484]", "[−0.0091, +0.0518]", "[−0.0091, +0.0494]"]
for r, w in enumerate(want_runs):
    x = iv["runs"][str(r)]["percentile_primary"]
    got = f"[{x['lower']:+.4f}, {x['upper']:+.4f}]".replace("-", "−")
    check("Bootstrap", f"Run {r} percentile interval", w, "corrected_interval_results.json", got, got == w and x["lower"] < 0 < x["upper"])
check("Bootstrap", "Design", "5 runs × 750 outer × 100 inner", "figure_values.json fig8", f"{len(iv['runs'])} × {fv['fig8']['B_out']} × {fv['fig8']['B_in']}",
      (len(iv["runs"]), fv["fig8"]["B_out"], fv["fig8"]["B_in"]) == (5, 750, 100))
check("Bootstrap", "Frozen verifier", "PASS 60/60", "figure_values.json fig8.verification", f"{fv['fig8']['verification']['verdict']} {fv['fig8']['verification']['status_counts']['PASS']}/60",
      fv["fig8"]["verification"]["verdict"] == "PASS" and fv["fig8"]["verification"]["status_counts"]["PASS"] == 60)
for m, want in (("M1", "0.833"), ("M2_primary", "0.862"), ("M2_sensitivity", "0.836")):
    got = f"{fv['fig9']['calibration'][m]['corrected']:.3f}"
    check("Calibration", f"{m} corrected slope", want, "figure_values.json fig9.calibration", got, got == want)
with open(ROOT / "reports/phase8c/score_vs_annotation.csv", newline="") as f:
    ann = {r["pathway"]: r for r in csv.DictReader(f) if r["annotation"] == "mmr.status"}
for pid, e2, q in (("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "0.0394", "1.3e-04"), ("HALLMARK_E2F_TARGETS", "0.0251", "0.0021"),
                   ("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "0.0187", "0.0070"), ("HALLMARK_P53_PATHWAY", "0.1716", "4.7e-16")):
    ge, gq = f"{float(ann[pid]['location_epsilon_squared']):.4f}", float(ann[pid]["location_q_BH"])
    gqs = f"{gq:.1e}" if gq < 1e-3 else f"{gq:.4f}"
    check("MMR association", f"{pid} ε² / q (n {ann[pid]['n']})", f"{e2} / {q}", "reports/phase8c/score_vs_annotation.csv", f"{ge} / {gqs}",
          ge == e2 and gqs == q and ann[pid]["n"] == "410" and ann[pid]["circular"] == "False")
d8c = text("docs/PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md")
for lab, pat in (("E2F 0.304", r"E2F_TARGETS \(r̄ = 0\.304\)"), ("EMT 0.241", r"EMT\s*\(0\.241\)"), ("Wnt 0.041", r"WNT \(0\.041\)"), ("p53 0.033", r"P53 \(0\.033\)")):
    check("Coherence", f"Implied inter-gene r̄ {lab}", lab, "PHASE8C §3.1 text", "found" if re.search(pat, d8c) else "not found", bool(re.search(pat, d8c)))
d6c = text("docs/PHASE6C_PATHWAY_PANEL_REVIEW.md")
check("Panel", "10-set condition number 77", 77, "PHASE6C §5 text", "found" if "Condition number 77" in d6c else "not found", "Condition number 77" in d6c)
check("Panel", "4-set / 5-set descriptive condition numbers", "2.95 / 6.95", "m2_panel JSON", f"{p4['collinearity_all461_descriptive_only']['condition_number']:.2f} / {p5['collinearity_all461_descriptive_only']['condition_number']:.2f}",
      f"{p4['collinearity_all461_descriptive_only']['condition_number']:.2f} / {p5['collinearity_all461_descriptive_only']['condition_number']:.2f}" == "2.95 / 6.95")
d5 = text("docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md")
check("External feasibility", "48 and 18 events vs ~338", "48 / 18 / ~338", "PHASE5 §7 item 9 text", "found" if "48 of ~338" in d5 and "18 of ~338" in d5 else "not found", "48 of ~338" in d5 and "18 of ~338" in d5)
d8b = text("docs/PHASE8B_PCA_BATCH_INVESTIGATION.md")
check("combat_mod", "Concordance with MMR, 410 samples", "0 discordant of 410", "PHASE8B text", "found" if re.search(r"0\s+discordant\*\* among the 410", d8b) else "not found", bool(re.search(r"0\s+discordant\*\* among the 410", d8b)))
check("Batch", "LOBO omitted batch sizes", "14–103", "execution JSON leave_one_batch_out", f"{V['lobo_sizes'][0]}–{V['lobo_sizes'][1]}", V["lobo_sizes"] == (14, 103))

# ---------------------------------------------------------------- claims and wording
required = {
    "cause-specific endpoint": r"cause-specific",
    "CIT internal, not external": r"CIT partition.{0,80}internal",
    "runs differ only in seed": r"differ only in random seed",
    "five-set sensitivity only": r"sensitivity( analysis)? only|sensitivity-only",
    "S-7 not replication": r"not (an independent )?replication",
    "combat_mod caveat": r"combat_mod.{0,200}(not an MMR variable|unresolved)",
    "cell composition not measured": r"(Cell composition|cell composition) (was )?(never measured|not measured|was not measured)",
    "batch characterised, not resolved": r"characterised, not resolved|characterised rather than resolved",
    "no proof of zero effect": r"does not (establish|show|prove).{0,40}(zero|no effect)",
}
for k, pat in required.items():
    for name in ("slides", "README", "research report", "web page"):
        hit = bool(re.search(pat, facing[name], re.S | re.I))
        check("Required statement", f"{k} ({name})", "present", name, "present" if hit else "missing", hit)
forbidden = {
    "tool-attribution term": r"\bAI\b|artificial intelligence|claude|coding assistant|AI-assisted|AI-generated",
    "significant improvement claim": r"significant(ly)? improv",
    "predictive success claim": r"\b(predicts|predicted) recurrence\b",
    "external validation claim": r"externally validated|was externally validated",
    "clinical utility claim": r"(demonstrat\w*|establish\w*|show\w*) clinical utility(?! was)",
    "cell composition confounding claim": r"cell composition (confounded|explains|explained)",
    "batch resolved claim": r"batch effects (were )?(removed|eliminated|resolved)(?!\W*(or|,))",
    "most aggressive subtype": r"most aggressive",
    "ε² vs HR comparison": r"tracked MMR more than outcome",
    "mechanistic discovery": r"\b(we|I) (discovered|identified) (a|the) mechanism",
}
for k, pat in forbidden.items():
    for name, t in list(facing.items()) + [("diagram sources", diagram_src)]:
        hits = [m.group(0) for m in re.finditer(pat, t, re.I)]
        # allow explicit negations such as "not resolved" or claims-to-avoid wording
        bad = [h for h in hits if not re.search(r"(not|never|no)\s", t[max(0, t.lower().find(h.lower()) - 40): t.lower().find(h.lower())], re.I)]
        check("Forbidden claim", f"{k} ({name})", "absent", name, "absent" if not bad else f"found: {bad[:2]}", not bad)

# ---------------------------------------------------------------- cross-deliverable number consistency
for v in ["+0.004316", "0.6016", "0.6059", "0.833", "0.862", "0.0394", "0.0251", "0.0187", "0.1716", "1.070 (0.951–1.202)", "0.906 (0.807–1.016)", "1.016 (0.906–1.138)", "0.995 (0.882–1.123)", "452", "136", "607", "21,815", "43,163", "54,675"]:
    for name in ("research report", "web page", "interview prep"):
        hit = v in facing[name]
        if not hit and name == "web page":
            m = re.match(r"([\d.]+) \(([\d.]+)–([\d.]+)\)", v)
            hit = bool(m) and f"[{m.group(1)},{m.group(2)},{m.group(3)}]" in web_src_compact
        check("Consistency", f"{v} in {name}", "present (rendered or page data)", name, "present" if hit else "missing", hit)
for w in want_runs:
    for name in ("research report",):
        check("Consistency", f"interval {w} in {name}", "present", name, "present" if w in facing[name] else "missing", w in facing[name])
# relative links in portfolio README resolve
for link in re.findall(r"\]\(([^)#]+)\)", readme):
    if link.startswith("http"):
        continue
    target = (P / link).resolve()
    check("Links", f"portfolio README link {link}", "resolves", "portfolio/README.md", "exists" if target.exists() else "missing", target.exists())
web_raw = text("portfolio/web/index.html")
check("Web page", "no external script or link hosts other than Google Fonts", "fonts only", "portfolio/web/index.html",
      ", ".join(sorted(set(re.findall(r"https?://([^/\"']+)", web_raw)))), set(re.findall(r"https?://([^/\"']+)", web_raw)) <= {"fonts.googleapis.com", "fonts.gstatic.com", "www.w3.org"})
check("Web page", "conceptual label on tissue illustration", "present", "portfolio/web/index.html", "present" if "not a reconstruction of GSE39582 samples" in web_raw else "missing", "not a reconstruction of GSE39582 samples" in web_raw)

# ---------------------------------------------------------------- artifact validation
from pptx import Presentation  # noqa: E402
from pypdf import PdfReader  # noqa: E402
from PIL import Image  # noqa: E402

pp = P / "presentation" / "GSE39582_recurrence_portfolio.pptx"
prs = Presentation(pp)
n_sl = len(prs.slides)
notes_ok = all(s.has_notes_slide and len(s.notes_slide.notes_text_frame.text) > 80 for s in prs.slides)
pics = sum(1 for s in prs.slides for sh in s.shapes if sh.shape_type == 13)
check("Artifacts", "PPTX opens; slides; speaker notes", f"21 slides (19 + 2 appendix), notes on all, {N_PICS} pictures", str(pp.relative_to(ROOT)), f"{n_sl} slides; notes {'ok' if notes_ok else 'missing'}; {pics} pictures", n_sl == 21 == N_SLIDES and notes_ok and pics == N_PICS)
titles_ok = all(s.shapes.title is not None and len(s.shapes.title.text.strip()) > 5 for s in prs.slides)
check("Accessibility", "PPTX: every slide has a real title placeholder", "21 of 21", str(pp.relative_to(ROOT)), sum(1 for s in prs.slides if s.shapes.title is not None and s.shapes.title.text.strip()), titles_ok)
alts = [sh._element.nvPicPr.cNvPr.get("descr") or "" for s in prs.slides for sh in s.shapes if sh.shape_type == 13]
alt_ok = len(alts) == N_PICS and all(len(a) > 60 and not a.endswith(".png") for a in alts)
check("Accessibility", "PPTX: every picture has descriptive alt text", f"{N_PICS} of {N_PICS}, not file names", str(pp.relative_to(ROOT)), f"{sum(1 for a in alts if len(a) > 60 and not a.endswith('.png'))} of {len(alts)}", alt_ok)
ptext = " ".join(sh.text_frame.text for s in prs.slides for sh in s.shapes if sh.has_text_frame)
check("Artifacts", "PPTX contains no tool-attribution term", "absent", "PPTX text", "absent" if not re.search(forbidden["tool-attribution term"], ptext, re.I) else "found", not re.search(forbidden["tool-attribution term"], ptext, re.I))
pdfp = P / "presentation" / "GSE39582_recurrence_portfolio.pdf"
rd = PdfReader(pdfp)
pdft = " ".join(pg.extract_text() or "" for pg in rd.pages)
check("Artifacts", "PDF valid; pages; key text", "21 pages; contains +0.004316", str(pdfp.relative_to(ROOT)), f"{len(rd.pages)} pages; {'+0.004316' in pdft}", len(rd.pages) == 21 and "0.004316" in pdft)
for d in ("D1_biological_rationale", "D2_processing_feature_engineering", "D3_pathway_interpretation", "D4_heterogeneity_limitations", "D5_model_comparison_uncertainty"):
    svg = (P / "diagrams" / f"{d}.svg").read_text()
    im = Image.open(P / "diagrams" / f"{d}.png")
    ok = svg.startswith("<?xml") and "<text" in svg and im.size == (2666, 1500)
    check("Artifacts", f"{d} SVG (editable text) + PNG", "renders; 2666×1500 PNG", "portfolio/diagrams", f"svg text={'<text' in svg}; png {im.size}", ok)
svg_keys = {
    "D2_processing_feature_engineering": ["54,675", "43,163", "21,815", "607", "2.95", "6.95"],
    "D3_pathway_interpretation": ["1.070 (0.951–1.202)", "0.906 (0.807–1.016)", "1.016 (0.906–1.138)", "0.995 (0.882–1.123)", "0.0394", "0.0251", "0.0187", "0.1716"],
    "D4_heterogeneity_limitations": ["0.0394", "0.1716", "combat_mod", "not tested".upper()],
    "D5_model_comparison_uncertainty": ["0.6016", "0.6059", "+0.004316", "0.833 / 0.862", "differ only in random seed"],
}
for d, keys in svg_keys.items():
    svg = (P / "diagrams" / f"{d}.svg").read_text()
    miss = [k for k in keys if k not in svg]
    check("Artifacts", f"{d} SVG text carries the plotted values (same figure object as the PNG)", ", ".join(keys), "portfolio/diagrams SVG", "all present" if not miss else f"missing {miss}", not miss)
order_ok = all(sl.notes_slide.notes_text_frame.text == sp["notes"] for sl, sp in zip(prs.slides, slides))
check("Artifacts", "Speaker notes attached to the intended slides, in order", "notes match slide specification 1–21", "PPTX notes vs make_deck.py", "match" if order_ok else "mismatch", order_ok and len(prs.slides) == len(slides))
titles_ok = [sp.get("title", "") for sp in slides][16] == "Decisions and audit trail"
check("Artifacts", "Slide 17 title", "Decisions and audit trail", "make_deck.py", slides[16].get("title"), titles_ok)
check("Artifacts", "Cover author and affiliation", "Prasad Akula; MS in Bioinformatics, Northeastern University", "PPTX slide 1",
      "present" if "Prasad Akula" in ptext and "MS in Bioinformatics, Northeastern University" in ptext else "missing",
      "Prasad Akula" in ptext and "MS in Bioinformatics, Northeastern University" in ptext)
for s in slides:
    for img in make_deck.images(s):
        check("Artifacts", f"Slide image exists: {Path(img).name}", "exists", "make_deck.py", "exists" if Path(img).is_file() else "missing", Path(img).is_file())
# ---------------------------------------------------------------- slide visuals (original vector drawings)
SV = P / "diagrams" / "slides"
for f in sorted(SV.glob("V*.svg")):
    svg = f.read_text()
    png = f.with_suffix(".png")
    ok = svg.startswith("<?xml") and "<text" in svg and png.is_file() and "<image" not in svg
    check("Slide visuals", f"{f.stem} SVG (editable text, no embedded raster) + PNG", "valid", "portfolio/diagrams/slides", f"text={'<text' in svg}; raster={'<image' in svg}; png={png.is_file()}", ok)
fmt = lambda t: f"{t[0]:.3f} ({t[1]:.3f}–{t[2]:.3f})"
pw = {x["label"]: x for x in V["pathways"]}
slide_keys = {
    "V04_wnt": ["ESTABLISHED MECHANISM", "CRC LITERATURE", "MEASURED IN THIS STUDY", "NOT TESTED", fmt(pw["Wnt/β-catenin"]["hr"]), "0.041"],
    "V05_e2f": ["ESTABLISHED MECHANISM", "CRC LITERATURE", "MEASURED IN THIS STUDY", "NOT TESTED", fmt(pw["E2F targets"]["hr"]), "sensitivity only"],
    "V06_p53": ["ESTABLISHED MECHANISM", "CRC LITERATURE", "MEASURED IN THIS STUDY", "NOT TESTED", fmt(pw["p53 pathway"]["hr"]), "0.1716"],
    "V07_emt": ["ESTABLISHED MECHANISM", "CRC LITERATURE", "MEASURED IN THIS STUDY", "NOT TESTED", "does not establish tumour-cell", "Fibroblasts"],
    "V08_bulk_tissue": ["NOT SEPARATED", "NOT MEASURED", "Conceptual illustration"],
    "V10_preprocessing_panel": ["54,675", "43,163", "21,815", "607", "2.95", "6.95"],
    "V12_nested_bootstrap": ["750", "Inner loop · 100", "paired", "5 runs", "differ only in random seed", "4.2 %", "60/60"],
    "V13_central_result": ["0.6016", "0.6059", "0.6242", "0.6382", "+0.004316", "did not demonstrate"],
    "V14_forest": [fmt(x["hr"]) for x in V["pathways"]] + [fmt(x["hr_s7"]) for x in V["pathways"]] + ["Sensitivity only", fmt(V["five_set_e2f"])],
    "V15_robustness": ["0.0394", "0.0251", "0.0187", "0.1716", "0.833", "0.862", "0.836", "characterised, not resolved"],
    "V16_project_established": ["What remains unresolved", "11/11", "not yet run", "4/4", "0.94", "0.95", "0.98", "gene-set baseline (grey band", "not tumour-cell EMT", "not a measure of p53", "explained by hypermutation", "+0.004316", "include zero", "provisionally months", "0 of 14"],
    "V17_decisions_audit": ["DC1", "DC2", "DC3", "DC4", "DC5", "DC6", "DC7", "60 / 60"],
    "V18_conclusion_questions": ["+0.004316", "includes zero", "pathway-audit", "Take-home", "NOT TESTED", "Next scientific questions", "≈337", "tumour cells or in stroma"],
    "V19_solution_tool": ["0.94", "0.95", "0.98", "failed quality gate", "pathway-audit", "Survival screen"],
}
for d, keys in slide_keys.items():
    svg = re.sub(r"\s+", " ", _html.unescape(" ".join(re.findall(r">([^<]*)<", (SV / f"{d}.svg").read_text()))))
    miss = [k for k in keys if k not in svg]
    check("Slide visuals", f"{d} text carries the required labels and frozen values", ", ".join(keys)[:80], "portfolio/diagrams/slides SVG", "all present" if not miss else f"missing {miss}", not miss)
# LOBO points plotted in V15 are the frozen per-batch coefficients (same values object)
lobo_ok = len(V["lobo_batches"]) == 9 and all(set(b["coef"]) == set(pw) for b in V["lobo_batches"])
check("Slide visuals", "V15 plots the 9 frozen leave-one-batch-out refits", "9 batches × 4 pathways", "execution JSON leave_one_batch_out", f"{len(V['lobo_batches'])} batches", lobo_ok)

# ---------------------------------------------------------------- post-hoc exploratory work (registered plan)
XD = ["docs/EXPLORATORY_ANALYSIS_PLAN.md", "docs/SCIENTIFIC_LIMITATION_RESOLUTION_REGISTER.md", "docs/ADDITIONAL_ANALYSES_METHODS_AND_RESULTS.md",
      "docs/EXTERNAL_COHORT_ELIGIBILITY_AUDIT.md", "docs/SCIENTIFIC_INTERPRETATION_AND_CONCLUSION.md", "docs/EXPLORATORY_INDEPENDENT_AUDIT_REPORT.md"]
for d in XD:
    check("Exploratory", f"{d} exists", "exists", d, "exists" if (ROOT / d).is_file() else "missing", (ROOT / d).is_file())
vj = jl("reports/exploratory/validation_checks.json")
check("Exploratory", "exploratory validation checks", "all PASS", "reports/exploratory/validation_checks.json", f"{vj['n'] - vj['n_fail']}/{vj['n']}", vj["n_fail"] == 0)
reg_txt = text("docs/SCIENTIFIC_LIMITATION_RESOLUTION_REGISTER.md")
bad_res = [l for l in reg_txt.splitlines() if l.startswith("| ") and re.search(r"\*\*RESOLVED\*\*", l)]
check("Exploratory", "register classifies no limitation as RESOLVED", "none", "register", f"{len(bad_res)} found", not bad_res)
XS = values.load_exploratory()
check("Exploratory", "slide 16 notes carry post-hoc values", "frailty θ and external counts", "make_deck.py",
      "present" if f"{XS['theta']:.2f}" in slides[15]["notes"] and str(XS["ext_assessed"]) in slides[15]["notes"] else "missing",
      f"{XS['theta']:.2f}" in slides[15]["notes"] and str(XS["ext_assessed"]) in slides[15]["notes"])
check("Exploratory", "slide 16 notes say post-hoc and not confirmatory", "present", "make_deck.py",
      "present" if "post-hoc" in slides[15]["notes"] and re.search(r"(not|none of this is) confirmatory", slides[15]["notes"]) else "missing",
      bool("post-hoc" in slides[15]["notes"] and re.search(r"(not|none of this is) confirmatory", slides[15]["notes"])))
for d in XD:
    t = text(d)
    for k in ("significant improvement claim", "external validation claim", "batch resolved claim", "tool-attribution term"):
        hits = [m.group(0) for m in re.finditer(forbidden[k], t, re.I)]
        bad = [h for h in hits if not re.search(r"(not|never|no)\s", t[max(0, t.lower().find(h.lower()) - 40): t.lower().find(h.lower())], re.I)]
        check("Exploratory", f"{k} ({d.split('/')[-1]})", "absent", d, "absent" if not bad else f"found: {bad[:2]}", not bad)

# ---------------------------------------------------------------- report
fails = [r for r in rows if r[5] == "FAIL"]
out = ["# Portfolio quality-control report", "",
       "**Generated by** `portfolio/src/qc_check.py` (8 November 2026). It re-reads the frozen sources, checks every portfolio number and required caveat, "
       "scans all portfolio-facing text for forbidden claims, and validates the presentation, PDF and diagram files. No model was fitted and no frozen file was changed.", "",
       f"**Result:** {len(rows)} checks, {len(rows) - len(fails)} PASS, {len(fails)} FAIL.", "",
       "Portfolio-facing text checked: slide text and speaker notes (`portfolio/src/make_deck.py`), `portfolio/SPEAKER_NOTES.md`, `portfolio/README.md`, "
       "`docs/PORTFOLIO_PRESENTATION_BLUEPRINT.md`, and the diagram text in `portfolio/src/make_diagrams.py`.", "",
       "| Area | Item | Reported in portfolio | Authoritative source | Verified value | Status |", "|---|---|---|---|---|---|"]
out += [f"| {a} | {i} | {r} | `{s}` | {v} | {st} |" if "/" in s or "." in s else f"| {a} | {i} | {r} | {s} | {v} | {st} |" for a, i, r, s, v, st in rows]
out += ["", "## Notes on interpretation-level checks", "",
        "- Forbidden-claim checks are pattern-based. A pattern preceded within 40 characters by an explicit negation (for example \"not resolved\") is accepted.",
        "- The blueprint's \"claims to avoid\" lists describe errors without quoting them and are therefore scanned like any other text.",
        "- The PPTX was validated by re-opening it programmatically (slide count, notes, pictures, text). Visual layout was checked on the matplotlib-rendered PDF pages, "
        "which use the same slide geometry; every PDF page was rasterised with macOS sips and inspected. "
        "PPTX visual inspection is BLOCKED: Keynote and PowerPoint are installed, but a scripted Keynote export timed out (AppleEvent -1712) and Quick Look produced no render, "
        "so the PPTX must be opened manually once (see the acceptance report).",
        "", "## Discrepancies with the research record", "", "None found." if not fails else "\n".join(f"- {r[0]} · {r[1]}: reported {r[2]}, verified {r[4]}" for r in fails), ""]
(P / "PORTFOLIO_QC_REPORT.md").write_text("\n".join(out), encoding="utf-8")
print(f"{len(rows)} checks, {len(fails)} FAIL")
for r in fails:
    print("FAIL", r)
sys.exit(1 if fails else 0)
