"""Stage F/H/I -- render every document from templates so that NO number is typed by hand.

Templates: documentation/*.template.md, presentation/*.template.html, PROJECT_COMPLETION_REPORT.template.md
Placeholders: {{ python-expression | format-spec }}  evaluated against the context below
              {{ TABLE:name }}                       a generated Markdown table
Rendered files are written next to their templates without the '.template' infix.
Also writes tables/*.md and reproducibility/environment.json.
Fails (exit 1) if any placeholder is left unrendered.

Usage: .venv/bin/python -B -u reports/phase9/final_project/analysis/s3_build_docs.py
"""
from __future__ import annotations

import base64, json, re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FP, ROOT, N_RUNS  # noqa: E402

R = json.loads((FP / "results/s1_interval_results.json").read_text())
S0 = json.loads((FP / "results/s0_starting_state.json").read_text())
IC = json.loads((FP / "validation/independent_check_results.json").read_text())
CV = json.loads((FP / "validation/cv_seed_state_evidence.json").read_text())
H5 = json.loads((ROOT / "reports/phase9/h5_spotcheck/h5_comparison.json").read_text())
QCP = FP / "audit/final_qc.json"
QC = json.loads(QCP.read_text()) if QCP.exists() else {"all_ok": None, "failed": [], "checks": {}}
RUNS = R["runs"]; BR = R["between_run"]; EP = R["EP1"]; TH = R["theta_hat"]; BD = R["bias_decomposition_observed"]


def f(x, spec=".5f"):
    return format(x, spec)


def rng(key, path, spec="+.4f"):
    vals = [eval(path, {}, {"x": x}) for x in RUNS]
    return f"{format(min(vals), spec)} to {format(max(vals), spec)}"


def md(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


TABLES = {}
TABLES["T1"] = md(["run", "median", "mean", "SD", "% < 0", "percentile 95 % (primary)", "basic 95 % (diag.)",
                   "BC z₀-only 95 % (diag.)", "z₀", "#{< θ̂}/500"],
                  [[x["run_id"], f(x["distribution"]["dC_corrected"]["median"], "+.5f"),
                    f(x["distribution"]["dC_corrected"]["mean"], "+.5f"), f(x["distribution"]["dC_corrected"]["sd"]),
                    f(100 * x["distribution"]["dC_corrected"]["frac_negative"], ".1f"),
                    f"[{f(x['percentile_primary']['lower'], '+.5f')}, {f(x['percentile_primary']['upper'], '+.5f')}]",
                    f"[{f(x['basic_diagnostic']['lower'], '+.5f')}, {f(x['basic_diagnostic']['upper'], '+.5f')}]",
                    f"[{f(x['bc_z0_diagnostic']['lower'], '+.5f')}, {f(x['bc_z0_diagnostic']['upper'], '+.5f')}]",
                    f(x["z0_diagnostic"]["z0"], "+.3f"), x["z0_diagnostic"]["n_less"]] for x in RUNS])
rows = []
for x in RUNS:
    k = x["kernel_MCSE_2_2"]
    for e in ("lower", "upper"):
        pe = k["percentile_endpoints"][e]
        rows.append([x["run_id"], e, f(pe["q"], "+.5f"), f(k["silverman_h"]),
                     *[f(pe[f"bw_x{m}"]["SE"], ".6f") for m in (0.5, 1.0, 2.0)],
                     *[f(pe[f"bw_x{m}"]["f_hat"], ".2f") for m in (0.5, 1.0, 2.0)],
                     f(pe["crosscheck_scipy_silverman_x1"]["SE"], ".6f")])
TABLES["T2"] = md(["run", "endpoint", "q", "Silverman h", "SE h×0.5", "SE h×1", "SE h×2", "f̂ h×0.5", "f̂ h×1",
                   "f̂ h×2", "SE (scipy 1.06σ rule, cross-check)"], rows)
lab = {"percentile_lower": "percentile lower", "percentile_upper": "percentile upper", "basic_lower": "basic lower",
       "basic_upper": "basic upper", "bc_z0_lower": "BC z₀ lower", "bc_z0_upper": "BC z₀ upper", "z0": "z₀",
       "percentile_width": "percentile width", "median_dC_corrected": "median ΔC_corr(k)",
       "mean_dC_corrected": "mean ΔC_corr(k)"}
TABLES["T3"] = md(["quantity", "run 0", "run 1", "run 2", "run 3", "run 4", "mean", "SD (ddof = 1, df = 4)", "range"],
                  [[lab[k], *[f(v, "+.5f") for v in b["values"]], f(b["mean"], "+.5f"), f(b["sd_ddof1"], ".6f"),
                    f(b["range"], ".6f")] for k, b in BR.items()])
TABLES["T4"] = md(["MC-SE estimator", "evaluations", "pass (≤ 0.05)", "max MC-SE / width"],
                  [[k, v["n"], v["n_pass"], f(v["max_ratio"], ".4f")] for k, v in EP["by_estimator"].items()])
TABLES["T4b"] = md(["run", "endpoint", "estimator", "MC-SE", "width", "MC-SE / width", "result"],
                   [[g["run"], g["endpoint"], g["estimator"], f(g["MCSE"], ".6f"), f(g["width"]), f(g["ratio"], ".4f"),
                     "PASS" if g["pass"] else "**FAIL**"] for g in EP["grid"]])
TABLES["T5"] = md(["run", "inner failures", "unreliable outer (6.5)", "H1", "bounding (6.4)", "H2", "W3 SD inflation",
                   "W5 spread lo / hi (bound)", "W6 A7 lo / hi", "W7 |z₀|"],
                  [[x["run_id"], x["accounting_6_2"]["inner_failures_total"], x["unreliable_6_5"]["n_unreliable"],
                    "breach" if x["unreliable_6_5"]["H1_breach"] else "no", x["bounding_6_4"]["status"].split(" ")[0],
                    "breach" if x["bounding_6_4"]["H2_breach"] else "no",
                    f(100 * x["warnings"]["W3_SD_inflation"], ".2f") + " %" + (" **breach**" if x["warnings"]["W3_breach"] else ""),
                    f"{f(x['warnings']['W5_spread']['lower'])} / {f(x['warnings']['W5_spread']['upper'])} "
                    f"({f(x['warnings']['W5_bound'])})" + (" **breach**" if x["warnings"]["W5_breach"] else ""),
                    f"{f(100*x['warnings']['W6_A7_fraction_of_width']['lower'], '.1f')} % / "
                    f"{f(100*x['warnings']['W6_A7_fraction_of_width']['upper'], '.1f')} %" + (" **breach**" if x["warnings"]["W6_breach"] else ""),
                    f(x["warnings"]["W7_A8_abs_z0"], ".3f") + (" **breach**" if x["warnings"]["W7_breach"] else "")]
                   for x in RUNS])
TABLES["T6"] = md(["run", "ΔC_app(k) mean", "ΔOpt(k) mean", "opt M1 mean", "opt M2 mean", "ΔC_corr skew",
                   "ΔC_corr excess kurtosis", "distinct patients (mean)", "M2 λ median", "M2 λ max"],
                  [[x["run_id"], f(x["distribution"]["dC_app"]["mean"], "+.5f"), f(x["distribution"]["dOpt"]["mean"], "+.5f"),
                    f(x["distribution"]["opt_M1_mean"]["mean"], "+.5f"), f(x["distribution"]["opt_M2_mean"]["mean"], "+.5f"),
                    f(x["distribution"]["dC_corrected"]["skewness"], ".2f"),
                    f(x["distribution"]["dC_corrected"]["excess_kurtosis"], ".2f"),
                    f(x["distribution"]["n_distinct_outer"]["mean"], ".1f"),
                    f(x["distribution"]["lambda_apparent_M2"]["median"], ".4g"),
                    f(x["distribution"]["lambda_apparent_M2"]["max"], ".4g")] for x in RUNS])
TABLES["T7"] = md(["quantity", "original data D (Phase 9 primary)", "mean over 2 500 outer replicates", "shift"],
                  [["ΔC_app", f(BD["dC_app_D"], "+.5f"), f(BD["mean_dC_app_k"], "+.5f"), f(BD["shift_dC_app"], "+.5f")],
                   ["ΔOpt", f(BD["dOpt_D"], "+.5f"), f(BD["mean_dOpt_k"], "+.5f"), f(BD["shift_dOpt"], "+.5f")],
                   ["ΔC_corr = ΔC_app − ΔOpt", f(BD["theta_hat"], "+.5f"), f(BD["mean_dC_corrected_k"], "+.5f"),
                    f(BD["shift_dC_corrected"], "+.5f")]])
TABLES["T8"] = md(["run", "z₀ at θ̂ − 2 MC-SE", "z₀ at θ̂", "z₀ at θ̂ + 2 MC-SE"],
                  [[r, *[f(v, "+.3f") for v in d.values()]] for r, d in R["z0_sensitivity_to_theta_MC_error"]["runs"].items()])
TABLES["H5"] = md(["k", "original ΔC_corr (float.hex)", "recomputed", "bit-identical", "multiplicity", "full record"],
                  [[c["k"], f"`{c['dC_corrected_original_hex']}`", f"`{c['dC_corrected_recomputed_hex']}`",
                    "✓" if c["dC_corrected_bit_identical"] else "✗", "✓" if c["multiplicity_sha_match"] else "✗",
                    "✓" if c["full_record_identical"] else "✗"] for c in H5["comparisons"]])
TABLES["CKPT"] = md(["run", "records", "unique k", "k = 0–499", "status ok", "paired-complete min", "inner failures",
                     "events / censored", "sha256/16"],
                    [[r, v["records"], v["unique_k"], "✓" if v["k_complete_0_to_499"] else "✗", v["status"].get("ok", 0),
                      v["paired_complete_min"], v["inner_failures"],
                      f"{v['events_outer'][0]} / {v['censored_outer'][0]}", f"`{v['sha256_16']}`"]
                     for r, v in S0["checkpoints"].items()])
TABLES["HASH17"] = md(["file", "sha256/16", "status"],
                      [[f"`{k}`", f"`{v['observed']}`", "OK" if v["ok"] else "**CHANGED**"] for k, v in S0["files_17"].items()])

for name, body in TABLES.items():
    (FP / "tables" / f"{name}.md").write_text(body + "\n")

# environment
sys.path.insert(0, str(ROOT))
from src.reproducibility import environment_report  # noqa: E402
ENV = environment_report()
(FP / "reproducibility" / "environment.json").write_text(json.dumps(ENV, indent=1))

CTX = {"R": R, "RUNS": RUNS, "BR": BR, "EP": EP, "TH": TH, "BD": BD, "S0": S0, "IC": IC, "CV": CV, "H5": H5,
       "ENV": ENV, "QC": QC, "rng": rng, "N_RUNS": N_RUNS, "min": min, "max": max, "sum": sum, "len": len, "abs": abs}


def img64(name):
    return "data:image/png;base64," + base64.b64encode((FP / "figures" / f"{name}.png").read_bytes()).decode()


CTX["img64"] = img64
PAT = re.compile(r"\{\{\s*(.+?)\s*\}\}")


def render(text):
    def sub(m):
        expr = m.group(1)
        if expr.startswith("TABLE:"):
            return TABLES[expr[6:].strip()]
        if " | " in expr:
            e, spec = expr.rsplit(" | ", 1)
            return format(eval(e, CTX), spec.strip())
        return str(eval(expr, CTX))
    return PAT.sub(sub, text)


templates = list(FP.rglob("*.template.*"))
left = []
for t in templates:
    out = t.with_name(t.name.replace(".template", ""))
    txt = render(t.read_text())
    if "{{" in txt:
        left.append(str(out))
    out.write_text(txt)
print("rendered:", [str(t.relative_to(FP)) for t in templates])
print("tables:", sorted(TABLES))
if left:
    print("UNRENDERED PLACEHOLDERS IN:", left); sys.exit(1)
