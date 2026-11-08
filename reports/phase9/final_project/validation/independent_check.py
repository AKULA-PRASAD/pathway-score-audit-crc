"""Independent re-computation of every reported section 9.1 number, by DIFFERENT code paths.

Does NOT import analysis/stats_core.py or analysis/common.py. Uses:
  * pure-Python sorted-list type-7 quantiles (no numpy.quantile)
  * statistics.NormalDist and statistics.stdev
  * scipy.stats.gaussian_kde for the density (a different KDE implementation)
  * raw json parsing of the checkpoints
then compares with results/s1_interval_results.json and with the CSV tables.
Writes ONLY validation/independent_check_results.json.

Usage: .venv/bin/python -B -u reports/phase9/final_project/validation/independent_check.py
"""
from __future__ import annotations

import csv, json, math, statistics
from pathlib import Path

from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[4]
FP = ROOT / "reports/phase9/final_project"
N = statistics.NormalDist()
TOL = 1e-12          # absolute; all quantities are O(1e-2) or larger in relative terms


def q7(xs, p):
    s = sorted(xs); h = (len(s) - 1) * p; lo = math.floor(h)
    return s[lo] + (h - lo) * (s[min(lo + 1, len(s) - 1)] - s[lo])


def load(r):
    rows = [json.loads(l) for l in (ROOT / f"reports/phase9/dc_interval/run{r}_outer.jsonl").read_text().splitlines() if l.strip()]
    rows.sort(key=lambda x: x["k"]); return rows


def main():
    res = json.loads((FP / "results/s1_interval_results.json").read_text())
    ph9 = json.loads((ROOT / "reports/phase9/execution/phase9_execution_results.json").read_text())
    theta = ph9["monte_carlo"]["M2_primary"]["corrected_C_mean"] - ph9["monte_carlo"]["M1"]["corrected_C_mean"]
    checks = []

    def chk(name, mine, theirs, tol=TOL):
        ok = (mine == theirs) if isinstance(mine, (bool, str)) else (
            (math.isinf(mine) and mine == theirs) or abs(mine - theirs) <= tol)
        checks.append({"check": name, "independent": mine, "reported": theirs, "ok": bool(ok)})

    chk("theta_hat", theta, res["theta_hat"]["value"])
    ends = {k: [] for k in ("pl", "pu")}
    for r in range(5):
        rows = load(r); x = [y["dC_corrected"] for y in rows]; R = res["runs"][r]
        lo, hi = q7(x, .025), q7(x, .975)
        chk(f"run{r} pct lower", lo, R["percentile_primary"]["lower"])
        chk(f"run{r} pct upper", hi, R["percentile_primary"]["upper"])
        chk(f"run{r} basic lower", 2 * theta - hi, R["basic_diagnostic"]["lower"])
        chk(f"run{r} basic upper", 2 * theta - lo, R["basic_diagnostic"]["upper"])
        nl = sum(v < theta for v in x); p = nl / len(x); z = N.inv_cdf(p)
        chk(f"run{r} z0", z, R["z0_diagnostic"]["z0"])
        chk(f"run{r} z0 ties", sum(v == theta for v in x), R["z0_diagnostic"]["n_tie"], tol=0)
        a1 = N.cdf(2 * z + N.inv_cdf(.025)); a2 = N.cdf(2 * z + N.inv_cdf(.975))
        chk(f"run{r} BC lower", q7(x, a1), R["bc_z0_diagnostic"]["lower"])
        chk(f"run{r} BC upper", q7(x, a2), R["bc_z0_diagnostic"]["upper"])
        s = statistics.stdev(x); iqr = q7(x, .75) - q7(x, .25)
        h = 0.9 * min(s, iqr / 1.34) * len(x) ** -0.2
        chk(f"run{r} silverman h", h, R["kernel_MCSE_2_2"]["silverman_h"])
        for e, pp, qq in (("lower", .025, lo), ("upper", .975, hi)):
            for m in (0.5, 1.0, 2.0):
                kde = sps.gaussian_kde(x, bw_method=(h * m) / s)
                se = math.sqrt(pp * (1 - pp) / len(x)) / float(kde(qq)[0])
                chk(f"run{r} kernel SE {e} bw x{m}", se,
                    R["kernel_MCSE_2_2"]["percentile_endpoints"][e][f"bw_x{m}"]["SE"], tol=1e-10)
        # MC-1 trace, final 40 %
        tr = {B: (q7(x[:B], .025), q7(x[:B], .975)) for B in (300, 400, 500)}
        for i, e in enumerate(("lower", "upper")):
            sp = max(t[i] for t in tr.values()) - min(t[i] for t in tr.values())
            chk(f"run{r} W5 spread {e}", sp, R["warnings"]["W5_spread"][e])
        # inner noise
        nv = []
        for y in rows:
            d = [(b2[0] - b2[1]) - (b1[0] - b1[1]) for b1, b2 in zip(y["inner_C"]["M1"], y["inner_C"]["M2"])]
            nv.append(statistics.variance(d) / len(d))
        noise = statistics.fmean(nv); obs = statistics.variance(x)
        chk(f"run{r} W3 SD inflation", math.sqrt(obs / (obs - noise)) - 1, R["warnings"]["W3_SD_inflation"], tol=1e-10)
        chk(f"run{r} A7 lower", abs(lo - (2 * theta - hi)) / (hi - lo), R["warnings"]["W6_A7_fraction_of_width"]["lower"])
        chk(f"run{r} failures", sum(len(y["failures"]) for y in rows), R["accounting_6_2"]["inner_failures_total"], tol=0)
        chk(f"run{r} frac negative", sum(v < 0 for v in x) / len(x), R["distribution"]["dC_corrected"]["frac_negative"])
        ends["pl"].append(lo); ends["pu"].append(hi)

    for key, e in (("pl", "lower"), ("pu", "upper")):
        chk(f"between-run SD {e}", statistics.stdev(ends[key]), res["between_run"][f"percentile_{e}"]["sd_ddof1"])
    # EP-1 grid recomputed from the independent values above
    n_pass, n = 0, 0
    for g in res["EP1"]["grid"]:
        r = g["run"]; w = ends["pu"][r] - ends["pl"][r]
        n += 1; n_pass += g["MCSE"] <= 0.05 * w
        chk(f"EP1 width run{r}", w, g["width"])
    chk("EP1 n_pass", n_pass, res["EP1"]["n_pass"], tol=0)

    # tables reproduce JSON
    with open(FP / "tables/T1_intervals_per_run.csv") as fh:
        for row in csv.DictReader(fh):
            r = int(row["run"]); R = res["runs"][r]
            chk(f"T1 run{r} pct_lower", float(row["pct_lower"]), R["percentile_primary"]["lower"], tol=0)
            chk(f"T1 run{r} basic_upper", float(row["basic_upper"]), R["basic_diagnostic"]["upper"], tol=0)
            chk(f"T1 run{r} z0", float(row["z0"]), R["z0_diagnostic"]["z0"], tol=0)
    with open(FP / "tables/T4_EP1_grid.csv") as fh:
        rows = list(csv.DictReader(fh))
        chk("T4 rows", len(rows), len(res["EP1"]["grid"]), tol=0)
        chk("T4 n_pass", sum(r_["pass"] == "True" for r_ in rows), res["EP1"]["n_pass"], tol=0)

    n_ok = sum(c["ok"] for c in checks)
    out = {"n_checks": len(checks), "n_ok": n_ok, "all_ok": n_ok == len(checks),
           "failed": [c for c in checks if not c["ok"]], "checks": checks}
    (FP / "validation/independent_check_results.json").write_text(json.dumps(out, indent=1))
    print(f"independent checks: {n_ok}/{len(checks)} agree; failed: {[c['check'] for c in out['failed']]}")
    return 0 if out["all_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
