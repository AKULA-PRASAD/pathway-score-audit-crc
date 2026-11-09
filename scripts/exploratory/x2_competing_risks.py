"""X2. Competing risks: recurrence vs death without recurrence. Registered plan §X2.

Python: event-coding gate, Aalen–Johansen cumulative incidence, cause-specific Cox (unpenalised).
R (x2_finegray.R): Fine–Gray subdistribution model, called from here.
Time unit of rfs.delay is PROVISIONAL (months); absolute quantities carry that label.
"""
import subprocess
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

import xcommon as X

warnings.filterwarnings("ignore")
HORIZONS = (36, 60)


def cs_cox(design, t, e, name):
    df = design.reset_index(drop=True).copy()
    df["_T"], df["_E"] = t, e
    m = CoxPHFitter().fit(df, "_T", "_E")
    s = m.summary
    return pd.DataFrame({"model": name, "term": s.index, "HR": s["exp(coef)"].values,
                         "lo": s["exp(coef) lower 95%"].values, "hi": s["exp(coef) upper 95%"].values, "p": s["p"].values})


def main():
    ctx = X.context()
    out = X.outdir("x2_competing_risks")
    ep = X.endpoint(ctx)
    rec = (ep.rfs_event_full == 1)
    death = (ep.rfs_event_full == 0) & (ep.os_event_source == 1) & (ep.os_time_source == ep.rfs_time_full)
    # ---- gate: must reproduce Phase 6A coding
    gate = {"recurrences": int(rec.sum()), "deaths_without_recurrence": int(death.sum()),
            "censored": int((~rec & ~death).sum()), "expected": {"recurrences": 136, "deaths_without_recurrence": 46}}
    assert gate["recurrences"] == 136 and gate["deaths_without_recurrence"] == 46, gate
    status = np.where(rec, 1, np.where(death, 2, 0))
    t = ep.rfs_time_full.to_numpy(float)
    # ---- Aalen–Johansen: deterministic, tie-aware estimator (R survival::survfit multistate).
    # (Deviation 5: lifelines' AalenJohansenFitter jitters tied times at random; replaced after independent review.)
    stage = ctx.cov.loc[ctx.ids, "tnm.stage"].astype(str).to_numpy()
    aj = pd.DataFrame({"time": t, "status": status, "group": np.where(stage == "2.0", "stage II", "stage III")})
    aj.to_csv(out / "_aj_input.csv", index=False)
    r0 = subprocess.run(["Rscript", str(X.ROOT / "scripts" / "exploratory" / "x2_aj.R"), str(out)], capture_output=True, text=True)
    if r0.returncode != 0:
        raise SystemExit(r0.stderr)
    cif = pd.read_csv(out / "_aj_table.csv")
    cif["cause"] = cif["cause"].map({"recurrence": "recurrence", "death": "death without recurrence"})
    cif.to_csv(out / "cumulative_incidence.csv", index=False)
    cv = pd.read_csv(out / "_aj_curves.csv")
    cv["cause"] = cv["cause"].map({"recurrence": "recurrence", "death": "death without recurrence"})
    cv.to_csv(out / "cumulative_incidence_curves.csv", index=False)
    for f_ in ("_aj_input.csv", "_aj_table.csv", "_aj_curves.csv"):
        (out / f_).unlink()
    # naive 1-KM (recurrence, deaths censored) for contrast
    from lifelines import KaplanMeierFitter
    km = KaplanMeierFitter().fit(t, rec.astype(int))
    naive = {h: float(1 - km.survival_function_at_times(h).iloc[0]) for h in HORIZONS}
    # ---- designs
    clin = X.clinical_design(ctx)
    z, _ = X.cohort_scores(ctx)
    full = pd.concat([clin, z.set_axis(clin.index)], axis=1)
    cs = pd.concat([cs_cox(clin, t, (status == 1).astype(int), "CS recurrence: clinical"),
                    cs_cox(full, t, (status == 1).astype(int), "CS recurrence: clinical + 4 scores"),
                    cs_cox(clin, t, (status == 2).astype(int), "CS death w/o recurrence: clinical"),
                    cs_cox(full, t, (status == 2).astype(int), "CS death w/o recurrence: clinical + 4 scores")])
    cs.to_csv(out / "cause_specific_cox.csv", index=False)
    # ---- Fine–Gray in R
    d = full.copy()
    d.columns = [c.replace("=", "_").replace(".", "_") for c in d.columns]
    d["time"], d["status"] = t, status
    d.to_csv(out / "_design_for_R.csv", index_label="sample_id")
    r = subprocess.run(["Rscript", str(X.ROOT / "scripts" / "exploratory" / "x2_finegray.R"), str(out)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(r.stderr)
    fg = pd.read_csv(out / "fine_gray.csv")
    X.write_json(out / "x2_summary.json", {
        "plan_section": "X2", "label": "exploratory sensitivity; time unit PROVISIONAL (months); CIF by R survfit (tie-aware, deterministic); death-without-recurrence models have 46 events for 10 terms (~4.6 per term): unstable",
        "event_coding_gate": gate, "n": len(t),
        "cif_at_horizons": cif.to_dict("records"),
        "naive_one_minus_KM_recurrence": naive,
        "R_version_output": r.stdout.strip().splitlines()[-3:],
    })
    (out / "_design_for_R.csv").unlink()
    print(gate)
    print(cif.round(4).to_string())
    print("naive 1-KM:", {k: round(v, 4) for k, v in naive.items()})
    print(cs.round(3).to_string())
    print(fg.round(3).to_string())


if __name__ == "__main__":
    main()
