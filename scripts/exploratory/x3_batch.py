"""X3. Batch estimability, LOBO instability diagnostics and frailty sensitivity. Registered plan §X3.

Reads the frozen leave-one-batch-out (LOBO) results; does not refit them. Frailty models are new,
unpenalised sensitivity fits in R (x3_frailty.R). The frozen ComBat matrix is not reprocessed.
"""
import json
import subprocess

import numpy as np
import pandas as pd
from scipy import stats

import xcommon as X


def main():
    ctx = X.context()
    out = X.outdir("x3_batch")
    ep = X.endpoint(ctx)
    part = pd.read_csv(X.ANALYSIS / "partition_and_batch.csv").set_index("sample_id").loc[ctx.ids]
    cov = ctx.cov.loc[ctx.ids]
    z, _ = X.cohort_scores(ctx)
    rec = ep.rfs_event_full == 1
    death = (ep.rfs_event_full == 0) & (ep.os_event_source == 1)
    # ---- per-batch description
    tab = pd.DataFrame({"batch": part.combat_batch, "rec": rec, "death": death,
                        "stage3": cov["tnm.stage"] == "3.0", "chemo": cov["chemotherapy.adjuvant"] == "Y",
                        "cit": part.dataset == "validation"}).join(z)
    g = tab.groupby("batch")
    desc = pd.DataFrame({"n": g.size(), "recurrences": g.rec.sum(), "deaths_wo_rec": g.death.sum(),
                         "stage3_pct": 100 * g.stage3.mean(), "chemo_pct": 100 * g.chemo.mean(),
                         "cit_partition_n": g.cit.sum()}).join(g[["EMT", "E2F", "Wnt", "p53"]].mean().add_suffix("_mean"))
    desc.round(3).to_csv(out / "batch_description.csv")
    # ---- frozen LOBO diagnostics
    exr = json.loads((X.ROOT / "reports/phase9/execution/phase9_execution_results.json").read_text())
    full = exr["apparent"]["M2_primary"]
    lobo = exr["batch_sensitivity"]["leave_one_batch_out"]
    rows = []
    for b in lobo:
        for pid, lab in X.PATHWAYS:
            d = b["coefficients"][pid] - full["coefficients"][pid]
            rows.append(dict(omitted_batch=b["omitted_batch"], n_omitted=b["n_omitted"], events_omitted=136 - b["events_fit"],
                             pathway=lab, beta_full=full["coefficients"][pid], se_full=full["se"][pid],
                             beta_lobo=b["coefficients"][pid], delta=d, delta_in_se=d / full["se"][pid],
                             sign_change=(b["coefficients"][pid] > 0) != (full["coefficients"][pid] > 0)))
    lb = pd.DataFrame(rows)
    lb.to_csv(out / "lobo_diagnostics.csv", index=False)
    rho = {}
    for lab in lb.pathway.unique():
        s = lb[lb.pathway == lab]
        r, p = stats.spearmanr(s.events_omitted, s.delta_in_se.abs())
        rho[lab] = {"spearman_abs_delta_vs_events_omitted": r, "p": p, "max_abs_delta_in_se": s.delta_in_se.abs().max(),
                    "n_sign_changes": int(s.sign_change.sum()), "all_within_1.96_se": bool((s.delta_in_se.abs() < 1.96).all()),
                    "beta_full": s.beta_full.iloc[0], "se_full": s.se_full.iloc[0]}
    # ---- frailty in R
    clin = X.clinical_design(ctx)
    dsg = pd.concat([clin, z.set_axis(clin.index)], axis=1)
    dsg.columns = [c.replace("=", "_").replace(".", "_") for c in dsg.columns]
    dsg["time"], dsg["event"], dsg["batch"] = ep.rfs_time_full.values, rec.astype(int).values, part.combat_batch.values
    dsg.to_csv(out / "_design_for_R.csv", index_label="sample_id")
    r = subprocess.run(["Rscript", str(X.ROOT / "scripts/exploratory/x3_frailty.R"), str(out)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(r.stderr)
    (out / "_design_for_R.csv").unlink()
    fr = json.loads((out / "frailty.json").read_text())
    strat = pd.read_csv(X.ROOT / "reports/phase8c/score_vs_annotation_stratified.csv")
    strat = strat[(strat.annotation == "combat_batch") & strat.pathway.isin([p for p, _ in X.PATHWAYS])]
    X.write_json(out / "x3_summary.json", {"plan_section": "X3", "label": "sensitivity / diagnostic",
                                          "n_batches": int(desc.shape[0]), "batch_sizes": desc.n.to_dict(),
                                          "events_per_batch": desc.recurrences.to_dict(),
                                          "lobo_by_pathway": rho, "frailty": fr,
                                          "residual_batch_score_association_stratified_by_combat_mod_p": dict(zip(strat.pathway.map(dict(X.PATHWAYS)), strat.location_p_strat))})
    print(desc.round(2).to_string())
    print(pd.DataFrame(rho).T.round(3).to_string())
    print(json.dumps(fr, indent=1)[:3000])


if __name__ == "__main__":
    main()
