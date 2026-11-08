"""Phase 9 corrected pipeline -- final interval computation for the scientific report
(plan docs/PHASE9_DELTAC_INTERVAL_PLAN.md v2.1, sections 2, 6, 7, 9.1 steps 5-6 and 11;
reporting format of the owner adjudication record, Decision 9).

Usage (repository root, BLAS/OpenMP thread variables unset):
    P=$(mktemp -d /tmp/p9c_pyc.XXXXXX)
    .venv/bin/python -I -B -X pycache_prefix=$P scripts/46_phase9_corrected_report.py \
        <owner-config.json> <output-directory>

Gate. Nothing is computed unless ALL of the following hold:
  * launched with -I -B -X pycache_prefix=<fresh, empty, external directory>;
  * the frozen code-identity check passes (this entry-point module excluded, as for
    scripts/45; its own SHA-256 is recorded in the output);
  * every run has the same number of records L, and its checkpoint is byte-identical to the
    one named by its latest end record;
  * <output_dir>/verification/L<L>.report.json exists with verdict PASS and
    reporting_permitted true, and names the current freeze, manifest and configuration.

Computation. Only frozen functions are used for the protocol quantities:
  orchestration.compute_sizing_decision (EP-1, section 6.5 exclusions, H1, W5, final action),
  statistics.run_summary (percentile primary; basic, z0, A7, A8 diagnostics; kernel MC-SE at
  0.5h, h, 2h), statistics.between_run, reliability.h2 (H2 bounding), reliability.w3_inner_noise,
  reliability.w4. The bias-corrected (z0-only) interval follows the historical analysis
  convention and is a diagnostic only. The pooled interval over all runs is descriptive only
  (Decision 9, item 5). theta-hat is the point estimate pre-specified by the plan (v2.1, line 6
  and section 3.1, line 139), taken from the Phase 9 primary execution; the percentile interval
  does not use it. Its alignment with the bootstrapped functional is documented in
  docs/PHASE9_THETA_ALIGNMENT_AUDIT.md. (The retention-and-disclosure proposal for theta-hat is
  Decision 3 of the adjudication record DRAFT; it is not one of the owner's recorded
  approvals, whose Decision 3 concerns B_in_max.)

Output: <output-directory>/corrected_interval_results.json and a Markdown table file, each
created exclusively (an existing file is never replaced). It never writes inside the pipeline
output directory and never touches a checkpoint, record, manifest or configuration.

Exit status: 0 written; 2 usage; 4 refused (gate not met); 5 output exists.
"""
import contextlib
import hashlib
import json
import math
import os
import sys
import types
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))                    # -I removes the script directory

from src.phase9_corrected import orchestration as orch  # noqa: E402
from src.phase9_corrected import reliability as rel  # noqa: E402
from src.phase9_corrected import statistics as st  # noqa: E402
from src.phase9_corrected.config import ConfigurationError, load_config  # noqa: E402

SELF_RELPATH = "scripts/46_phase9_corrected_report.py"
THETA_SOURCE = "reports/phase9/execution/phase9_execution_results.json"
CONDITIONAL = ("Conditional on the observed 136 recurrence events and 316 censored observations "
               "in the 452-patient A-19 population.")


def fail(code: int, msg: str):
    print(f"p9c-report: REFUSED (code {code}): {msg}", file=sys.stderr)
    raise SystemExit(code)


@contextlib.contextmanager
def _entry_point_excluded():
    main = sys.modules.get("__main__")
    if main is None or Path(getattr(main, "__file__", "")).resolve() != Path(__file__).resolve():
        fail(4, "unexpected entry-point module")
    sys.modules["__main__"] = types.ModuleType("__main__")
    try:
        yield
    finally:
        sys.modules["__main__"] = main


def theta_hat_from_execution(path: Path) -> dict:
    """Plan v2.1 section 3.1: theta-hat = DeltaC_corrected(D) = corrected_C_mean(M2 primary) -
    corrected_C_mean(M1) from the Phase 9 primary execution results (pinned by the input manifest)."""
    d = json.loads(path.read_bytes())
    mc = d["monte_carlo"]
    m2, m1 = mc["M2_primary"]["corrected_C_mean"], mc["M1"]["corrected_C_mean"]
    return {"value": m2 - m1, "C_corrected_M2_primary": m2, "C_corrected_M1": m1,
            "source": THETA_SOURCE, "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "authority": "plan v2.1 line 6 and section 3.1 line 139 (theta_hat = DeltaC_corrected(D))",
            "alignment_audit": "docs/PHASE9_THETA_ALIGNMENT_AUDIT.md",
            "disclosures": [
                "Computed on D by scripts/38 with 5 x 500 event-stratified optimism replicates; each "
                "outer replicate of the interval uses B_in = 100 (same functional, different inner "
                "Monte Carlo precision; plan section 1).",
                "Its Monte Carlo SE (about 0.00019) excludes the cross-validation randomness on D "
                "that selected lambda; the apparent DeltaC on D is one CV realisation "
                "(default_rng(20261102) per model, scripts/38 lines 90-92).",
                "Historical seeding convention: CV draws come from one advancing generator after "
                "each resample (scripts/38 lines 117-129); the apparent-CV generator on D starts "
                "from the same state as the first resample of seed 20261102 (not evaluated; bounded "
                "in the alignment audit)."]}


def bc_z0_interval(v, theta) -> dict:
    """Bias-corrected (z0-only, a = 0) percentile interval -- diagnostic only."""
    from scipy.stats import norm
    z = st.z0(v, theta)
    if not math.isfinite(z):
        return {"lower": None, "upper": None, "alpha1": None, "alpha2": None, "boundary": True}
    a1 = float(norm.cdf(2 * z + norm.ppf(st.ALPHA / 2)))
    a2 = float(norm.cdf(2 * z + norm.ppf(1 - st.ALPHA / 2)))
    return {"lower": st.quantile(v, a1), "upper": st.quantile(v, a2), "alpha1": a1, "alpha2": a2,
            "boundary": False}


def describe(v) -> dict:
    v = np.asarray(v, dtype=np.float64)
    return {"n": int(v.size), "mean": float(v.mean()), "sd": float(v.std(ddof=1)),
            "median": float(np.median(v)), "min": float(v.min()), "max": float(v.max()),
            "share_below_zero": float(np.mean(v < 0))}


def analyse(config, records_by_run: dict, L: int, theta: float) -> dict:
    """All protocol quantities at level L from frozen functions. Pure: reads records only."""
    dec = orch.compute_sizing_decision(config, L, records_by_run)
    runs, endpoints = {}, {}
    pooled = []
    for r in sorted(records_by_run):
        recs = [records_by_run[r][k] for k in range(L)]
        pr = dec["per_run"][str(r)] if str(r) in dec["per_run"] else dec["per_run"][r]
        mcse_min = pr["endpoint_mcse_min"]
        s65 = rel.section_6_5(recs, min(mcse_min.values()))
        usable = s65["usable_k"]
        byk = {x["k"]: x for x in recs}
        v = np.array([byk[k]["dC_corrected"] for k in usable], dtype=np.float64)
        pooled.append(v)
        summ = st.run_summary(v, theta)
        lo, hi = summ["percentile"]
        endpoints[r] = (lo, hi)
        med = float(np.median(v))
        mid = 0.5 * (lo + hi)
        h2 = rel.h2(recs, usable, mcse_min)
        w3 = rel.w3_inner_noise(recs, usable)
        n_fail = sum(len(x.get("failures") or []) for x in recs)
        h6 = sum(1 for x in recs if x.get("status") == "ok"
                 and (x.get("n_used") or {}).get("M1") != (x.get("n_used") or {}).get("M2"))
        cm1 = [x["C_corrected"]["M1"] for x in recs if x.get("status") == "ok"]
        cm2 = [x["C_corrected"]["M2"] for x in recs if x.get("status") == "ok"]
        runs[str(r)] = {
            "n_records": len(recs), "n_usable": len(usable), "section_6_5": s65, "H1": s65["H1"],
            "percentile_primary": {"lower": lo, "upper": hi, "width": hi - lo},
            "basic_diagnostic": dict(zip(("lower", "upper"), summ["basic"])),
            "z0": summ["z0"], "z0_counts": summ["z0_counts"], "bc_z0_diagnostic": bc_z0_interval(v, theta),
            "A7": summ["A7"], "A8_abs_z0": summ["A8"],
            "A7_components": {"tail_asymmetry_mid_minus_median": mid - med,
                              "median_offset_median_minus_theta": med - theta,
                              "note": "A7 = 2|m - theta|/W, with m - theta = (m - median) + (median - theta)"},
            "kernel_mcse": summ["mcse"], "endpoint_mcse_min_decision_5B": mcse_min,
            "H2": h2, "W3": w3, "W5": pr["w5"], "W5_decision": pr["w5_decision"],
            "inner_fit_failures": n_fail, "H6_unequal_paired_counts": h6,
            "distribution_dC_corrected": describe(v),
            "mean_C_corrected": {"M1": float(np.mean(cm1)), "M2": float(np.mean(cm2))},
            "mean_dC_app": float(np.mean([x["dC_app"] for x in recs])),
            "mean_dC_app_as_run": float(np.mean([x["dC_app_as_run"] for x in recs
                                                 if x.get("dC_app_as_run") is not None])),
            "mean_dOpt": float(np.mean([x["dOpt"] for x in recs])),
        }
    between = st.between_run(endpoints)
    w4 = []
    for e in ("lower", "upper"):
        b = between[e]["sd_ddof1"]
        for r in sorted(records_by_run):
            k = runs[str(r)]["kernel_mcse"][e]["SE"][1.0] if 1.0 in runs[str(r)]["kernel_mcse"][e]["SE"] \
                else runs[str(r)]["kernel_mcse"][e]["SE"]["1.0"]
            w4.append({"run": r, "endpoint": e, "kernel_SE_bw_x1": k, "between_run_SD": b, **rel.w4(k, b)})
    allv = np.concatenate(pooled)
    plo, phi = st.percentile_ci(allv)
    return orch._jsonable({
        "level_B_out": L, "B_in": config.B_in, "theta_hat": theta, "alpha": st.ALPHA,
        "conditional_statement": CONDITIONAL,
        "sizing_decision_at_level": dec,
        "runs": runs, "between_run": between, "W4": w4,
        "pooled_descriptive": {"note": "DESCRIPTIVE ONLY -- not a primary or separately established "
                                       "interval (adjudication Decision 9, item 5)",
                               "n": int(allv.size), "percentile": [plo, phi],
                               "distribution": describe(allv)},
        "hard_failures": {
            "H1": [runs[k]["H1"] for k in sorted(runs)],
            "H2": [runs[k]["H2"]["H2"] for k in sorted(runs)],
            "H3": "triggered" if dec["action"] == "stop_label_imprecise" else "not triggered",
            "H6_unequal_paired_counts_total": sum(runs[k]["H6_unequal_paired_counts"] for k in runs)},
        "warnings": {
            "W1_inner_failures_total": sum(runs[k]["inner_fit_failures"] for k in runs),
            "W3": {k: runs[k]["W3"]["W3"] for k in sorted(runs)},
            "W4_any": any(x["W4"] for x in w4),
            "W5": {k: runs[k]["W5"].get("breach") for k in sorted(runs)},
            "W6_A7_gt_0.25": {k: runs[k]["A7"] > 0.25 for k in sorted(runs)},
            "W7_A8_gt_0.25": {k: runs[k]["A8_abs_z0"] > 0.25 for k in sorted(runs)}},
    })


def gate(config, out: Path, freeze_info: dict, msha: str):
    """Returns (L, records_by_run, verification_report). Refuses unless reporting is permitted."""
    records, levels, ck = {}, set(), {}
    for r in config.run_ids:
        raw = (out / "runs" / f"run{r}.jsonl").read_bytes()
        st_ = orch.invocation_state(out, r)
        invs = st_["invocations"]
        if not invs:
            fail(4, f"run {r} has no invocation records")
        last = invs[max(invs)]
        end = last.get("end")
        if not end or end.get("outcome") != "completed":
            fail(4, f"run {r}: latest invocation has no completed end record")
        if hashlib.sha256(raw).hexdigest() != end.get("checkpoint_sha256"):
            fail(4, f"run {r}: checkpoint differs from its latest end record")
        byk = orch._records_by_k(raw)
        records[r], ck[str(r)] = byk, hashlib.sha256(raw).hexdigest()
        levels.add(len(byk))
    if len(levels) != 1:
        fail(4, f"runs are at different levels: {sorted(levels)}")
    L = levels.pop()
    vp = out / "verification" / f"L{L:05d}.report.json"
    try:
        rep = json.loads(vp.read_bytes())
    except (OSError, ValueError):
        fail(4, f"verification report {vp} missing or unreadable")
    if rep.get("verdict") != "PASS" or rep.get("reporting_permitted") is not True:
        fail(4, f"verification at level {L} does not permit reporting (verdict {rep.get('verdict')})")
    if rep.get("level") != L:
        fail(4, "verification report level differs from the checkpoint level")
    if rep.get("freeze_record_sha256") != freeze_info["freeze_record_sha256"] \
            or rep.get("manifest_sha256") != msha or rep.get("config_sha256") != config.sha256():
        fail(4, "verification report identity differs from the current freeze, manifest or configuration")
    return L, records, rep, ck, hashlib.sha256(vp.read_bytes()).hexdigest()


def main(argv) -> int:
    if len(argv) != 2:
        fail(2, "usage: 46_phase9_corrected_report.py <owner-config.json> <output-directory>")
    fl = orch._runtime_flags()
    if not (fl.get("isolated") is True and fl.get("dont_write_bytecode") is True
            and isinstance(fl.get("pycache_prefix"), str) and fl.get("pycache_prefix")):
        fail(4, "not launched with -I -B -X pycache_prefix=<fresh, empty, external directory>")
    dest = Path(argv[1]).resolve()
    jpath, mpath = dest / "corrected_interval_results.json", dest / "corrected_interval_tables.md"
    if jpath.exists() or mpath.exists():
        fail(5, f"output already exists in {dest}; it is never replaced")
    try:
        config = load_config(argv[0], ROOT)
        config.validate_for_execution()
    except ConfigurationError as ex:
        fail(4, f"configuration: {ex}")
    out = Path(config.output_dir)
    if dest == out.resolve() or out.resolve() in dest.parents:
        fail(4, "the report must not be written inside the pipeline output directory")
    try:
        with _entry_point_excluded():
            freeze_info = orch.verify_code_identity(config)
        mraw = (out / "manifest" / "input_manifest.json").read_bytes()
        with _entry_point_excluded():
            orch.verify_code_identity(config, manifest=json.loads(mraw))
    except (orch.CodeIdentityRefusal, OSError, ValueError) as ex:
        fail(4, f"identity: {ex}")
    msha = hashlib.sha256(mraw).hexdigest()
    L, records, rep, ck, vsha = gate(config, out, freeze_info, msha)
    th = theta_hat_from_execution(ROOT / THETA_SOURCE)
    res = analyse(config, records, L, th["value"])
    res = {"generated_utc": datetime.now(timezone.utc).isoformat(),
           "identity": {"freeze_record_sha256": freeze_info["freeze_record_sha256"],
                        "manifest_sha256": msha, "config_sha256": config.sha256(),
                        "checkpoint_sha256": ck, "verification_report_sha256": vsha,
                        "writer": {"path": SELF_RELPATH,
                                   "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}},
           "verification": {"verdict": rep["verdict"], "status_counts": rep.get("status_counts")},
           "theta_hat_source": th, **res}
    dest.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(res, indent=1, sort_keys=True) + "\n").encode()
    fd = os.open(jpath, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        os.write(fd, data)
        os.fsync(fd)
    finally:
        os.close(fd)
    fd = os.open(mpath, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        os.write(fd, tables_markdown(res).encode())
        os.fsync(fd)
    finally:
        os.close(fd)
    print(json.dumps({"written": [str(jpath), str(mpath)], "level": L,
                      "results_sha256": hashlib.sha256(data).hexdigest()}, indent=1))
    return 0


def tables_markdown(res: dict) -> str:
    f = lambda x: f"{x:+.5f}"  # noqa: E731
    L = res["level_B_out"]
    lines = [f"# Corrected nested-bootstrap interval tables (B_out = {L}, B_in = {res['B_in']})", "",
             f"Generated {res['generated_utc']} by `{SELF_RELPATH}`. {res['conditional_statement']}", "",
             "## T1. Per-run 95 % percentile intervals (primary) and diagnostics", "",
             "| Run | n usable | Percentile (primary) | Width | Basic (diagnostic) | BC z0 (diagnostic) | z0 | A7 | A8 |",
             "|---|---:|---|---:|---|---|---:|---:|---:|"]
    for r, x in sorted(res["runs"].items()):
        p, b, c = x["percentile_primary"], x["basic_diagnostic"], x["bc_z0_diagnostic"]
        bc = f"[{f(c['lower'])}, {f(c['upper'])}]" if not c["boundary"] else "boundary"
        lines.append(f"| {r} | {x['n_usable']} | [{f(p['lower'])}, {f(p['upper'])}] | {p['width']:.5f} | "
                     f"[{f(b['lower'])}, {f(b['upper'])}] | {bc} | {x['z0']:+.3f} | {x['A7']:.3f} | {x['A8_abs_z0']:.3f} |")
    br = res["between_run"]
    lines += ["", "## T2. Between-run endpoint variability (SD ddof = 1, df = 4)", "",
              "| Endpoint | Five values | SD | Range |", "|---|---|---:|---:|"]
    for e in ("lower", "upper"):
        lines.append(f"| {e} | {', '.join(f(v) for v in br[e]['values'])} | {br[e]['sd_ddof1']:.5f} | {br[e]['range']:.5f} |")
    lines += ["", "## T3. Kernel endpoint MC-SE (bandwidth multipliers 0.5, 1, 2) and EP-1 ratio", "",
              "| Run | Endpoint | SE ×0.5 | SE ×1 | SE ×2 | EP-1 ratio (Option B) |", "|---|---|---:|---:|---:|---:|"]
    rows = {(x["run"], x["endpoint"]): x for x in res["sizing_decision_at_level"]["ep1"]["rows"]}
    for r, x in sorted(res["runs"].items()):
        for e in ("lower", "upper"):
            se = x["kernel_mcse"][e]["SE"]
            g = lambda m: se.get(m, se.get(str(m)))  # noqa: E731
            row = rows.get((int(r), e)) or rows.get((r, e))
            lines.append(f"| {r} | {e} | {g(0.5):.5f} | {g(1.0):.5f} | {g(2.0):.5f} | {row['ratio']:.4f} |")
    pd = res["pooled_descriptive"]
    lines += ["", f"EP-1 rho_max = {res['sizing_decision_at_level']['ep1']['rho_max']:.4f}; "
                  f"sizing action at this level: `{res['sizing_decision_at_level']['action']}`.", "",
              "## T4. Pooled draws (descriptive only, not a primary interval)", "",
              f"n = {pd['n']}; percentile [{f(pd['percentile'][0])}, {f(pd['percentile'][1])}]; "
              f"mean {f(pd['distribution']['mean'])}, median {f(pd['distribution']['median'])}, "
              f"share below zero {pd['distribution']['share_below_zero']:.3f}.", ""]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
