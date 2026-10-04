"""End-to-end run: data -> cleaning -> cells -> models -> sensitivity -> power -> report."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from . import simulate, clean, spatial, model, power, report


def _need(curve):
    rows = []
    for d, g in curve.groupby("delta"):
        ok = g[g.power >= 0.8]
        rows.append(dict(true_difference_t_ha=d,
                         replicates_per_treatment=int(ok.total_reps_per_trt.min()) if len(ok) else ">32"))
    return pd.DataFrame(rows)


def run(out="outputs", seed=42):
    out = Path(out)
    out.mkdir(exist_ok=True)
    raw, strips, soil = simulate.simulate_trial(seed=seed)
    raw.drop(columns="geometry").to_csv(out / "raw_harvester.csv", index=False)

    cl, log, lag, scores = clean.clean_harvest(raw, strips)
    log.to_csv(out / "cleaning_log.csv", index=False)
    cells = spatial.to_cells(cl, soil)
    fit = model.fit_all(cells)
    c = fit["cells"]
    fit["table"].to_csv(out / "treatment_effects.csv", index=False)

    coords = c[["x", "y"]].to_numpy()
    r1, r2 = fit["m1"].resid.to_numpy(), fit["m2"].resid.to_numpy()
    mor1, mor2 = spatial.morans_i(r1, coords), spatial.morans_i(r2, coords)
    v1, v2 = spatial.empirical_variogram(r1, coords), spatial.empirical_variogram(r2, coords)

    # sensitivity: no lag correction, no edge buffer, leave-one-field-out
    sens = [fit["table"][fit["table"].model == "Mixed + soil/trend covariates"].assign(model="Baseline")]
    for label, kw in [("No lag correction", dict(apply_lag=False)), ("No edge buffer", dict(do_buffer=False))]:
        cl2, *_ = clean.clean_harvest(raw, strips, **kw)
        c2 = model.prepare(spatial.to_cells(cl2, soil))
        sens.append(model.effects_table(model._fit_mixed(c2, "yield_t_ha ~ C(treatment) + som_c + yc + I(yc**2)"), label))
    sens = pd.concat(sens, ignore_index=True)
    loo = model.leave_one_field_out(cells)

    vc = model.variance_components(fit["m2"])
    curve = power.power_curve(vc)
    need = _need(curve)
    curve.to_csv(out / "power_curve.csv", index=False)

    report.fig_cleaning(raw, cl, out / "fig1_cleaning.png")
    report.fig_lag(scores, lag, out / "fig2_flow_lag.png")
    report.fig_variogram(v1, v2, out / "fig3_variogram.png")
    report.fig_effects(fit["table"], simulate.TRUE_EFFECTS, out / "fig4_effects.png")
    report.fig_power(curve, out / "fig5_power.png")
    report.write_report(out / "report.md", log, lag, fit["table"], simulate.TRUE_EFFECTS,
                        mor1, mor2, vc, loo, sens, need, len(c))
    return dict(lag=lag, table=fit["table"], log=log, mor=(mor1, mor2), vc=vc, need=need, sens=sens, loo=loo)
