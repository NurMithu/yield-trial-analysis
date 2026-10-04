"""Figures and the plain-language report."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fig_cleaning(raw, clean, path, field="F1"):
    fig, ax = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    r = raw[raw["field_true"] == field]
    c = clean[clean["field"] == field]
    ax[0].scatter(r.x, r.y, c=r.yield_wet, s=3, cmap="viridis", vmin=0, vmax=14)
    ax[0].set_title(f"Raw harvester points ({field})")
    sc = ax[1].scatter(c.x, c.y, c=c.yield_15, s=3, cmap="viridis", vmin=4, vmax=12)
    ax[1].set_title("After cleaning (15% moisture)")
    fig.colorbar(sc, ax=ax, label="Yield, t/ha")
    for a in ax:
        a.set_xlabel("Easting (m)")
    ax[0].set_ylabel("Northing (m)")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_lag(scores, best, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(scores.index, scores.values, "o-")
    ax.axvline(best, color="red", ls="--", label=f"selected lag = {best} s")
    ax.set_xlabel("Candidate flow lag (s)")
    ax.set_ylabel("Adjacent-pass disagreement")
    ax.legend()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_variogram(v1, v2, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(v1.lag_m, v1.semivariance, "o-", label="Mixed: field + strip")
    ax.plot(v2.lag_m, v2.semivariance, "s-", label="+ soil/trend covariates")
    ax.set_xlabel("Distance (m)")
    ax.set_ylabel("Semivariance of residuals")
    ax.legend()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_effects(table, truth, path):
    fig, ax = plt.subplots(figsize=(8, 4))
    pos, labels, y = [], [], 0.0
    for contrast in ("B - A", "C - A"):
        for _, r in table[table.contrast == contrast].iterrows():
            ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.lo95], [r.hi95 - r.estimate]], fmt="o", color="C0")
            pos.append(y)
            labels.append(f"{contrast} | {r.model}")
            y += 1
        ax.axvline(truth[contrast[0]], color="red", ls=":", lw=1)
        y += 0.6
    ax.set_yticks(pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Yield difference vs control, t/ha (95% CI); red = simulated truth")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_power(curve, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    for d, g in curve.groupby("delta"):
        ax.plot(g.total_reps_per_trt, g.power, "o-", label=f"{d} t/ha")
    ax.axhline(0.8, color="grey", ls="--")
    ax.set_xlabel("Replicate strips per treatment (total)")
    ax.set_ylabel("Power (alpha = 0.05)")
    ax.legend(title="True difference")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def md_table(df, nd=3):
    d = df.round(nd)
    head = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n"
    return head + "\n".join("| " + " | ".join(map(str, r)) + " |" for r in d.values)


def write_report(path, log, lag, table, truth, mor1, mor2, vc, loo, sens, need, n_cells):
    md = ["# Treatment effect report\n", "## Summary for non-statisticians\n",
          f"- {n_cells} analysis cells built from {int(log.iloc[0].after):,} raw GPS records "
          f"({int(log.iloc[-1].after):,} kept after quality control).",
          f"- The yield sensor lags the GPS position by about **{lag} seconds**; corrected before analysis."]
    for _, r in table[table.model == "Mixed + soil/trend covariates"].iterrows():
        md.append(f"- **{r.contrast}**: {r.estimate:+.2f} t/ha (95% CI {r.lo95:+.2f} to {r.hi95:+.2f}); "
                  f"simulated truth {truth[r.contrast[0]]:+.2f}.")
    md += ["\n## Data-quality log\n", md_table(log),
           "\n## Treatment effects by model\n", md_table(table),
           "\n## Spatial autocorrelation of residuals (Moran's I)\n",
           f"- Mixed (field + strip): I = {mor1[0]:.3f} (p = {mor1[1]:.3f})",
           f"- With soil/trend covariates: I = {mor2[0]:.3f} (p = {mor2[1]:.3f})",
           "\n## Variance components (model with covariates)\n",
           f"field = {vc['field']:.3f}, strip = {vc['strip']:.3f}, residual = {vc['resid']:.3f}",
           "\n## Sensitivity checks\n", md_table(sens),
           "\n### Leave-one-field-out\n", md_table(loo),
           "\n## Replicates needed for 80% power\n", md_table(need)]
    Path(path).write_text("\n".join(md))
