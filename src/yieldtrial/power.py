"""Simulation-based power analysis from fitted variance components."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats


def simulate_power(delta, var_field, var_strip, var_resid, reps_per_trt, n_fields,
                   cells_per_strip=12.0, nsim=1000, alpha=0.05, seed=1):
    """Power to detect a B-vs-A difference of `delta` t/ha.

    Layout: each field holds `reps_per_trt` strips of each of 3 treatments.
    Strip mean = treatment + field effect + strip effect + averaged cell noise;
    analysed by a t-test on the treatment coefficient with field as a fixed block.
    """
    rng = np.random.default_rng(seed)
    trt = np.tile(np.repeat([0, 1, 2], reps_per_trt), n_fields)
    fld = np.repeat(np.arange(n_fields), 3 * reps_per_trt)
    X = np.column_stack([np.ones(len(trt)), trt == 1, trt == 2] +
                        [fld == j for j in range(1, n_fields)]).astype(float)
    dof = len(trt) - X.shape[1]
    if dof < 2:
        return float("nan")
    XtXi = np.linalg.inv(X.T @ X)
    sd = np.sqrt(var_strip + var_resid / cells_per_strip)
    hits = 0
    for _ in range(nsim):
        y = delta * (trt == 1) + np.sqrt(var_field) * rng.normal(size=n_fields)[fld] + sd * rng.normal(size=len(trt))
        b = XtXi @ X.T @ y
        r = y - X @ b
        t = b[1] / np.sqrt((r @ r / dof) * XtXi[1, 1])
        hits += 2 * (1 - stats.t.cdf(abs(t), dof)) < alpha
    return hits / nsim


def power_curve(vc, deltas=(0.15, 0.25, 0.40), reps=range(1, 9), n_fields=4, inflate=1.0):
    rows = [dict(delta=d, reps_per_trt_per_field=r, total_reps_per_trt=r * n_fields,
                 power=simulate_power(d, vc["field"], vc["strip"] * inflate, vc["resid"], r, n_fields))
            for d in deltas for r in reps]
    return pd.DataFrame(rows)
