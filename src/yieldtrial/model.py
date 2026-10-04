"""Treatment-effect models: naive OLS vs mixed models, plus sensitivity fits."""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

VC = {"strip": "0 + C(strip_id)"}


def _fit_mixed(df, formula):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return smf.mixedlm(formula, df, groups=df["field"], vc_formula=VC, re_formula="1").fit(
            reml=True, method="lbfgs")


def effects_table(res, label):
    ci = res.conf_int()
    rows = []
    for t in ("B", "C"):
        k = f"C(treatment)[T.{t}]"
        rows.append(dict(model=label, contrast=f"{t} - A", estimate=res.params[k], se=res.bse[k],
                         lo95=ci.loc[k, 0], hi95=ci.loc[k, 1]))
    return pd.DataFrame(rows)


def prepare(cells):
    c = cells.copy()
    c["som_c"] = c["som_pct"] - c["som_pct"].mean()
    c["yc"] = (c["y"] - c["y"].mean()) / 100.0
    return c


def fit_all(cells):
    c = prepare(cells)
    ols = smf.ols("yield_t_ha ~ C(treatment) + C(field)", c).fit()
    m1 = _fit_mixed(c, "yield_t_ha ~ C(treatment)")
    m2 = _fit_mixed(c, "yield_t_ha ~ C(treatment) + som_c + yc + I(yc**2)")
    tab = pd.concat([effects_table(ols, "Naive OLS (cells independent)"),
                     effects_table(m1, "Mixed: field + strip"),
                     effects_table(m2, "Mixed + soil/trend covariates")], ignore_index=True)
    return dict(ols=ols, m1=m1, m2=m2, table=tab, cells=c)


def variance_components(res):
    return dict(field=float(res.cov_re.iloc[0, 0]), strip=float(res.vcomp[0]), resid=float(res.scale))


def leave_one_field_out(cells):
    c = prepare(cells)
    out = []
    for f in sorted(c["field"].unique()):
        r = _fit_mixed(c[c["field"] != f], "yield_t_ha ~ C(treatment)")
        out.append(effects_table(r, f"drop {f}"))
    return pd.concat(out, ignore_index=True)
