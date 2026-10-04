import sys
sys.path.insert(0, "src")
import numpy as np
from yieldtrial import simulate, clean


def test_flow_lag_recovered():
    raw, strips, _ = simulate.simulate_trial(seed=7, n_fields=2)
    df = clean.assign_passes(raw.drop(columns="geometry"))
    df = df[~clean.flag_gps_jumps(df)]
    lag, _ = clean.estimate_flow_lag(df)
    assert lag == simulate.TRUE_DELAY_S


def test_gps_jumps_removed():
    raw, strips, _ = simulate.simulate_trial(seed=7, n_fields=2)
    cl, log, *_ = clean.clean_harvest(raw, strips)
    assert (log.removed >= 0).all()
    assert cl["yield_15"].between(2, 20).all()


def test_robust_outliers_flags_extremes():
    import pandas as pd
    s = pd.Series([8.0] * 50 + [30.0, 0.0])
    assert clean.robust_outliers(s).sum() >= 2
