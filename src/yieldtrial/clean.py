"""Cleaning and quality control of GPS yield-monitor data. Every step is logged."""
from __future__ import annotations
import numpy as np
import pandas as pd
import geopandas as gpd

TARGET_MOISTURE = 15.0


class CleaningLog:
    def __init__(self):
        self.rows = []

    def add(self, step, before, after, note=""):
        self.rows.append(dict(step=step, before=before, after=after, removed=before - after, note=note))

    def frame(self):
        return pd.DataFrame(self.rows)


def assign_passes(df, gap_s=5.0):
    df = df.sort_values("time").copy()
    df["pass_id"] = (df["time"].diff().dt.total_seconds().fillna(0) > gap_s).cumsum()
    return df


def flag_gps_jumps(df, max_step_m=10.0):
    """True where the step to a neighbouring point is implausible for the machine speed."""
    g = df.groupby("pass_id")
    prev = np.hypot(g["x"].diff(), g["y"].diff())
    nxt = np.hypot(g["x"].diff(-1), g["y"].diff(-1))
    return ((prev > max_step_m) & (nxt > max_step_m)) | ((prev > max_step_m) & nxt.isna()) | \
           ((nxt > max_step_m) & prev.isna())


def _apply_lag(df, lag):
    return df.groupby("pass_id")["yield_wet"].shift(-lag)


def estimate_flow_lag(df, candidates=range(0, 13), bin_m=3.0, trim=12):
    """Pick the time shift that makes adjacent, opposite-direction passes agree.

    Passes alternate direction, so a wrong lag displaces neighbouring passes in opposite
    directions along the track. The correct lag minimises their disagreement.
    """
    scores = {}
    for lag in candidates:
        d = df.copy()
        d["y_lag"] = _apply_lag(d, lag)
        d["idx"] = d.groupby("pass_id").cumcount()
        d["n"] = d.groupby("pass_id")["idx"].transform("max")
        d = d[(d["idx"] > trim) & (d["idx"] < d["n"] - trim) & d["y_lag"].notna()]
        d["ybin"] = (d["y"] // bin_m).astype(int)
        prof = d.groupby(["pass_id", "ybin"])["y_lag"].mean().unstack(0)
        diffs = [(prof[p] - prof[p + 1]).dropna() for p in prof.columns if p + 1 in prof.columns]
        d_all = pd.concat(diffs)
        scores[lag] = float(d_all.abs().median())  # robust to spikes and strip boundaries
    best = min(scores, key=scores.get)
    return best, pd.Series(scores, name="mean_sq_diff")


def robust_outliers(s, k=4.0):
    med = s.median()
    mad = 1.4826 * (s - med).abs().median()
    return (s - med).abs() > k * mad


def clean_harvest(raw, strips, buffer_m=1.5, apply_lag=True, do_buffer=True):
    log = CleaningLog()
    df = assign_passes(pd.DataFrame(raw.drop(columns="geometry")))
    log.add("raw records", len(df), len(df))

    n = len(df)
    df = df[~flag_gps_jumps(df)]
    log.add("GPS jumps", n, len(df), "step >10 m")

    lag, scores = 0, None
    if apply_lag:
        lag, scores = estimate_flow_lag(df)
        df["yield_wet"] = _apply_lag(df, lag)
        n = len(df)
        df = df[df["yield_wet"].notna()]
        log.add("flow-lag correction", n, len(df), f"estimated lag = {lag} s")

    n = len(df)
    df = df[df["header_down"]]
    log.add("header up", n, len(df))

    df["idx"] = df.groupby("pass_id").cumcount()
    df["n"] = df.groupby("pass_id")["idx"].transform("max")
    n = len(df)
    df = df[(df["idx"] >= 3) & (df["idx"] <= df["n"] - 3)]
    log.add("start/end-of-pass trim", n, len(df), "3 s each end")

    n = len(df)
    df = df[(df["speed"] >= 1.0) & (df["speed"] <= 4.0)]
    log.add("speed limits", n, len(df), "1.0-4.0 m/s")

    n = len(df)
    df["yield_15"] = df["yield_wet"] * (100 - df["moisture"]) / (100 - TARGET_MOISTURE)
    df["field"] = df["field_true"]
    bad = df.groupby("field")["yield_15"].transform(lambda s: robust_outliers(s))
    df = df[~bad]
    log.add("yield outliers", n, len(df), "median/MAD, k=4, by field")

    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x, df.y), crs=raw.crs)
    st = strips.copy()
    if do_buffer:
        st["geometry"] = st.geometry.buffer(-buffer_m)
    n = len(gdf)
    gdf = gpd.sjoin(gdf, st[["strip_id", "treatment", "geometry"]], how="inner", predicate="within")
    gdf = gdf.drop(columns="index_right")
    log.add("clip to strips", n, len(gdf), f"inner buffer {buffer_m if do_buffer else 0} m")
    return gdf, log.frame(), lag, scores
