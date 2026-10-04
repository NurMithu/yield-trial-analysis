"""Generate harvester, layout and soil data with a known ground truth.

The generator reproduces the faults found in real yield-monitor exports
(flow delay, header-up records, GPS jumps, outliers, moisture variation) so the
cleaning and modelling steps can be validated against the true treatment effects.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box

CRS = "EPSG:32631"
TRUE_EFFECTS = {"A": 0.0, "B": 0.25, "C": 0.50}   # t/ha relative to control A
TRUE_DELAY_S = 6
STRIP_W, STRIP_L, SWATH = 18.0, 240.0, 6.0


def _trend(x, y, p):
    a, b, c, d = p
    return 0.5 * np.sin(y / 70 + a) + 0.4 * np.cos(x / 55 + b) + 0.3 * np.sin((x + y) / 90 + c) + d


def _fine(y, ph):
    """Short-range patchiness along the field (headland compaction, soil patches)."""
    return 0.45 * np.sin(y / 10 + ph) + 0.25 * np.sin(y / 23 + 2 * ph)


def simulate_trial(seed: int = 42, n_fields: int = 4, reps: int = 2):
    rng = np.random.default_rng(seed)
    treatments = list(TRUE_EFFECTS)
    strips, pts, soil = [], [], []
    t0 = pd.Timestamp("2025-08-20 08:00:00")
    clock = 0.0
    for f in range(n_fields):
        fid = f"F{f+1}"
        x0 = f * 600.0
        field_eff = rng.normal(0, 0.4)
        tp = (rng.uniform(0, 6), rng.uniform(0, 6), rng.uniform(0, 6), 0.0)
        ph = rng.uniform(0, 6)
        layout = rng.permutation(treatments * reps)
        strip_eff = rng.normal(0, 0.15, len(layout))
        for s, trt in enumerate(layout):
            sid = f"{fid}-S{s+1}"
            sx = x0 + s * STRIP_W
            strips.append(dict(field=fid, strip_id=sid, treatment=trt,
                               geometry=box(sx, 0, sx + STRIP_W, STRIP_L)))
            for k in range(int(STRIP_W // SWATH)):
                pass_no = s * 3 + k
                cx = sx + SWATH * (k + 0.5)
                up = pass_no % 2 == 0
                speed = np.clip(rng.normal(2.3, 0.15, 400), 1.6, 3.0)
                ypos = np.cumsum(speed)
                ypos = ypos[ypos < STRIP_L]
                n = len(ypos)
                y = ypos if up else STRIP_L - ypos
                x = cx + rng.normal(0, 0.35, n)
                sp = speed[:n]
                true15 = (8.0 + field_eff + strip_eff[s] + TRUE_EFFECTS[trt]
                          + _trend(x, y, tp) + _fine(y, ph) + rng.normal(0, 0.12, n))
                moist = np.clip(16 + 1.5 * np.sin(y / 60) + rng.normal(0, 0.6, n), 12, 22)
                wet_true = true15 * (100 - 15) / (100 - moist)
                obs = np.empty(n)                       # flow delay
                obs[TRUE_DELAY_S:] = wet_true[: n - TRUE_DELAY_S]
                obs[:TRUE_DELAY_S] = wet_true[0] * np.linspace(0.1, 0.8, TRUE_DELAY_S)
                obs = obs + rng.normal(0, 0.35, n)
                header = np.ones(n, bool)
                header[:2] = False
                header[-2:] = False
                obs[-3:] *= rng.uniform(0.2, 0.6, 3)
                out = rng.random(n) < 0.005             # sensor spikes / dropouts
                obs[out] *= rng.choice([0.0, 3.0], out.sum())
                jump = rng.random(n) < 0.003            # GPS jumps
                x = x + jump * rng.choice([-1, 1], n) * rng.uniform(30, 80, n)
                y = y + jump * rng.choice([-1, 1], n) * rng.uniform(30, 80, n)
                tt = t0 + pd.to_timedelta(clock + np.arange(n), unit="s")
                clock += n + 40
                pts.append(pd.DataFrame(dict(time=tt, x=x, y=y, speed=sp, yield_wet=obs,
                                             moisture=moist, header_down=header, field_true=fid)))
        gx = np.linspace(x0 + 5, x0 + len(layout) * STRIP_W - 5, 8)
        gy = np.linspace(10, STRIP_L - 10, 6)
        for xx in gx:
            for yy in gy:
                som = 3.0 + 1.2 * _trend(xx, yy, tp) + rng.normal(0, 0.1)
                soil.append(dict(field=fid, x=xx, y=yy, som_pct=som, ph=6.5 + rng.normal(0, 0.2)))
    harvest = pd.concat(pts, ignore_index=True)
    harvest = gpd.GeoDataFrame(harvest, geometry=gpd.points_from_xy(harvest.x, harvest.y), crs=CRS)
    strips = gpd.GeoDataFrame(strips, crs=CRS)
    soil = pd.DataFrame(soil)
    soil = gpd.GeoDataFrame(soil, geometry=gpd.points_from_xy(soil.x, soil.y), crs=CRS)
    return harvest, strips, soil
