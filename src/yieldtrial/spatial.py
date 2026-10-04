"""Aggregation to analysis cells, soil join, spatial diagnostics."""
from __future__ import annotations
import numpy as np
import pandas as pd
import geopandas as gpd
from libpysal.weights import KNN
import warnings
from esda.moran import Moran


def to_cells(clean, soil, cell_len=15.0):
    d = clean.copy()
    d["ybin"] = (d["y"] // cell_len).astype(int)
    cells = (d.groupby(["field", "strip_id", "treatment", "ybin"])
               .agg(yield_t_ha=("yield_15", "mean"), x=("x", "mean"), y=("y", "mean"),
                    n_pts=("yield_15", "size")).reset_index())
    cells = cells[cells["n_pts"] >= 5]
    g = gpd.GeoDataFrame(cells, geometry=gpd.points_from_xy(cells.x, cells.y), crs=clean.crs)
    g = gpd.sjoin_nearest(g, soil[["som_pct", "ph", "geometry"]], how="left").drop(columns="index_right")
    g = g[~g.index.duplicated()]
    return g.reset_index(drop=True)


def morans_i(values, coords, k=8):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        w = KNN.from_array(coords, k=k)
    w.transform = "r"
    m = Moran(values, w)
    return float(m.I), float(m.p_sim)


def empirical_variogram(values, coords, n_bins=12, max_dist=150.0):
    dist = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1))
    sq = 0.5 * (values[:, None] - values[None, :]) ** 2
    iu = np.triu_indices(len(values), 1)
    dist, sq = dist[iu], sq[iu]
    edges = np.linspace(0, max_dist, n_bins + 1)
    idx = np.digitize(dist, edges) - 1
    out = [((edges[b] + edges[b + 1]) / 2, sq[idx == b].mean(), int((idx == b).sum()))
           for b in range(n_bins) if (idx == b).sum() > 20]
    return pd.DataFrame(out, columns=["lag_m", "semivariance", "pairs"])
