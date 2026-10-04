# yield-trial-analysis

A reproducible Python pipeline that takes raw GPS harvester yield-monitor records from replicated
strip trials and produces cleaned, spatially joined data, mixed-model treatment effects with
confidence intervals, sensitivity checks, a power analysis, and a plain-language report.

> **About the data.** The repository generates its own multi-field strip-trial data
> (`src/yieldtrial/simulate.py`) with a *known ground truth* (true treatment effects of +0.25 and
> +0.50 t/ha, a 6 s flow delay, GPS jumps, header-up records, sensor spikes, moisture variation).
> No farm data is included. This allows every step to be validated against the truth, and the
> same code runs on real exports by replacing the loader.

## Quick start

```bash
pip install -r requirements.txt
python run_pipeline.py        # writes figures, tables and report.md to outputs/
pytest                        # data-quality and cleaning tests
```

## Pipeline

| Stage | Module | What it does |
|---|---|---|
| 1. Cleaning | `clean.py` | Pass detection, GPS-jump removal, **flow-lag estimation**, header-up and start/end trimming, speed limits, median/MAD yield outliers, moisture standardisation to 15%, clip to strip polygons with an inner buffer. Every step logs before/after counts. |
| 2. Spatial join | `spatial.py` | Aggregates points to cells within strips, joins nearest soil sample (organic matter, pH), residual diagnostics (Moran's I, empirical variogram). |
| 3. Modelling | `model.py` | Naive OLS vs mixed models (`statsmodels` MixedLM): random intercept for field and variance component for strip within field; second model adds soil and trend covariates. |
| 4. Sensitivity | `pipeline.py` | No lag correction, no edge buffer, leave-one-field-out. |
| 5. Power | `power.py` | Simulation-based power from the fitted variance components. |
| 6. Report | `report.py` | Figures and `outputs/report.md` for non-statisticians. |

## Key methods

**Flow-lag estimation.** Grain takes seconds to reach the sensor, so yield is recorded against the
wrong position. Passes alternate direction, so a wrong lag displaces neighbouring passes in opposite
directions. The lag minimising the median absolute difference between adjacent passes (binned along
the track) is selected. It recovers the true 6 s delay across seeds (tested in `tests/`).

**Why mixed models.** Cells in one strip are not independent. The naive model treats them as
independent, so its confidence intervals are too narrow (SE 0.059 vs 0.101 for the mixed model for
B vs A), and its interval for B - A excludes the true effect. The mixed model carries the correct
uncertainty (see `outputs/fig4_effects.png`).

**Spatial structure.** Residual Moran's I falls from 0.46 (field + strip model) to 0.25 once soil
organic matter and a trend are added. Spatial dependence remains, so for production use a
spatial correlation structure (e.g. a Gaussian-process term) would be the next step; this is a
known limitation of `statsmodels` MixedLM, which has no built-in spatial covariance.

## Results (seed 42)

| Contrast | Estimate (t/ha) | 95% CI | True effect |
|---|---|---|---|
| B - A | +0.19 | +0.06 to +0.32 | +0.25 |
| C - A | +0.37 | +0.24 to +0.50 | +0.50 |

Estimates come from the mixed model with covariates. Full tables, the cleaning log and
sensitivity results are in `outputs/report.md`.

**Replicates for 80% power** (4 fields, variance components from the fitted model):
about 16 strips per treatment for a 0.15 t/ha difference, 8 for 0.25, 4 for 0.40.

## Figures

| | |
|---|---|
| `fig1_cleaning.png` | Raw vs cleaned harvester points |
| `fig2_flow_lag.png` | Flow-lag objective curve |
| `fig3_variogram.png` | Residual semivariance by model |
| `fig4_effects.png` | Treatment effects by model against the true effect |
| `fig5_power.png` | Power vs replicates |

## Limitations

- Estimates from a single trial of 4 fields and 24 strips are noisy; the intervals reflect that.
- Power uses strip means with a field block and a t-test, a simplification of the full mixed model.
- The inner buffer removed no points here because simulated passes are centred in strips; real
  data with swath overlap and wander will lose more.
- Spatial dependence is handled through covariates and trend, not an explicit covariance model.

## Layout

```
src/yieldtrial/   simulate, clean, spatial, model, power, report, pipeline
tests/            cleaning and lag-recovery tests
outputs/          figures, tables, report.md (generated)
run_pipeline.py   entry point
```
