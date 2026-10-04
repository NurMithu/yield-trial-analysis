# Treatment effect report

## Summary for non-statisticians

- 336 analysis cells built from 7,485 raw GPS records (6,418 kept after quality control).
- The yield sensor lags the GPS position by about **6 seconds**; corrected before analysis.
- **B - A**: +0.19 t/ha (95% CI +0.06 to +0.32); simulated truth +0.25.
- **C - A**: +0.37 t/ha (95% CI +0.24 to +0.50); simulated truth +0.50.

## Data-quality log

| step | before | after | removed | note |
|---|---|---|---|---|
| raw records | 7485 | 7485 | 0 |  |
| GPS jumps | 7485 | 7460 | 25 | step >10 m |
| flow-lag correction | 7460 | 7028 | 432 | estimated lag = 6 s |
| header up | 7028 | 6884 | 144 |  |
| start/end-of-pass trim | 6884 | 6452 | 432 | 3 s each end |
| speed limits | 6452 | 6452 | 0 | 1.0-4.0 m/s |
| yield outliers | 6452 | 6418 | 34 | median/MAD, k=4, by field |
| clip to strips | 6418 | 6418 | 0 | inner buffer 1.5 m |

## Treatment effects by model

| model | contrast | estimate | se | lo95 | hi95 |
|---|---|---|---|---|---|
| Naive OLS (cells independent) | B - A | 0.156 | 0.059 | 0.039 | 0.273 |
| Naive OLS (cells independent) | C - A | 0.239 | 0.059 | 0.122 | 0.356 |
| Mixed: field + strip | B - A | 0.156 | 0.101 | -0.043 | 0.354 |
| Mixed: field + strip | C - A | 0.239 | 0.101 | 0.041 | 0.438 |
| Mixed + soil/trend covariates | B - A | 0.189 | 0.067 | 0.059 | 0.32 |
| Mixed + soil/trend covariates | C - A | 0.371 | 0.068 | 0.238 | 0.504 |

## Spatial autocorrelation of residuals (Moran's I)

- Mixed (field + strip): I = 0.462 (p = 0.001)
- With soil/trend covariates: I = 0.250 (p = 0.001)

## Variance components (model with covariates)

field = 0.233, strip = 0.009, residual = 0.127

## Sensitivity checks

| model | contrast | estimate | se | lo95 | hi95 |
|---|---|---|---|---|---|
| Baseline | B - A | 0.189 | 0.067 | 0.059 | 0.32 |
| Baseline | C - A | 0.371 | 0.068 | 0.238 | 0.504 |
| No lag correction | B - A | 0.179 | 0.07 | 0.041 | 0.317 |
| No lag correction | C - A | 0.376 | 0.071 | 0.237 | 0.514 |
| No edge buffer | B - A | 0.189 | 0.067 | 0.059 | 0.32 |
| No edge buffer | C - A | 0.371 | 0.068 | 0.238 | 0.504 |

### Leave-one-field-out

| model | contrast | estimate | se | lo95 | hi95 |
|---|---|---|---|---|---|
| drop F1 | B - A | 0.105 | 0.128 | -0.145 | 0.355 |
| drop F1 | C - A | 0.276 | 0.128 | 0.026 | 0.526 |
| drop F2 | B - A | 0.214 | 0.109 | 0.0 | 0.428 |
| drop F2 | C - A | 0.206 | 0.109 | -0.008 | 0.42 |
| drop F3 | B - A | 0.223 | 0.12 | -0.012 | 0.458 |
| drop F3 | C - A | 0.248 | 0.12 | 0.013 | 0.483 |
| drop F4 | B - A | 0.08 | 0.11 | -0.135 | 0.295 |
| drop F4 | C - A | 0.227 | 0.11 | 0.011 | 0.442 |

## Replicates needed for 80% power

| true_difference_t_ha | replicates_per_treatment |
|---|---|
| 0.15 | 16.0 |
| 0.25 | 8.0 |
| 0.4 | 4.0 |