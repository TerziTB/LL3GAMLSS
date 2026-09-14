# LL3GAMLSS publication stress test

Run date: 2026-08-17  
Master seed: 20260817  
Replicates: 500 per design cell (4,500 datasets)  
R: 4.4.1 on Windows 10 x64

## Design

Data were generated from a three-parameter log-logistic distribution with a
constant threshold (`mu = -50`), constant shape (`nu = 1.8`), and a linear
scale predictor `log(sigma) = log(70) + 0.35 x`, where `x` spans `[-1, 1]`.
Sample sizes were 30, 60, and 100. The three scenarios were independent
uniform probabilities (iid), Gaussian AR(1) dependence with `rho = 0.55`, and
iid observations with 5% positive contamination of `4 * sigma`.

Each dataset was fitted with stationary and `sigma ~ x` LL3 GAMLSS models.
A replicate was counted as complete only when both fits converged, the selected
sigma model passed the support/boundary inference diagnostic, and all reported
statistics were finite. AIC and BIC selection rates are therefore calculated
only among complete replicates. Their reported Monte Carlo standard errors use
`sqrt(p * (1 - p) / m)`.

## Results

| Scenario | n | Complete | Boundary | Slope bias | Slope RMSE | AIC selects trend | BIC selects trend |
|---|---:|---:|---:|---:|---:|---:|---:|
| iid | 30 | 95.2% | 4.8% | 0.015 | 0.329 | 40.8% (2.3%) | 26.1% (2.0%) |
| iid | 60 | 99.6% | 0.2% | 0.009 | 0.216 | 57.4% (2.2%) | 32.9% (2.1%) |
| iid | 100 | 100.0% | 0.0% | 0.013 | 0.159 | 76.6% (1.9%) | 51.4% (2.2%) |
| AR(1) | 30 | 94.8% | 5.2% | 0.039 | 0.528 | 55.5% (2.3%) | 42.8% (2.3%) |
| AR(1) | 60 | 99.4% | 0.6% | -0.031 | 0.421 | 55.5% (2.2%) | 40.2% (2.2%) |
| AR(1) | 100 | 100.0% | 0.0% | -0.001 | 0.315 | 69.4% (2.1%) | 52.6% (2.2%) |
| Contaminated | 30 | 93.8% | 6.2% | 0.031 | 0.368 | 40.7% (2.3%) | 23.7% (2.0%) |
| Contaminated | 60 | 99.4% | 0.6% | 0.050 | 0.246 | 61.0% (2.2%) | 32.4% (2.1%) |
| Contaminated | 100 | 100.0% | 0.0% | 0.025 | 0.194 | 70.6% (2.0%) | 43.6% (2.2%) |

Overall, 4,411 of 4,500 replicates (98.0%) were inference-ready and complete.
Of the remaining 89 replicates, 88 were correctly identified threshold-boundary
solutions and one iid `n = 60` replicate failed to return a fit. Every returned
sigma-trend fit reported convergence. Boundary contact disappeared in all three
scenarios at `n = 100`.

## Interpretation

The results support numerical robustness and approximately unbiased recovery of
the scale trend, with precision improving as record length increases. They also
show an important practical limitation: a real but moderate trend is frequently
missed in records of 30 to 60 observations, especially by BIC. Serial dependence
substantially increases slope RMSE and residual autocorrelation, so an iid
likelihood fit should not be treated as providing autocorrelation-adjusted
uncertainty. Positive contamination modestly increases bias at larger sample
sizes. Short-record boundary solutions must remain excluded from ordinary
likelihood inference, as enforced by the package diagnostic.

This stress experiment complements, rather than replaces, the existing null,
parameter-recovery, joint-nonstationarity, direct-optimizer, and interval-
coverage simulations in this repository.
