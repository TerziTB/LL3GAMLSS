# LL3GAMLSS 0.1.1

- Added an end-to-end seasonal model-selection simulation at 1-, 3-, 6-, and
  12-month accumulation scales under iid and AR(1) dependence.
- Added direct validation against the maximum-likelihood engine distributed
  with SPEI 1.8.1, including explicit generalized-logistic parameter mapping
  and compatibility diagnostics.
- Added empirical normalized-PIT and moving-block bootstrap analyses.
- Added machine-readable validation, runtime, and release-readiness records.
- `LL3_model_table()` now reports inference readiness and withholds deviance,
  AIC, BIC, and delta values from boundary-contact or otherwise non-regular
  fits.

# LL3GAMLSS 0.1.0

- Added the support-safe LL3 family for classic `gamlss`.
- Added constant-threshold, covariate-dependent scale and shape models.
- Added independent optimization, diagnostic, and uncertainty helpers.
- Added calendar-month, multi-timescale LL3-SPEI fitting.
- Added a circular moving-block residual bootstrap.
- Material threshold-link violations now stop instead of being silently clipped.
- Boundary-contact fits are explicitly marked unsuitable for ordinary
  Hessian, AIC/BIC, or Wald inference.
