# Validation

This directory contains the numerical, simulation, workflow, and release
checks used to validate LL3GAMLSS 0.1.1. The archived results are compact
research outputs; restricted observations and large fit caches are not
included.

## Quick checks

Run commands from the repository root after installing the package and its
suggested dependencies.

```sh
Rscript validation/validate_probability_functions.R
Rscript validation/validate_derivatives.R
Rscript validation/validate_stationary_MLE.R
Rscript validation/validate_optimizer_and_lower_sensitivity.R
Rscript validation/validate_spei_maxlik_equivalence.R
Rscript validation/audit_internal_consistency.R
Rscript validation/audit_reported_results.R
Rscript -e "testthat::test_local('.', reporter='summary', stop_on_failure=TRUE)"
Rscript examples/README_smoke_test.R
```

The consolidated local check runner is:

```sh
Rscript scripts/run_clean_release_validation.R
```

## Monte Carlo and seasonal workflow checks

The simulation scripts cover stationary-null selection, scale and joint
non-stationarity, interval coverage, short records, serial dependence,
contamination, and the full twelve-calendar-month model-selection workflow.

```sh
Rscript validation/simulate_stationary_null.R
Rscript validation/simulate_sigma_nonstationarity.R
Rscript validation/simulate_joint_nonstationarity.R
Rscript validation/simulate_interval_coverage.R
Rscript validation/simulate_publication_stress.R
Rscript validation/simulate_seasonal_model_selection.R
Rscript validation/summarize_seasonal_model_selection.R
```

These simulations can be computationally expensive. Their archived seeds,
settings, summaries, raw seasonal-selection results, Monte Carlo standard
errors, failure counts, and session information are retained under
`validation/results/`.

## Output map

- `results/core_numerical/`: probability, derivative, FAdist, and independent
  stationary-likelihood comparisons.
- `results/seasonal_selection/`: end-to-end monthly and multiscale selection
  simulation, including seeds and settings.
- `results/spei_maxlik_equivalence/`: stationary maximum-likelihood benchmark
  against the implementation supplied with SPEI.
- `results/internal_consistency/`: checks connecting formulas, parameter
  counts, fit denominators, tables, and figures.
- `results/reported_results/`: audit of values reported in the manuscript.
- `results/release_checks/`: package-test, smoke-test, clean-install, and
  `R CMD check` logs.
- `VALIDATION_SUMMARY.json`: machine-readable consolidation of the principal
  validation results.

## Validated scope

The validated non-stationary candidate set keeps the LL3 threshold constant
and permits parametric covariates in scale, shape, or both. Boundary-contact
fits are not treated as inference-ready. Covariate-dependent thresholds,
smooth-term prediction in the operational SPEI wrapper, spatial pooling, and
claims that LL3 is universally superior to other distributions are outside the
validated scope. AICc differences are model-ranking evidence, not formal
hypothesis tests.
