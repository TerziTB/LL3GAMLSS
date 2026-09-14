#!/usr/bin/env Rscript

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("The jsonlite package is required to build VALIDATION_SUMMARY.json.")
}
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read_csv <- function(...) utils::read.csv(file.path(root, ...), check.names = FALSE)
read_optional <- function(...) {
  path <- file.path(root, ...)
  if (file.exists(path)) utils::read.csv(path, check.names = FALSE) else data.frame()
}

description <- read.dcf(file.path(root, "DESCRIPTION"))
probability <- read_csv("validation", "results", "core_numerical", "probability_function_checks.csv")
fadist <- read_csv("validation", "results", "core_numerical", "FAdist_probability_agreement.csv")
derivatives <- read_csv("validation", "results", "core_numerical", "derivative_checks.csv")
optimizer <- read_csv("validation", "results", "core_numerical", "stationary_MLE_logLik_agreement.csv")

summary <- list(
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  package = list(name = description[1, "Package"], version = description[1, "Version"]),
  package_checks = read_optional("validation", "results", "release_checks", "r_cmd_check_summary.csv"),
  clean_repository_validation = read_optional("validation", "results", "release_checks", "clean_release_validation.csv"),
  numerical = list(
    maximum_cdf_quantile_inversion_error = max(probability$inversion_error),
    maximum_density_integral_error = max(probability$integral_error),
    maximum_FAdist_density_difference = fadist$density[[1L]],
    maximum_FAdist_CDF_difference = fadist$CDF[[1L]],
    maximum_FAdist_quantile_difference = fadist$quantile[[1L]],
    maximum_gradient_error = max(derivatives$maximum_gradient_error),
    maximum_hessian_error = max(derivatives$maximum_hessian_error),
    stationary_logLik_range = diff(range(optimizer$logLik))
  ),
  internal_consistency = read_csv("validation", "results", "internal_consistency", "internal_consistency_checks.csv"),
  monte_carlo = list(
    stationary_null = read_csv("validation", "results", "stationary_null_summary.csv"),
    scale_trend = read_csv("validation", "results", "sigma_nonstationarity_summary.csv"),
    joint_trend = read_csv("validation", "results", "joint_nonstationarity_summary.csv"),
    interval_coverage = read_csv("validation", "results", "interval_coverage_recovered_summary.csv"),
    stress = read_csv("validation", "results", "publication_stress_summary.csv"),
    seasonal_end_to_end = read_csv("validation", "results", "seasonal_selection", "seasonal_selection_summary.csv")
  ),
  SPEI_maximum_likelihood_equivalence = read_csv(
    "validation", "results", "spei_maxlik_equivalence", "spei_maxlik_equivalence_summary.csv"
  ),
  empirical_dependence = list(
    normalized_PIT = read_optional("analysis", "seyhan_results", "dependence_bootstrap", "strong_case_PIT_summary.csv"),
    moving_block_bootstrap = read_optional("analysis", "seyhan_results", "dependence_bootstrap", "representative_moving_block_bootstrap.csv")
  ),
  ERA5_Land = list(
    runtime = read_optional("analysis", "era5_results", "final_full_195001_202512", "runtime_hardware.csv"),
    package_versions = read_optional("analysis", "era5_results", "final_full_195001_202512", "package_versions.csv"),
    fit_audit = read_optional("analysis", "era5_results", "final_full_195001_202512", "fit_audit.csv")
  ),
  known_external_items = c(
    "SPEI 1.8.1 public fit='max-lik' dispatch returned non-finite coefficients; the shipped internal maximum-likelihood engine was benchmarked separately.",
    "A public repository URL and persistent archive DOI remain to be supplied by the author before manuscript submission."
  )
)

jsonlite::write_json(
  summary,
  file.path(root, "VALIDATION_SUMMARY.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16,
  na = "null"
)
cat("Wrote VALIDATION_SUMMARY.json\n")
