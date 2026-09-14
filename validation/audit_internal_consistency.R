project_dir <- local({
  source(file.path("validation", "validation_helpers.R"), local = TRUE)
  resolve_project_dir()
})
setwd(project_dir)
source(file.path("validation", "validation_helpers.R"))
load_LL3_source(project_dir)

output_dir <- file.path(project_dir, "validation", "results", "internal_consistency")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

formula_checks <- data.frame(
  check = c(
    "cdf_quantile_inversion",
    "density_formula",
    "support_safe_mu_inverse",
    "support_lower_below_data",
    "seasonal_parameter_counts",
    "station_AICc_recalculation",
    "ERA5_AICc_recalculation"
  ),
  passed = FALSE,
  maximum_absolute_difference = NA_real_,
  detail = NA_character_,
  stringsAsFactors = FALSE
)

parameters <- data.frame(mu = c(-50, -10, 0), sigma = c(70, 15, 3), nu = c(1.8, 2.4, 1.3))
probabilities <- seq(0.001, 0.999, length.out = 501)
inversion_error <- max(vapply(seq_len(nrow(parameters)), function(i) {
  max(abs(pLL3(
    qLL3(probabilities, parameters$mu[i], parameters$sigma[i], parameters$nu[i]),
    parameters$mu[i], parameters$sigma[i], parameters$nu[i]
  ) - probabilities))
}, numeric(1)))
formula_checks$passed[formula_checks$check == "cdf_quantile_inversion"] <- inversion_error < 1e-12
formula_checks$maximum_absolute_difference[formula_checks$check == "cdf_quantile_inversion"] <- inversion_error
formula_checks$detail[formula_checks$check == "cdf_quantile_inversion"] <- "pLL3(qLL3(p)) on p=0.001,...,0.999"

x <- seq(-49.5, 300, length.out = 1000)
mu <- -50
sigma <- 70
nu <- 1.8
z <- (x - mu) / sigma
manual_density <- (nu / sigma) * z^(nu - 1) / (1 + z^nu)^2
density_error <- max(abs(dLL3(x, mu, sigma, nu) - manual_density))
formula_checks$passed[formula_checks$check == "density_formula"] <- density_error < 1e-12
formula_checks$maximum_absolute_difference[formula_checks$check == "density_formula"] <- density_error
formula_checks$detail[formula_checks$check == "density_formula"] <- "Implemented density versus manuscript expression"

family <- LL3(lower = -100)
eta <- seq(-10, 10, length.out = 101)
mu_values <- family$mu.linkinv(eta)
link_error <- max(abs(family$mu.linkfun(mu_values) - eta))
formula_checks$passed[formula_checks$check == "support_safe_mu_inverse"] <-
  link_error < 1e-10 && all(mu_values < -100)
formula_checks$maximum_absolute_difference[formula_checks$check == "support_safe_mu_inverse"] <- link_error
formula_checks$detail[formula_checks$check == "support_safe_mu_inverse"] <- "mu=L-exp(eta) remains below L"

set.seed(1123)
y <- rLL3(75, -50, 70, 1.8)
lower <- LL3_support_lower(y)
formula_checks$passed[formula_checks$check == "support_lower_below_data"] <- lower < min(y)
formula_checks$detail[formula_checks$check == "support_lower_below_data"] <-
  sprintf("lower=%.12g; min(y)=%.12g", lower, min(y))

expected_df <- c(stationary = 36, sigma_time = 48, nu_time = 48, sigma_nu_time = 60)
station_path <- file.path(project_dir, "analysis", "seyhan_results", "model_comparison.csv")
era_candidates <- file.path(
  project_dir,
  "analysis",
  "era5_results",
  c("final_full_195001_202512", "full_195001_202512"),
  "cell_candidate_metrics.csv"
)
era_path <- era_candidates[file.exists(era_candidates)][1L]

audit_aicc <- function(path) {
  if (!file.exists(path)) return(c(parameter_counts = FALSE, max_difference = NA_real_))
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  observed_df <- tapply(values$df, values$model, unique)
  counts_ok <- all(vapply(names(expected_df), function(name) {
    length(observed_df[[name]]) == 1L && observed_df[[name]] == expected_df[[name]]
  }, logical(1)))
  recalculated <- with(values, deviance + 2 * df + 2 * df * (df + 1) / (n - df - 1))
  c(parameter_counts = counts_ok, max_difference = max(abs(values$AICc - recalculated), na.rm = TRUE))
}

station_audit <- audit_aicc(station_path)
era_audit <- audit_aicc(era_path)
formula_checks$passed[formula_checks$check == "seasonal_parameter_counts"] <-
  isTRUE(as.logical(station_audit[["parameter_counts"]])) &&
  isTRUE(as.logical(era_audit[["parameter_counts"]]))
formula_checks$detail[formula_checks$check == "seasonal_parameter_counts"] <-
  "Expected K: M0=36, Msigma=48, Mnu=48, Msigma+nu=60"

formula_checks$passed[formula_checks$check == "station_AICc_recalculation"] <-
  is.finite(station_audit[["max_difference"]]) && station_audit[["max_difference"]] < 1e-8
formula_checks$maximum_absolute_difference[formula_checks$check == "station_AICc_recalculation"] <-
  station_audit[["max_difference"]]
formula_checks$detail[formula_checks$check == "station_AICc_recalculation"] <- basename(station_path)

formula_checks$passed[formula_checks$check == "ERA5_AICc_recalculation"] <-
  is.finite(era_audit[["max_difference"]]) && era_audit[["max_difference"]] < 1e-8
formula_checks$maximum_absolute_difference[formula_checks$check == "ERA5_AICc_recalculation"] <-
  era_audit[["max_difference"]]
formula_checks$detail[formula_checks$check == "ERA5_AICc_recalculation"] <- basename(era_path)

utils::write.csv(formula_checks, file.path(output_dir, "internal_consistency_checks.csv"), row.names = FALSE)
writeLines(c(
  "# Internal consistency audit",
  "",
  sprintf("Audit date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  sprintf("Checks passed: %d of %d.", sum(formula_checks$passed), nrow(formula_checks)),
  "",
  "The implemented LL3 density, CDF, quantile, support-safe threshold link, seasonal parameter counts, and archived AICc values agree with the manuscript definitions and generated result files.",
  "",
  "Monthly record length means the number of years contributing to one calendar-month fit. The end-to-end seasonal simulation records both total monthly observations and per-calendar-month sample size explicitly."
), file.path(output_dir, "INTERNAL_CONSISTENCY_AUDIT.md"))
write_session_info(file.path(output_dir, "sessionInfo.txt"))

if (!all(formula_checks$passed)) {
  stop("One or more internal-consistency checks failed.")
}
print(formula_checks)
