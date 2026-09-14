args <- commandArgs(trailingOnly = TRUE)
argument_value <- function(prefix, default) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit)) sub(prefix, "", hit[[1L]]) else default
}

project_dir <- local({
  source(file.path("validation", "validation_helpers.R"), local = TRUE)
  resolve_project_dir()
})
setwd(project_dir)
source(file.path("validation", "validation_helpers.R"))
load_LL3_source(project_dir)

B <- as.integer(argument_value("--B=", "499"))
excluded_stations <- strsplit(argument_value("--exclude-stations=", "17934"), ",", fixed = TRUE)[[1L]]
input_values <- argument_value(
  "--values=",
  file.path(project_dir, "analysis", "seyhan_results", "spei_values.csv")
)
input_selection <- argument_value(
  "--selection=",
  file.path(project_dir, "analysis", "seyhan_results", "selected_model_summary.csv")
)
if (!file.exists(input_values) || !file.exists(input_selection)) {
  stop(
    "Restricted station-derived inputs are unavailable. Supply --values= and --selection= ",
    "using outputs from analysis/seyhan_real_data_analysis.R."
  )
}

output_dir <- file.path(project_dir, "analysis", "seyhan_results", "dependence_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
values <- utils::read.csv(input_values, stringsAsFactors = FALSE)
selection <- utils::read.csv(input_selection, stringsAsFactors = FALSE)
values <- values[!as.character(values$station) %in% excluded_stations, , drop = FALSE]
selection <- selection[!as.character(selection$station) %in% excluded_stations, , drop = FALSE]
values$date <- as.Date(values$date)
strong <- selection[selection$inference_status == "strong_support", , drop = FALSE]
control <- empirical_LL3_control()

acf1 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  as.numeric(stats::acf(x, lag.max = 1L, plot = FALSE)$acf[[2L]])
}

fit_selected_months <- function(station, scale_months, model) {
  current <- values[
    as.character(values$station) == as.character(station) &
      values$scale_months == scale_months,
    , drop = FALSE
  ]
  current <- current[order(current$date), ]
  current$time_scaled <- as.numeric(scale(seq_len(nrow(current))))
  current$month <- as.integer(format(current$date, "%m"))
  fits <- vector("list", 12L)
  data_by_month <- vector("list", 12L)
  residual_rows <- list()
  for (month in seq_len(12L)) {
    rows <- current$month == month & is.finite(current$accumulated_balance_mm)
    fit_data <- data.frame(
      response = current$accumulated_balance_mm[rows],
      time_scaled = current$time_scaled[rows]
    )
    formulas <- switch(
      model,
      sigma_time = list(sigma = ~time_scaled, nu = ~1),
      nu_time = list(sigma = ~1, nu = ~time_scaled),
      sigma_nu_time = list(sigma = ~time_scaled, nu = ~time_scaled),
      stop("Unknown selected model: ", model)
    )
    fit <- suppressWarnings(fit_LL3_gamlss(
      response ~ 1,
      sigma.formula = formulas$sigma,
      nu.formula = formulas$nu,
      data = fit_data,
      control = control,
      boundary_action = "error"
    ))
    fits[[month]] <- fit
    data_by_month[[month]] <- fit_data
    residual <- LL3_residual_diagnostics(fit, fit_data$response, max_lag = 1L)
    residual_rows[[month]] <- data.frame(
      station = station,
      scale_months = scale_months,
      selected_model = model,
      month = month,
      n = nrow(fit_data),
      normalized_PIT_acf1 = if (length(residual$acf)) residual$acf[[1L]] else NA_real_,
      normalized_PIT_mean = residual$summary[["mean"]],
      normalized_PIT_sd = residual$summary[["sd"]],
      sigma_slope = if (model %in% c("sigma_time", "sigma_nu_time")) {
        stats::coef(fit, what = "sigma")[["time_scaled"]]
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  list(fits = fits, data = data_by_month, residuals = do.call(rbind, residual_rows), current = current)
}

case_fits <- vector("list", nrow(strong))
residuals <- list()
chronological <- list()
for (i in seq_len(nrow(strong))) {
  fitted <- fit_selected_months(
    strong$station[i],
    strong$scale_months[i],
    strong$best_nonstationary_model[i]
  )
  case_fits[[i]] <- fitted
  residuals[[i]] <- fitted$residuals
  chronological[[i]] <- data.frame(
    station = strong$station[i],
    scale_months = strong$scale_months[i],
    selected_model = strong$best_nonstationary_model[i],
    n = sum(is.finite(fitted$current$LL3_selected_nonstationary)),
    chronological_normalized_PIT_acf1 = acf1(fitted$current$LL3_selected_nonstationary),
    maximum_absolute_monthly_acf1 = max(abs(fitted$residuals$normalized_PIT_acf1), na.rm = TRUE),
    median_monthly_acf1 = stats::median(fitted$residuals$normalized_PIT_acf1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
monthly_residuals <- do.call(rbind, residuals)
case_residuals <- do.call(rbind, chronological)

representative <- do.call(rbind, lapply(split(strong, strong$scale_months), function(d) {
  d[which.max(d$delta_AICc_favoring_nonstationary), , drop = FALSE]
}))
representative <- representative[order(representative$scale_months), ]
bootstrap_rows <- list()
for (i in seq_len(nrow(representative))) {
  case_index <- which(
    strong$station == representative$station[i] &
      strong$scale_months == representative$scale_months[i]
  )[[1L]]
  fitted <- case_fits[[case_index]]
  month <- fitted$residuals$month[
    which.max(abs(fitted$residuals$sigma_slope))
  ]
  fit <- fitted$fits[[month]]
  fit_data <- fitted$data[[month]]
  result <- LL3_moving_block_bootstrap(
    fit = fit,
    data = fit_data,
    response = "response",
    mu.formula = response ~ 1,
    sigma.formula = ~time_scaled,
    nu.formula = ~1,
    B = B,
    seed = 973000L + representative$scale_months[i] * 100L + month,
    control = control,
    show_progress = FALSE,
    boundary_action = "error"
  )
  interval <- result$percentile_95[, "sigma:time_scaled"]
  bootstrap_rows[[i]] <- data.frame(
    station = representative$station[i],
    scale_months = representative$scale_months[i],
    selected_model = representative$best_nonstationary_model[i],
    delta_AICc = representative$delta_AICc_favoring_nonstationary[i],
    representative_month = month,
    original_sigma_slope = stats::coef(fit, what = "sigma")[["time_scaled"]],
    bootstrap_lower_95 = interval[[1L]],
    bootstrap_upper_95 = interval[[2L]],
    B_requested = B,
    B_valid = nrow(result$estimates),
    invalid_refit_rate = 1 - result$convergence_rate,
    boundary_refit_rate = result$boundary_failure_rate,
    block_length = result$block_length,
    normalized_PIT_acf1 = fitted$residuals$normalized_PIT_acf1[fitted$residuals$month == month],
    stringsAsFactors = FALSE
  )
}
bootstrap_summary <- do.call(rbind, bootstrap_rows)

utils::write.csv(monthly_residuals, file.path(output_dir, "strong_case_monthly_PIT_diagnostics.csv"), row.names = FALSE)
utils::write.csv(case_residuals, file.path(output_dir, "strong_case_PIT_summary.csv"), row.names = FALSE)
utils::write.csv(bootstrap_summary, file.path(output_dir, "representative_moving_block_bootstrap.csv"), row.names = FALSE)
writeLines(c(
  sprintf("strong_station_scale_cases=%d", nrow(strong)),
  sprintf("bootstrap_representative_cases=%d", nrow(representative)),
  sprintf("bootstrap_replicates=%d", B),
  sprintf("excluded_stations=%s", paste(excluded_stations, collapse = ",")),
  "exclusion_reason=Pozanti long deterministic climatological infill identified by the station QC workflow",
  "representative_case_rule=largest delta AICc within each accumulation scale",
  "representative_month_rule=largest absolute fitted sigma slope within the selected case",
  "block_length=max(2,ceiling(n^(1/3)))",
  "interval=percentile 95 percent moving-block residual bootstrap",
  "iid Wald uncertainty is not interpreted as dependence-robust"
), file.path(output_dir, "dependence_bootstrap_settings.txt"))
write_session_info(file.path(output_dir, "sessionInfo.txt"))
print(case_residuals)
print(bootstrap_summary)
