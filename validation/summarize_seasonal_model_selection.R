project_dir <- local({
  source(file.path("validation", "validation_helpers.R"), local = TRUE)
  resolve_project_dir()
})
setwd(project_dir)
source(file.path("validation", "validation_helpers.R"))

output_dir <- file.path(project_dir, "validation", "results", "seasonal_selection")
raw_path <- file.path(output_dir, "seasonal_selection_raw.csv")
if (!file.exists(raw_path)) {
  stop("Run validation/simulate_seasonal_model_selection.R first.")
}
raw <- utils::read.csv(raw_path, stringsAsFactors = FALSE)
sigma_slope_truth <- 0.25

summarize_rate <- function(condition) rate_with_mcse(condition)

summarize_group <- function(d) {
  stationary_truth <- unique(d$truth) == "M0"
  sigma_truth <- unique(d$truth) == "M_sigma"
  false_2 <- summarize_rate(d$selected_delta2 != "M0")
  false_10 <- summarize_rate(d$selected_delta10 != "M0")
  correct_2 <- summarize_rate(d$selected_delta2 == "M_sigma")
  completion <- summarize_rate(d$all_candidates_complete)
  selection_completion <- summarize_rate(!is.na(d$selected_delta2))
  boundary <- summarize_rate(d$any_boundary_contact)
  slope_target <- if (sigma_truth && unique(d$scale_months) == 1L) sigma_slope_truth else NA_real_
  errors <- if (is.finite(slope_target)) d$M_sigma_slope_median - slope_target else numeric()
  errors <- errors[is.finite(errors)]
  bias <- if (length(errors)) mean(errors) else NA_real_
  bias_mcse <- if (length(errors) > 1L) stats::sd(errors) / sqrt(length(errors)) else NA_real_
  rmse <- if (length(errors)) sqrt(mean(errors^2)) else NA_real_
  rmse_mcse <- if (length(errors) > 1L && rmse > 0) {
    stats::sd(errors^2) / sqrt(length(errors)) / (2 * rmse)
  } else {
    NA_real_
  }
  data.frame(
    truth = unique(d$truth),
    dependence = unique(d$dependence),
    phi = unique(d$phi),
    scale_months = unique(d$scale_months),
    replicates = nrow(d),
    total_monthly_observations = unique(d$total_monthly_observations),
    calendar_month_n_min = min(d$minimum_calendar_month_n),
    calendar_month_n_max = max(d$minimum_calendar_month_n),
    false_selection_delta2 = if (stationary_truth) false_2[["estimate"]] else NA_real_,
    false_selection_delta2_mcse = if (stationary_truth) false_2[["mcse"]] else NA_real_,
    false_selection_delta2_n = if (stationary_truth) false_2[["denominator"]] else NA_real_,
    false_selection_delta10 = if (stationary_truth) false_10[["estimate"]] else NA_real_,
    false_selection_delta10_mcse = if (stationary_truth) false_10[["mcse"]] else NA_real_,
    false_selection_delta10_n = if (stationary_truth) false_10[["denominator"]] else NA_real_,
    detection_power_delta2 = if (sigma_truth) false_2[["estimate"]] else NA_real_,
    detection_power_delta2_mcse = if (sigma_truth) false_2[["mcse"]] else NA_real_,
    detection_power_delta2_n = if (sigma_truth) false_2[["denominator"]] else NA_real_,
    detection_power_delta10 = if (sigma_truth) false_10[["estimate"]] else NA_real_,
    detection_power_delta10_mcse = if (sigma_truth) false_10[["mcse"]] else NA_real_,
    detection_power_delta10_n = if (sigma_truth) false_10[["denominator"]] else NA_real_,
    correct_Msigma_delta2 = if (sigma_truth) correct_2[["estimate"]] else NA_real_,
    correct_Msigma_delta2_mcse = if (sigma_truth) correct_2[["mcse"]] else NA_real_,
    correct_Msigma_delta2_n = if (sigma_truth) correct_2[["denominator"]] else NA_real_,
    slope_target = slope_target,
    slope_bias = bias,
    slope_bias_mcse = bias_mcse,
    slope_rmse = rmse,
    slope_rmse_mcse = rmse_mcse,
    slope_n = length(errors),
    boundary_contact_rate = boundary[["estimate"]],
    boundary_contact_mcse = boundary[["mcse"]],
    boundary_contact_n = boundary[["denominator"]],
    fit_completion_rate = completion[["estimate"]],
    fit_completion_mcse = completion[["mcse"]],
    fit_completion_n = completion[["denominator"]],
    selection_completion_rate = selection_completion[["estimate"]],
    selection_completion_mcse = selection_completion[["mcse"]],
    selection_completion_n = selection_completion[["denominator"]],
    stringsAsFactors = FALSE
  )
}

groups <- split(raw, interaction(raw$truth, raw$dependence, raw$scale_months, drop = TRUE))
summary <- do.call(rbind, lapply(groups, summarize_group))
summary <- summary[order(summary$truth, summary$dependence, summary$scale_months), ]
rownames(summary) <- NULL
utils::write.csv(summary, file.path(output_dir, "seasonal_selection_summary.csv"), row.names = FALSE)
print(summary)
