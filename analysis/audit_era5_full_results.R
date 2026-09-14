#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args)) args[[1L]] else file.path(
  "analysis", "era5_results", "final_full_195001_202512"
)

cell_selection <- read.csv(file.path(result_dir, "cell_model_selection.csv"))
cell_candidates <- read.csv(file.path(result_dir, "cell_candidate_metrics.csv"))
basin_selection <- read.csv(file.path(result_dir, "basin_model_selection.csv"))
basin_candidates <- read.csv(file.path(result_dir, "basin_candidate_metrics.csv"))
basin_nspei <- read.csv(file.path(result_dir, "basin_nspei.csv"))
area_fractions <- read.csv(file.path(result_dir, "area_weighted_model_fractions.csv"))

audit <- data.frame(
  metric = c(
    "unique_cells",
    "cell_scale_selections",
    "cell_candidate_model_groups",
    "cell_candidate_monthly_fits",
    "cell_candidate_monthly_fits_converged",
    "cell_candidate_monthly_fits_inference_ready",
    "cell_candidate_boundary_contacts",
    "cell_scale_fit_failures",
    "basin_candidate_model_groups",
    "basin_candidate_monthly_fits",
    "basin_candidate_monthly_fits_converged",
    "basin_candidate_monthly_fits_inference_ready",
    "basin_candidate_boundary_contacts"
  ),
  value = c(
    length(unique(cell_selection$cell_id)),
    nrow(cell_selection),
    nrow(cell_candidates),
    sum(cell_candidates$fitted_months),
    sum(cell_candidates$converged_months),
    sum(cell_candidates$inference_ready_months),
    sum(cell_candidates$boundary_contact_months),
    sum(cell_selection$selected_model == "fit_failure"),
    nrow(basin_candidates),
    sum(basin_candidates$fitted_months),
    sum(basin_candidates$converged_months),
    sum(basin_candidates$inference_ready_months),
    sum(basin_candidates$boundary_contact_months)
  )
)

status_counts <- as.data.frame(with(
  cell_selection,
  table(scale_months, selection_status)
))
names(status_counts) <- c("scale_months", "selection_status", "cell_count")
status_counts <- status_counts[status_counts$cell_count > 0, ]

cell_calibration <- do.call(rbind, lapply(
  split(cell_selection, cell_selection$scale_months),
  function(d) data.frame(
    scale_months = unique(d$scale_months),
    cells = nrow(d),
    mean_of_cell_nspei_means = mean(d$nspei_mean, na.rm = TRUE),
    minimum_cell_nspei_mean = min(d$nspei_mean, na.rm = TRUE),
    maximum_cell_nspei_mean = max(d$nspei_mean, na.rm = TRUE),
    mean_of_cell_nspei_sds = mean(d$nspei_sd, na.rm = TRUE),
    minimum_cell_nspei_sd = min(d$nspei_sd, na.rm = TRUE),
    maximum_cell_nspei_sd = max(d$nspei_sd, na.rm = TRUE)
  )
))

basin_calibration <- do.call(rbind, lapply(
  split(basin_nspei, basin_nspei$scale_months),
  function(d) data.frame(
    scale_months = unique(d$scale_months),
    selected_model = unique(d$selected_model),
    finite_months = sum(is.finite(d$nSPEI)),
    nspei_mean = mean(d$nSPEI, na.rm = TRUE),
    nspei_sd = sd(d$nSPEI, na.rm = TRUE),
    nspei_minimum = min(d$nSPEI, na.rm = TRUE),
    nspei_maximum = max(d$nSPEI, na.rm = TRUE),
    extreme_dry_months = sum(d$nSPEI <= -2, na.rm = TRUE)
  )
))

slope_rows <- cell_selection$selected_model != "stationary" &
  cell_selection$selected_model != "fit_failure"
slope_data <- cell_selection[slope_rows, ]
slope_groups <- interaction(
  slope_data$scale_months,
  slope_data$selected_model,
  drop = TRUE
)
slope_summary <- do.call(rbind, lapply(split(slope_data, slope_groups), function(d) {
  summarize_slope <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(NA_real_, NA_real_, NA_real_))
    c(stats::median(x), stats::quantile(x, 0.1), stats::quantile(x, 0.9))
  }
  sigma <- summarize_slope(d$selected_sigma_slope_median)
  nu <- summarize_slope(d$selected_nu_slope_median)
  data.frame(
    scale_months = unique(d$scale_months),
    selected_model = unique(d$selected_model),
    cells = nrow(d),
    median_delta_AICc = median(d$delta_AICc_favoring_best_nonstationary),
    sigma_slope_median = sigma[1],
    sigma_slope_p10 = sigma[2],
    sigma_slope_p90 = sigma[3],
    nu_slope_median = nu[1],
    nu_slope_p10 = nu[2],
    nu_slope_p90 = nu[3]
  )
}))

write.csv(audit, file.path(result_dir, "fit_audit.csv"), row.names = FALSE)
write.csv(status_counts, file.path(result_dir, "selection_status_counts.csv"), row.names = FALSE)
write.csv(cell_calibration, file.path(result_dir, "cell_nspei_calibration.csv"), row.names = FALSE)
write.csv(basin_calibration, file.path(result_dir, "basin_nspei_calibration.csv"), row.names = FALSE)
write.csv(slope_summary, file.path(result_dir, "selected_slope_summary.csv"), row.names = FALSE)

cat("Fit audit:\n")
print(audit, row.names = FALSE)
cat("\nArea-weighted selected-model fractions (%):\n")
print(area_fractions[, c("scale_months", "selected_model", "cell_count", "basin_area_percent")], row.names = FALSE)
cat("\nBasin selections:\n")
print(basin_selection[, c(
  "scale_months", "selected_model", "selection_status",
  "delta_AICc_favoring_best_nonstationary"
)], row.names = FALSE)
cat("\nBasin nSPEI calibration:\n")
print(basin_calibration, row.names = FALSE)
cat("\nCell selection status counts:\n")
print(status_counts, row.names = FALSE)
