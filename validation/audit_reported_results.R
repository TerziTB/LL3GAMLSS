#!/usr/bin/env Rscript

source(file.path("validation", "validation_helpers.R"))
project_dir <- resolve_project_dir()
setwd(project_dir)

checks <- list()
add_check <- function(name, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    observed = paste(observed, collapse = ";"),
    expected = paste(expected, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

station <- utils::read.csv(file.path("analysis", "seyhan_results", "selected_model_summary.csv"))
station <- station[as.character(station$station) != "17934", , drop = FALSE]
strong <- station[station$inference_status == "strong_support", , drop = FALSE]
add_check("station_scale_rows", nrow(station) == 28L, nrow(station), 28L)
add_check("retained_stations", length(unique(station$station)) == 7L, length(unique(station$station)), 7L)
add_check("strong_station_scale_cases", nrow(strong) == 17L, nrow(strong), 17L)
add_check(
  "strong_cases_select_scale",
  all(strong$best_nonstationary_model == "sigma_time"),
  paste(unique(strong$best_nonstationary_model), collapse = ","),
  "sigma_time"
)

dependence_dir <- file.path("analysis", "seyhan_results", "dependence_bootstrap")
pit <- utils::read.csv(file.path(dependence_dir, "strong_case_PIT_summary.csv"))
monthly_pit <- utils::read.csv(file.path(dependence_dir, "strong_case_monthly_PIT_diagnostics.csv"))
bootstrap <- utils::read.csv(file.path(dependence_dir, "representative_moving_block_bootstrap.csv"))
add_check("PIT_case_rows", nrow(pit) == 17L, nrow(pit), 17L)
add_check("PIT_monthly_rows", nrow(monthly_pit) == 204L, nrow(monthly_pit), 204L)
add_check("bootstrap_cases", nrow(bootstrap) == 4L, nrow(bootstrap), 4L)
add_check("bootstrap_refits", all(bootstrap$B_valid == 499L), bootstrap$B_valid, rep(499L, 4L))
add_check(
  "excluded_station_absent",
  !any(as.character(c(pit$station, monthly_pit$station, bootstrap$station)) == "17934"),
  any(as.character(c(pit$station, monthly_pit$station, bootstrap$station)) == "17934"),
  FALSE
)

seasonal_raw <- utils::read.csv(file.path(
  "validation", "results", "seasonal_selection", "seasonal_selection_raw.csv"
))
seasonal_summary <- utils::read.csv(file.path(
  "validation", "results", "seasonal_selection", "seasonal_selection_summary.csv"
))
group_counts <- table(seasonal_raw$truth, seasonal_raw$dependence, seasonal_raw$scale_months)
add_check("seasonal_raw_rows", nrow(seasonal_raw) == 480L, nrow(seasonal_raw), 480L)
add_check("seasonal_group_denominators", all(group_counts == 30L), range(group_counts), "30 per cell")
add_check(
  "seasonal_monthly_record_length",
  all(seasonal_raw$total_monthly_observations == 792L),
  unique(seasonal_raw$total_monthly_observations),
  792L
)
for (i in seq_len(nrow(seasonal_summary))) {
  s <- seasonal_summary[i, ]
  d <- seasonal_raw[
    seasonal_raw$truth == s$truth & seasonal_raw$dependence == s$dependence &
      seasonal_raw$scale_months == s$scale_months,
    , drop = FALSE
  ]
  if (s$truth == "M0") {
    observed <- mean(d$selected_delta2 != "M0", na.rm = TRUE)
    expected <- s$false_selection_delta2
  } else {
    observed <- mean(d$selected_delta2 != "M0", na.rm = TRUE)
    expected <- s$detection_power_delta2
  }
  add_check(
    sprintf("seasonal_rate_%s_%s_%02d", s$truth, s$dependence, s$scale_months),
    isTRUE(all.equal(observed, expected, tolerance = 1e-12)),
    format(observed, digits = 16),
    format(expected, digits = 16)
  )
}

spei <- utils::read.csv(file.path(
  "validation", "results", "spei_maxlik_equivalence", "spei_maxlik_equivalence_raw.csv"
))
add_check("SPEI_equivalence_rows", nrow(spei) == 39L, nrow(spei), 39L)
add_check("SPEI_simulated_rows", sum(spei$source == "simulated") == 36L, sum(spei$source == "simulated"), 36L)
add_check("SPEI_empirical_rows", sum(spei$source == "empirical") == 3L, sum(spei$source == "empirical"), 3L)

era_audit <- utils::read.csv(file.path(
  "analysis", "era5_results", "final_full_195001_202512", "fit_audit.csv"
))
era_value <- setNames(era_audit$value, era_audit$metric)
era_expected <- c(
  unique_cells = 274L,
  cell_scale_selections = 1096L,
  cell_candidate_monthly_fits = 52608L,
  cell_candidate_monthly_fits_converged = 52608L,
  cell_candidate_monthly_fits_inference_ready = 52608L,
  cell_candidate_boundary_contacts = 0L,
  cell_scale_fit_failures = 0L,
  basin_candidate_monthly_fits = 192L,
  basin_candidate_monthly_fits_converged = 192L,
  basin_candidate_monthly_fits_inference_ready = 192L,
  basin_candidate_boundary_contacts = 0L
)
for (name in names(era_expected)) {
  add_check(
    paste0("ERA5_", name),
    era_value[[name]] == era_expected[[name]],
    era_value[[name]],
    era_expected[[name]]
  )
}
era_reproduction <- utils::read.csv(file.path(
  "analysis", "era5_results", "final_full_195001_202512", "fresh_run_reproducibility.csv"
))
add_check(
  "ERA5_fresh_run_reproduced",
  all(era_reproduction$reproduced),
  sum(era_reproduction$reproduced),
  nrow(era_reproduction)
)

figure_paths <- c(
  file.path("manuscript", "figures", c(
    "Fig3_simulation_validation.png",
    "Fig4_model_support_heatmap.png",
    "Fig5_spei12_comparison.png"
  )),
  file.path("analysis", "era5_results", "final_full_195001_202512", "era5_all_scales_selected_models.png")
)
add_check("manuscript_figures_exist", all(file.exists(figure_paths)), sum(file.exists(figure_paths)), length(figure_paths))

result <- do.call(rbind, checks)
output_dir <- file.path("validation", "results", "reported_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, file.path(output_dir, "reported_results_audit.csv"), row.names = FALSE)
writeLines(c(
  "# Reported-results audit",
  "",
  sprintf("Checks passed: %d of %d.", sum(result$passed), nrow(result)),
  "",
  "The audit reconciles manuscript station, seasonal-simulation, SPEI-equivalence, bootstrap, ERA5-Land, and figure counts with their archived generated files."
), file.path(output_dir, "REPORTED_RESULTS_AUDIT.md"))
write_session_info(file.path(output_dir, "sessionInfo.txt"))
print(result)
if (!all(result$passed)) stop("Reported-results audit failed.")
