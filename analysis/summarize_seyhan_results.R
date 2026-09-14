#!/usr/bin/env Rscript

output_dir <- file.path("analysis", "seyhan_results")
summary <- read.csv(file.path(output_dir, "selected_model_summary.csv"), stringsAsFactors = FALSE)
models <- read.csv(file.path(output_dir, "model_comparison.csv"), stringsAsFactors = FALSE)
values <- read.csv(file.path(output_dir, "spei_values.csv"), stringsAsFactors = FALSE)
quality <- read.csv(file.path(output_dir, "data_quality.csv"), stringsAsFactors = FALSE)
metadata <- read.csv(file.path(output_dir, "station_metadata.csv"), stringsAsFactors = FALSE)

screen_path <- file.path(output_dir, "annual_repetition_screen.csv")
excluded_stations <- character()
if (file.exists(screen_path)) {
  repetition_screen <- read.csv(screen_path, stringsAsFactors = FALSE)
  excluded_flag <- tolower(as.character(repetition_screen$excluded)) %in% c("true", "1")
  excluded_stations <- repetition_screen$station[excluded_flag]
  summary <- summary[!summary$station %in% excluded_stations, , drop = FALSE]
  models <- models[!models$station %in% excluded_stations, , drop = FALSE]
  values <- values[!values$station %in% excluded_stations, , drop = FALSE]
  quality <- quality[!quality$station %in% excluded_stations, , drop = FALSE]
  metadata <- metadata[!metadata$station %in% excluded_stations, , drop = FALSE]
}

benchmark_rows <- list()
for (station in unique(values$station)) {
  for (scale_months in sort(unique(values$scale_months))) {
    x <- values[values$station == station & values$scale_months == scale_months, ]
    keep <- is.finite(x$SPEI_PWM_stationary) & is.finite(x$LL3_stationary)
    benchmark_rows[[length(benchmark_rows) + 1L]] <- data.frame(
      station = station,
      scale_months = scale_months,
      n = sum(keep),
      correlation_PWM_vs_LL3_MLE = cor(x$SPEI_PWM_stationary[keep], x$LL3_stationary[keep]),
      mean_absolute_difference = mean(abs(x$SPEI_PWM_stationary[keep] - x$LL3_stationary[keep])),
      max_absolute_difference = max(abs(x$SPEI_PWM_stationary[keep] - x$LL3_stationary[keep])),
      PWM_extreme_dry_count = sum(x$SPEI_PWM_stationary <= -2, na.rm = TRUE),
      LL3_MLE_extreme_dry_count = sum(x$LL3_stationary <= -2, na.rm = TRUE)
    )
  }
}
benchmark <- do.call(rbind, benchmark_rows)
write.csv(benchmark, file.path(output_dir, "stationary_benchmark.csv"), row.names = FALSE)

station_summary_rows <- lapply(unique(summary$station), function(station) {
  x <- summary[summary$station == station, ]
  data.frame(
    station = station,
    station_name = metadata$station_name[match(station, metadata$station)],
    scales_with_strong_support = sum(x$inference_status == "strong_support"),
    scales_without_material_support = sum(x$inference_status == "no_material_support"),
    median_delta_AICc = median(x$delta_AICc_favoring_nonstationary),
    maximum_delta_AICc = max(x$delta_AICc_favoring_nonstationary),
    median_stationary_nonstationary_MAE = median(x$mean_absolute_difference),
    maximum_index_difference = max(x$max_absolute_difference),
    total_stationary_extreme_dry = sum(x$stationary_extreme_dry_count),
    total_nonstationary_extreme_dry = sum(x$nonstationary_extreme_dry_count)
  )
})
station_summary <- do.call(rbind, station_summary_rows)
write.csv(station_summary, file.path(output_dir, "station_summary.csv"), row.names = FALSE)

scale_summary_rows <- lapply(sort(unique(summary$scale_months)), function(scale_months) {
  x <- summary[summary$scale_months == scale_months, ]
  data.frame(
    scale_months = scale_months,
    station_count = nrow(x),
    strong_support_count = sum(x$inference_status == "strong_support"),
    no_material_support_count = sum(x$inference_status == "no_material_support"),
    median_delta_AICc = median(x$delta_AICc_favoring_nonstationary),
    maximum_delta_AICc = max(x$delta_AICc_favoring_nonstationary),
    median_index_MAE = median(x$mean_absolute_difference),
    maximum_index_difference = max(x$max_absolute_difference)
  )
})
scale_summary <- do.call(rbind, scale_summary_rows)
write.csv(scale_summary, file.path(output_dir, "scale_summary.csv"), row.names = FALSE)

raw <- jsonlite::fromJSON(
  file.path("artifact_work", "inspection", "workbook_values.json"),
  simplifyVector = FALSE
)
issue_rows <- list()
for (station in names(raw)) {
  rows <- raw[[station]]
  header <- unlist(rows[[1]], use.names = FALSE)
  matrix <- do.call(rbind, lapply(rows[-1], function(x) unlist(x, use.names = FALSE)))
  d <- as.data.frame(matrix, stringsAsFactors = FALSE)
  names(d) <- header
  d$Date <- as.Date(as.numeric(d$Date), origin = "1899-12-30")
  for (nm in setdiff(header, "Date")) d[[nm]] <- as.numeric(d[[nm]])
  bad <- d$Min_Temp > d$Max_Temp |
    d$Average_Temp < d$Min_Temp |
    d$Average_Temp > d$Max_Temp
  if (any(bad)) {
    issue <- d[bad, ]
    issue$station <- station
    issue$issue <- ifelse(
      issue$Min_Temp > issue$Max_Temp,
      "minimum_temperature_above_maximum",
      "mean_temperature_outside_minimum_maximum"
    )
    issue$inside_common_analysis_period <-
      issue$Date >= as.Date(quality$analysis_common_start[1]) &
      issue$Date <= as.Date(quality$analysis_common_end[1])
    issue_rows[[length(issue_rows) + 1L]] <- issue[c(
      "station", "Date", "Precipitation", "Average_Temp", "Min_Temp",
      "Max_Temp", "issue", "inside_common_analysis_period"
    )]
  }
}
issues <- do.call(rbind, issue_rows)
issues <- issues[!issues$station %in% excluded_stations, , drop = FALSE]
write.csv(issues, file.path(output_dir, "data_quality_issues.csv"), row.names = FALSE)

selected_model_rows <- merge(
  summary[c("station", "scale_months", "best_nonstationary_model")],
  models,
  by.x = c("station", "scale_months", "best_nonstationary_model"),
  by.y = c("station", "scale_months", "model"),
  all.x = TRUE
)

overview <- data.frame(
  metric = c(
    "Stations analysed",
    "Accumulation scales",
    "Station-scale comparisons",
    "Strong nonstationary support",
    "No material nonstationary support",
    "Selected models with all 12 calendar months converged",
    "Selected models with threshold-boundary contact",
    "Rejected candidate model rows with threshold-boundary contact",
    "Median stationary/nonstationary correlation",
    "Median stationary/nonstationary mean absolute difference",
    "Maximum observed stationary/nonstationary difference",
    "Median stationary LL3-MLE / SPEI-PWM correlation"
  ),
  value = c(
    length(unique(summary$station)),
    paste(sort(unique(summary$scale_months)), collapse = ", "),
    nrow(summary),
    sum(summary$inference_status == "strong_support"),
    sum(summary$inference_status == "no_material_support"),
    sum(selected_model_rows$converged_months == 12L),
    sum(selected_model_rows$boundary_contact_months > 0L),
    sum(models$boundary_contact_months > 0L),
    format(median(summary$correlation_stationary_nonstationary), digits = 4),
    format(median(summary$mean_absolute_difference), digits = 4),
    format(max(summary$max_absolute_difference), digits = 4),
    format(median(benchmark$correlation_PWM_vs_LL3_MLE), digits = 4)
  ),
  stringsAsFactors = FALSE
)
write.csv(overview, file.path(output_dir, "overview.csv"), row.names = FALSE)

diagnostics <- read.csv(file.path(output_dir, "index_diagnostics.csv"), stringsAsFactors = FALSE)
diagnostics <- diagnostics[!diagnostics$station %in% excluded_stations, , drop = FALSE]
workbook_data <- list(
  overview = overview,
  station_summary = station_summary,
  scale_summary = scale_summary,
  station_scale_results = summary,
  model_comparison = models,
  index_diagnostics = diagnostics,
  monthly_spei = values,
  stationary_benchmark = benchmark,
  data_quality = quality,
  quality_issues = issues,
  station_metadata = metadata
)
jsonlite::write_json(
  workbook_data,
  file.path(output_dir, "workbook_data.json"),
  dataframe = "rows",
  na = "null",
  digits = NA,
  auto_unbox = TRUE
)

cat("Created result summaries in", normalizePath(output_dir, winslash = "/"), "\n")
