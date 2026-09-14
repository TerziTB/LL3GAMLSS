#!/usr/bin/env Rscript

result_dir <- if (length(commandArgs(trailingOnly = TRUE))) {
  commandArgs(trailingOnly = TRUE)[[1L]]
} else {
  file.path("analysis", "era5_results", "final_full_195001_202512")
}

validation <- read.csv(file.path("analysis", "era5_results", "era5_basin_validation_monthly.csv"))
validation$date <- as.Date(validation$date)
basin_nspei <- read.csv(file.path(result_dir, "basin_nspei.csv"))
basin_nspei$date <- as.Date(basin_nspei$date)

calendar_month_standardize <- function(values, dates) {
  month <- as.integer(format(dates, "%m"))
  result <- rep(NA_real_, length(values))
  for (current_month in 1:12) {
    rows <- month == current_month & is.finite(values)
    result[rows] <- as.numeric(scale(values[rows]))
  }
  result
}

validation$root_zone_soil_water_anomaly <- calendar_month_standardize(
  validation$root_zone_soil_water_m3_m3,
  validation$date
)

merged <- merge(
  basin_nspei,
  validation[, c("date", "root_zone_soil_water_anomaly")],
  by = "date",
  all.x = TRUE,
  sort = TRUE
)

lag_results <- do.call(rbind, lapply(sort(unique(merged$scale_months)), function(scale_months) {
  data <- merged[merged$scale_months == scale_months, ]
  do.call(rbind, lapply(0:6, function(lag_months) {
    response <- if (lag_months == 0L) {
      data$root_zone_soil_water_anomaly
    } else {
      c(data$root_zone_soil_water_anomaly[-seq_len(lag_months)], rep(NA_real_, lag_months))
    }
    keep <- is.finite(data$nSPEI) & is.finite(response)
    data.frame(
      scale_months = scale_months,
      soil_response_lag_months = lag_months,
      n = sum(keep),
      pearson_correlation = stats::cor(data$nSPEI[keep], response[keep], method = "pearson"),
      spearman_correlation = stats::cor(data$nSPEI[keep], response[keep], method = "spearman")
    )
  }))
}))

event_results <- do.call(rbind, lapply(sort(unique(merged$scale_months)), function(scale_months) {
  data <- merged[merged$scale_months == scale_months, ]
  dry <- is.finite(data$nSPEI) & data$nSPEI <= -1 &
    is.finite(data$root_zone_soil_water_anomaly)
  data.frame(
    scale_months = scale_months,
    dry_nspei_months = sum(dry),
    mean_soil_anomaly_during_dry_nspei = mean(data$root_zone_soil_water_anomaly[dry]),
    percent_with_negative_soil_anomaly = 100 * mean(data$root_zone_soil_water_anomaly[dry] < 0)
  )
}))

best_lag <- do.call(rbind, lapply(split(lag_results, lag_results$scale_months), function(data) {
  data[which.max(data$pearson_correlation), , drop = FALSE]
}))

write.csv(lag_results, file.path(result_dir, "soil_moisture_nspei_lag_correlations.csv"), row.names = FALSE)
write.csv(event_results, file.path(result_dir, "soil_moisture_nspei_dry_event_consistency.csv"), row.names = FALSE)
write.csv(best_lag, file.path(result_dir, "soil_moisture_nspei_best_lag.csv"), row.names = FALSE)

cat("Best root-zone soil-moisture response lag by nSPEI scale:\n")
print(best_lag)
cat("\nDry-event consistency:\n")
print(event_results)
