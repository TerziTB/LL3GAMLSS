#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ncdf4))

main_file <- file.path("analysis", "era5_input", "final_main", "data_stream-moda.nc")
prepared_file <- file.path("analysis", "era5_results", "era5_basin_monthly.rds")
output_dir <- file.path("analysis", "era5_results")

prepared <- readRDS(prepared_file)
nc <- nc_open(main_file)
on.exit(nc_close(nc), add = TRUE)
main_tp <- ncvar_get(nc, "tp")
seconds <- ncvar_get(nc, "valid_time")
dates <- as.Date(as.POSIXct(seconds, origin = "1970-01-01", tz = "UTC"))

first <- as.Date(format(dates, "%Y-%m-01"))
next_first <- as.Date(format(first + 35, "%Y-%m-01"))
month_days <- as.integer(next_first - first)
ncell <- dim(main_tp)[[1L]] * dim(main_tp)[[2L]]
main_monthly_mm <- t(matrix(main_tp, nrow = ncell, ncol = length(dates))) *
  1000 * month_days

affected <- dates >= as.Date("2022-09-01") & dates <= as.Date("2024-02-01")
cell_ids <- prepared$cell_metadata$cell_id
weights <- prepared$cell_metadata$basin_area_weight_km2
original <- main_monthly_mm[affected, cell_ids, drop = FALSE]
corrected <- prepared$precipitation_mm[affected, , drop = FALSE]

weighted_mean <- function(values) {
  apply(values, 1L, stats::weighted.mean, w = weights)
}

comparison <- data.frame(
  month = format(dates[affected], "%Y-%m"),
  original_basin_precipitation_mm = weighted_mean(original),
  corrected_basin_precipitation_mm = weighted_mean(corrected)
)
comparison$difference_mm <-
  comparison$corrected_basin_precipitation_mm - comparison$original_basin_precipitation_mm

write.csv(
  comparison,
  file.path(output_dir, "era5_precipitation_correction_comparison.csv"),
  row.names = FALSE
)

difference <- corrected - original
writeLines(c(
  paste("Affected months:", nrow(comparison)),
  paste("Affected basin-cell values:", length(difference)),
  paste("Values differing by more than 0.01 mm:", sum(abs(difference) > 0.01)),
  sprintf("Cell-month correction range: %.3f to %.3f mm", min(difference), max(difference)),
  sprintf("Median cell-month correction: %.3f mm", stats::median(difference)),
  sprintf(
    "Basin-month correction range: %.3f to %.3f mm",
    min(comparison$difference_mm),
    max(comparison$difference_mm)
  )
), file.path(output_dir, "era5_precipitation_correction_qc.txt"))

print(comparison)
