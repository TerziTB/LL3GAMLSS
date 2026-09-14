#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
reference_dir <- if (length(args) >= 1L) args[[1L]] else file.path(
  "analysis", "era5_results", "full_195001_202512"
)
fresh_dir <- if (length(args) >= 2L) args[[2L]] else file.path(
  "analysis", "era5_results", "final_full_195001_202512"
)
files <- c(
  "area_weighted_model_fractions.csv",
  "basin_candidate_metrics.csv",
  "basin_model_selection.csv",
  "basin_nspei_calibration.csv",
  "basin_nspei.csv",
  "cell_candidate_metrics.csv",
  "cell_model_selection.csv",
  "cell_nspei_calibration.csv",
  "fit_audit.csv",
  "selected_slope_summary.csv",
  "selection_status_counts.csv",
  "soil_moisture_nspei_best_lag.csv",
  "soil_moisture_nspei_dry_event_consistency.csv",
  "soil_moisture_nspei_lag_correlations.csv"
)

compare_file <- function(file) {
  reference <- utils::read.csv(file.path(reference_dir, file), check.names = FALSE)
  fresh <- utils::read.csv(file.path(fresh_dir, file), check.names = FALSE)
  structure_equal <- identical(dim(reference), dim(fresh)) && identical(names(reference), names(fresh))
  numeric_columns <- intersect(names(reference)[vapply(reference, is.numeric, logical(1))], names(fresh))
  nonnumeric_columns <- setdiff(names(reference), numeric_columns)
  numeric_na_equal <- structure_equal && all(vapply(numeric_columns, function(column) {
    identical(is.na(reference[[column]]), is.na(fresh[[column]]))
  }, logical(1)))
  numeric_difference <- if (structure_equal && length(numeric_columns)) {
    values <- unlist(Map(function(x, y) abs(x - y), reference[numeric_columns], fresh[numeric_columns]))
    if (all(is.na(values))) 0 else max(values, na.rm = TRUE)
  } else if (structure_equal) 0 else Inf
  text_equal <- structure_equal && all(vapply(nonnumeric_columns, function(column) {
    identical(as.character(reference[[column]]), as.character(fresh[[column]]))
  }, logical(1)))
  data.frame(
    file = file,
    structure_equal = structure_equal,
    numeric_missingness_equal = numeric_na_equal,
    nonnumeric_equal = text_equal,
    maximum_absolute_numeric_difference = numeric_difference,
    reproduced = structure_equal && numeric_na_equal && text_equal &&
      is.finite(numeric_difference) && numeric_difference <= 1e-10,
    stringsAsFactors = FALSE
  )
}

missing <- files[
  !file.exists(file.path(reference_dir, files)) |
    !file.exists(file.path(fresh_dir, files))
]
if (length(missing)) stop("Missing comparison file(s): ", paste(missing, collapse = ", "))
result <- do.call(rbind, lapply(files, compare_file))
utils::write.csv(result, file.path(fresh_dir, "fresh_run_reproducibility.csv"), row.names = FALSE)
print(result)
if (!all(result$reproduced)) stop("Fresh ERA5-Land outputs differ from the archived reference.")
