project_dir <- local({
  source(file.path("validation", "validation_helpers.R"), local = TRUE)
  resolve_project_dir()
})
setwd(project_dir)

output_dir <- file.path(project_dir, "analysis", "manuscript_inputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  station_metadata = file.path("analysis", "seyhan_results", "station_metadata.csv"),
  station_summary = file.path("analysis", "seyhan_results", "station_summary.csv"),
  scale_summary = file.path("analysis", "seyhan_results", "scale_summary.csv"),
  station_selection = file.path("analysis", "seyhan_results", "selected_model_summary.csv"),
  stationary_benchmark = file.path("analysis", "seyhan_results", "stationary_benchmark.csv"),
  station_values = file.path("analysis", "seyhan_results", "spei_values.csv"),
  era_area = file.path("analysis", "era5_results", "final_full_195001_202512", "area_weighted_model_fractions.csv"),
  era_basin = file.path("analysis", "era5_results", "final_full_195001_202512", "basin_model_selection.csv"),
  era_fit = file.path("analysis", "era5_results", "final_full_195001_202512", "fit_audit.csv"),
  era_soil = file.path("analysis", "era5_results", "final_full_195001_202512", "soil_moisture_nspei_best_lag.csv")
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop(
    "Required archived analysis outputs are missing: ",
    paste(unname(missing), collapse = ", "),
    ". Regenerate restricted station outputs locally before rebuilding these inputs."
  )
}

excluded_stations <- "17934"

copy_csv <- function(name) {
  value <- utils::read.csv(required[[name]], stringsAsFactors = FALSE)
  if ("station" %in% names(value)) {
    value <- value[!as.character(value$station) %in% excluded_stations, , drop = FALSE]
  }
  utils::write.csv(value, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  value
}

invisible(lapply(c(
  "station_metadata", "station_summary", "scale_summary", "station_selection",
  "stationary_benchmark", "era_area", "era_basin", "era_fit", "era_soil"
), copy_csv))

station_values <- utils::read.csv(required[["station_values"]], stringsAsFactors = FALSE)
station_values <- station_values[
  !as.character(station_values$station) %in% excluded_stations,
  , drop = FALSE
]
figure_values <- station_values[
  as.character(station_values$station) %in% c("17351", "17802", "17837", "17840") &
    station_values$scale_months == 12L,
  c("station", "date", "scale_months", "LL3_stationary", "LL3_selected_nonstationary")
]
utils::write.csv(
  figure_values,
  file.path(output_dir, "figure5_stationary_nonstationary.csv"),
  row.names = FALSE
)

files <- list.files(output_dir, pattern = "\\.csv$", full.names = FALSE)
files <- setdiff(files, "manifest.csv")
source_lookup <- c(
  station_metadata.csv = required[["station_metadata"]],
  station_summary.csv = required[["station_summary"]],
  scale_summary.csv = required[["scale_summary"]],
  station_selection.csv = required[["station_selection"]],
  stationary_benchmark.csv = required[["stationary_benchmark"]],
  era_area.csv = required[["era_area"]],
  era_basin.csv = required[["era_basin"]],
  era_fit.csv = required[["era_fit"]],
  era_soil.csv = required[["era_soil"]],
  figure5_stationary_nonstationary.csv =
    "analysis/seyhan_results/spei_values.csv (date and index columns only)"
)
utils::write.csv(
  data.frame(file = files, source = unname(source_lookup[files]), stringsAsFactors = FALSE),
  file.path(output_dir, "manifest.csv"),
  row.names = FALSE
)

cat("Created compact manuscript inputs in", normalizePath(output_dir, winslash = "/"), "\n")
