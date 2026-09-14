#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
main_file <- if (length(args) >= 1L) args[[1L]] else file.path(
  "analysis", "era5_input", "era5_land_moda_1950_2025.nc"
)
patch_file <- if (length(args) >= 2L && nzchar(args[[2L]])) args[[2L]] else NA_character_
output_dir <- file.path("analysis", "era5_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(main_file)) {
  stop(
    "ERA5-Land NetCDF input is not included in the repository. Supply the main ",
    "monthly-means file as the first command-line argument; optionally supply the ",
    "2022-2024 hourly-monthly precipitation correction as the second argument."
  )
}
if (!is.na(patch_file) && !file.exists(patch_file)) {
  stop("The supplied ERA5-Land precipitation-correction file does not exist: ", patch_file)
}

.libPaths(c(file.path(getwd(), ".rlib-realdata"), .libPaths()))
suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(SPEI)
})

read_monthly_file <- function(path, variables) {
  nc <- nc_open(path)
  on.exit(nc_close(nc), add = TRUE)
  available <- names(nc$var)
  missing_variables <- setdiff(variables, available)
  if (length(missing_variables)) {
    stop("Missing NetCDF variables: ", paste(missing_variables, collapse = ", "))
  }
  seconds <- ncvar_get(nc, "valid_time")
  dates <- as.Date(as.POSIXct(seconds, origin = "1970-01-01", tz = "UTC"))
  list(
    dates = dates,
    longitude = ncvar_get(nc, "longitude"),
    latitude = ncvar_get(nc, "latitude"),
    values = setNames(lapply(variables, function(v) ncvar_get(nc, v)), variables),
    units = setNames(lapply(variables, function(v) ncatt_get(nc, v, "units")$value), variables),
    expver = if ("expver" %in% available) ncvar_get(nc, "expver") else NA
  )
}

netcdf_variable_names <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc), add = TRUE)
  names(nc$var)
}

days_in_month <- function(dates) {
  first <- as.Date(format(dates, "%Y-%m-01"))
  next_first <- as.Date(format(first + 35, "%Y-%m-01"))
  as.integer(next_first - first)
}

weighted_row_mean <- function(x, weights) {
  apply(x, 1L, function(row) {
    keep <- is.finite(row) & is.finite(weights) & weights > 0
    if (!any(keep)) return(NA_real_)
    stats::weighted.mean(row[keep], weights[keep])
  })
}

optional_variables <- intersect(
  c("sd", "swvl1", "swvl2", "swvl3"),
  netcdf_variable_names(main_file)
)
main <- read_monthly_file(main_file, c("t2m", "tp", optional_variables))
expected <- seq(as.Date("1950-01-01"), as.Date("2025-12-01"), by = "month")
if (length(main$dates) != length(expected) || any(unclass(main$dates) != unclass(expected))) {
  stop("The main NetCDF does not contain the expected complete 1950-2025 monthly sequence.")
}
if (!identical(main$units$t2m, "K") || !identical(main$units$tp, "m")) {
  stop("Unexpected units: t2m=", main$units$t2m, ", tp=", main$units$tp)
}

nlongitude <- length(main$longitude)
nlatitude <- length(main$latitude)
ntime <- length(main$dates)
ncell <- nlongitude * nlatitude

temperature_all <- t(matrix(main$values$t2m, nrow = ncell, ncol = ntime)) - 273.15
precipitation_daily_all <- t(matrix(main$values$tp, nrow = ncell, ncol = ntime))
precipitation_all <- precipitation_daily_all * 1000 * days_in_month(main$dates)

patch_status <- "not supplied; affected precipitation months set to missing"
affected <- main$dates >= as.Date("2022-09-01") & main$dates <= as.Date("2024-02-01")
if (!is.na(patch_file) && nzchar(patch_file) && file.exists(patch_file)) {
  patch <- read_monthly_file(patch_file, "tp")
  if (!isTRUE(all.equal(main$longitude, patch$longitude)) ||
      !isTRUE(all.equal(main$latitude, patch$latitude))) {
    stop("The correction file grid does not match the main NetCDF grid.")
  }
  patch_cells <- length(patch$longitude) * length(patch$latitude)
  patch_matrix <- t(matrix(patch$values$tp, nrow = patch_cells, ncol = length(patch$dates)))
  main_month_keys <- format(main$dates[affected], "%Y-%m")
  patch_month_keys <- format(patch$dates, "%Y-%m")
  if (anyDuplicated(patch_month_keys)) stop("The correction file contains duplicate calendar months.")
  patch_rows <- match(main_month_keys, patch_month_keys)
  if (anyNA(patch_rows)) stop("The correction file does not cover every affected month.")
  precipitation_all[affected, ] <-
    patch_matrix[patch_rows, , drop = FALSE] * 1000 * days_in_month(main$dates[affected])
  patch_status <- paste("corrected with", normalizePath(patch_file, winslash = "/"))
} else {
  precipitation_all[affected, ] <- NA_real_
}

grid <- rast(
  ncols = nlongitude,
  nrows = nlatitude,
  xmin = min(main$longitude) - 0.05,
  xmax = max(main$longitude) + 0.05,
  ymin = min(main$latitude) - 0.05,
  ymax = max(main$latitude) + 0.05,
  crs = "EPSG:4326"
)
values(grid) <- as.vector(main$values$t2m[, , 1L])

basin_path <- file.path("analysis", "era5_input", "seyhan_hydrobasins_lev08.gpkg")
basin <- vect(basin_path)
overlap <- extract(grid, basin, cells = TRUE, weights = TRUE, exact = TRUE)
names(overlap)[names(overlap) == names(grid)] <- "first_temperature_k"
overlap <- overlap[is.finite(overlap$first_temperature_k) & overlap$weight > 0, ]
overlap <- overlap[!duplicated(overlap$cell), ]

cell_area <- values(cellSize(grid, unit = "km"), mat = FALSE)
centres <- xyFromCell(grid, overlap$cell)
cell_metadata <- data.frame(
  cell_id = overlap$cell,
  longitude = centres[, 1L],
  latitude = centres[, 2L],
  basin_fraction = overlap$weight,
  cell_area_km2 = cell_area[overlap$cell],
  basin_area_weight_km2 = overlap$weight * cell_area[overlap$cell]
)
cell_metadata$cell_key <- sprintf("cell_%04d", cell_metadata$cell_id)

temperature <- temperature_all[, cell_metadata$cell_id, drop = FALSE]
precipitation <- precipitation_all[, cell_metadata$cell_id, drop = FALSE]

pet <- matrix(NA_real_, nrow = ntime, ncol = nrow(cell_metadata))
start <- c(as.integer(format(main$dates[[1L]], "%Y")), as.integer(format(main$dates[[1L]], "%m")))
for (i in seq_len(nrow(cell_metadata))) {
  pet[, i] <- as.numeric(SPEI::thornthwaite(
    Tave = stats::ts(temperature[, i], start = start, frequency = 12),
    lat = cell_metadata$latitude[[i]],
    na.rm = FALSE,
    verbose = FALSE
  ))
}
water_balance <- precipitation - pet

validation <- list()
if (all(c("swvl1", "swvl2", "swvl3") %in% names(main$values))) {
  soil_layers <- lapply(c("swvl1", "swvl2", "swvl3"), function(variable) {
    all_cells <- t(matrix(main$values[[variable]], nrow = ncell, ncol = ntime))
    all_cells[, cell_metadata$cell_id, drop = FALSE]
  })
  validation$root_zone_soil_water_m3_m3 <-
    0.07 * soil_layers[[1L]] + 0.21 * soil_layers[[2L]] + 0.72 * soil_layers[[3L]]
}
if ("sd" %in% names(main$values)) {
  snow_all <- t(matrix(main$values$sd, nrow = ncell, ncol = ntime))
  validation$snow_water_equivalent_mm <- pmax(
    snow_all[, cell_metadata$cell_id, drop = FALSE] * 1000,
    0
  )
}

weights <- cell_metadata$basin_area_weight_km2
basin_monthly <- data.frame(
  date = main$dates,
  precipitation_mm = weighted_row_mean(precipitation, weights),
  mean_temperature_c = weighted_row_mean(temperature, weights),
  PET_Thornthwaite_mm = weighted_row_mean(pet, weights),
  water_balance_mm = weighted_row_mean(water_balance, weights)
)

basin_validation <- data.frame(date = main$dates)
if ("root_zone_soil_water_m3_m3" %in% names(validation)) {
  basin_validation$root_zone_soil_water_m3_m3 <- weighted_row_mean(
    validation$root_zone_soil_water_m3_m3,
    weights
  )
}
if ("snow_water_equivalent_mm" %in% names(validation)) {
  basin_validation$snow_water_equivalent_mm <- weighted_row_mean(
    validation$snow_water_equivalent_mm,
    weights
  )
}

grid_footprint <- as.polygons(ext(grid), crs = crs(grid))
covered_basin <- intersect(basin, grid_footprint)
basin_area_km2 <- sum(expanse(basin, unit = "km"))
covered_area_km2 <- sum(expanse(covered_basin, unit = "km"))
coverage_percent <- 100 * covered_area_km2 / basin_area_km2

result <- list(
  source_file = normalizePath(main_file, winslash = "/"),
  patch_status = patch_status,
  dates = main$dates,
  cell_metadata = cell_metadata,
  temperature_c = temperature,
  precipitation_mm = precipitation,
  pet_thornthwaite_mm = pet,
  water_balance_mm = water_balance,
  validation = validation,
  basin_monthly = basin_monthly,
  basin_validation = basin_validation,
  basin_boundary = list(
    source = "HydroBASINS v1c level 8",
    main_bas = 2080001450,
    coverage_percent = coverage_percent
  )
)

saveRDS(result, file.path(output_dir, "era5_basin_monthly.rds"), compress = "xz")
write.csv(cell_metadata, file.path(output_dir, "era5_basin_cell_metadata.csv"), row.names = FALSE)
write.csv(basin_monthly, file.path(output_dir, "era5_basin_area_weighted_monthly.csv"), row.names = FALSE)
write.csv(basin_validation, file.path(output_dir, "era5_basin_validation_monthly.csv"), row.names = FALSE)
writeLines(c(
  paste("Source:", result$source_file),
  paste("Boundary: HydroBASINS v1c level 8 MAIN_BAS 2080001450"),
  paste("Intersecting valid cells:", nrow(cell_metadata)),
  sprintf("Basin area represented by current NetCDF: %.3f%%", coverage_percent),
  paste("Precipitation correction:", patch_status),
  paste("Affected dates:", min(main$dates[affected]), "to", max(main$dates[affected])),
  paste("Optional validation variables:", if (length(validation)) paste(names(validation), collapse = ", ") else "none"),
  paste("PET: SPEI::thornthwaite version", as.character(packageVersion("SPEI")))
), file.path(output_dir, "era5_preprocessing_qc.txt"))

cat("Prepared", nrow(cell_metadata), "basin-intersecting valid cells.\n")
cat(sprintf("Current NetCDF represents %.3f%% of the HydroBASINS Seyhan area.\n", coverage_percent))
cat("Precipitation correction:", patch_status, "\n")
