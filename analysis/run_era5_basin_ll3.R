#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
run_started <- Sys.time()
mode <- if ("--full" %in% args) {
  "full"
} else if ("--pilot-spatial" %in% args) {
  "pilot_spatial"
} else {
  "quick"
}
argument_value <- function(prefix, default) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit)) sub(prefix, "", hit[[1L]]) else default
}
end_date <- as.Date(argument_value("--end=", "2025-12-01"))
workers <- as.integer(argument_value("--workers=", "4"))
workers <- max(1L, workers)
output_tag <- argument_value("--output-tag=", "")
refresh_cache <- "--refresh-cache" %in% args
if (nzchar(output_tag) && !grepl("^[A-Za-z0-9_.-]+$", output_tag)) {
  stop("--output-tag may contain only letters, numbers, dots, underscores, and hyphens.")
}

.libPaths(c(file.path(getwd(), ".rlib-realdata"), .libPaths()))
suppressPackageStartupMessages(library(LL3GAMLSS))
source(file.path("analysis", "era5_ll3_functions.R"))

input <- readRDS(file.path("analysis", "era5_results", "era5_basin_monthly.rds"))
keep <- input$dates <= end_date
dates <- input$dates[keep]
balance_matrix <- input$water_balance_mm[keep, , drop = FALSE]
cell_metadata <- input$cell_metadata
scales <- c(1L, 3L, 6L, 12L)

period_tag <- paste0(format(min(dates), "%Y%m"), "_", format(max(dates), "%Y%m"))
run_tag <- if (nzchar(output_tag)) output_tag else paste(mode, period_tag, sep = "_")
output_dir <- file.path("analysis", "era5_results", run_tag)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
task_cache_dir <- file.path(output_dir, "task_cache")
dir.create(task_cache_dir, recursive = TRUE, showWarnings = FALSE)

message("Fitting area-weighted basin series at four accumulation scales.")
basin_results <- lapply(scales, function(scale_months) {
  fit_LL3_selected_nspei(
    input$basin_monthly$water_balance_mm[keep],
    dates,
    scale_months,
    "basin_area_weighted"
  )
})
basin_summary <- do.call(rbind, lapply(basin_results, `[[`, "summary"))
basin_metrics <- do.call(rbind, lapply(basin_results, function(x) {
  cbind(series_id = "basin_area_weighted", scale_months = x$summary$scale_months, x$model_metrics)
}))
basin_values <- do.call(rbind, lapply(basin_results, function(x) {
  data.frame(
    date = x$date,
    scale_months = x$summary$scale_months,
    accumulated_balance_mm = x$accumulated_balance_mm,
    selected_model = x$summary$selected_model,
    nSPEI = x$nspei
  )
}))
write.csv(basin_summary, file.path(output_dir, "basin_model_selection.csv"), row.names = FALSE)
write.csv(basin_metrics, file.path(output_dir, "basin_candidate_metrics.csv"), row.names = FALSE)
write.csv(basin_values, file.path(output_dir, "basin_nspei.csv"), row.names = FALSE)

if (mode == "quick") {
  cell_indices <- unique(round(seq(1, ncol(balance_matrix), length.out = min(8L, ncol(balance_matrix)))))
  cell_scales <- 12L
} else if (mode == "pilot_spatial") {
  cell_indices <- seq_len(ncol(balance_matrix))
  cell_scales <- 12L
} else {
  cell_indices <- seq_len(ncol(balance_matrix))
  cell_scales <- scales
}

tasks <- expand.grid(
  cell_index = cell_indices,
  scale_months = cell_scales,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
task_list <- split(tasks, seq_len(nrow(tasks)))
message("Fitting ", nrow(tasks), " cell-scale tasks using ", workers, " worker(s).")

run_task <- function(task) {
  i <- task$cell_index[[1L]]
  scale_months <- task$scale_months[[1L]]
  cache_file <- file.path(
    task_cache_dir,
    sprintf("%s_scale_%02d.rds", cell_metadata$cell_key[[i]], scale_months)
  )
  if (!refresh_cache && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(cached)) return(cached)
  }
  result <- fit_LL3_selected_nspei(
    balance_matrix[, i],
    dates,
    scale_months,
    cell_metadata$cell_key[[i]]
  )
  temporary_file <- tempfile(pattern = "ll3_task_", tmpdir = task_cache_dir, fileext = ".rds")
  saveRDS(result, temporary_file, compress = FALSE)
  if (!file.rename(temporary_file, cache_file)) {
    unlink(temporary_file)
    stop("Could not finalize task checkpoint: ", cache_file)
  }
  result
}

if (workers == 1L || length(task_list) == 1L) {
  cell_results <- lapply(task_list, run_task)
} else {
  cluster <- parallel::makePSOCKcluster(min(workers, length(task_list)))
  project_dir <- normalizePath(getwd(), winslash = "/")
  parallel::clusterExport(cluster, "project_dir")
  parallel::clusterEvalQ(cluster, {
    setwd(project_dir)
    .libPaths(c(file.path(project_dir, ".rlib-realdata"), .libPaths()))
    suppressPackageStartupMessages(library(LL3GAMLSS))
    NULL
  })
  parallel::clusterExport(
    cluster,
    c(
      "LL3_joint_model_metrics", "LL3_index_from_parameters",
      "fit_LL3_selected_nspei", "balance_matrix", "dates", "cell_metadata",
      "task_cache_dir", "refresh_cache"
    ),
    envir = environment()
  )
  cell_results <- parallel::parLapplyLB(cluster, task_list, run_task)
  parallel::stopCluster(cluster)
}

cell_summary <- do.call(rbind, lapply(cell_results, `[[`, "summary"))
cell_summary$cell_id <- as.integer(sub("cell_", "", cell_summary$series_id))
cell_summary <- merge(cell_summary, cell_metadata, by = "cell_id", all.x = TRUE, sort = FALSE)
cell_metrics <- do.call(rbind, lapply(cell_results, function(x) {
  cbind(series_id = x$summary$series_id, scale_months = x$summary$scale_months, x$model_metrics)
}))

write.csv(cell_summary, file.path(output_dir, "cell_model_selection.csv"), row.names = FALSE)
write.csv(cell_metrics, file.path(output_dir, "cell_candidate_metrics.csv"), row.names = FALSE)
saveRDS(cell_results, file.path(output_dir, "cell_selected_nspei.rds"), compress = "xz")
saveRDS(basin_results, file.path(output_dir, "basin_selected_nspei.rds"), compress = "xz")

area_summary <- do.call(rbind, lapply(split(cell_summary, cell_summary$scale_months), function(d) {
  models <- c("stationary", "sigma_time", "nu_time", "sigma_nu_time", NA_character_)
  labels <- c("stationary", "sigma_time", "nu_time", "sigma_nu_time", "fit_failure")
  do.call(rbind, lapply(seq_along(models), function(i) {
    selected <- if (is.na(models[[i]])) is.na(d$selected_model) else d$selected_model == models[[i]]
    data.frame(
      scale_months = unique(d$scale_months),
      selected_model = labels[[i]],
      cell_count = sum(selected),
      basin_area_weight_km2 = sum(d$basin_area_weight_km2[selected], na.rm = TRUE),
      basin_area_percent = 100 * sum(d$basin_area_weight_km2[selected], na.rm = TRUE) /
        sum(d$basin_area_weight_km2, na.rm = TRUE)
    )
  }))
}))
write.csv(area_summary, file.path(output_dir, "area_weighted_model_fractions.csv"), row.names = FALSE)

writeLines(c(
  paste("Mode:", mode),
  paste("Analysis period:", min(dates), "to", max(dates)),
  paste("Source coverage percent:", input$basin_boundary$coverage_percent),
  paste("Precipitation correction:", input$patch_status),
  "Selection order: fit M0/Msigma/Mnu/Msigma+nu, reject non-ready fits, compare joint calendar-month AICc, require delta AICc >= 2, then compute nSPEI from selected conditional LL3 CDF.",
  paste("Cell-scale tasks:", nrow(tasks)),
  paste("Workers:", workers),
  paste("Cache policy:", if (refresh_cache) "fresh fits" else "reuse valid task caches")
), file.path(output_dir, "run_methodology.txt"))

run_finished <- Sys.time()
package_names <- c(
  "LL3GAMLSS", "gamlss", "gamlss.dist", "SPEI", "terra", "ncdf4",
  "jsonlite", "FAdist", "numDeriv", "testthat"
)
package_versions <- vapply(package_names, function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
}, character(1))
ram_bytes <- tryCatch({
  output <- system2(
    "powershell",
    c(
      "-NoProfile", "-Command",
      "Add-Type -AssemblyName Microsoft.VisualBasic; ([Microsoft.VisualBasic.Devices.ComputerInfo]::new()).TotalPhysicalMemory"
    ),
    stdout = TRUE,
    stderr = FALSE
  )
  as.numeric(output[[1L]])
}, error = function(e) NA_real_)
runtime <- data.frame(
  started = format(run_started, "%Y-%m-%d %H:%M:%S %Z"),
  finished = format(run_finished, "%Y-%m-%d %H:%M:%S %Z"),
  wall_clock_seconds = as.numeric(difftime(run_finished, run_started, units = "secs")),
  workers = workers,
  logical_cpus = parallel::detectCores(logical = TRUE),
  processor = Sys.getenv("PROCESSOR_IDENTIFIER", unset = NA_character_),
  available_memory_bytes = ram_bytes,
  operating_system = paste(Sys.info()[c("sysname", "release", "version", "machine")], collapse = "; "),
  R_version = R.version.string,
  cache_policy = if (refresh_cache) "fresh fits" else "reuse valid task caches",
  cell_scale_tasks = nrow(tasks),
  stringsAsFactors = FALSE
)
utils::write.csv(runtime, file.path(output_dir, "runtime_hardware.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(package = package_names, version = package_versions),
  file.path(output_dir, "package_versions.csv"),
  row.names = FALSE
)
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message("Completed ERA5-Land LL3 run: ", normalizePath(output_dir, winslash = "/"))
