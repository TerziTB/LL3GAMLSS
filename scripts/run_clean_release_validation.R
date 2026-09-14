#!/usr/bin/env Rscript

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run this script from the cleaned LL3GAMLSS repository root.")
}

output_dir <- file.path(root, "validation", "results", "release_checks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
arguments <- commandArgs(trailingOnly = TRUE)
library_argument <- grep("^--library=", arguments, value = TRUE)
dependency_argument <- grep("^--dependency-library=", arguments, value = TRUE)
release_library <- if (length(library_argument)) {
  sub("^--library=", "", library_argument[[1L]])
} else {
  Sys.getenv("LL3_RELEASE_LIBRARY", unset = "")
}
dependency_library <- if (length(dependency_argument)) {
  sub("^--dependency-library=", "", dependency_argument[[1L]])
} else {
  Sys.getenv("LL3_RELEASE_DEPENDENCY_LIBRARY", unset = "")
}
if (nzchar(release_library)) {
  normalized_library <- normalizePath(release_library, winslash = "/", mustWork = TRUE)
  normalized_dependency_library <- if (nzchar(dependency_library)) {
    normalizePath(dependency_library, winslash = "/", mustWork = TRUE)
  } else character()
  child_libraries <- paste(
    unique(c(normalized_library, normalized_dependency_library, .libPaths())),
    collapse = .Platform$path.sep
  )
  Sys.setenv(R_LIBS = child_libraries, R_LIBS_USER = "")
}
library_prefix <- ""

jobs <- list(
  installed_version = c(
    "-e",
    paste0(
      library_prefix,
      "expected <- unname(read.dcf('DESCRIPTION')[1,'Version']); actual <- as.character(packageVersion('LL3GAMLSS')); stopifnot(identical(actual, expected)); cat(actual, 'clean installation loaded', '\\n')"
    )
  ),
  probability_functions = c("validation/validate_probability_functions.R"),
  analytical_derivatives = c("validation/validate_derivatives.R"),
  stationary_optimizer_agreement = c("validation/validate_stationary_MLE.R"),
  package_tests = c(
    "-e",
    paste0(library_prefix, "testthat::test_local('.', reporter='summary', stop_on_failure=TRUE)")
  ),
  readme_smoke = c("-e", paste0(library_prefix, "source('examples/README_smoke_test.R'); cat('README_SMOKE_OK\\n')"))
)

rows <- lapply(names(jobs), function(name) {
  started <- Sys.time()
  log_path <- file.path(output_dir, paste0(name, ".log"))
  job_arguments <- jobs[[name]]
  if (length(job_arguments) >= 2L && identical(job_arguments[[1L]], "-e")) {
    job_arguments[[2L]] <- shQuote(job_arguments[[2L]])
  }
  status <- system2(rscript, job_arguments, stdout = log_path, stderr = log_path)
  data.frame(
    check = name,
    status = as.integer(status),
    passed = identical(as.integer(status), 0L),
    wall_clock_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    log = file.path("validation", "results", "release_checks", basename(log_path)),
    stringsAsFactors = FALSE
  )
})
guard_jobs <- list(
  restricted_station_input_guard = list(
    arguments = c("analysis/seyhan_real_data_analysis.R"),
    pattern = "access-controlled station input is not included"
  ),
  missing_ERA5_input_guard = list(
    arguments = c("analysis/prepare_era5_basin_data.R"),
    pattern = "ERA5-Land NetCDF input is not included"
  )
)
guard_rows <- lapply(names(guard_jobs), function(name) {
  started <- Sys.time()
  log_path <- file.path(output_dir, paste0(name, ".log"))
  status <- system2(
    rscript,
    guard_jobs[[name]]$arguments,
    stdout = log_path,
    stderr = log_path
  )
  output <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
  passed <- status != 0L && grepl(guard_jobs[[name]]$pattern, output, fixed = TRUE)
  data.frame(
    check = name,
    status = as.integer(status),
    passed = passed,
    wall_clock_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    log = file.path("validation", "results", "release_checks", basename(log_path)),
    stringsAsFactors = FALSE
  )
})
rows <- c(rows, guard_rows)
summary <- do.call(rbind, rows)
utils::write.csv(summary, file.path(output_dir, "clean_release_validation.csv"), row.names = FALSE)
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
print(summary)
if (!all(summary$passed)) stop("One or more clean-release validation jobs failed.")
