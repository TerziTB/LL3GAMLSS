#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: run_r_cmd_check.R <no-manual|as-cran> <source-tarball> [dependency-library]")
}
mode <- match.arg(args[[1L]], c("no-manual", "as-cran"))
tarball <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
Sys.setlocale("LC_ALL", "English_United States.utf8")
Sys.setenv(
  LC_ALL = "English_United States.utf8",
  LANG = "English_United States.utf8",
  LANGUAGE = ""
)
if (length(args) >= 3L) {
  dependency_library <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
  check_libraries <- paste(unique(c(dependency_library, .libPaths())), collapse = .Platform$path.sep)
  Sys.setenv(R_LIBS = check_libraries, R_LIBS_USER = check_libraries)
}
check_args <- c(
  "CMD", "check",
  if (mode == "as-cran") "--as-cran",
  "--no-manual",
  tarball
)
status <- system2(file.path(R.home("bin"), "R.exe"), check_args)
quit(status = status)
