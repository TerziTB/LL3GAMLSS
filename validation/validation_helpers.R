resolve_project_dir <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  script_argument <- grep("^--file=", command, value = TRUE)
  candidates <- c(
    if (length(script_argument)) {
      dirname(dirname(normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/")))
    },
    normalizePath(getwd(), winslash = "/"),
    normalizePath(file.path(getwd(), ".."), winslash = "/")
  )
  candidates <- unique(candidates)
  valid <- vapply(candidates, function(path) {
    file.exists(file.path(path, "DESCRIPTION")) && dir.exists(file.path(path, "R"))
  }, logical(1))
  if (!any(valid)) {
    stop("Run this script from the LL3GAMLSS project or validation directory.")
  }
  candidates[which(valid)[[1L]]]
}

load_LL3_source <- function(project_dir = resolve_project_dir()) {
  source_order <- c(
    "LL3_distribution.R",
    "LL3_gamlss_family.R",
    "LL3_diagnostics.R",
    "LL3_fitting.R",
    "LL3_spei.R"
  )
  invisible(lapply(file.path(project_dir, "R", source_order), source, local = .GlobalEnv))
}

empirical_LL3_control <- function() {
  LL3_default_control(
    n_cycles = 2000,
    mu_step = 0.03,
    sigma_step = 0.03,
    nu_step = 0.03,
    trace = FALSE
  )
}

seasonal_candidate_metrics <- function(fits, model_name) {
  valid <- vapply(fits, inherits, logical(1), what = "gamlss")
  valid_fits <- fits[valid]
  if (!length(valid_fits)) {
    return(data.frame(
      model = model_name,
      fitted_months = 0L,
      converged_months = 0L,
      inference_ready_months = 0L,
      boundary_contact_months = NA_integer_,
      n = NA_integer_,
      df = NA_real_,
      deviance = NA_real_,
      AICc = NA_real_,
      sigma_slope_median = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  diagnostics <- lapply(valid_fits, function(fit) {
    tryCatch(LL3_boundary_diagnostic(fit), error = function(e) NULL)
  })
  ready <- vapply(diagnostics, function(x) {
    !is.null(x) && isTRUE(x$inference_ready[[1L]])
  }, logical(1))
  boundary <- vapply(diagnostics, function(x) {
    !is.null(x) && isTRUE(x$exact_boundary_contact[[1L]])
  }, logical(1))
  deviance <- sum(vapply(valid_fits, stats::deviance, numeric(1)))
  parameters <- sum(vapply(valid_fits, function(x) x$df.fit, numeric(1)))
  observations <- sum(vapply(valid_fits, function(x) length(x$y), integer(1)))
  aic <- deviance + 2 * parameters
  aicc <- if (observations > parameters + 1) {
    aic + 2 * parameters * (parameters + 1) /
      (observations - parameters - 1)
  } else {
    NA_real_
  }
  slope <- vapply(valid_fits, function(fit) {
    coefficient <- try(stats::coef(fit, what = "sigma")[["time_scaled"]], silent = TRUE)
    if (inherits(coefficient, "try-error") || !length(coefficient)) NA_real_ else coefficient
  }, numeric(1))

  data.frame(
    model = model_name,
    fitted_months = length(valid_fits),
    converged_months = sum(vapply(valid_fits, LL3_fit_converged, logical(1))),
    inference_ready_months = sum(ready),
    boundary_contact_months = sum(boundary),
    n = observations,
    df = parameters,
    deviance = deviance,
    AICc = aicc,
    sigma_slope_median = if (any(is.finite(slope))) stats::median(slope, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}

rate_with_mcse <- function(indicator) {
  indicator <- indicator[!is.na(indicator)]
  if (!length(indicator)) {
    return(c(estimate = NA_real_, mcse = NA_real_, denominator = 0))
  }
  estimate <- mean(indicator)
  c(
    estimate = estimate,
    mcse = sqrt(estimate * (1 - estimate) / length(indicator)),
    denominator = length(indicator)
  )
}

write_session_info <- function(path) {
  output <- capture.output(utils::sessionInfo())
  writeLines(output, path, useBytes = TRUE)
  invisible(path)
}
