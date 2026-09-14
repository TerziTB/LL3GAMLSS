LL3_joint_model_metrics <- function(fits, model_name) {
  valid <- vapply(fits, inherits, logical(1), what = "gamlss")
  valid_fits <- fits[valid]
  if (length(valid_fits) == 0L) {
    return(data.frame(
      model = model_name, fitted_months = 0L, converged_months = 0L,
      inference_ready_months = 0L, boundary_contact_months = NA_integer_,
      n = NA_integer_, df = NA_real_, deviance = NA_real_, AIC = NA_real_,
      AICc = NA_real_, BIC = NA_real_, sigma_slope_median = NA_real_,
      nu_slope_median = NA_real_
    ))
  }

  diagnostics <- lapply(valid_fits, function(fit) {
    tryCatch(LL3_boundary_diagnostic(fit), error = function(e) NULL)
  })
  inference_ready <- vapply(diagnostics, function(x) {
    !is.null(x) && isTRUE(x$inference_ready[[1L]])
  }, logical(1))
  boundary_contact <- vapply(diagnostics, function(x) {
    !is.null(x) && isTRUE(x$exact_boundary_contact[[1L]])
  }, logical(1))

  deviance <- sum(vapply(valid_fits, stats::deviance, numeric(1)))
  degrees_freedom <- sum(vapply(valid_fits, function(x) x$df.fit, numeric(1)))
  observations <- sum(vapply(valid_fits, function(x) length(x$y), integer(1)))
  aic <- deviance + 2 * degrees_freedom
  aicc <- if (observations > degrees_freedom + 1) {
    aic + 2 * degrees_freedom * (degrees_freedom + 1) /
      (observations - degrees_freedom - 1)
  } else {
    NA_real_
  }

  slope <- function(fit, what) {
    value <- try(stats::coef(fit, what = what)[["time_scaled"]], silent = TRUE)
    if (inherits(value, "try-error") || length(value) == 0L) NA_real_ else value
  }
  sigma_slopes <- vapply(valid_fits, slope, numeric(1), what = "sigma")
  nu_slopes <- vapply(valid_fits, slope, numeric(1), what = "nu")

  data.frame(
    model = model_name,
    fitted_months = length(valid_fits),
    converged_months = sum(vapply(valid_fits, LL3_fit_converged, logical(1))),
    inference_ready_months = sum(inference_ready),
    boundary_contact_months = sum(boundary_contact),
    n = observations,
    df = degrees_freedom,
    deviance = deviance,
    AIC = aic,
    AICc = aicc,
    BIC = deviance + log(observations) * degrees_freedom,
    sigma_slope_median = if (any(is.finite(sigma_slopes))) {
      stats::median(sigma_slopes, na.rm = TRUE)
    } else NA_real_,
    nu_slope_median = if (any(is.finite(nu_slopes))) {
      stats::median(nu_slopes, na.rm = TRUE)
    } else NA_real_
  )
}

LL3_index_from_parameters <- function(y, parameters, epsilon = 1e-10) {
  complete <- is.finite(y) & stats::complete.cases(parameters)
  probability <- rep(NA_real_, length(y))
  probability[complete] <- pLL3(
    y[complete], parameters$mu[complete], parameters$sigma[complete], parameters$nu[complete]
  )
  index <- rep(NA_real_, length(y))
  index[complete] <- stats::qnorm(pmin(pmax(probability[complete], epsilon), 1 - epsilon))
  index
}

fit_LL3_selected_nspei <- function(
    water_balance,
    dates,
    scale_months,
    series_id,
    minimum_delta_aicc = 2,
    control = LL3_default_control(
      n_cycles = 2000,
      mu_step = 0.03,
      sigma_step = 0.03,
      nu_step = 0.03,
      trace = FALSE
    )) {
  model_names <- c("stationary", "sigma_time", "nu_time", "sigma_nu_time")
  accumulated <- LL3_accumulate(water_balance, scale = scale_months)
  month <- as.integer(format(dates, "%m"))
  time_scaled <- as.numeric(scale(seq_along(dates)))

  parameter_store <- setNames(lapply(model_names, function(unused) {
    data.frame(
      mu = rep(NA_real_, length(dates)),
      sigma = rep(NA_real_, length(dates)),
      nu = rep(NA_real_, length(dates))
    )
  }), model_names)
  monthly_fits <- setNames(lapply(model_names, function(unused) vector("list", 12L)), model_names)
  failures <- character()

  for (current_month in seq_len(12L)) {
    rows <- month == current_month & is.finite(accumulated)
    if (sum(rows) < 30L) {
      failures <- c(failures, paste0("month_", current_month, "_insufficient_data"))
      next
    }
    fit_data <- data.frame(
      response = accumulated[rows],
      time_scaled = time_scaled[rows]
    )
    candidates <- tryCatch(
      suppressWarnings(fit_LL3_candidate_models(
        fit_data,
        response = "response",
        covariate = "time_scaled",
        control = control
      )),
      error = function(e) e
    )
    if (inherits(candidates, "error")) {
      failures <- c(failures, paste0("month_", current_month, "_", conditionMessage(candidates)))
      next
    }

    mapped <- list(
      stationary = candidates$stationary,
      sigma_time = candidates$sigma,
      nu_time = candidates$nu,
      sigma_nu_time = candidates$joint
    )
    for (model_name in model_names) {
      fit <- mapped[[model_name]]
      monthly_fits[[model_name]][[current_month]] <- fit
      parameters <- extract_LL3_parameters(fit)
      parameter_store[[model_name]][rows, ] <- data.frame(
        mu = parameters$mu,
        sigma = parameters$sigma,
        nu = parameters$nu
      )
    }
  }

  metrics <- do.call(rbind, lapply(model_names, function(model_name) {
    LL3_joint_model_metrics(monthly_fits[[model_name]], model_name)
  }))
  eligible <- with(
    metrics,
    fitted_months == 12L & converged_months == 12L &
      inference_ready_months == 12L & is.finite(AICc)
  )
  stationary_row <- which(metrics$model == "stationary")
  stationary_ready <- length(stationary_row) == 1L && eligible[[stationary_row]]

  nonstationary_rows <- which(metrics$model != "stationary" & eligible)
  best_nonstationary_row <- if (length(nonstationary_rows)) {
    nonstationary_rows[[which.min(metrics$AICc[nonstationary_rows])]]
  } else integer()
  delta <- if (stationary_ready && length(best_nonstationary_row)) {
    metrics$AICc[[stationary_row]] - metrics$AICc[[best_nonstationary_row]]
  } else NA_real_

  if (!stationary_ready) {
    selected_model <- NA_character_
    selection_status <- "stationary_fit_not_inference_ready"
  } else if (length(best_nonstationary_row) && is.finite(delta) && delta >= minimum_delta_aicc) {
    selected_model <- metrics$model[[best_nonstationary_row]]
    selection_status <- if (delta < 4) {
      "weak_nonstationarity"
    } else if (delta < 10) {
      "moderate_nonstationarity"
    } else {
      "strong_nonstationarity"
    }
  } else {
    selected_model <- "stationary"
    selection_status <- "stationary_selected"
  }

  selected_parameters <- if (is.na(selected_model)) {
    data.frame(
      mu = rep(NA_real_, length(dates)),
      sigma = rep(NA_real_, length(dates)),
      nu = rep(NA_real_, length(dates))
    )
  } else {
    parameter_store[[selected_model]]
  }
  nspei <- LL3_index_from_parameters(accumulated, selected_parameters)
  selected_metric <- if (is.na(selected_model)) NULL else {
    metrics[metrics$model == selected_model, , drop = FALSE]
  }

  summary <- data.frame(
    series_id = series_id,
    scale_months = scale_months,
    selected_model = selected_model,
    selection_status = selection_status,
    delta_AICc_favoring_best_nonstationary = delta,
    stationary_AICc = metrics$AICc[metrics$model == "stationary"],
    sigma_time_AICc = metrics$AICc[metrics$model == "sigma_time"],
    nu_time_AICc = metrics$AICc[metrics$model == "nu_time"],
    sigma_nu_time_AICc = metrics$AICc[metrics$model == "sigma_nu_time"],
    selected_sigma_slope_median = if (is.null(selected_metric)) NA_real_ else selected_metric$sigma_slope_median,
    selected_nu_slope_median = if (is.null(selected_metric)) NA_real_ else selected_metric$nu_slope_median,
    finite_nspei = sum(is.finite(nspei)),
    nspei_mean = mean(nspei, na.rm = TRUE),
    nspei_sd = stats::sd(nspei, na.rm = TRUE),
    nspei_minimum = suppressWarnings(min(nspei, na.rm = TRUE)),
    nspei_maximum = suppressWarnings(max(nspei, na.rm = TRUE)),
    extreme_dry_count = sum(nspei <= -2, na.rm = TRUE),
    fit_failures = paste(failures, collapse = " | "),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    model_metrics = metrics,
    date = dates,
    accumulated_balance_mm = accumulated,
    nspei = nspei,
    selected_parameters = selected_parameters
  )
}
