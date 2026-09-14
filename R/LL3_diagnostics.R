# Diagnostics and model-comparison helpers for LL3 fits.

.extract_LL3_parameter <- function(object, what) {
  value <- try(stats::fitted(object, what = what, type = "response"), silent = TRUE)
  if (!inherits(value, "try-error") && length(value) > 0L) {
    return(as.numeric(value))
  }

  fallback <- paste0(what, ".fv")
  if (!is.null(object[[fallback]])) return(as.numeric(object[[fallback]]))
  stop("Could not extract fitted parameter ", what, ".")
}

extract_LL3_parameters <- function(object) {
  data.frame(
    mu = .extract_LL3_parameter(object, "mu"),
    sigma = .extract_LL3_parameter(object, "sigma"),
    nu = .extract_LL3_parameter(object, "nu")
  )
}

LL3_fit_converged <- function(fit) {
  if (!is.null(fit$converged)) return(isTRUE(fit$converged))
  value <- try(stats::deviance(fit), silent = TRUE)
  !inherits(value, "try-error") && length(value) == 1L && is.finite(value)
}

check_LL3_fit <- function(
    object,
    y = NULL,
    boundary_tolerance_multiplier = 100
) {
  parameters <- extract_LL3_parameters(
    object
  )

  if (is.null(y)) {
    y <- object$y
  }

  if (is.null(y)) {
    stop(
      "Supply y because it could not be extracted from the fitted object."
    )
  }

  y <- as.numeric(y)

  if (length(y) != nrow(parameters)) {
    stop(
      "The response length does not match the fitted parameter vectors."
    )
  }

  mu <- as.numeric(
    parameters$mu
  )

  sigma <- as.numeric(
    parameters$sigma
  )

  nu <- as.numeric(
    parameters$nu
  )

  lower <- attr(
    object,
    "LL3_lower"
  )

  if (
    is.null(lower) &&
      !is.null(object$LL3_lower)
  ) {
    lower <- object$LL3_lower
  }

  finite_parameters <- all(
    is.finite(mu) &
      is.finite(sigma) &
      is.finite(nu)
  )

  converged <- isTRUE(
    LL3_fit_converged(object)
  )

  support_ok <-
    finite_parameters &&
    all(
      y > mu
    )

  sigma_positive <-
    finite_parameters &&
    all(
      sigma > 0
    )

  nu_positive <-
    finite_parameters &&
    all(
      nu > 0
    )

  lower_available <-
    length(lower) == 1L &&
    is.finite(lower)

  if (lower_available) {
    numerical_scale <- max(
      c(
        1,
        abs(lower),
        abs(y),
        abs(mu)
      ),
      na.rm = TRUE
    )

    boundary_tolerance <-
      boundary_tolerance_multiplier *
      .Machine$double.eps *
      numerical_scale

    lower_gap <- lower - mu

    lower_below_data <-
      lower < min(y)

    mu_below_lower_exact <-
      finite_parameters &&
      all(
        lower_gap > 0
      )

    mu_not_materially_above_lower <-
      finite_parameters &&
      all(
        mu <=
          lower + boundary_tolerance
      )

    mu_boundary_contact <-
      finite_parameters &&
      any(
        lower_gap <= boundary_tolerance
      )

    material_boundary_violation <-
      finite_parameters &&
      any(
        mu >
          lower + boundary_tolerance
      )
  } else {
    lower_below_data <- NA
    mu_below_lower_exact <- NA
    mu_not_materially_above_lower <- NA
    mu_boundary_contact <- NA
    material_boundary_violation <- NA
  }

  # Return only logical values.
  #
  # The compatibility field `mu_below_lower` is retained for validation
  # scripts that use the original check name.

  c(
    converged =
      converged,

    finite_parameters =
      finite_parameters,

    support_ok =
      support_ok,

    sigma_positive =
      sigma_positive,

    nu_positive =
      nu_positive,

    lower_below_data =
      lower_below_data,

    mu_below_lower =
      mu_below_lower_exact,

    mu_below_lower_exact =
      mu_below_lower_exact,

    mu_not_materially_above_lower =
      mu_not_materially_above_lower,

    mu_boundary_contact =
      mu_boundary_contact,

    material_boundary_violation =
      material_boundary_violation
  )
}
LL3_pit <- function(y, mu, sigma, nu, epsilon = 1e-10) {
  p <- pLL3(y, mu = mu, sigma = sigma, nu = nu)
  pmin(pmax(p, epsilon), 1 - epsilon)
}

LL3_index <- function(y, mu, sigma, nu, epsilon = 1e-10) {
  stats::qnorm(LL3_pit(y, mu, sigma, nu, epsilon = epsilon))
}

LL3_residual_diagnostics <- function(object, y = NULL, max_lag = 24) {
  pars <- extract_LL3_parameters(object)
  if (is.null(y)) y <- object$y
  y <- as.numeric(y)
  z <- LL3_index(y, pars$mu, pars$sigma, pars$nu)
  lag_max <- min(max_lag, length(z) - 1L)
  acf_values <- if (lag_max >= 1L) {
    as.numeric(stats::acf(z, lag.max = lag_max, plot = FALSE)$acf[-1])
  } else numeric(0)

  list(
    residuals = z,
    summary = c(
      mean = mean(z),
      sd = stats::sd(z),
      tail_rate_95 = mean(abs(z) > stats::qnorm(0.975)),
      minimum = min(z),
      maximum = max(z)
    ),
    acf = acf_values,
    max_absolute_acf = if (length(acf_values)) max(abs(acf_values)) else NA_real_
  )
}

LL3_boundary_diagnostic <- function(
    object,
    y = NULL,
    boundary_tolerance_multiplier = 100
) {
  pars <- extract_LL3_parameters(object)

  if (is.null(y)) {
    y <- object$y
  }

  if (is.null(y)) {
    stop(
      "Supply y because it could not be extracted from the fit."
    )
  }

  y <- as.numeric(y)

  lower <- attr(
    object,
    "LL3_lower"
  )

  if (
    is.null(lower) &&
      !is.null(object$LL3_lower)
  ) {
    lower <- object$LL3_lower
  }

  if (
    is.null(lower) ||
      length(lower) != 1L ||
      !is.finite(lower)
  ) {
    return(
      data.frame(
        minimum_y_minus_mu =
          min(y - pars$mu),

        minimum_lower_minus_mu =
          NA_real_,

        minimum_y_minus_lower =
          NA_real_,

        boundary_tolerance =
          NA_real_,

        exact_boundary_contact =
          NA,

        material_boundary_violation =
          NA,

        inference_ready =
          NA,

        recommendation =
          "LL3 lower boundary is unavailable; inference readiness cannot be assessed."
      )
    )
  }

  numerical_scale <- max(
    1,
    abs(lower),
    abs(y),
    abs(pars$mu),
    na.rm = TRUE
  )

  tolerance <-
    boundary_tolerance_multiplier *
    .Machine$double.eps *
    numerical_scale

  lower_gap <- lower - pars$mu

  exact_boundary_contact <-
    any(lower_gap <= tolerance)

  material_boundary_violation <-
    any(pars$mu > lower + tolerance)

  inference_ready <-
    !exact_boundary_contact &&
    !material_boundary_violation &&
    all(is.finite(pars$mu)) &&
    all(is.finite(pars$sigma) & pars$sigma > 0) &&
    all(is.finite(pars$nu) & pars$nu > 0) &&
    all(y > pars$mu)

  data.frame(
    minimum_y_minus_mu =
      min(y - pars$mu),

    minimum_lower_minus_mu =
      min(lower_gap),

    minimum_y_minus_lower =
      min(y - lower),

    boundary_tolerance =
      tolerance,

    exact_boundary_contact =
      exact_boundary_contact,

    material_boundary_violation =
      material_boundary_violation,

    inference_ready =
      inference_ready,

    recommendation =
      if (inference_ready) {
        "Interior solution: ordinary likelihood diagnostics may be used."
      } else if (material_boundary_violation) {
        "Invalid threshold estimate: do not report this fit."
      } else {
        paste(
          "Threshold is on the numerical support boundary; ordinary Hessian,",
          "AIC/BIC, and Wald inference are non-regular. Refit or use a",
          "pre-specified threshold and sensitivity analysis."
        )
      }
  )
}


LL3_assert_inference_ready <- function(
    object,
    y = NULL,
    action = c("error", "warning")
) {
  action <- match.arg(action)

  checks <- check_LL3_fit(
    object,
    y = y
  )

  diagnostic <- LL3_boundary_diagnostic(
    object,
    y = y
  )

  required_true <- c(
    "converged",
    "finite_parameters",
    "support_ok",
    "sigma_positive",
    "nu_positive",
    "lower_below_data",
    "mu_not_materially_above_lower"
  )

  failed <- required_true[
    is.na(checks[required_true]) |
      !checks[required_true]
  ]

  boundary_problem <-
    is.na(diagnostic$inference_ready[1]) ||
    !isTRUE(diagnostic$inference_ready[1])

  if (length(failed) || boundary_problem) {
    message_text <- paste0(
      "LL3 fit is not ready for ordinary likelihood inference. ",
      if (length(failed)) {
        paste0("Failed checks: ", paste(failed, collapse = ", "), ". ")
      } else {
        ""
      },
      diagnostic$recommendation[1]
    )

    if (identical(action, "error")) {
      stop(message_text, call. = FALSE)
    }

    warning(message_text, call. = FALSE)
  }

  invisible(diagnostic)
}

LL3_model_table <- function(models, n = NULL) {
  if (!is.list(models) || !length(models)) stop("models must be a non-empty list.")
  if (is.null(names(models))) names(models) <- paste0("model", seq_along(models))

  rows <- lapply(seq_along(models), function(i) {
    model <- models[[i]]
    if (inherits(model, "try-error") || inherits(model, "error")) {
      return(data.frame(
        model = names(models)[i], converged = FALSE,
        inference_ready = FALSE, df = NA_real_,
        deviance = NA_real_, AIC = NA_real_, BIC = NA_real_
      ))
    }

    inference_ready <- tryCatch(
      {
        LL3_assert_inference_ready(model, y = model$y, action = "error")
        TRUE
      },
      error = function(e) FALSE
    )
    n_i <- if (is.null(n)) length(model$y) else n
    df <- model$df.fit
    dev <- if (inference_ready) stats::deviance(model) else NA_real_
    aic <- if (inference_ready) dev + 2 * df else NA_real_
    bic <- if (inference_ready) dev + log(n_i) * df else NA_real_
    data.frame(
      model = names(models)[i],
      converged = LL3_fit_converged(model),
      inference_ready = inference_ready,
      df = df,
      deviance = dev,
      AIC = aic,
      BIC = bic
    )
  })
  result <- do.call(rbind, rows)
  result$delta_AIC <- NA_real_
  finite_aic <- is.finite(result$AIC)
  if (any(finite_aic)) {
    result$delta_AIC[finite_aic] <-
      result$AIC[finite_aic] - min(result$AIC[finite_aic])
  }
  result$delta_BIC <- NA_real_
  finite_bic <- is.finite(result$BIC)
  if (any(finite_bic)) {
    result$delta_BIC[finite_bic] <-
      result$BIC[finite_bic] - min(result$BIC[finite_bic])
  }
  rownames(result) <- NULL
  result[order(result$AIC, na.last = TRUE), ]
}

plot_LL3_diagnostics <- function(object, y = NULL) {
  diagnostic <- LL3_residual_diagnostics(object, y = y)
  z <- diagnostic$residuals
  graphics::hist(z, breaks = "FD", xlab = "LL3 normalized residual",
                 main = "Normalized residuals")
  stats::qqnorm(z, main = "Normal Q-Q plot")
  stats::qqline(z)
  stats::acf(z, main = "Residual autocorrelation")
  invisible(diagnostic)
}
