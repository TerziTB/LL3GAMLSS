# validation/simulate_interval_coverage.R
#
# Monte Carlo coverage validation for numerical-Hessian Wald intervals
# under joint sigma-nu nonstationarity.
#
# Data-generating model:
#
#   mu = -50
#   log(sigma_t) = log(70) + 0.55 * x_t
#   log(nu_t)    = log(1.8) - 0.20 * x_t
#
# The script:
#
#   1. fits the joint LL3 GAMLSS model;
#   2. polishes the solution using an independent formula-only likelihood;
#   3. calculates the numerical observed Hessian;
#   4. constructs 95% Wald intervals;
#   5. estimates interval coverage for the sigma and nu slopes.
#
# The independent likelihood does not call dLL3().


# ============================================================================
# 1. Project setup
# ============================================================================

resolve_project_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0L) {
    script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE))
  }

  candidates <- unique(c(normalizePath(getwd(), winslash = "/"), normalizePath(file.path(getwd(), ".."), winslash = "/")))
  valid <- candidates[file.exists(file.path(candidates, "DESCRIPTION")) & dir.exists(file.path(candidates, "R"))]
  if (length(valid) == 0L) stop("Run this script with Rscript or from the LL3GAMLSS project root.")
  valid[[1L]]
}

project_dir <- resolve_project_dir()

results_dir <- file.path(
  project_dir,
  "validation",
  "results"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

source(file.path(project_dir, "R", "LL3_distribution.R"))
source(file.path(project_dir, "R", "LL3_gamlss_family.R"))
source(file.path(project_dir, "R", "LL3_fitting.R"))
source(file.path(project_dir, "R", "LL3_diagnostics.R"))

if (!requireNamespace("gamlss", quietly = TRUE)) {
  stop("The gamlss package must be installed.")
}

if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop(
    "The numDeriv package must be installed. Run:\n",
    "install.packages('numDeriv')"
  )
}


# ============================================================================
# 2. Simulation settings
# ============================================================================

repetitions <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "LL3_COVERAGE_REPS",
      "500"
    )
  )
)

if (
  length(repetitions) != 1L ||
    is.na(repetitions) ||
    repetitions < 1L
) {
  stop(
    "LL3_COVERAGE_REPS must be a positive integer."
  )
}

sample_sizes <- c(
  100L,
  250L
)

truth <- list(
  mu = -50,

  sigma_intercept = log(70),
  sigma_slope = 0.55,

  nu_intercept = log(1.8),
  nu_slope = -0.20
)

confidence_level <- 0.95

critical_value <- stats::qnorm(
  1 - (1 - confidence_level) / 2
)

master_seed <- 2031L

set.seed(master_seed)

simulation_seeds <- sample.int(
  .Machine$integer.max,
  size = repetitions * length(sample_sizes),
  replace = FALSE
)


# ============================================================================
# 3. GAMLSS controls
# ============================================================================

primary_control <- gamlss::gamlss.control(
  c.crit = 1e-3,
  n.cyc = 1000,

  mu.step = 0.03,
  sigma.step = 0.04,
  nu.step = 0.04,

  gd.tol = Inf,
  autostep = TRUE,
  trace = FALSE
)

retry_control <- gamlss::gamlss.control(
  c.crit = 1e-3,
  n.cyc = 3000,

  mu.step = 0.03,
  sigma.step = 0.04,
  nu.step = 0.04,

  gd.tol = Inf,
  autostep = TRUE,
  trace = FALSE
)


# ============================================================================
# 4. Numerical utilities
# ============================================================================

stable_log1pexp <- function(x) {
  result <- numeric(length(x))

  positive <- x > 0

  result[positive] <-
    x[positive] +
    log1p(
      exp(
        -x[positive]
      )
    )

  result[!positive] <-
    log1p(
      exp(
        x[!positive]
      )
    )

  result
}


data_scale <- function(y) {
  candidates <- c(
    diff(range(y)),
    stats::IQR(y),
    stats::mad(y),
    stats::sd(y),
    1
  )

  candidates <- candidates[
    is.finite(candidates) &
      candidates > 0
  ]

  max(candidates)
}


safe_mean <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(x)
}


safe_sd <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(NA_real_)
  }

  stats::sd(x)
}


safe_median <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  stats::median(x)
}


safe_rate <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(as.logical(x))
}


collapse_messages <- function(x) {
  x <- unique(
    x[
      !is.na(x) &
        nzchar(x)
    ]
  )

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(
    x,
    collapse = " | "
  )
}


recover_linear_coefficients <- function(
    values,
    x
) {
  design_matrix <- cbind(
    intercept = 1,
    slope = x
  )

  as.numeric(
    qr.coef(
      qr(design_matrix),
      values
    )
  )
}


# ============================================================================
# 5. Independent LL3 likelihood
#
# Parameter order:
#
#   theta[1] = alpha_mu
#   theta[2] = sigma intercept
#   theta[3] = sigma slope
#   theta[4] = nu intercept
#   theta[5] = nu slope
#
# with:
#
#   mu = support_lower - exp(alpha_mu)
# ============================================================================

direct_LL3_negative_loglikelihood <- function(
    theta,
    y,
    x,
    support_lower
) {
  theta <- as.numeric(theta)

  if (
    length(theta) != 5L ||
      any(!is.finite(theta))
  ) {
    return(1e100)
  }

  alpha_mu <- theta[1L]

  eta_sigma <-
    theta[2L] +
    theta[3L] * x

  eta_nu <-
    theta[4L] +
    theta[5L] * x

  if (
    alpha_mu < -50 ||
      alpha_mu > 50 ||
      any(!is.finite(eta_sigma)) ||
      any(!is.finite(eta_nu)) ||
      any(abs(eta_sigma) > 30) ||
      any(abs(eta_nu) > 20)
  ) {
    return(1e100)
  }

  mu <- support_lower -
    exp(alpha_mu)

  distance <- y - mu

  if (
    any(!is.finite(distance)) ||
      any(distance <= 0)
  ) {
    return(1e100)
  }

  nu <- exp(eta_nu)

  if (
    any(!is.finite(nu)) ||
      any(nu <= 0)
  ) {
    return(1e100)
  }

  log_z <-
    log(distance) -
    eta_sigma

  log_density <-
    eta_nu -
    eta_sigma +
    (nu - 1) * log_z -
    2 * stable_log1pexp(
      nu * log_z
    )

  if (any(!is.finite(log_density))) {
    return(1e100)
  }

  objective <- -sum(log_density)

  if (!is.finite(objective)) {
    return(1e100)
  }

  objective
}


# ============================================================================
# 6. Fit-capture helpers
# ============================================================================

capture_fit <- function(arguments) {
  warning_messages <- character()

  fitted_object <- withCallingHandlers(
    tryCatch(
      do.call(
        fit_LL3_gamlss,
        arguments
      ),
      error = function(e) {
        e
      }
    ),

    warning = function(w) {
      warning_messages <<- c(
        warning_messages,
        conditionMessage(w)
      )

      invokeRestart("muffleWarning")
    }
  )

  list(
    fit = fitted_object,

    warnings = collapse_messages(
      warning_messages
    )
  )
}


validate_fit <- function(
    fit,
    y
) {
  output <- list(
    success = FALSE,
    parameters = NULL,
    message = NA_character_
  )

  if (inherits(fit, "error")) {
    output$message <- conditionMessage(fit)
    return(output)
  }

  if (inherits(fit, "try-error")) {
    output$message <- as.character(fit)
    return(output)
  }

  parameters <- try(
    extract_LL3_parameters(fit),
    silent = TRUE
  )

  if (inherits(parameters, "try-error")) {
    output$message <- as.character(parameters)
    return(output)
  }

  mu <- as.numeric(parameters$mu)
  sigma <- as.numeric(parameters$sigma)
  nu <- as.numeric(parameters$nu)

  checks <- c(
    converged =
      isTRUE(fit$converged),

    correct_lengths =
      length(mu) == length(y) &&
      length(sigma) == length(y) &&
      length(nu) == length(y),

    finite_parameters =
      all(
        is.finite(
          c(
            mu,
            sigma,
            nu
          )
        )
      ),

    support_ok =
      all(y > mu),

    sigma_positive =
      all(sigma > 0),

    nu_positive =
      all(nu > 0)
  )

  output$parameters <- list(
    mu = mu,
    sigma = sigma,
    nu = nu
  )

  output$success <- all(checks)

  if (!output$success) {
    output$message <- paste0(
      "Failed checks: ",
      paste(
        names(checks)[!checks],
        collapse = ", "
      ),
      "."
    )
  }

  output
}


fit_one_model <- function(
    mu.formula,
    sigma.formula,
    nu.formula,
    data,
    lower,
    mu.start = NULL,
    sigma.start = NULL,
    nu.start = NULL
) {
  construct_arguments <- function(control) {
    arguments <- list(
      mu.formula = mu.formula,
      sigma.formula = sigma.formula,
      nu.formula = nu.formula,

      data = data,
      lower = lower,

      information = "opg",
      allow_nonstationary_mu = FALSE,

      control = control
    )

    if (!is.null(mu.start)) {
      safe_margin <- max(
        1e-8,
        sqrt(.Machine$double.eps) *
          max(
            1,
            abs(lower),
            abs(data$y)
          )
      )

      arguments$mu.start <- pmin(
        as.numeric(mu.start),
        lower - safe_margin
      )
    }

    if (!is.null(sigma.start)) {
      arguments$sigma.start <-
        as.numeric(sigma.start)
    }

    if (!is.null(nu.start)) {
      arguments$nu.start <-
        as.numeric(nu.start)
    }

    arguments
  }

  first_attempt <- capture_fit(
    construct_arguments(
      primary_control
    )
  )

  first_validation <- validate_fit(
    first_attempt$fit,
    data$y
  )

  if (first_validation$success) {
    return(
      list(
        success = TRUE,
        fit = first_attempt$fit,
        parameters = first_validation$parameters,
        retried = FALSE,
        warning = first_attempt$warnings,
        message = NA_character_
      )
    )
  }

  second_attempt <- capture_fit(
    construct_arguments(
      retry_control
    )
  )

  second_validation <- validate_fit(
    second_attempt$fit,
    data$y
  )

  list(
    success = second_validation$success,
    fit = second_attempt$fit,
    parameters = second_validation$parameters,
    retried = TRUE,

    warning = collapse_messages(
      c(
        first_attempt$warnings,
        second_attempt$warnings
      )
    ),

    message = collapse_messages(
      c(
        first_validation$message,
        second_validation$message
      )
    )
  )
}


# ============================================================================
# 7. Fit the joint model with auxiliary starts
# ============================================================================

fit_joint_model <- function(
    data,
    lower
) {
  stationary <- fit_one_model(
    mu.formula = y ~ 1,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,

    data = data,
    lower = lower
  )

  if (!stationary$success) {
    return(
      list(
        success = FALSE,
        stage = "stationary fit",
        message = stationary$message
      )
    )
  }

  stationary_parameters <- stationary$parameters

  sigma_only <- fit_one_model(
    mu.formula = y ~ 1,
    sigma.formula = ~ x,
    nu.formula = ~ 1,

    data = data,
    lower = lower,

    mu.start =
      stationary_parameters$mu,

    sigma.start =
      stationary_parameters$sigma,

    nu.start =
      stationary_parameters$nu
  )

  nu_only <- fit_one_model(
    mu.formula = y ~ 1,
    sigma.formula = ~ 1,
    nu.formula = ~ x,

    data = data,
    lower = lower,

    mu.start =
      stationary_parameters$mu,

    sigma.start =
      stationary_parameters$sigma,

    nu.start =
      stationary_parameters$nu
  )

  sigma_start <-
    if (sigma_only$success) {
      sigma_only$parameters$sigma
    } else {
      stationary_parameters$sigma
    }

  nu_start <-
    if (nu_only$success) {
      nu_only$parameters$nu
    } else {
      stationary_parameters$nu
    }

  joint <- fit_one_model(
    mu.formula = y ~ 1,
    sigma.formula = ~ x,
    nu.formula = ~ x,

    data = data,
    lower = lower,

    mu.start =
      stationary_parameters$mu,

    sigma.start =
      sigma_start,

    nu.start =
      nu_start
  )

  if (!joint$success) {
    return(
      list(
        success = FALSE,
        stage = "joint fit",
        message = joint$message
      )
    )
  }

  parameters <- joint$parameters

  fitted_mu <- mean(
    parameters$mu
  )

  sigma_coefficients <- recover_linear_coefficients(
    log(parameters$sigma),
    data$x
  )

  nu_coefficients <- recover_linear_coefficients(
    log(parameters$nu),
    data$x
  )

  theta <- c(
    alpha_mu =
      log(
        lower -
          fitted_mu
      ),

    sigma_intercept =
      sigma_coefficients[1L],

    sigma_slope =
      sigma_coefficients[2L],

    nu_intercept =
      nu_coefficients[1L],

    nu_slope =
      nu_coefficients[2L]
  )

  list(
    success = TRUE,
    fit = joint$fit,
    theta = theta,

    warning = collapse_messages(
      c(
        stationary$warning,
        sigma_only$warning,
        nu_only$warning,
        joint$warning
      )
    ),

    retried = any(
      c(
        stationary$retried,
        sigma_only$retried,
        nu_only$retried,
        joint$retried
      )
    )
  )
}


# ============================================================================
# 8. Independent likelihood polishing
# ============================================================================

polish_likelihood <- function(
    theta_start,
    y,
    x,
    support_lower
) {
  scale_y <- data_scale(y)

  parameter_lower_bounds <- c(
    log(
      max(
        1e-12 * scale_y,
        1e-12
      )
    ),
    -20,
    -10,
    -10,
    -10
  )

  parameter_upper_bounds <- c(
    log(
      max(
        1e4 * scale_y,
        1e4
      )
    ),
    20,
    10,
    10,
    10
  )

  clip_start <- function(start) {
    pmin(
      pmax(
        start,
        parameter_lower_bounds + 1e-8
      ),
      parameter_upper_bounds - 1e-8
    )
  }

  generic_start <- c(
    log(scale_y),
    log(max(stats::sd(y), 1)),
    0,
    log(2),
    0
  )

  starts <- rbind(
    theta_start,

    theta_start +
      c(
        0.25,
        0,
        0,
        0,
        0
      ),

    theta_start -
      c(
        0.25,
        0,
        0,
        0,
        0
      ),

    theta_start +
      c(
        0,
        0.10,
        0.10,
        0.10,
        0.10
      ),

    generic_start
  )

  starts <- t(
    apply(
      starts,
      1L,
      clip_start
    )
  )

  candidates <- list()

  for (start_index in seq_len(nrow(starts))) {
    current_start <- starts[
      start_index,
      ,
      drop = TRUE
    ]

    current_fit <- try(
      stats::nlminb(
        start = current_start,

        objective =
          direct_LL3_negative_loglikelihood,

        y = y,
        x = x,
        support_lower = support_lower,

        lower = parameter_lower_bounds,
        upper = parameter_upper_bounds,

        control = list(
          eval.max = 30000,
          iter.max = 15000,
          rel.tol = 1e-11,
          x.tol = 1e-9,
          trace = 0
        )
      ),
      silent = TRUE
    )

    if (
      !inherits(current_fit, "try-error") &&
        current_fit$convergence == 0 &&
        length(current_fit$par) == 5L &&
        all(is.finite(current_fit$par)) &&
        is.finite(current_fit$objective) &&
        current_fit$objective < 1e99
    ) {
      candidates[[length(candidates) + 1L]] <-
        list(
          theta = as.numeric(current_fit$par),
          objective = as.numeric(current_fit$objective)
        )
    }
  }

  if (length(candidates) == 0L) {
    return(
      list(
        success = FALSE,
        message =
          "Every likelihood-polishing attempt failed."
      )
    )
  }

  objectives <- vapply(
    candidates,
    function(candidate) {
      candidate$objective
    },
    numeric(1)
  )

  best <- candidates[[which.min(objectives)]]

  final_fit <- try(
    stats::optim(
      par = best$theta,

      fn =
        direct_LL3_negative_loglikelihood,

      y = y,
      x = x,
      support_lower = support_lower,

      method = "L-BFGS-B",

      lower = parameter_lower_bounds,
      upper = parameter_upper_bounds,

      control = list(
        maxit = 30000,
        factr = 10,
        pgtol = 1e-12,
        trace = 0
      )
    ),
    silent = TRUE
  )

  if (
    !inherits(final_fit, "try-error") &&
      final_fit$convergence == 0 &&
      is.finite(final_fit$value) &&
      final_fit$value <= best$objective
  ) {
    best <- list(
      theta = as.numeric(final_fit$par),
      objective = as.numeric(final_fit$value)
    )
  }

  list(
    success = TRUE,
    theta = best$theta,
    logLik = -best$objective
  )
}


# ============================================================================
# 9. Numerical Hessian and covariance
# ============================================================================

calculate_hessian_covariance <- function(
    theta,
    y,
    x,
    support_lower
) {
  objective_at_theta <- function(theta_value) {
    direct_LL3_negative_loglikelihood(
      theta = theta_value,
      y = y,
      x = x,
      support_lower = support_lower
    )
  }

  gradient <- try(
    numDeriv::grad(
      func = objective_at_theta,
      x = theta,
      method = "Richardson"
    ),
    silent = TRUE
  )

  if (
    inherits(gradient, "try-error") ||
      any(!is.finite(gradient))
  ) {
    return(
      list(
        success = FALSE,
        message =
          "Numerical gradient calculation failed."
      )
    )
  }

  hessian <- try(
    numDeriv::hessian(
      func = objective_at_theta,
      x = theta,
      method = "Richardson"
    ),
    silent = TRUE
  )

  if (
    inherits(hessian, "try-error") ||
      any(!is.finite(hessian))
  ) {
    return(
      list(
        success = FALSE,
        message =
          "Numerical Hessian calculation failed."
      )
    )
  }

  hessian <- (
    hessian +
      t(hessian)
  ) / 2

  eigenvalues <- try(
    eigen(
      hessian,
      symmetric = TRUE,
      only.values = TRUE
    )$values,
    silent = TRUE
  )

  if (
    inherits(eigenvalues, "try-error") ||
      any(!is.finite(eigenvalues)) ||
      any(eigenvalues <= 0)
  ) {
    return(
      list(
        success = FALSE,
        message =
          "The numerical Hessian was not positive definite."
      )
    )
  }

  condition_number <-
    max(eigenvalues) /
    min(eigenvalues)

  if (
    !is.finite(condition_number) ||
      condition_number > 1e12
  ) {
    return(
      list(
        success = FALSE,
        message =
          "The numerical Hessian was excessively ill-conditioned."
      )
    )
  }

  covariance <- try(
    chol2inv(
      chol(hessian)
    ),
    silent = TRUE
  )

  if (
    inherits(covariance, "try-error") ||
      any(!is.finite(covariance))
  ) {
    return(
      list(
        success = FALSE,
        message =
          "The numerical Hessian could not be inverted."
      )
    )
  }

  standard_errors <- sqrt(
    diag(covariance)
  )

  if (
    any(!is.finite(standard_errors)) ||
      any(standard_errors <= 0)
  ) {
    return(
      list(
        success = FALSE,
        message =
          "The numerical-Hessian standard errors were invalid."
      )
    )
  }

  list(
    success = TRUE,
    covariance = covariance,
    standard_errors = standard_errors,

    maximum_absolute_gradient =
      max(abs(gradient)),

    minimum_hessian_eigenvalue =
      min(eigenvalues),

    maximum_hessian_eigenvalue =
      max(eigenvalues),

    condition_number =
      condition_number
  )
}


# ============================================================================
# 10. Empty result row
# ============================================================================

new_result_row <- function(
    n,
    seed
) {
  data.frame(
    n = as.integer(n),
    seed = as.integer(seed),

    complete = FALSE,

    fit_success = FALSE,
    optimizer_success = FALSE,
    hessian_success = FALSE,

    fit_retried = FALSE,

    logLik = NA_real_,
    likelihood_improvement = NA_real_,

    sigma_slope_estimate = NA_real_,
    sigma_slope_standard_error = NA_real_,
    sigma_slope_lower = NA_real_,
    sigma_slope_upper = NA_real_,
    sigma_slope_covered = NA,

    nu_slope_estimate = NA_real_,
    nu_slope_standard_error = NA_real_,
    nu_slope_lower = NA_real_,
    nu_slope_upper = NA_real_,
    nu_slope_covered = NA,

    maximum_absolute_gradient = NA_real_,
    hessian_minimum_eigenvalue = NA_real_,
    hessian_condition_number = NA_real_,

    warning = NA_character_,
    error_stage = NA_character_,
    error_message = NA_character_,

    stringsAsFactors = FALSE
  )
}


# ============================================================================
# 11. Run simulations
# ============================================================================

number_of_rows <-
  repetitions *
  length(sample_sizes)

rows <- vector(
  "list",
  number_of_rows
)

row_index <- 0L

for (current_n in sample_sizes) {
  for (replication in seq_len(repetitions)) {
    row_index <- row_index + 1L

    current_seed <- simulation_seeds[
      row_index
    ]

    result <- new_result_row(
      n = current_n,
      seed = current_seed
    )

    set.seed(current_seed)

    x <- seq(
      -1,
      1,
      length.out = current_n
    )

    sigma_true <- exp(
      truth$sigma_intercept +
        truth$sigma_slope * x
    )

    nu_true <- exp(
      truth$nu_intercept +
        truth$nu_slope * x
    )

    y <- rLL3(
      n = current_n,
      mu = truth$mu,
      sigma = sigma_true,
      nu = nu_true
    )

    simulation_data <- data.frame(
      y = y,
      x = x
    )

    support_lower <- LL3_support_lower(y)

    fitted <- fit_joint_model(
      data = simulation_data,
      lower = support_lower
    )

    if (!fitted$success) {
      result$error_stage <- fitted$stage
      result$error_message <- fitted$message

      rows[[row_index]] <- result
      next
    }

    result$fit_success <- TRUE
    result$fit_retried <- fitted$retried
    result$warning <- fitted$warning

    initial_logLik <-
      -direct_LL3_negative_loglikelihood(
        theta = fitted$theta,
        y = y,
        x = x,
        support_lower = support_lower
      )

    polished <- polish_likelihood(
      theta_start = fitted$theta,
      y = y,
      x = x,
      support_lower = support_lower
    )

    if (!polished$success) {
      result$error_stage <-
        "likelihood polishing"

      result$error_message <-
        polished$message

      rows[[row_index]] <- result
      next
    }

    result$optimizer_success <- TRUE
    result$logLik <- polished$logLik

    result$likelihood_improvement <-
      polished$logLik -
      initial_logLik

    hessian_result <- calculate_hessian_covariance(
      theta = polished$theta,
      y = y,
      x = x,
      support_lower = support_lower
    )

    if (!hessian_result$success) {
      result$error_stage <-
        "numerical Hessian"

      result$error_message <-
        hessian_result$message

      rows[[row_index]] <- result
      next
    }

    result$hessian_success <- TRUE

    result$maximum_absolute_gradient <-
      hessian_result$maximum_absolute_gradient

    result$hessian_minimum_eigenvalue <-
      hessian_result$minimum_hessian_eigenvalue

    result$hessian_condition_number <-
      hessian_result$condition_number

    standard_errors <-
      hessian_result$standard_errors

    sigma_slope_estimate <-
      polished$theta[3L]

    sigma_slope_se <-
      standard_errors[3L]

    sigma_slope_lower <-
      sigma_slope_estimate -
      critical_value *
      sigma_slope_se

    sigma_slope_upper <-
      sigma_slope_estimate +
      critical_value *
      sigma_slope_se

    nu_slope_estimate <-
      polished$theta[5L]

    nu_slope_se <-
      standard_errors[5L]

    nu_slope_lower <-
      nu_slope_estimate -
      critical_value *
      nu_slope_se

    nu_slope_upper <-
      nu_slope_estimate +
      critical_value *
      nu_slope_se

    result$sigma_slope_estimate <-
      sigma_slope_estimate

    result$sigma_slope_standard_error <-
      sigma_slope_se

    result$sigma_slope_lower <-
      sigma_slope_lower

    result$sigma_slope_upper <-
      sigma_slope_upper

    result$sigma_slope_covered <-
      sigma_slope_lower <=
      truth$sigma_slope &&
      sigma_slope_upper >=
      truth$sigma_slope

    result$nu_slope_estimate <-
      nu_slope_estimate

    result$nu_slope_standard_error <-
      nu_slope_se

    result$nu_slope_lower <-
      nu_slope_lower

    result$nu_slope_upper <-
      nu_slope_upper

    result$nu_slope_covered <-
      nu_slope_lower <=
      truth$nu_slope &&
      nu_slope_upper >=
      truth$nu_slope

    result$complete <- TRUE

    rows[[row_index]] <- result
  }

  message(
    "Completed n = ",
    current_n
  )
}


# ============================================================================
# 12. Combine raw results
# ============================================================================

coverage_results <- do.call(
  rbind,
  rows
)

rownames(coverage_results) <- NULL


# ============================================================================
# 13. Summarize coverage
# ============================================================================

summary_rows <- lapply(
  sample_sizes,
  function(current_n) {
    all_runs <- coverage_results[
      coverage_results$n == current_n,
      ,
      drop = FALSE
    ]

    complete_runs <- all_runs[
      all_runs$complete %in% TRUE,
      ,
      drop = FALSE
    ]

    complete_count <- nrow(
      complete_runs
    )

    sigma_coverage <- safe_rate(
      complete_runs$sigma_slope_covered
    )

    nu_coverage <- safe_rate(
      complete_runs$nu_slope_covered
    )

    data.frame(
      n = current_n,

      attempted =
        nrow(all_runs),

      complete =
        complete_count,

      completion_rate =
        safe_rate(all_runs$complete),

      fit_success_rate =
        safe_rate(all_runs$fit_success),

      optimizer_success_rate =
        safe_rate(all_runs$optimizer_success),

      hessian_success_rate =
        safe_rate(all_runs$hessian_success),

      retry_rate =
        safe_rate(all_runs$fit_retried),

      sigma_slope_bias =
        safe_mean(
          complete_runs$
            sigma_slope_estimate -
            truth$sigma_slope
        ),

      sigma_slope_rmse =
        sqrt(
          safe_mean(
            (
              complete_runs$
                sigma_slope_estimate -
                truth$sigma_slope
            )^2
          )
        ),

      sigma_slope_empirical_sd =
        safe_sd(
          complete_runs$
            sigma_slope_estimate
        ),

      sigma_slope_mean_se =
        safe_mean(
          complete_runs$
            sigma_slope_standard_error
        ),

      sigma_slope_coverage =
        sigma_coverage,

      sigma_coverage_MCSE =
        if (
          is.finite(sigma_coverage) &&
            complete_count > 0L
        ) {
          sqrt(
            sigma_coverage *
              (1 - sigma_coverage) /
              complete_count
          )
        } else {
          NA_real_
        },

      sigma_slope_mean_width =
        safe_mean(
          complete_runs$
            sigma_slope_upper -
            complete_runs$
              sigma_slope_lower
        ),

      nu_slope_bias =
        safe_mean(
          complete_runs$
            nu_slope_estimate -
            truth$nu_slope
        ),

      nu_slope_rmse =
        sqrt(
          safe_mean(
            (
              complete_runs$
                nu_slope_estimate -
                truth$nu_slope
            )^2
          )
        ),

      nu_slope_empirical_sd =
        safe_sd(
          complete_runs$
            nu_slope_estimate
        ),

      nu_slope_mean_se =
        safe_mean(
          complete_runs$
            nu_slope_standard_error
        ),

      nu_slope_coverage =
        nu_coverage,

      nu_coverage_MCSE =
        if (
          is.finite(nu_coverage) &&
            complete_count > 0L
        ) {
          sqrt(
            nu_coverage *
              (1 - nu_coverage) /
              complete_count
          )
        } else {
          NA_real_
        },

      nu_slope_mean_width =
        safe_mean(
          complete_runs$
            nu_slope_upper -
            complete_runs$
              nu_slope_lower
        ),

      median_likelihood_improvement =
        safe_median(
          complete_runs$
            likelihood_improvement
        ),

      maximum_likelihood_improvement =
        if (
          nrow(complete_runs) > 0L
        ) {
          max(
            complete_runs$
              likelihood_improvement,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },

      median_maximum_gradient =
        safe_median(
          complete_runs$
            maximum_absolute_gradient
        ),

      median_hessian_condition_number =
        safe_median(
          complete_runs$
            hessian_condition_number
        ),

      stringsAsFactors = FALSE
    )
  }
)

coverage_summary <- do.call(
  rbind,
  summary_rows
)

rownames(coverage_summary) <- NULL


# ============================================================================
# 14. Failures
# ============================================================================

coverage_failures <- coverage_results[
  !coverage_results$complete,
  c(
    "n",
    "seed",

    "fit_success",
    "optimizer_success",
    "hessian_success",

    "fit_retried",

    "warning",
    "error_stage",
    "error_message"
  ),
  drop = FALSE
]


# ============================================================================
# 15. Save output
# ============================================================================

raw_output_file <- file.path(
  results_dir,
  "interval_coverage_numerical_hessian_raw.csv"
)

summary_output_file <- file.path(
  results_dir,
  "interval_coverage_numerical_hessian_summary.csv"
)

failure_output_file <- file.path(
  results_dir,
  "interval_coverage_numerical_hessian_failures.csv"
)

utils::write.csv(
  coverage_results,
  raw_output_file,
  row.names = FALSE
)

utils::write.csv(
  coverage_summary,
  summary_output_file,
  row.names = FALSE
)

utils::write.csv(
  coverage_failures,
  failure_output_file,
  row.names = FALSE
)


# ============================================================================
# 16. Print results
# ============================================================================

cat(
  "\n============================================================\n",
  "Numerical-Hessian interval-coverage summary\n",
  "============================================================\n",
  sep = ""
)

print(
  coverage_summary,
  digits = 12,
  row.names = FALSE
)

cat(
  "\nIncomplete simulations: ",
  nrow(coverage_failures),
  "\n",
  sep = ""
)

if (nrow(coverage_failures) > 0L) {
  print(
    coverage_failures,
    row.names = FALSE
  )
}

cat(
  "\nFiles written:\n",
  raw_output_file,
  "\n",
  summary_output_file,
  "\n",
  failure_output_file,
  "\n",
  sep = ""
)
