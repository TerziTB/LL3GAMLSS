# validation/validate_optimizer_and_lower_sensitivity.R
#
# Reviewer-facing LL3 validation:
#
#   1. Identify direct-optimizer likelihood discrepancies from the completed
#      joint nonstationarity simulation.
#
#   2. Recreate those datasets and compare GAMLSS against a completely
#      separate formula-only LL3 likelihood using many optimizer starts.
#
#   3. Test sensitivity to the numerical lower support bound using:
#
#          project default
#          min(y) - 1e-6 * data scale
#          min(y) - 1e-5 * data scale
#          min(y) - 1e-4 * data scale
#          min(y) - 1e-3 * data scale
#
#   4. Refit stationary, sigma-only, nu-only, and joint models at every
#      support-bound setting and compare fitted curves, quantiles, likelihoods,
#      and AIC/BIC selections.
#
# This script does not call dLL3() inside the independent direct likelihood.


# ============================================================================
# 1. User settings
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


# Joint simulation truth.
#
# These values match the joint nonstationarity experiment used previously:
#
#   mu_t = -50
#   log(sigma_t) = log(70) + 0.55 * x_t
#   log(nu_t)    = log(1.8) - 0.20 * x_t
#
# Change these values only if the original joint simulation used different
# values.

truth_mu <- -50

truth_sigma_intercept <- log(70)
truth_sigma_slope <- 0.55

truth_nu_intercept <- log(1.8)
truth_nu_slope <- -0.20


# Original discrepancy criterion.

direct_discrepancy_threshold <- 1e-3


# Number of randomized starts used by the independent direct optimizer.

number_of_direct_random_starts <- 40L


# Number of ordinary non-boundary cases added per sample size to the
# lower-bound sensitivity experiment.

representative_cases_per_n <- 2L


# Maximum number of closest-to-boundary cases added per sample size.

closest_boundary_cases_per_n <- 3L


# Explicit lower-bound offsets.

lower_offset_multipliers <- c(
  1e-6,
  1e-5,
  1e-4,
  1e-3
)


# Conditional quantiles compared in the sensitivity analysis.

quantile_probabilities <- c(
  0.02,
  0.05,
  0.10,
  0.50,
  0.90,
  0.95,
  0.98
)


# ============================================================================
# 2. Load the LL3 project
# ============================================================================

source(
  file.path(
    project_dir,
    "R",
    "LL3_distribution.R"
  )
)

source(
  file.path(
    project_dir,
    "R",
    "LL3_gamlss_family.R"
  )
)

source(
  file.path(
    project_dir,
    "R",
    "LL3_fitting.R"
  )
)

source(
  file.path(
    project_dir,
    "R",
    "LL3_diagnostics.R"
  )
)


if (!requireNamespace("gamlss", quietly = TRUE)) {
  stop(
    "The gamlss package must be installed."
  )
}


required_project_functions <- c(
  "rLL3",
  "qLL3",
  "LL3_support_lower",
  "fit_LL3_gamlss",
  "extract_LL3_parameters"
)


missing_project_functions <- required_project_functions[
  !vapply(
    required_project_functions,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  )
]


if (length(missing_project_functions) > 0L) {
  stop(
    "Missing project functions: ",
    paste(
      missing_project_functions,
      collapse = ", "
    )
  )
}


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


raw_file <- file.path(
  results_dir,
  "joint_nonstationarity_raw.csv"
)


if (!file.exists(raw_file)) {
  stop(
    "Could not find: ",
    raw_file
  )
}


raw_results <- utils::read.csv(
  raw_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (!all(c("n", "seed") %in% names(raw_results))) {
  stop(
    "joint_nonstationarity_raw.csv must contain n and seed columns."
  )
}


# ============================================================================
# 3. General utility functions
# ============================================================================

first_existing_column <- function(
    data,
    candidates,
    required = FALSE
) {
  result <- candidates[
    candidates %in% names(data)
  ]

  if (length(result) > 0L) {
    return(result[1L])
  }

  if (required) {
    stop(
      "None of the required columns was found: ",
      paste(
        candidates,
        collapse = ", "
      ),
      "\n\nAvailable columns are:\n",
      paste(
        names(data),
        collapse = ", "
      )
    )
  }

  NA_character_
}


as_logical_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  if (is.numeric(x)) {
    return(
      !is.na(x) &
        x != 0
    )
  }

  normalized <- toupper(
    trimws(
      as.character(x)
    )
  )

  normalized %in% c(
    "TRUE",
    "T",
    "YES",
    "Y",
    "1"
  )
}


safe_mean <- function(x) {
  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(x)
}


safe_median <- function(x) {
  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  stats::median(x)
}


safe_max <- function(x) {
  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  max(x)
}


safe_rate <- function(x) {
  x <- x[
    !is.na(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(
    as.logical(x)
  )
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


data_scale_function <- function(y) {
  spread_values <- c(
    diff(range(y)),
    stats::IQR(y),
    stats::mad(y),
    stats::sd(y),
    1
  )

  max(
    spread_values[
      is.finite(spread_values)
    ]
  )
}


stable_log1pexp <- function(x) {
  result <- numeric(
    length(x)
  )

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


extract_numeric_value <- function(
    row,
    column_name
) {
  if (
    is.na(column_name) ||
      !column_name %in% names(row)
  ) {
    return(NA_real_)
  }

  value <- suppressWarnings(
    as.numeric(
      row[[column_name]][1L]
    )
  )

  if (
    length(value) != 1L ||
      !is.finite(value)
  ) {
    return(NA_real_)
  }

  value
}


# ============================================================================
# 4. Recreate joint nonstationary datasets
# ============================================================================

simulate_joint_dataset <- function(
    n,
    seed
) {
  n <- as.integer(n)
  seed <- as.integer(seed)

  if (
    is.na(n) ||
      n < 3L ||
      is.na(seed)
  ) {
    stop(
      "Invalid n or seed."
    )
  }

  set.seed(seed)

  x <- seq(
    -1,
    1,
    length.out = n
  )

  mu_true <- rep(
    truth_mu,
    n
  )

  sigma_true <- exp(
    truth_sigma_intercept +
      truth_sigma_slope * x
  )

  nu_true <- exp(
    truth_nu_intercept +
      truth_nu_slope * x
  )

  y <- rLL3(
    n = n,
    mu = mu_true,
    sigma = sigma_true,
    nu = nu_true
  )

  list(
    data = data.frame(
      y = y,
      x = x,
      time = x
    ),

    y = y,
    x = x,

    mu_true = mu_true,
    sigma_true = sigma_true,
    nu_true = nu_true
  )
}


# ============================================================================
# 5. GAMLSS controls
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
# 6. Fit-capture and validation helpers
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

      invokeRestart(
        "muffleWarning"
      )
    }
  )

  list(
    fit = fitted_object,

    warnings = collapse_messages(
      warning_messages
    )
  )
}


validate_fitted_model <- function(
    fit,
    y,
    lower
) {
  result <- list(
    success = FALSE,
    parameters = NULL,
    message = NA_character_
  )

  if (inherits(fit, "error")) {
    result$message <- conditionMessage(fit)
    return(result)
  }

  if (inherits(fit, "try-error")) {
    result$message <- as.character(fit)
    return(result)
  }

  parameters <- try(
    extract_LL3_parameters(fit),
    silent = TRUE
  )

  if (inherits(parameters, "try-error")) {
    result$message <- as.character(parameters)
    return(result)
  }

  mu <- as.numeric(parameters$mu)
  sigma <- as.numeric(parameters$sigma)
  nu <- as.numeric(parameters$nu)

  expected_length <- length(y)

  if (
    length(mu) != expected_length ||
      length(sigma) != expected_length ||
      length(nu) != expected_length
  ) {
    result$message <-
      "The fitted parameter vectors have incorrect lengths."

    return(result)
  }

  numerical_scale <- max(
    1,
    abs(lower),
    abs(y),
    abs(mu)
  )

  boundary_tolerance <-
    1000 *
    .Machine$double.eps *
    numerical_scale

  checks <- c(
    converged =
      isTRUE(fit$converged),

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
      all(
        y > mu
      ),

    sigma_positive =
      all(
        sigma > 0
      ),

    nu_positive =
      all(
        nu > 0
      ),

    lower_not_materially_violated =
      all(
        mu <=
          lower +
          boundary_tolerance
      )
  )

  failed_checks <- names(checks)[
    !checks
  ]

  result$parameters <- list(
    mu = mu,
    sigma = sigma,
    nu = nu
  )

  result$success <-
    length(failed_checks) == 0L

  if (!result$success) {
    result$message <- paste0(
      "Failed checks: ",
      paste(
        failed_checks,
        collapse = ", "
      ),
      "."
    )
  }

  result
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
      arguments$mu.start <- mu.start
    }

    if (!is.null(sigma.start)) {
      arguments$sigma.start <- sigma.start
    }

    if (!is.null(nu.start)) {
      arguments$nu.start <- nu.start
    }

    arguments
  }

  first_attempt <- capture_fit(
    construct_arguments(
      primary_control
    )
  )

  first_validation <- validate_fitted_model(
    fit = first_attempt$fit,
    y = data$y,
    lower = lower
  )

  if (first_validation$success) {
    return(
      list(
        success = TRUE,
        fit = first_attempt$fit,
        parameters = first_validation$parameters,
        retried = FALSE,
        warnings = first_attempt$warnings,
        message = NA_character_
      )
    )
  }

  second_attempt <- capture_fit(
    construct_arguments(
      retry_control
    )
  )

  second_validation <- validate_fitted_model(
    fit = second_attempt$fit,
    y = data$y,
    lower = lower
  )

  list(
    success = second_validation$success,
    fit = second_attempt$fit,
    parameters = second_validation$parameters,
    retried = TRUE,

    warnings = collapse_messages(
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
# 7. Fit all four candidate models
# ============================================================================

fit_candidate_set <- function(
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
        complete = FALSE,
        stage = "stationary",
        message = stationary$message,
        models = list(
          stationary = stationary
        )
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

    mu.start = stationary_parameters$mu,
    sigma.start = stationary_parameters$sigma,
    nu.start = stationary_parameters$nu
  )

  nu_only <- fit_one_model(
    mu.formula = y ~ 1,
    sigma.formula = ~ 1,
    nu.formula = ~ x,

    data = data,
    lower = lower,

    mu.start = stationary_parameters$mu,
    sigma.start = stationary_parameters$sigma,
    nu.start = stationary_parameters$nu
  )

  joint_sigma_start <-
    if (sigma_only$success) {
      sigma_only$parameters$sigma
    } else {
      stationary_parameters$sigma
    }

  joint_nu_start <-
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

    mu.start = stationary_parameters$mu,
    sigma.start = joint_sigma_start,
    nu.start = joint_nu_start
  )

  model_results <- list(
    stationary = stationary,
    sigma = sigma_only,
    nu = nu_only,
    joint = joint
  )

  model_success <- vapply(
    model_results,
    function(x) {
      isTRUE(x$success)
    },
    logical(1)
  )

  if (!all(model_success)) {
    failed_models <- names(model_success)[
      !model_success
    ]

    failure_messages <- vapply(
      model_results[
        failed_models
      ],
      function(x) {
        if (
          is.null(x$message) ||
            is.na(x$message)
        ) {
          "unknown failure"
        } else {
          x$message
        }
      },
      character(1)
    )

    return(
      list(
        complete = FALSE,
        stage = "candidate models",

        message = paste0(
          paste(
            failed_models,
            failure_messages,
            sep = ": ",
            collapse = " | "
          )
        ),

        models = model_results
      )
    )
  }

  fitted_models <- lapply(
    model_results,
    function(x) {
      x$fit
    }
  )

  AIC_values <- try(
    vapply(
      fitted_models,
      stats::AIC,
      numeric(1)
    ),
    silent = TRUE
  )

  deviances <- try(
    vapply(
      fitted_models,
      stats::deviance,
      numeric(1)
    ),
    silent = TRUE
  )

  degrees_of_freedom <- try(
    vapply(
      fitted_models,
      function(model) {
        as.numeric(
          model$df.fit
        )
      },
      numeric(1)
    ),
    silent = TRUE
  )

  if (
    inherits(AIC_values, "try-error") ||
      inherits(deviances, "try-error") ||
      inherits(degrees_of_freedom, "try-error")
  ) {
    return(
      list(
        complete = FALSE,
        stage = "information criteria",
        message = "AIC or BIC calculation failed.",
        models = model_results
      )
    )
  }

  BIC_values <-
    deviances +
    log(nrow(data)) *
    degrees_of_freedom

  if (
    any(!is.finite(AIC_values)) ||
      any(!is.finite(BIC_values))
  ) {
    return(
      list(
        complete = FALSE,
        stage = "information criteria",
        message = "AIC or BIC contained non-finite values.",
        models = model_results
      )
    )
  }

  list(
    complete = TRUE,
    stage = NA_character_,
    message = NA_character_,

    models = model_results,

    fits = fitted_models,

    AIC = AIC_values,
    BIC = BIC_values,

    AIC_selected = names(
      which.min(AIC_values)
    ),

    BIC_selected = names(
      which.min(BIC_values)
    )
  )
}


# ============================================================================
# 8. Extract joint-model coefficient curves
# ============================================================================

extract_joint_curves <- function(
    fitted_model,
    x
) {
  parameters <- extract_LL3_parameters(
    fitted_model
  )

  mu <- as.numeric(
    parameters$mu
  )

  sigma <- as.numeric(
    parameters$sigma
  )

  nu <- as.numeric(
    parameters$nu
  )

  design_matrix <- cbind(
    intercept = 1,
    slope = x
  )

  sigma_coefficients <- as.numeric(
    qr.coef(
      qr(design_matrix),
      log(sigma)
    )
  )

  nu_coefficients <- as.numeric(
    qr.coef(
      qr(design_matrix),
      log(nu)
    )
  )

  list(
    mu = mean(mu),

    sigma_intercept =
      sigma_coefficients[1L],

    sigma_slope =
      sigma_coefficients[2L],

    nu_intercept =
      nu_coefficients[1L],

    nu_slope =
      nu_coefficients[2L],

    sigma = sigma,
    nu = nu
  )
}


LL3_quantile_matrix_formula <- function(
    probabilities,
    mu,
    sigma,
    nu
) {
  result <- vapply(
    probabilities,
    function(probability) {
      mu +
        sigma *
        exp(
          log(
            probability /
              (1 - probability)
          ) /
            nu
        )
    },
    numeric(length(sigma))
  )

  if (is.null(dim(result))) {
    result <- matrix(
      result,
      ncol = 1L
    )
  }

  result
}


# ============================================================================
# 9. Completely independent formula-only LL3 likelihood
# ============================================================================

direct_LL3_negative_loglikelihood <- function(
    theta,
    y,
    x,
    lower
) {
  theta <- as.numeric(theta)

  if (
    length(theta) != 5L ||
      any(!is.finite(theta))
  ) {
    return(1e100)
  }

  alpha_mu <- theta[1L]

  if (
    alpha_mu < -50 ||
      alpha_mu > 50
  ) {
    return(1e100)
  }

  threshold_gap <- exp(
    alpha_mu
  )

  mu <- lower -
    threshold_gap

  eta_sigma <-
    theta[2L] +
    theta[3L] * x

  eta_nu <-
    theta[4L] +
    theta[5L] * x

  if (
    any(!is.finite(eta_sigma)) ||
      any(!is.finite(eta_nu)) ||
      any(abs(eta_sigma) > 30) ||
      any(abs(eta_nu) > 20)
  ) {
    return(1e100)
  }

  distance <- y - mu

  if (
    any(!is.finite(distance)) ||
      any(distance <= 0)
  ) {
    return(1e100)
  }

  nu <- exp(
    eta_nu
  )

  log_z <-
    log(distance) -
    eta_sigma

  power_argument <-
    nu *
    log_z

  log_density <-
    eta_nu -
    eta_sigma +
    (nu - 1) *
    log_z -
    2 *
    stable_log1pexp(
      power_argument
    )

  if (
    any(!is.finite(log_density))
  ) {
    return(1e100)
  }

  negative_loglikelihood <-
    -sum(
      log_density
    )

  if (!is.finite(negative_loglikelihood)) {
    return(1e100)
  }

  negative_loglikelihood
}


direct_parameter_curves <- function(
    theta,
    x,
    lower
) {
  mu <- lower -
    exp(theta[1L])

  sigma <- exp(
    theta[2L] +
      theta[3L] * x
  )

  nu <- exp(
    theta[4L] +
      theta[5L] * x
  )

  list(
    mu = mu,

    sigma_intercept =
      theta[2L],

    sigma_slope =
      theta[3L],

    nu_intercept =
      theta[4L],

    nu_slope =
      theta[5L],

    sigma = sigma,
    nu = nu
  )
}


# ============================================================================
# 10. Strong multistart direct optimizer
# ============================================================================

run_direct_optimizer <- function(
    y,
    x,
    lower,
    gamlss_curves,
    seed,
    random_starts =
      number_of_direct_random_starts
) {
  data_scale <- data_scale_function(y)

  initial_gap <-
    lower -
    gamlss_curves$mu

  if (
    !is.finite(initial_gap) ||
      initial_gap <= 0
  ) {
    initial_gap <- max(
      data_scale,
      1
    )
  }

  base_start <- c(
    log(initial_gap),

    gamlss_curves$sigma_intercept,
    gamlss_curves$sigma_slope,

    gamlss_curves$nu_intercept,
    gamlss_curves$nu_slope
  )

  generic_start <- c(
    log(
      max(
        data_scale,
        1
      )
    ),

    log(
      max(
        stats::sd(y),
        1
      )
    ),

    0,

    log(2),

    0
  )

  lower_bounds <- c(
    log(
      max(
        1e-12 * data_scale,
        1e-12
      )
    ),

    -20,
    -10,
    -10,
    -10
  )

  upper_bounds <- c(
    log(
      max(
        1e4 * data_scale,
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
        lower_bounds +
          1e-8
      ),
      upper_bounds -
        1e-8
    )
  }

  deterministic_starts <- rbind(
    base_start,
    generic_start,

    base_start +
      c(
        0.25,
        0,
        0,
        0,
        0
      ),

    base_start -
      c(
        0.25,
        0,
        0,
        0,
        0
      ),

    base_start +
      c(
        0,
        0.20,
        0.10,
        0,
        0
      ),

    base_start -
      c(
        0,
        0.20,
        0.10,
        0,
        0
      ),

    base_start +
      c(
        0,
        0,
        0,
        0.20,
        0.10
      ),

    base_start -
      c(
        0,
        0,
        0,
        0.20,
        0.10
      )
  )

  set.seed(
    as.integer(seed) +
      91027L
  )

  random_start_matrix <- matrix(
    stats::rnorm(
      random_starts * 5L
    ),
    ncol = 5L
  )

  perturbation_scales <- c(
    0.75,
    0.50,
    0.40,
    0.50,
    0.40
  )

  random_starts_matrix <- sweep(
    random_start_matrix,
    2L,
    perturbation_scales,
    `*`
  )

  random_starts_matrix <- sweep(
    random_starts_matrix,
    2L,
    base_start,
    `+`
  )

  starts <- rbind(
    deterministic_starts,
    random_starts_matrix
  )

  starts <- t(
    apply(
      starts,
      1L,
      clip_start
    )
  )

  optimizer_candidates <- list()

  add_candidate <- function(
      theta,
      objective,
      converged,
      method
  ) {
    if (
      length(theta) == 5L &&
        all(is.finite(theta)) &&
        is.finite(objective) &&
        objective < 1e99
    ) {
      optimizer_candidates[[length(optimizer_candidates) +
            1L]] <<- list(
        theta = as.numeric(theta),
        objective = as.numeric(objective),
        converged = isTRUE(converged),
        method = method
      )
    }

    invisible(NULL)
  }

  for (
    start_index in
      seq_len(
        nrow(starts)
      )
  ) {
    current_start <-
      starts[
        start_index,
        ,
        drop = TRUE
      ]

    optim_result <- try(
      stats::optim(
        par = current_start,

        fn =
          direct_LL3_negative_loglikelihood,

        y = y,
        x = x,
        lower = lower,

        method = "L-BFGS-B",

        lower = lower_bounds,
        upper = upper_bounds,

        control = list(
          maxit = 20000,
          factr = 1e3,
          pgtol = 1e-10,
          trace = 0
        )
      ),
      silent = TRUE
    )

    if (!inherits(optim_result, "try-error")) {
      add_candidate(
        theta =
          optim_result$par,

        objective =
          optim_result$value,

        converged =
          optim_result$convergence == 0,

        method =
          paste0(
            "L-BFGS-B start ",
            start_index
          )
      )
    }

    nlminb_result <- try(
      stats::nlminb(
        start = current_start,

        objective =
          direct_LL3_negative_loglikelihood,

        y = y,
        x = x,
        lower = lower,

        lower = lower_bounds,
        upper = upper_bounds,

        control = list(
          eval.max = 50000,
          iter.max = 20000,
          rel.tol = 1e-12,
          x.tol = 1e-10,
          trace = 0
        )
      ),
      silent = TRUE
    )

    if (!inherits(nlminb_result, "try-error")) {
      add_candidate(
        theta =
          nlminb_result$par,

        objective =
          nlminb_result$objective,

        converged =
          nlminb_result$convergence == 0,

        method =
          paste0(
            "nlminb start ",
            start_index
          )
      )
    }
  }

  if (length(optimizer_candidates) == 0L) {
    return(
      list(
        success = FALSE,
        message =
          "Every direct-optimizer attempt failed."
      )
    )
  }

  convergence_flags <- vapply(
    optimizer_candidates,
    function(candidate) {
      candidate$converged
    },
    logical(1)
  )

  candidate_pool <-
    if (any(convergence_flags)) {
      which(convergence_flags)
    } else {
      seq_along(
        optimizer_candidates
      )
    }

  candidate_objectives <- vapply(
    optimizer_candidates[
      candidate_pool
    ],
    function(candidate) {
      candidate$objective
    },
    numeric(1)
  )

  best_candidate_index <-
    candidate_pool[
      which.min(
        candidate_objectives
      )
    ]

  best_candidate <-
    optimizer_candidates[[best_candidate_index]]

  final_polish <- try(
    stats::optim(
      par =
        best_candidate$theta,

      fn =
        direct_LL3_negative_loglikelihood,

      y = y,
      x = x,
      lower = lower,

      method = "L-BFGS-B",

      lower = lower_bounds,
      upper = upper_bounds,

      control = list(
        maxit = 50000,
        factr = 10,
        pgtol = 1e-12,
        trace = 0
      )
    ),
    silent = TRUE
  )

  if (
    !inherits(
      final_polish,
      "try-error"
    ) &&
      is.finite(
        final_polish$value
      ) &&
      final_polish$value <
        best_candidate$objective
  ) {
    best_candidate <- list(
      theta =
        final_polish$par,

      objective =
        final_polish$value,

      converged =
        final_polish$convergence == 0,

      method =
        paste0(
          best_candidate$method,
          " + final polish"
        )
    )
  }

  list(
    success = TRUE,

    theta =
      best_candidate$theta,

    logLik =
      -best_candidate$objective,

    converged =
      best_candidate$converged,

    method =
      best_candidate$method,

    attempts =
      length(optimizer_candidates),

    converged_attempts =
      sum(convergence_flags)
  )
}


# ============================================================================
# 11. Identify the original discrepancy rows
# ============================================================================

direct_success_column <- first_existing_column(
  raw_results,
  c(
    "direct_optimizer_success",
    "direct_success"
  )
)


direct_difference_column <- first_existing_column(
  raw_results,
  c(
    "direct_logLik_difference",
    "direct_loglik_difference",
    "direct_abs_logLik_difference",
    "absolute_direct_logLik_difference"
  ),
  required = TRUE
)


raw_gamlss_loglik_column <- first_existing_column(
  raw_results,
  c(
    "direct_gamlss_logLik",
    "gamlss_logLik",
    "joint_logLik"
  )
)


raw_direct_loglik_column <- first_existing_column(
  raw_results,
  c(
    "direct_optimizer_logLik",
    "direct_logLik"
  )
)


direct_difference_values <- suppressWarnings(
  as.numeric(
    raw_results[[direct_difference_column]]
  )
)


direct_success_flags <-
  if (!is.na(direct_success_column)) {
    as_logical_flag(
      raw_results[[direct_success_column]]
    )
  } else {
    is.finite(
      direct_difference_values
    )
  }


discrepancy_rows <- raw_results[
  direct_success_flags &
    is.finite(
      direct_difference_values
    ) &
    abs(
      direct_difference_values
    ) >=
      direct_discrepancy_threshold,
  ,
  drop = FALSE
]


discrepancy_rows <- discrepancy_rows[
  order(
    discrepancy_rows$n,
    discrepancy_rows$seed
  ),
  ,
  drop = FALSE
]


utils::write.csv(
  discrepancy_rows,
  file.path(
    results_dir,
    "joint_direct_discrepancies_original.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n",
  "Original direct-optimizer discrepancies\n",
  "============================================================\n",
  sep = ""
)


if (nrow(discrepancy_rows) == 0L) {
  cat(
    "No rows exceeded the discrepancy threshold of ",
    direct_discrepancy_threshold,
    ".\n",
    sep = ""
  )
} else {
  columns_to_print <- unique(
    c(
      "n",
      "seed",
      raw_gamlss_loglik_column,
      raw_direct_loglik_column,
      direct_difference_column
    )
  )

  columns_to_print <- columns_to_print[
    !is.na(columns_to_print) &
      columns_to_print %in%
      names(discrepancy_rows)
  ]

  print(
    discrepancy_rows[
      ,
      columns_to_print,
      drop = FALSE
    ],
    digits = 15,
    row.names = FALSE
  )
}


# ============================================================================
# 12. Rerun every discrepant case
# ============================================================================

direct_rerun_rows <- list()


if (nrow(discrepancy_rows) > 0L) {
  for (
    discrepancy_index in
      seq_len(
        nrow(discrepancy_rows)
      )
  ) {
    source_row <- discrepancy_rows[
      discrepancy_index,
      ,
      drop = FALSE
    ]

    current_n <- as.integer(
      source_row$n[1L]
    )

    current_seed <- as.integer(
      source_row$seed[1L]
    )

    simulated <- simulate_joint_dataset(
      n = current_n,
      seed = current_seed
    )

    default_lower <- LL3_support_lower(
      simulated$y
    )

    candidate_set <- fit_candidate_set(
      data = simulated$data,
      lower = default_lower
    )

    empty_result <- data.frame(
      n = current_n,
      seed = current_seed,

      gamlss_success = FALSE,
      direct_success = FALSE,

      raw_gamlss_logLik =
        extract_numeric_value(
          source_row,
          raw_gamlss_loglik_column
        ),

      raw_direct_logLik =
        extract_numeric_value(
          source_row,
          raw_direct_loglik_column
        ),

      raw_absolute_difference =
        extract_numeric_value(
          source_row,
          direct_difference_column
        ),

      refitted_gamlss_logLik = NA_real_,
      strengthened_direct_logLik = NA_real_,

      direct_minus_gamlss = NA_real_,
      absolute_logLik_difference = NA_real_,

      likelihood_winner = NA_character_,

      gamlss_mu = NA_real_,
      direct_mu = NA_real_,
      absolute_mu_difference = NA_real_,

      gamlss_sigma_intercept = NA_real_,
      direct_sigma_intercept = NA_real_,

      gamlss_sigma_slope = NA_real_,
      direct_sigma_slope = NA_real_,

      gamlss_nu_intercept = NA_real_,
      direct_nu_intercept = NA_real_,

      gamlss_nu_slope = NA_real_,
      direct_nu_slope = NA_real_,

      maximum_absolute_log_sigma_difference =
        NA_real_,

      maximum_absolute_log_nu_difference =
        NA_real_,

      maximum_scaled_quantile_difference =
        NA_real_,

      agreement_below_1e_3 = NA,
      agreement_below_1e_6 = NA,

      direct_method = NA_character_,
      direct_attempts = NA_integer_,
      direct_converged_attempts = NA_integer_,

      message = NA_character_,

      stringsAsFactors = FALSE
    )

    if (!candidate_set$complete) {
      empty_result$message <- paste0(
        "GAMLSS refit failed at stage ",
        candidate_set$stage,
        ": ",
        candidate_set$message
      )

      direct_rerun_rows[[length(direct_rerun_rows) +
            1L]] <- empty_result

      next
    }

    joint_fit <-
      candidate_set$fits$joint

    gamlss_curves <- extract_joint_curves(
      fitted_model = joint_fit,
      x = simulated$x
    )

    gamlss_logLik <- as.numeric(
      stats::logLik(
        joint_fit
      )
    )

    direct_fit <- run_direct_optimizer(
      y = simulated$y,
      x = simulated$x,
      lower = default_lower,

      gamlss_curves =
        gamlss_curves,

      seed =
        current_seed
    )

    empty_result$gamlss_success <- TRUE
    empty_result$refitted_gamlss_logLik <-
      gamlss_logLik

    empty_result$gamlss_mu <-
      gamlss_curves$mu

    empty_result$gamlss_sigma_intercept <-
      gamlss_curves$sigma_intercept

    empty_result$gamlss_sigma_slope <-
      gamlss_curves$sigma_slope

    empty_result$gamlss_nu_intercept <-
      gamlss_curves$nu_intercept

    empty_result$gamlss_nu_slope <-
      gamlss_curves$nu_slope

    if (!direct_fit$success) {
      empty_result$message <-
        direct_fit$message

      direct_rerun_rows[[length(direct_rerun_rows) +
            1L]] <- empty_result

      next
    }

    direct_curves <- direct_parameter_curves(
      theta = direct_fit$theta,
      x = simulated$x,
      lower = default_lower
    )

    direct_minus_gamlss <-
      direct_fit$logLik -
      gamlss_logLik

    absolute_difference <-
      abs(
        direct_minus_gamlss
      )

    gamlss_quantiles <-
      LL3_quantile_matrix_formula(
        probabilities =
          quantile_probabilities,

        mu =
          gamlss_curves$mu,

        sigma =
          gamlss_curves$sigma,

        nu =
          gamlss_curves$nu
      )

    direct_quantiles <-
      LL3_quantile_matrix_formula(
        probabilities =
          quantile_probabilities,

        mu =
          direct_curves$mu,

        sigma =
          direct_curves$sigma,

        nu =
          direct_curves$nu
      )

    current_data_scale <-
      data_scale_function(
        simulated$y
      )

    likelihood_winner <-
      if (
        direct_minus_gamlss >
          1e-7
      ) {
        "direct optimizer"
      } else if (
        direct_minus_gamlss <
          -1e-7
      ) {
        "GAMLSS"
      } else {
        "numerical tie"
      }

    empty_result$direct_success <- TRUE

    empty_result$strengthened_direct_logLik <-
      direct_fit$logLik

    empty_result$direct_minus_gamlss <-
      direct_minus_gamlss

    empty_result$absolute_logLik_difference <-
      absolute_difference

    empty_result$likelihood_winner <-
      likelihood_winner

    empty_result$direct_mu <-
      direct_curves$mu

    empty_result$absolute_mu_difference <-
      abs(
        direct_curves$mu -
          gamlss_curves$mu
      )

    empty_result$direct_sigma_intercept <-
      direct_curves$sigma_intercept

    empty_result$direct_sigma_slope <-
      direct_curves$sigma_slope

    empty_result$direct_nu_intercept <-
      direct_curves$nu_intercept

    empty_result$direct_nu_slope <-
      direct_curves$nu_slope

    empty_result$
      maximum_absolute_log_sigma_difference <-
      max(
        abs(
          log(
            direct_curves$sigma
          ) -
            log(
              gamlss_curves$sigma
            )
        )
      )

    empty_result$
      maximum_absolute_log_nu_difference <-
      max(
        abs(
          log(
            direct_curves$nu
          ) -
            log(
              gamlss_curves$nu
            )
        )
      )

    empty_result$
      maximum_scaled_quantile_difference <-
      max(
        abs(
          direct_quantiles -
            gamlss_quantiles
        )
      ) /
      current_data_scale

    empty_result$agreement_below_1e_3 <-
      absolute_difference <
      1e-3

    empty_result$agreement_below_1e_6 <-
      absolute_difference <
      1e-6

    empty_result$direct_method <-
      direct_fit$method

    empty_result$direct_attempts <-
      direct_fit$attempts

    empty_result$direct_converged_attempts <-
      direct_fit$converged_attempts

    direct_rerun_rows[[length(direct_rerun_rows) +
          1L]] <- empty_result
  }
}


if (length(direct_rerun_rows) > 0L) {
  direct_rerun_results <- do.call(
    rbind,
    direct_rerun_rows
  )
} else {
  direct_rerun_results <- data.frame(
    n = integer(),
    seed = integer(),
    stringsAsFactors = FALSE
  )
}


utils::write.csv(
  direct_rerun_results,
  file.path(
    results_dir,
    "joint_direct_discrepancies_rerun.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n",
  "Strengthened direct-optimizer reruns\n",
  "============================================================\n",
  sep = ""
)


if (nrow(direct_rerun_results) == 0L) {
  cat(
    "No discrepancy rows required rerunning.\n"
  )
} else {
  print(
    direct_rerun_results,
    digits = 15,
    row.names = FALSE
  )
}


# ============================================================================
# 13. Select lower-bound sensitivity cases
# ============================================================================

boundary_contact_column <- first_existing_column(
  raw_results,
  c(
    "boundary_contact",
    "joint_boundary_contact",
    "mu_boundary_contact"
  )
)


boundary_gap_column <- first_existing_column(
  raw_results,
  c(
    "minimum_lower_minus_mu",
    "min_lower_minus_mu",
    "joint_minimum_lower_minus_mu",
    "minimum_lower_mu_gap"
  )
)


case_parts <- list()


add_case_part <- function(
    data,
    reason
) {
  if (
    is.null(data) ||
      nrow(data) == 0L
  ) {
    return(
      invisible(NULL)
    )
  }

  part <- data.frame(
    n = as.integer(data$n),
    seed = as.integer(data$seed),
    reason = reason,
    stringsAsFactors = FALSE
  )

  case_parts[[length(case_parts) +
        1L]] <<- part

  invisible(NULL)
}


add_case_part(
  discrepancy_rows,
  "direct discrepancy"
)


if (!is.na(boundary_contact_column)) {
  boundary_flags <- as_logical_flag(
    raw_results[[boundary_contact_column]]
  )

  add_case_part(
    raw_results[
      boundary_flags,
      ,
      drop = FALSE
    ],
    "recorded boundary contact"
  )
}


if (!is.na(boundary_gap_column)) {
  boundary_gap_values <- suppressWarnings(
    as.numeric(
      raw_results[[boundary_gap_column]]
    )
  )

  for (
    current_n in
      sort(
        unique(
          raw_results$n
        )
      )
  ) {
    eligible_indices <- which(
      raw_results$n ==
        current_n &
        is.finite(
          boundary_gap_values
        )
    )

    if (length(eligible_indices) > 0L) {
      ordered_indices <- eligible_indices[
        order(
          boundary_gap_values[
            eligible_indices
          ]
        )
      ]

      selected_indices <- head(
        ordered_indices,
        closest_boundary_cases_per_n
      )

      add_case_part(
        raw_results[
          selected_indices,
          ,
          drop = FALSE
        ],
        "closest fitted threshold to lower bound"
      )
    }
  }
}


set.seed(2028L)


for (
  current_n in
    sort(
      unique(
        raw_results$n
      )
    )
) {
  eligible <- raw_results[
    raw_results$n ==
      current_n,
    ,
    drop = FALSE
  ]

  if (nrow(eligible) == 0L) {
    next
  }

  number_to_select <- min(
    representative_cases_per_n,
    nrow(eligible)
  )

  selected_indices <- sample(
    seq_len(
      nrow(eligible)
    ),
    size = number_to_select,
    replace = FALSE
  )

  add_case_part(
    eligible[
      selected_indices,
      ,
      drop = FALSE
    ],
    "representative case"
  )
}


if (length(case_parts) == 0L) {
  stop(
    "No lower-bound sensitivity cases could be selected."
  )
}


case_rows_raw <- do.call(
  rbind,
  case_parts
)


sensitivity_cases <- stats::aggregate(
  reason ~ n + seed,
  data = case_rows_raw,

  FUN = function(x) {
    paste(
      unique(x),
      collapse = "; "
    )
  }
)


sensitivity_cases <- sensitivity_cases[
  order(
    sensitivity_cases$n,
    sensitivity_cases$seed
  ),
  ,
  drop = FALSE
]


rownames(
  sensitivity_cases
) <- NULL


utils::write.csv(
  sensitivity_cases,
  file.path(
    results_dir,
    "lower_bound_sensitivity_cases.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n",
  "Lower-bound sensitivity cases\n",
  "============================================================\n",
  sep = ""
)


print(
  sensitivity_cases,
  row.names = FALSE
)


# ============================================================================
# 14. Run lower-bound sensitivity experiment
# ============================================================================

sensitivity_rows <- list()


for (
  case_index in
    seq_len(
      nrow(
        sensitivity_cases
      )
    )
) {
  current_case <- sensitivity_cases[
    case_index,
    ,
    drop = FALSE
  ]

  current_n <- as.integer(
    current_case$n[1L]
  )

  current_seed <- as.integer(
    current_case$seed[1L]
  )

  current_reason <-
    current_case$reason[1L]

  simulated <- simulate_joint_dataset(
    n = current_n,
    seed = current_seed
  )

  current_scale <- data_scale_function(
    simulated$y
  )

  default_lower <- LL3_support_lower(
    simulated$y
  )

  lower_definitions <- data.frame(
    lower_label = c(
      "project_default",
      paste0(
        "offset_",
        format(
          lower_offset_multipliers,
          scientific = TRUE,
          trim = TRUE
        )
      )
    ),

    lower_value = c(
      default_lower,

      min(
        simulated$y
      ) -
        lower_offset_multipliers *
        current_scale
    ),

    requested_multiplier = c(
      NA_real_,
      lower_offset_multipliers
    ),

    stringsAsFactors = FALSE
  )

  lower_definitions$effective_multiplier <-
    (
      min(
        simulated$y
      ) -
        lower_definitions$lower_value
    ) /
    current_scale

  local_rows <- vector(
    "list",
    nrow(
      lower_definitions
    )
  )

  local_curves <- vector(
    "list",
    nrow(
      lower_definitions
    )
  )

  names(local_curves) <-
    lower_definitions$lower_label

  for (
    lower_index in
      seq_len(
        nrow(
          lower_definitions
        )
      )
  ) {
    lower_label <-
      lower_definitions$
        lower_label[
          lower_index
        ]

    lower_value <-
      lower_definitions$
        lower_value[
          lower_index
        ]

    cat(
      "Sensitivity case ",
      case_index,
      "/",
      nrow(sensitivity_cases),
      ": n = ",
      current_n,
      ", seed = ",
      current_seed,
      ", lower = ",
      lower_label,
      "\n",
      sep = ""
    )

    candidate_set <- fit_candidate_set(
      data = simulated$data,
      lower = lower_value
    )

    result_row <- data.frame(
      n = current_n,
      seed = current_seed,
      case_reason = current_reason,

      lower_label = lower_label,
      lower_value = lower_value,

      requested_multiplier =
        lower_definitions$
          requested_multiplier[
            lower_index
          ],

      effective_multiplier =
        lower_definitions$
          effective_multiplier[
            lower_index
          ],

      complete = FALSE,

      AIC_selected = NA_character_,
      BIC_selected = NA_character_,

      stationary_logLik = NA_real_,
      sigma_logLik = NA_real_,
      nu_logLik = NA_real_,
      joint_logLik = NA_real_,

      joint_mu = NA_real_,

      sigma_intercept = NA_real_,
      sigma_slope = NA_real_,

      nu_intercept = NA_real_,
      nu_slope = NA_real_,

      minimum_lower_minus_mu = NA_real_,
      boundary_tolerance = NA_real_,
      boundary_contact = NA,
      material_boundary_violation = NA,

      delta_joint_logLik_from_default = NA_real_,
      absolute_mu_difference_scaled = NA_real_,

      maximum_absolute_log_sigma_difference =
        NA_real_,

      maximum_absolute_log_nu_difference =
        NA_real_,

      maximum_scaled_quantile_difference =
        NA_real_,

      AIC_selection_same_as_default = NA,
      BIC_selection_same_as_default = NA,

      error_stage =
        candidate_set$stage,

      error_message =
        candidate_set$message,

      stringsAsFactors = FALSE
    )

    if (!candidate_set$complete) {
      local_rows[[lower_index]] <- result_row

      next
    }

    joint_curves <- extract_joint_curves(
      fitted_model =
        candidate_set$fits$joint,

      x =
        simulated$x
    )

    joint_quantiles <-
      LL3_quantile_matrix_formula(
        probabilities =
          quantile_probabilities,

        mu =
          joint_curves$mu,

        sigma =
          joint_curves$sigma,

        nu =
          joint_curves$nu
      )

    lower_minus_mu <-
      lower_value -
      joint_curves$mu

    boundary_tolerance <-
      1000 *
      .Machine$double.eps *
      max(
        1,
        abs(lower_value),
        abs(joint_curves$mu),
        current_scale
      )

    result_row$complete <- TRUE

    result_row$AIC_selected <-
      candidate_set$AIC_selected

    result_row$BIC_selected <-
      candidate_set$BIC_selected

    result_row$stationary_logLik <-
      as.numeric(
        stats::logLik(
          candidate_set$
            fits$
            stationary
        )
      )

    result_row$sigma_logLik <-
      as.numeric(
        stats::logLik(
          candidate_set$
            fits$
            sigma
        )
      )

    result_row$nu_logLik <-
      as.numeric(
        stats::logLik(
          candidate_set$
            fits$
            nu
        )
      )

    result_row$joint_logLik <-
      as.numeric(
        stats::logLik(
          candidate_set$
            fits$
            joint
        )
      )

    result_row$joint_mu <-
      joint_curves$mu

    result_row$sigma_intercept <-
      joint_curves$sigma_intercept

    result_row$sigma_slope <-
      joint_curves$sigma_slope

    result_row$nu_intercept <-
      joint_curves$nu_intercept

    result_row$nu_slope <-
      joint_curves$nu_slope

    result_row$minimum_lower_minus_mu <-
      lower_minus_mu

    result_row$boundary_tolerance <-
      boundary_tolerance

    result_row$boundary_contact <-
      lower_minus_mu <=
      boundary_tolerance

    result_row$material_boundary_violation <-
      lower_minus_mu <
      -boundary_tolerance

    result_row$error_stage <-
      NA_character_

    result_row$error_message <-
      NA_character_

    local_rows[[lower_index]] <- result_row

    local_curves[[lower_label]] <- list(
      mu = joint_curves$mu,
      sigma = joint_curves$sigma,
      nu = joint_curves$nu,
      quantiles = joint_quantiles
    )
  }

  local_results <- do.call(
    rbind,
    local_rows
  )

  default_index <- which(
    local_results$lower_label ==
      "project_default"
  )

  if (
    length(default_index) == 1L &&
      isTRUE(
        local_results$
          complete[
            default_index
          ]
      )
  ) {
    default_curves <-
      local_curves[["project_default"]]

    default_joint_logLik <-
      local_results$
        joint_logLik[
          default_index
        ]

    default_mu <-
      local_results$
        joint_mu[
          default_index
        ]

    default_AIC_selection <-
      local_results$
        AIC_selected[
          default_index
        ]

    default_BIC_selection <-
      local_results$
        BIC_selected[
          default_index
        ]

    for (
      local_index in
        seq_len(
          nrow(local_results)
        )
    ) {
      if (
        !isTRUE(
          local_results$
            complete[
              local_index
            ]
        )
      ) {
        next
      }

      current_label <-
        local_results$
          lower_label[
            local_index
          ]

      current_curves <-
        local_curves[[current_label]]

      local_results$
        delta_joint_logLik_from_default[
          local_index
        ] <-
        local_results$
          joint_logLik[
            local_index
          ] -
        default_joint_logLik

      local_results$
        absolute_mu_difference_scaled[
          local_index
        ] <-
        abs(
          local_results$
            joint_mu[
              local_index
            ] -
            default_mu
        ) /
        current_scale

      local_results$
        maximum_absolute_log_sigma_difference[
          local_index
        ] <-
        max(
          abs(
            log(
              current_curves$sigma
            ) -
              log(
                default_curves$sigma
              )
          )
        )

      local_results$
        maximum_absolute_log_nu_difference[
          local_index
        ] <-
        max(
          abs(
            log(
              current_curves$nu
            ) -
              log(
                default_curves$nu
              )
          )
        )

      local_results$
        maximum_scaled_quantile_difference[
          local_index
        ] <-
        max(
          abs(
            current_curves$quantiles -
              default_curves$quantiles
          )
        ) /
        current_scale

      local_results$
        AIC_selection_same_as_default[
          local_index
        ] <-
        identical(
          local_results$
            AIC_selected[
              local_index
            ],
          default_AIC_selection
        )

      local_results$
        BIC_selection_same_as_default[
          local_index
        ] <-
        identical(
          local_results$
            BIC_selected[
              local_index
            ],
          default_BIC_selection
        )
    }
  }

  for (
    local_index in
      seq_len(
        nrow(
          local_results
        )
      )
  ) {
    sensitivity_rows[[length(sensitivity_rows) +
          1L]] <-
      local_results[
        local_index,
        ,
        drop = FALSE
      ]
  }
}


lower_sensitivity_results <- do.call(
  rbind,
  sensitivity_rows
)


rownames(
  lower_sensitivity_results
) <- NULL


utils::write.csv(
  lower_sensitivity_results,
  file.path(
    results_dir,
    "lower_bound_sensitivity_raw.csv"
  ),
  row.names = FALSE
)


# ============================================================================
# 15. Summarize lower-bound sensitivity
# ============================================================================

sensitivity_summary_rows <- lapply(
  unique(
    lower_sensitivity_results$
      lower_label
  ),

  function(current_label) {
    current_results <-
      lower_sensitivity_results[
        lower_sensitivity_results$
          lower_label ==
          current_label,
        ,
        drop = FALSE
      ]

    complete_results <-
      current_results[
        current_results$complete %in% TRUE,
        ,
        drop = FALSE
      ]

    data.frame(
      lower_label =
        current_label,

      cases =
        nrow(
          current_results
        ),

      complete =
        nrow(
          complete_results
        ),

      completion_rate =
        safe_rate(
          current_results$complete
        ),

      effective_multiplier_median =
        safe_median(
          current_results$
            effective_multiplier
        ),

      AIC_selection_agreement_rate =
        safe_rate(
          complete_results$
            AIC_selection_same_as_default
        ),

      BIC_selection_agreement_rate =
        safe_rate(
          complete_results$
            BIC_selection_same_as_default
        ),

      boundary_contact_rate =
        safe_rate(
          complete_results$
            boundary_contact
        ),

      material_boundary_violation_rate =
        safe_rate(
          complete_results$
            material_boundary_violation
        ),

      maximum_absolute_delta_joint_logLik =
        safe_max(
          abs(
            complete_results$
              delta_joint_logLik_from_default
          )
        ),

      median_absolute_delta_joint_logLik =
        safe_median(
          abs(
            complete_results$
              delta_joint_logLik_from_default
          )
        ),

      maximum_scaled_mu_difference =
        safe_max(
          complete_results$
            absolute_mu_difference_scaled
        ),

      maximum_log_sigma_curve_difference =
        safe_max(
          complete_results$
            maximum_absolute_log_sigma_difference
        ),

      maximum_log_nu_curve_difference =
        safe_max(
          complete_results$
            maximum_absolute_log_nu_difference
        ),

      maximum_scaled_quantile_difference =
        safe_max(
          complete_results$
            maximum_scaled_quantile_difference
        ),

      median_scaled_quantile_difference =
        safe_median(
          complete_results$
            maximum_scaled_quantile_difference
        ),

      stringsAsFactors = FALSE
    )
  }
)


lower_sensitivity_summary <- do.call(
  rbind,
  sensitivity_summary_rows
)


rownames(
  lower_sensitivity_summary
) <- NULL


utils::write.csv(
  lower_sensitivity_summary,
  file.path(
    results_dir,
    "lower_bound_sensitivity_summary.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n",
  "Lower-bound sensitivity summary\n",
  "============================================================\n",
  sep = ""
)


print(
  lower_sensitivity_summary,
  digits = 12,
  row.names = FALSE
)


# ============================================================================
# 16. Basic result checks
# ============================================================================

if (nrow(direct_rerun_results) > 0L) {
  unresolved_direct_cases <-
    direct_rerun_results[
      direct_rerun_results$direct_success %in% TRUE &
        is.finite(
          direct_rerun_results$
            absolute_logLik_difference
        ) &
        direct_rerun_results$
          absolute_logLik_difference >=
          direct_discrepancy_threshold,
      ,
      drop = FALSE
    ]

  cat(
    "\nDirect cases still exceeding ",
    direct_discrepancy_threshold,
    ": ",
    nrow(unresolved_direct_cases),
    "\n",
    sep = ""
  )

  if (nrow(unresolved_direct_cases) > 0L) {
    print(
      unresolved_direct_cases[
        ,
        c(
          "n",
          "seed",
          "refitted_gamlss_logLik",
          "strengthened_direct_logLik",
          "direct_minus_gamlss",
          "absolute_logLik_difference",
          "likelihood_winner",
          "maximum_scaled_quantile_difference"
        ),
        drop = FALSE
      ],
      digits = 15,
      row.names = FALSE
    )
  }
}


incomplete_sensitivity_cases <-
  lower_sensitivity_results[
    !lower_sensitivity_results$complete,
    ,
    drop = FALSE
  ]


cat(
  "\nIncomplete lower-bound fits: ",
  nrow(incomplete_sensitivity_cases),
  "\n",
  sep = ""
)


if (nrow(incomplete_sensitivity_cases) > 0L) {
  print(
    incomplete_sensitivity_cases[
      ,
      c(
        "n",
        "seed",
        "lower_label",
        "error_stage",
        "error_message"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}


# ============================================================================
# 17. Output locations
# ============================================================================

cat(
  "\n============================================================\n",
  "Files written\n",
  "============================================================\n",
  sep = ""
)


output_files <- c(
  file.path(
    results_dir,
    "joint_direct_discrepancies_original.csv"
  ),

  file.path(
    results_dir,
    "joint_direct_discrepancies_rerun.csv"
  ),

  file.path(
    results_dir,
    "lower_bound_sensitivity_cases.csv"
  ),

  file.path(
    results_dir,
    "lower_bound_sensitivity_raw.csv"
  ),

  file.path(
    results_dir,
    "lower_bound_sensitivity_summary.csv"
  )
)


cat(
  paste0(
    output_files,
    collapse = "\n"
  ),
  "\n"
)


cat(
  "\nValidation script completed.\n"
)


# Run from the project root with:
#
# source("validation/validate_optimizer_and_lower_sensitivity.R")
