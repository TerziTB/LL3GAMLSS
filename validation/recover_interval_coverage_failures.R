# validation/recover_interval_coverage_failures.R
#
# Targeted recovery of incomplete interval-coverage simulations.
#
# This script:
#
#   1. reads the existing interval-coverage output;
#   2. reruns only incomplete seeds;
#   3. uses a multistart formula-only LL3 optimizer;
#   4. does not require an optimizer convergence code of zero;
#   5. instead requires a small numerical gradient, a positive-definite
#      Hessian, acceptable conditioning, and finite standard errors;
#   6. warm-starts GAMLSS from the recovered optimum when the original
#      GAMLSS fit failed;
#   7. merges recovered results with the original completed simulations;
#   8. writes updated coverage summaries without overwriting the originals.


# ============================================================================
# 1. Project setup
# ============================================================================

resolve_project_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0L) {
    script_path <- normalizePath(
      sub("^--file=", "", script_arg[[1L]]),
      winslash = "/",
      mustWork = TRUE
    )
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
  stop("The gamlss package must be installed.")
}

if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop(
    "The numDeriv package must be installed. Run:\n",
    "install.packages('numDeriv')"
  )
}


# ============================================================================
# 2. Input files
# ============================================================================

raw_input_file <- file.path(
  results_dir,
  "interval_coverage_numerical_hessian_raw.csv"
)

if (!file.exists(raw_input_file)) {
  stop(
    "Missing file:\n",
    raw_input_file
  )
}

coverage_results <- utils::read.csv(
  raw_input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "n",
  "seed",
  "complete",
  "fit_success",
  "optimizer_success",
  "hessian_success",
  "fit_retried",
  "logLik",
  "likelihood_improvement",
  "sigma_slope_estimate",
  "sigma_slope_standard_error",
  "sigma_slope_lower",
  "sigma_slope_upper",
  "sigma_slope_covered",
  "nu_slope_estimate",
  "nu_slope_standard_error",
  "nu_slope_lower",
  "nu_slope_upper",
  "nu_slope_covered",
  "maximum_absolute_gradient",
  "hessian_minimum_eigenvalue",
  "hessian_condition_number",
  "warning",
  "error_stage",
  "error_message"
)

missing_columns <- setdiff(
  required_columns,
  names(coverage_results)
)

if (length(missing_columns) > 0L) {
  stop(
    "The raw coverage file is missing these columns:\n",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ============================================================================
# 3. Simulation specification
# ============================================================================

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

maximum_gradient_tolerance <- 1e-3
maximum_hessian_condition_number <- 1e12

number_of_random_starts <- 24L
number_of_candidates_to_diagnose <- 20L
number_of_candidates_to_polish <- 6L


# ============================================================================
# 4. General utilities
# ============================================================================

as_logical_safe <- function(x) {
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


safe_sd <- function(x) {
  x <- x[
    is.finite(x)
  ]

  if (length(x) < 2L) {
    return(NA_real_)
  }

  stats::sd(x)
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


# ============================================================================
# 5. Recreate one simulation
# ============================================================================

simulate_coverage_dataset <- function(
    n,
    seed
) {
  n <- as.integer(n)
  seed <- as.integer(seed)

  set.seed(seed)

  x <- seq(
    -1,
    1,
    length.out = n
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
    n = n,
    mu = truth$mu,
    sigma = sigma_true,
    nu = nu_true
  )

  list(
    y = y,
    x = x,

    data = data.frame(
      y = y,
      x = x
    )
  )
}


# ============================================================================
# 6. Independent formula-only LL3 likelihood
#
# Parameter order:
#
#   theta[1] = alpha_mu
#   theta[2] = sigma intercept
#   theta[3] = sigma slope
#   theta[4] = nu intercept
#   theta[5] = nu slope
#
# where:
#
#   mu = support_lower - exp(alpha_mu)
# ============================================================================

make_LL3_objective <- function(
    y,
    x,
    support_lower
) {
  force(y)
  force(x)
  force(support_lower)

  function(theta) {
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

    mu <-
      support_lower -
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

    objective <- -sum(
      log_density
    )

    if (!is.finite(objective)) {
      return(1e100)
    }

    objective
  }
}


# ============================================================================
# 7. Starting values and parameter bounds
# ============================================================================

construct_optimizer_setup <- function(
    y,
    x,
    support_lower,
    seed
) {
  scale_y <- data_scale(y)

  starting_values <- .LL3_start_values(y)

  initial_gap <-
    support_lower -
    as.numeric(
      starting_values["mu"]
    )

  if (
    !is.finite(initial_gap) ||
      initial_gap <= 0
  ) {
    initial_gap <- scale_y
  }

  sigma_start <- as.numeric(
    starting_values["sigma"]
  )

  nu_start <- as.numeric(
    starting_values["nu"]
  )

  if (
    !is.finite(sigma_start) ||
      sigma_start <= 0
  ) {
    sigma_start <- max(
      stats::sd(y),
      1
    )
  }

  if (
    !is.finite(nu_start) ||
      nu_start <= 0
  ) {
    nu_start <- 2
  }

  base_start <- c(
    log(initial_gap),
    log(sigma_start),
    0,
    log(nu_start),
    0
  )

  alternative_start <- c(
    log(scale_y),
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
        as.numeric(start),
        parameter_lower_bounds + 1e-8
      ),
      parameter_upper_bounds - 1e-8
    )
  }

  deterministic_starts <- rbind(
    base_start,
    alternative_start,

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
        0.25,
        0,
        0,
        0
      ),

    base_start -
      c(
        0,
        0.25,
        0,
        0,
        0
      ),

    base_start +
      c(
        0,
        0,
        0.25,
        0,
        0
      ),

    base_start -
      c(
        0,
        0,
        0.25,
        0,
        0
      ),

    base_start +
      c(
        0,
        0,
        0,
        0.25,
        0
      ),

    base_start -
      c(
        0,
        0,
        0,
        0.25,
        0
      ),

    base_start +
      c(
        0,
        0,
        0,
        0,
        0.25
      ),

    base_start -
      c(
        0,
        0,
        0,
        0,
        0.25
      )
  )

  random_seed <- (
    as.double(seed) +
      731041
  ) %% .Machine$integer.max

  set.seed(
    as.integer(random_seed)
  )

  random_perturbations <- matrix(
    stats::rnorm(
      number_of_random_starts * 5L
    ),
    ncol = 5L
  )

  random_perturbations <- sweep(
    random_perturbations,
    2L,
    c(
      0.8,
      0.5,
      0.4,
      0.5,
      0.4
    ),
    `*`
  )

  random_starts <- sweep(
    random_perturbations,
    2L,
    base_start,
    `+`
  )

  starts <- rbind(
    deterministic_starts,
    random_starts
  )

  starts <- t(
    apply(
      starts,
      1L,
      clip_start
    )
  )

  list(
    starts = starts,

    parameter_lower_bounds =
      parameter_lower_bounds,

    parameter_upper_bounds =
      parameter_upper_bounds
  )
}


# ============================================================================
# 8. Candidate storage
# ============================================================================

add_optimizer_candidate <- function(
    candidates,
    theta,
    objective,
    method,
    reported_converged
) {
  theta <- as.numeric(theta)
  objective <- as.numeric(objective)

  if (
    length(theta) == 5L &&
      all(is.finite(theta)) &&
      length(objective) == 1L &&
      is.finite(objective) &&
      objective < 1e99
  ) {
    candidates[[length(candidates) + 1L]] <- list(
      theta = theta,
      objective = objective,
      method = method,
      reported_converged =
        isTRUE(reported_converged)
    )
  }

  candidates
}


deduplicate_candidates <- function(candidates) {
  if (length(candidates) == 0L) {
    return(candidates)
  }

  keys <- vapply(
    candidates,
    function(candidate) {
      paste(
        round(
          candidate$theta,
          digits = 10
        ),
        collapse = ","
      )
    },
    character(1)
  )

  candidates[
    !duplicated(keys)
  ]
}


sort_candidates <- function(candidates) {
  if (length(candidates) == 0L) {
    return(candidates)
  }

  objectives <- vapply(
    candidates,
    function(candidate) {
      candidate$objective
    },
    numeric(1)
  )

  candidates[
    order(objectives)
  ]
}


# ============================================================================
# 9. Multistart optimization
# ============================================================================

generate_optimizer_candidates <- function(
    objective,
    setup
) {
  starts <- setup$starts

  parameter_lower_bounds <-
    setup$parameter_lower_bounds

  parameter_upper_bounds <-
    setup$parameter_upper_bounds

  candidates <- list()

  for (start_index in seq_len(nrow(starts))) {
    current_start <- starts[
      start_index,
      ,
      drop = TRUE
    ]

    starting_objective <- objective(
      current_start
    )

    candidates <- add_optimizer_candidate(
      candidates = candidates,
      theta = current_start,
      objective = starting_objective,
      method = paste0(
        "unoptimized start ",
        start_index
      ),
      reported_converged = FALSE
    )

    nlminb_result <- try(
      stats::nlminb(
        start = current_start,
        objective = objective,

        lower = parameter_lower_bounds,
        upper = parameter_upper_bounds,

        control = list(
          eval.max = 50000,
          iter.max = 30000,
          rel.tol = 1e-12,
          x.tol = 1e-10,
          trace = 0
        )
      ),
      silent = TRUE
    )

    if (!inherits(nlminb_result, "try-error")) {
      candidates <- add_optimizer_candidate(
        candidates = candidates,
        theta = nlminb_result$par,
        objective = nlminb_result$objective,
        method = paste0(
          "nlminb start ",
          start_index
        ),
        reported_converged =
          nlminb_result$convergence == 0
      )
    }

    lbfgsb_result <- try(
      stats::optim(
        par = current_start,
        fn = objective,
        method = "L-BFGS-B",

        lower = parameter_lower_bounds,
        upper = parameter_upper_bounds,

        control = list(
          maxit = 30000,
          factr = 100,
          pgtol = 1e-10,
          trace = 0
        )
      ),
      silent = TRUE
    )

    if (!inherits(lbfgsb_result, "try-error")) {
      candidates <- add_optimizer_candidate(
        candidates = candidates,
        theta = lbfgsb_result$par,
        objective = lbfgsb_result$value,
        method = paste0(
          "L-BFGS-B start ",
          start_index
        ),
        reported_converged =
          lbfgsb_result$convergence == 0
      )
    }
  }

  candidates <- deduplicate_candidates(
    candidates
  )

  candidates <- sort_candidates(
    candidates
  )

  if (length(candidates) == 0L) {
    return(candidates)
  }

  number_to_polish <- min(
    number_of_candidates_to_polish,
    length(candidates)
  )

  initial_best_candidates <- candidates[
    seq_len(number_to_polish)
  ]

  for (candidate_index in seq_along(initial_best_candidates)) {
    candidate <- initial_best_candidates[[candidate_index]]

    bfgs_result <- try(
      stats::optim(
        par = candidate$theta,
        fn = objective,
        method = "BFGS",

        control = list(
          maxit = 50000,
          reltol = 1e-12,
          trace = 0
        )
      ),
      silent = TRUE
    )

    if (!inherits(bfgs_result, "try-error")) {
      candidates <- add_optimizer_candidate(
        candidates = candidates,
        theta = bfgs_result$par,
        objective = bfgs_result$value,
        method = paste0(
          candidate$method,
          " + BFGS"
        ),
        reported_converged =
          bfgs_result$convergence == 0
      )
    }

    final_lbfgsb <- try(
      stats::optim(
        par = candidate$theta,
        fn = objective,
        method = "L-BFGS-B",

        lower =
          setup$parameter_lower_bounds,

        upper =
          setup$parameter_upper_bounds,

        control = list(
          maxit = 50000,
          factr = 10,
          pgtol = 1e-12,
          trace = 0
        )
      ),
      silent = TRUE
    )

    if (!inherits(final_lbfgsb, "try-error")) {
      candidates <- add_optimizer_candidate(
        candidates = candidates,
        theta = final_lbfgsb$par,
        objective = final_lbfgsb$value,
        method = paste0(
          candidate$method,
          " + strict L-BFGS-B"
        ),
        reported_converged =
          final_lbfgsb$convergence == 0
      )
    }
  }

  candidates <- deduplicate_candidates(
    candidates
  )

  sort_candidates(
    candidates
  )
}


# ============================================================================
# 10. Independent gradient and Hessian diagnostics
# ============================================================================

diagnose_candidate <- function(
    candidate,
    objective
) {
  result <- list(
    success = FALSE,

    theta = candidate$theta,
    objective = candidate$objective,

    method = candidate$method,

    reported_converged =
      candidate$reported_converged,

    maximum_absolute_gradient =
      NA_real_,

    minimum_hessian_eigenvalue =
      NA_real_,

    maximum_hessian_eigenvalue =
      NA_real_,

    condition_number =
      NA_real_,

    covariance =
      NULL,

    standard_errors =
      NULL,

    message =
      NA_character_
  )

  gradient <- try(
    numDeriv::grad(
      func = objective,
      x = candidate$theta,
      method = "Richardson"
    ),
    silent = TRUE
  )

  if (
    inherits(gradient, "try-error") ||
      any(!is.finite(gradient))
  ) {
    result$message <-
      "Numerical gradient calculation failed."

    return(result)
  }

  result$maximum_absolute_gradient <-
    max(
      abs(gradient)
    )

  hessian <- try(
    numDeriv::hessian(
      func = objective,
      x = candidate$theta,
      method = "Richardson"
    ),
    silent = TRUE
  )

  if (
    inherits(hessian, "try-error") ||
      any(!is.finite(hessian))
  ) {
    result$message <-
      "Numerical Hessian calculation failed."

    return(result)
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
      any(!is.finite(eigenvalues))
  ) {
    result$message <-
      "Hessian eigenvalue calculation failed."

    return(result)
  }

  result$minimum_hessian_eigenvalue <-
    min(eigenvalues)

  result$maximum_hessian_eigenvalue <-
    max(eigenvalues)

  if (any(eigenvalues <= 0)) {
    result$message <-
      "The numerical Hessian was not positive definite."

    return(result)
  }

  result$condition_number <-
    max(eigenvalues) /
    min(eigenvalues)

  if (
    !is.finite(result$condition_number) ||
      result$condition_number >
        maximum_hessian_condition_number
  ) {
    result$message <-
      "The numerical Hessian was excessively ill-conditioned."

    return(result)
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
    result$message <-
      "The Hessian could not be inverted."

    return(result)
  }

  standard_errors <- sqrt(
    diag(covariance)
  )

  if (
    any(!is.finite(standard_errors)) ||
      any(standard_errors <= 0)
  ) {
    result$message <-
      "The Hessian standard errors were invalid."

    return(result)
  }

  result$covariance <- covariance
  result$standard_errors <- standard_errors

  if (
    result$maximum_absolute_gradient >
      maximum_gradient_tolerance
  ) {
    result$message <- paste0(
      "Maximum absolute gradient exceeded ",
      maximum_gradient_tolerance,
      "."
    )

    return(result)
  }

  result$success <- TRUE
  result$message <- NA_character_

  result
}


# ============================================================================
# 11. Recover one failed simulation
# ============================================================================

recover_one_simulation <- function(
    n,
    seed
) {
  simulated <- simulate_coverage_dataset(
    n = n,
    seed = seed
  )

  y <- simulated$y
  x <- simulated$x

  support_lower <- LL3_support_lower(y)

  objective <- make_LL3_objective(
    y = y,
    x = x,
    support_lower = support_lower
  )

  setup <- construct_optimizer_setup(
    y = y,
    x = x,
    support_lower = support_lower,
    seed = seed
  )

  candidates <- generate_optimizer_candidates(
    objective = objective,
    setup = setup
  )

  if (length(candidates) == 0L) {
    return(
      list(
        success = FALSE,
        message =
          "No finite optimizer candidates were produced."
      )
    )
  }

  number_to_diagnose <- min(
    number_of_candidates_to_diagnose,
    length(candidates)
  )

  diagnostics <- vector(
    "list",
    number_to_diagnose
  )

  for (candidate_index in seq_len(number_to_diagnose)) {
    diagnostics[[candidate_index]] <- diagnose_candidate(
      candidate = candidates[[candidate_index]],
      objective = objective
    )
  }

  acceptable <- which(
    vapply(
      diagnostics,
      function(diagnostic) {
        isTRUE(diagnostic$success)
      },
      logical(1)
    )
  )

  if (length(acceptable) == 0L) {
    best_diagnostic <- diagnostics[[1L]]

    return(
      list(
        success = FALSE,

        message = paste0(
          "No candidate passed the gradient/Hessian checks. ",
          "Best candidate: ",
          best_diagnostic$message
        ),

        best_objective =
          best_diagnostic$objective,

        best_gradient =
          best_diagnostic$
            maximum_absolute_gradient,

        best_minimum_eigenvalue =
          best_diagnostic$
            minimum_hessian_eigenvalue,

        best_condition_number =
          best_diagnostic$
            condition_number,

        candidate_count =
          length(candidates)
      )
    )
  }

  acceptable_objectives <- vapply(
    diagnostics[acceptable],
    function(diagnostic) {
      diagnostic$objective
    },
    numeric(1)
  )

  best_index <- acceptable[
    which.min(
      acceptable_objectives
    )
  ]

  best <- diagnostics[[best_index]]

  list(
    success = TRUE,

    theta = best$theta,
    logLik = -best$objective,

    covariance = best$covariance,
    standard_errors = best$standard_errors,

    method = best$method,

    optimizer_reported_converged =
      best$reported_converged,

    maximum_absolute_gradient =
      best$maximum_absolute_gradient,

    minimum_hessian_eigenvalue =
      best$minimum_hessian_eigenvalue,

    maximum_hessian_eigenvalue =
      best$maximum_hessian_eigenvalue,

    condition_number =
      best$condition_number,

    candidate_count =
      length(candidates),

    diagnosed_candidate_count =
      number_to_diagnose,

    y = y,
    x = x,
    support_lower = support_lower
  )
}


# ============================================================================
# 12. Warm-start GAMLSS at a recovered direct optimum
#
# This is used only when the original GAMLSS fit failed.
# ============================================================================

warm_start_gamlss <- function(
    recovery
) {
  theta <- recovery$theta
  y <- recovery$y
  x <- recovery$x

  support_lower <-
    recovery$support_lower

  mu_start <- rep(
    support_lower -
      exp(theta[1L]),
    length(y)
  )

  sigma_start <- exp(
    theta[2L] +
      theta[3L] * x
  )

  nu_start <- exp(
    theta[4L] +
      theta[5L] * x
  )

  data <- data.frame(
    y = y,
    x = x
  )

  step_values <- c(
    0.005,
    0.010,
    0.020,
    0.040
  )

  valid_fits <- list()
  messages <- character()

  for (step_value in step_values) {
    current_control <- gamlss::gamlss.control(
      c.crit = 1e-7,
      n.cyc = 10000,

      mu.step = step_value,
      sigma.step = step_value,
      nu.step = step_value,

      gd.tol = Inf,
      autostep = TRUE,
      trace = FALSE
    )

    warning_messages <- character()

    current_fit <- withCallingHandlers(
      tryCatch(
        fit_LL3_gamlss(
          mu.formula = y ~ 1,
          sigma.formula = ~ x,
          nu.formula = ~ x,

          data = data,
          lower = support_lower,

          information = "opg",
          allow_nonstationary_mu = FALSE,

          mu.start = mu_start,
          sigma.start = sigma_start,
          nu.start = nu_start,

          control = current_control
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

    if (inherits(current_fit, "error")) {
      messages <- c(
        messages,
        paste0(
          "step ",
          step_value,
          ": ",
          conditionMessage(current_fit)
        )
      )

      next
    }

    parameters <- try(
      extract_LL3_parameters(
        current_fit
      ),
      silent = TRUE
    )

    if (inherits(parameters, "try-error")) {
      messages <- c(
        messages,
        paste0(
          "step ",
          step_value,
          ": parameter extraction failed"
        )
      )

      next
    }

    fitted_mu <- as.numeric(
      parameters$mu
    )

    fitted_sigma <- as.numeric(
      parameters$sigma
    )

    fitted_nu <- as.numeric(
      parameters$nu
    )

    valid <- isTRUE(
      current_fit$converged
    ) &&
      all(is.finite(fitted_mu)) &&
      all(is.finite(fitted_sigma)) &&
      all(is.finite(fitted_nu)) &&
      all(y > fitted_mu) &&
      all(fitted_sigma > 0) &&
      all(fitted_nu > 0)

    if (!valid) {
      messages <- c(
        messages,
        paste0(
          "step ",
          step_value,
          ": failed convergence or validity checks"
        )
      )

      next
    }

    current_logLik <- as.numeric(
      stats::logLik(
        current_fit
      )
    )

    valid_fits[[length(valid_fits) + 1L]] <- list(
      fit = current_fit,
      step = step_value,
      logLik = current_logLik,
      warning = collapse_messages(
        warning_messages
      )
    )
  }

  if (length(valid_fits) == 0L) {
    return(
      list(
        success = FALSE,
        message = collapse_messages(
          messages
        )
      )
    )
  }

  log_likelihoods <- vapply(
    valid_fits,
    function(result) {
      result$logLik
    },
    numeric(1)
  )

  best <- valid_fits[[which.max(log_likelihoods)]]

  difference <-
    recovery$logLik -
    best$logLik

  list(
    success =
      abs(difference) < 1e-3,

    logLik =
      best$logLik,

    direct_minus_gamlss =
      difference,

    step =
      best$step,

    warning =
      best$warning,

    message =
      if (
        abs(difference) < 1e-3
      ) {
        NA_character_
      } else {
        paste0(
          "Warm-start GAMLSS remained ",
          format(
            abs(difference),
            digits = 8
          ),
          " log-likelihood units from the direct optimum."
        )
      }
  )
}


# ============================================================================
# 13. Add recovery columns
# ============================================================================

coverage_results$direct_inference_complete <-
  as_logical_safe(
    coverage_results$optimizer_success
  ) &
  as_logical_safe(
    coverage_results$hessian_success
  )

coverage_results$recovery_attempted <- FALSE
coverage_results$recovery_success <- NA
coverage_results$recovery_method <- NA_character_

coverage_results$recovery_optimizer_reported_converged <-
  NA

coverage_results$recovery_candidate_count <- NA_integer_

coverage_results$recovery_diagnosed_candidate_count <-
  NA_integer_

coverage_results$recovery_gamlss_warmstart_attempted <-
  FALSE

coverage_results$recovery_gamlss_warmstart_success <-
  NA

coverage_results$recovery_gamlss_logLik <- NA_real_

coverage_results$recovery_direct_minus_gamlss <-
  NA_real_

coverage_results$recovery_message <- NA_character_


# ============================================================================
# 14. Recover incomplete rows
# ============================================================================

incomplete_indices <- which(
  !as_logical_safe(
    coverage_results$complete
  )
)

cat(
  "\nIncomplete simulations to recover: ",
  length(incomplete_indices),
  "\n",
  sep = ""
)

recovery_detail_rows <- vector(
  "list",
  length(incomplete_indices)
)

for (recovery_index in seq_along(incomplete_indices)) {
  row_index <- incomplete_indices[
    recovery_index
  ]

  current_n <- as.integer(
    coverage_results$n[
      row_index
    ]
  )

  current_seed <- as.integer(
    coverage_results$seed[
      row_index
    ]
  )

  cat(
    "Recovery ",
    recovery_index,
    "/",
    length(incomplete_indices),
    ": n = ",
    current_n,
    ", seed = ",
    current_seed,
    "\n",
    sep = ""
  )

  coverage_results$recovery_attempted[
    row_index
  ] <- TRUE

  recovery <- recover_one_simulation(
    n = current_n,
    seed = current_seed
  )

  if (!recovery$success) {
    coverage_results$recovery_success[
      row_index
    ] <- FALSE

    coverage_results$recovery_message[
      row_index
    ] <- recovery$message

    recovery_detail_rows[[recovery_index]] <- data.frame(
      n = current_n,
      seed = current_seed,

      recovery_success = FALSE,

      method = NA_character_,

      optimizer_reported_converged = NA,

      candidate_count =
        if (
          !is.null(
            recovery$candidate_count
          )
        ) {
          recovery$candidate_count
        } else {
          NA_integer_
        },

      maximum_absolute_gradient =
        if (
          !is.null(
            recovery$best_gradient
          )
        ) {
          recovery$best_gradient
        } else {
          NA_real_
        },

      minimum_hessian_eigenvalue =
        if (
          !is.null(
            recovery$best_minimum_eigenvalue
          )
        ) {
          recovery$
            best_minimum_eigenvalue
        } else {
          NA_real_
        },

      hessian_condition_number =
        if (
          !is.null(
            recovery$best_condition_number
          )
        ) {
          recovery$
            best_condition_number
        } else {
          NA_real_
        },

      logLik = NA_real_,

      sigma_slope_estimate = NA_real_,
      sigma_slope_standard_error = NA_real_,
      sigma_slope_covered = NA,

      nu_slope_estimate = NA_real_,
      nu_slope_standard_error = NA_real_,
      nu_slope_covered = NA,

      gamlss_warmstart_attempted = FALSE,
      gamlss_warmstart_success = NA,
      gamlss_direct_difference = NA_real_,

      message = recovery$message,

      stringsAsFactors = FALSE
    )

    next
  }

  theta <- recovery$theta
  standard_errors <- recovery$standard_errors

  sigma_slope_estimate <- theta[3L]
  sigma_slope_se <- standard_errors[3L]

  sigma_slope_lower <-
    sigma_slope_estimate -
    critical_value *
    sigma_slope_se

  sigma_slope_upper <-
    sigma_slope_estimate +
    critical_value *
    sigma_slope_se

  sigma_slope_covered <-
    sigma_slope_lower <=
    truth$sigma_slope &&
    sigma_slope_upper >=
    truth$sigma_slope

  nu_slope_estimate <- theta[5L]
  nu_slope_se <- standard_errors[5L]

  nu_slope_lower <-
    nu_slope_estimate -
    critical_value *
    nu_slope_se

  nu_slope_upper <-
    nu_slope_estimate +
    critical_value *
    nu_slope_se

  nu_slope_covered <-
    nu_slope_lower <=
    truth$nu_slope &&
    nu_slope_upper >=
    truth$nu_slope

  coverage_results$optimizer_success[
    row_index
  ] <- TRUE

  coverage_results$hessian_success[
    row_index
  ] <- TRUE

  coverage_results$direct_inference_complete[
    row_index
  ] <- TRUE

  coverage_results$logLik[
    row_index
  ] <- recovery$logLik

  coverage_results$sigma_slope_estimate[
    row_index
  ] <- sigma_slope_estimate

  coverage_results$sigma_slope_standard_error[
    row_index
  ] <- sigma_slope_se

  coverage_results$sigma_slope_lower[
    row_index
  ] <- sigma_slope_lower

  coverage_results$sigma_slope_upper[
    row_index
  ] <- sigma_slope_upper

  coverage_results$sigma_slope_covered[
    row_index
  ] <- sigma_slope_covered

  coverage_results$nu_slope_estimate[
    row_index
  ] <- nu_slope_estimate

  coverage_results$nu_slope_standard_error[
    row_index
  ] <- nu_slope_se

  coverage_results$nu_slope_lower[
    row_index
  ] <- nu_slope_lower

  coverage_results$nu_slope_upper[
    row_index
  ] <- nu_slope_upper

  coverage_results$nu_slope_covered[
    row_index
  ] <- nu_slope_covered

  coverage_results$maximum_absolute_gradient[
    row_index
  ] <- recovery$maximum_absolute_gradient

  coverage_results$hessian_minimum_eigenvalue[
    row_index
  ] <- recovery$minimum_hessian_eigenvalue

  coverage_results$hessian_condition_number[
    row_index
  ] <- recovery$condition_number

  coverage_results$recovery_success[
    row_index
  ] <- TRUE

  coverage_results$recovery_method[
    row_index
  ] <- recovery$method

  coverage_results$recovery_optimizer_reported_converged[
    row_index
  ] <- recovery$optimizer_reported_converged

  coverage_results$recovery_candidate_count[
    row_index
  ] <- recovery$candidate_count

  coverage_results$recovery_diagnosed_candidate_count[
    row_index
  ] <- recovery$diagnosed_candidate_count

  original_fit_success <- isTRUE(
    as_logical_safe(
      coverage_results$fit_success[
        row_index
      ]
    )
  )

  gamlss_warmstart_attempted <- FALSE
  gamlss_warmstart_success <- NA
  gamlss_direct_difference <- NA_real_
  warmstart_message <- NA_character_

  if (!original_fit_success) {
    gamlss_warmstart_attempted <- TRUE

    warmstart <- warm_start_gamlss(
      recovery
    )

    gamlss_warmstart_success <-
      warmstart$success

    if (!is.null(warmstart$logLik)) {
      coverage_results$recovery_gamlss_logLik[
        row_index
      ] <- warmstart$logLik
    }

    if (
      !is.null(
        warmstart$direct_minus_gamlss
      )
    ) {
      gamlss_direct_difference <-
        warmstart$direct_minus_gamlss

      coverage_results$recovery_direct_minus_gamlss[
        row_index
      ] <- gamlss_direct_difference
    }

    if (isTRUE(warmstart$success)) {
      coverage_results$fit_success[
        row_index
      ] <- TRUE

      coverage_results$fit_retried[
        row_index
      ] <- TRUE

      coverage_results$warning[
        row_index
      ] <- collapse_messages(
        c(
          coverage_results$warning[
            row_index
          ],
          warmstart$warning
        )
      )
    } else {
      warmstart_message <-
        warmstart$message
    }
  }

  coverage_results$recovery_gamlss_warmstart_attempted[
    row_index
  ] <- gamlss_warmstart_attempted

  coverage_results$recovery_gamlss_warmstart_success[
    row_index
  ] <- gamlss_warmstart_success

  current_fit_success <- isTRUE(
    as_logical_safe(
      coverage_results$fit_success[
        row_index
      ]
    )
  )

  coverage_results$complete[
    row_index
  ] <- current_fit_success &&
    isTRUE(
      coverage_results$optimizer_success[
        row_index
      ]
    ) &&
    isTRUE(
      coverage_results$hessian_success[
        row_index
      ]
    )

  if (
    isTRUE(
      coverage_results$complete[
        row_index
      ]
    )
  ) {
    coverage_results$error_stage[
      row_index
    ] <- NA_character_

    coverage_results$error_message[
      row_index
    ] <- NA_character_
  } else {
    coverage_results$error_stage[
      row_index
    ] <- "GAMLSS warm-start recovery"

    coverage_results$error_message[
      row_index
    ] <- warmstart_message
  }

  coverage_results$recovery_message[
    row_index
  ] <- warmstart_message

  recovery_detail_rows[[recovery_index]] <- data.frame(
    n = current_n,
    seed = current_seed,

    recovery_success = TRUE,

    method = recovery$method,

    optimizer_reported_converged =
      recovery$
        optimizer_reported_converged,

    candidate_count =
      recovery$candidate_count,

    maximum_absolute_gradient =
      recovery$
        maximum_absolute_gradient,

    minimum_hessian_eigenvalue =
      recovery$
        minimum_hessian_eigenvalue,

    hessian_condition_number =
      recovery$condition_number,

    logLik =
      recovery$logLik,

    sigma_slope_estimate =
      sigma_slope_estimate,

    sigma_slope_standard_error =
      sigma_slope_se,

    sigma_slope_covered =
      sigma_slope_covered,

    nu_slope_estimate =
      nu_slope_estimate,

    nu_slope_standard_error =
      nu_slope_se,

    nu_slope_covered =
      nu_slope_covered,

    gamlss_warmstart_attempted =
      gamlss_warmstart_attempted,

    gamlss_warmstart_success =
      gamlss_warmstart_success,

    gamlss_direct_difference =
      gamlss_direct_difference,

    message =
      collapse_messages(
        c(
          warmstart_message
        )
      ),

    stringsAsFactors = FALSE
  )
}


# ============================================================================
# 15. Combine recovery details
# ============================================================================

recovery_details <- do.call(
  rbind,
  recovery_detail_rows
)

rownames(
  recovery_details
) <- NULL


# ============================================================================
# 16. Recalculate coverage summaries
# ============================================================================

sample_sizes <- sort(
  unique(
    coverage_results$n
  )
)

summary_rows <- lapply(
  sample_sizes,
  function(current_n) {
    all_runs <- coverage_results[
      coverage_results$n ==
        current_n,
      ,
      drop = FALSE
    ]

    pipeline_runs <- all_runs[
      as_logical_safe(
        all_runs$complete
      ),
      ,
      drop = FALSE
    ]

    direct_runs <- all_runs[
      as_logical_safe(
        all_runs$
          direct_inference_complete
      ),
      ,
      drop = FALSE
    ]

    pipeline_sigma_coverage <- safe_rate(
      pipeline_runs$
        sigma_slope_covered
    )

    pipeline_nu_coverage <- safe_rate(
      pipeline_runs$
        nu_slope_covered
    )

    direct_sigma_coverage <- safe_rate(
      direct_runs$
        sigma_slope_covered
    )

    direct_nu_coverage <- safe_rate(
      direct_runs$
        nu_slope_covered
    )

    pipeline_count <- nrow(
      pipeline_runs
    )

    direct_count <- nrow(
      direct_runs
    )

    data.frame(
      n = current_n,

      attempted =
        nrow(all_runs),

      pipeline_complete =
        pipeline_count,

      pipeline_completion_rate =
        pipeline_count /
        nrow(all_runs),

      direct_inference_complete =
        direct_count,

      direct_inference_completion_rate =
        direct_count /
        nrow(all_runs),

      fit_success_rate =
        safe_rate(
          all_runs$fit_success
        ),

      optimizer_success_rate =
        safe_rate(
          all_runs$optimizer_success
        ),

      hessian_success_rate =
        safe_rate(
          all_runs$hessian_success
        ),

      recovery_attempt_rate =
        safe_rate(
          all_runs$
            recovery_attempted
        ),

      recovery_success_rate_among_attempted =
        safe_rate(
          all_runs$
            recovery_success[
              as_logical_safe(
                all_runs$
                  recovery_attempted
              )
            ]
        ),

      sigma_slope_bias =
        safe_mean(
          pipeline_runs$
            sigma_slope_estimate -
            truth$sigma_slope
        ),

      sigma_slope_rmse =
        sqrt(
          safe_mean(
            (
              pipeline_runs$
                sigma_slope_estimate -
                truth$sigma_slope
            )^2
          )
        ),

      sigma_slope_empirical_sd =
        safe_sd(
          pipeline_runs$
            sigma_slope_estimate
        ),

      sigma_slope_mean_se =
        safe_mean(
          pipeline_runs$
            sigma_slope_standard_error
        ),

      sigma_slope_coverage =
        pipeline_sigma_coverage,

      sigma_coverage_MCSE =
        if (
          is.finite(
            pipeline_sigma_coverage
          ) &&
            pipeline_count > 0L
        ) {
          sqrt(
            pipeline_sigma_coverage *
              (
                1 -
                  pipeline_sigma_coverage
              ) /
              pipeline_count
          )
        } else {
          NA_real_
        },

      direct_sigma_slope_coverage =
        direct_sigma_coverage,

      direct_sigma_coverage_MCSE =
        if (
          is.finite(
            direct_sigma_coverage
          ) &&
            direct_count > 0L
        ) {
          sqrt(
            direct_sigma_coverage *
              (
                1 -
                  direct_sigma_coverage
              ) /
              direct_count
          )
        } else {
          NA_real_
        },

      nu_slope_bias =
        safe_mean(
          pipeline_runs$
            nu_slope_estimate -
            truth$nu_slope
        ),

      nu_slope_rmse =
        sqrt(
          safe_mean(
            (
              pipeline_runs$
                nu_slope_estimate -
                truth$nu_slope
            )^2
          )
        ),

      nu_slope_empirical_sd =
        safe_sd(
          pipeline_runs$
            nu_slope_estimate
        ),

      nu_slope_mean_se =
        safe_mean(
          pipeline_runs$
            nu_slope_standard_error
        ),

      nu_slope_coverage =
        pipeline_nu_coverage,

      nu_coverage_MCSE =
        if (
          is.finite(
            pipeline_nu_coverage
          ) &&
            pipeline_count > 0L
        ) {
          sqrt(
            pipeline_nu_coverage *
              (
                1 -
                  pipeline_nu_coverage
              ) /
              pipeline_count
          )
        } else {
          NA_real_
        },

      direct_nu_slope_coverage =
        direct_nu_coverage,

      direct_nu_coverage_MCSE =
        if (
          is.finite(
            direct_nu_coverage
          ) &&
            direct_count > 0L
        ) {
          sqrt(
            direct_nu_coverage *
              (
                1 -
                  direct_nu_coverage
              ) /
              direct_count
          )
        } else {
          NA_real_
        },

      median_maximum_gradient =
        safe_median(
          pipeline_runs$
            maximum_absolute_gradient
        ),

      maximum_maximum_gradient =
        safe_max(
          pipeline_runs$
            maximum_absolute_gradient
        ),

      median_hessian_condition_number =
        safe_median(
          pipeline_runs$
            hessian_condition_number
        ),

      maximum_hessian_condition_number =
        safe_max(
          pipeline_runs$
            hessian_condition_number
        ),

      stringsAsFactors = FALSE
    )
  }
)

recovered_summary <- do.call(
  rbind,
  summary_rows
)

rownames(
  recovered_summary
) <- NULL


# ============================================================================
# 17. Remaining failures
# ============================================================================

remaining_failures <- coverage_results[
  !as_logical_safe(
    coverage_results$complete
  ),
  ,
  drop = FALSE
]

direct_inference_failures <- coverage_results[
  !as_logical_safe(
    coverage_results$
      direct_inference_complete
  ),
  ,
  drop = FALSE
]


# ============================================================================
# 18. Save outputs
# ============================================================================

recovered_raw_file <- file.path(
  results_dir,
  "interval_coverage_recovered_raw.csv"
)

recovered_summary_file <- file.path(
  results_dir,
  "interval_coverage_recovered_summary.csv"
)

recovery_details_file <- file.path(
  results_dir,
  "interval_coverage_recovery_details.csv"
)

remaining_failures_file <- file.path(
  results_dir,
  "interval_coverage_recovered_failures.csv"
)

direct_failures_file <- file.path(
  results_dir,
  "interval_coverage_direct_inference_failures.csv"
)

utils::write.csv(
  coverage_results,
  recovered_raw_file,
  row.names = FALSE
)

utils::write.csv(
  recovered_summary,
  recovered_summary_file,
  row.names = FALSE
)

utils::write.csv(
  recovery_details,
  recovery_details_file,
  row.names = FALSE
)

utils::write.csv(
  remaining_failures,
  remaining_failures_file,
  row.names = FALSE
)

utils::write.csv(
  direct_inference_failures,
  direct_failures_file,
  row.names = FALSE
)


# ============================================================================
# 19. Print results
# ============================================================================

cat(
  "\n============================================================\n",
  "Targeted interval-coverage recovery summary\n",
  "============================================================\n",
  sep = ""
)

print(
  recovered_summary,
  digits = 12,
  row.names = FALSE
)

cat(
  "\n============================================================\n",
  "Recovery details\n",
  "============================================================\n",
  sep = ""
)

print(
  recovery_details,
  digits = 12,
  row.names = FALSE
)

cat(
  "\nRemaining GAMLSS-pipeline failures: ",
  nrow(remaining_failures),
  "\n",
  sep = ""
)

cat(
  "Remaining direct-inference failures: ",
  nrow(direct_inference_failures),
  "\n",
  sep = ""
)

if (nrow(remaining_failures) > 0L) {
  print(
    remaining_failures[
      ,
      intersect(
        c(
          "n",
          "seed",
          "fit_success",
          "optimizer_success",
          "hessian_success",
          "recovery_success",
          "recovery_message",
          "error_stage",
          "error_message"
        ),
        names(remaining_failures)
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}

cat(
  "\nFiles written:\n",
  recovered_raw_file,
  "\n",
  recovered_summary_file,
  "\n",
  recovery_details_file,
  "\n",
  remaining_failures_file,
  "\n",
  direct_failures_file,
  "\n",
  sep = ""
)
