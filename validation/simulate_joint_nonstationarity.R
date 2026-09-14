# validation/simulate_joint_nonstationarity.R
#
# Monte Carlo validation for the joint nonstationary LL3 model:
#
#   mu_t = constant
#   log(sigma_t) = beta_sigma_0 + beta_sigma_1 * time_t
#   log(nu_t)    = beta_nu_0    + beta_nu_1    * time_t
#
# Candidate models:
#
#   stationary: sigma ~ 1,    nu ~ 1
#   sigma:      sigma ~ time, nu ~ 1
#   nu:         sigma ~ 1,    nu ~ time
#   joint:      sigma ~ time, nu ~ time
#
# Important distinction:
#
#   joint_fit_success:
#     The joint model converged and satisfies distribution-support,
#     positivity, and numerical-boundary checks.
#
#   comparison_complete:
#     The stationary, sigma-only, nu-only, and joint models all succeeded,
#     allowing a complete AIC/BIC comparison.
#
# A failure of an auxiliary comparison model does not invalidate a
# successful joint-model fit.


# =========================================================================
# 1. Locate the project directory
# =========================================================================

locate_LL3_project_root <- function() {
  source_file <- tryCatch(
    normalizePath(
      sys.frame(1)$ofile,
      winslash = "/",
      mustWork = TRUE
    ),
    error = function(e) {
      NA_character_
    }
  )

  if (!is.na(source_file)) {
    candidate <- dirname(
      dirname(source_file)
    )

    if (
      file.exists(
        file.path(
          candidate,
          "R",
          "LL3_distribution.R"
        )
      )
    ) {
      return(candidate)
    }
  }

  if (
    file.exists(
      file.path(
        "R",
        "LL3_distribution.R"
      )
    )
  ) {
    return(
      normalizePath(
        ".",
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  candidate <- file.path(
    getwd(),
    "LL3_GAMLSS_submission"
  )

  if (
    file.exists(
      file.path(
        candidate,
        "R",
        "LL3_distribution.R"
      )
    )
  ) {
    return(
      normalizePath(
        candidate,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  stop(
    "Could not locate the LL3_GAMLSS_submission project directory."
  )
}


LL3_PROJECT_ROOT <- locate_LL3_project_root()


# =========================================================================
# 2. Load project code
# =========================================================================

source(
  file.path(
    LL3_PROJECT_ROOT,
    "R",
    "LL3_distribution.R"
  )
)

source(
  file.path(
    LL3_PROJECT_ROOT,
    "R",
    "LL3_gamlss_family.R"
  )
)

source(
  file.path(
    LL3_PROJECT_ROOT,
    "R",
    "LL3_fitting.R"
  )
)

source(
  file.path(
    LL3_PROJECT_ROOT,
    "R",
    "LL3_diagnostics.R"
  )
)


if (!requireNamespace("gamlss", quietly = TRUE)) {
  stop(
    "Install the gamlss package before running this simulation."
  )
}


# =========================================================================
# 3. Simulation settings
# =========================================================================

repetitions <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "LL3_REPS",
      "200"
    )
  )
)

if (
  length(repetitions) != 1L ||
    is.na(repetitions) ||
    repetitions < 1L
) {
  stop(
    "LL3_REPS must be a positive integer."
  )
}


direct_checks_per_sample_size <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "LL3_DIRECT_CHECKS_PER_N",
      "10"
    )
  )
)

if (
  length(direct_checks_per_sample_size) != 1L ||
    is.na(direct_checks_per_sample_size) ||
    direct_checks_per_sample_size < 0L
) {
  stop(
    "LL3_DIRECT_CHECKS_PER_N must be a nonnegative integer."
  )
}


direct_optimizer_starts <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "LL3_DIRECT_STARTS",
      "10"
    )
  )
)

if (
  length(direct_optimizer_starts) != 1L ||
    is.na(direct_optimizer_starts) ||
    direct_optimizer_starts < 1L
) {
  stop(
    "LL3_DIRECT_STARTS must be a positive integer."
  )
}


sample_sizes <- c(
  60L,
  100L,
  250L,
  1000L
)


truth <- list(
  mu = -50,

  sigma_intercept = log(70),
  sigma_slope = 0.55,

  nu_intercept = log(1.8),
  nu_slope = -0.20
)


master_seed <- 2026L

set.seed(master_seed)

simulation_seeds <- sample.int(
  .Machine$integer.max,
  size = repetitions * length(sample_sizes),
  replace = FALSE
)


# =========================================================================
# 4. Utility functions
# =========================================================================

extract_try_error_message <- function(x) {
  if (!inherits(x, "try-error")) {
    return(NA_character_)
  }

  condition <- attr(
    x,
    "condition"
  )

  if (!is.null(condition)) {
    return(
      conditionMessage(condition)
    )
  }

  as.character(x)
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


safe_max <- function(x) {
  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  max(x)
}


safe_rmse <- function(
    estimate,
    target
) {
  keep <- is.finite(estimate)

  estimate <- estimate[
    keep
  ]

  if (length(estimate) == 0L) {
    return(NA_real_)
  }

  sqrt(
    mean(
      (estimate - target)^2
    )
  )
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


safe_correlation <- function(
    x,
    y
) {
  keep <- is.finite(x) & is.finite(y)

  x <- x[
    keep
  ]

  y <- y[
    keep
  ]

  if (
    length(x) < 2L ||
      stats::sd(x) == 0 ||
      stats::sd(y) == 0
  ) {
    return(NA_real_)
  }

  stats::cor(
    x,
    y
  )
}


get_model_convergence <- function(model) {
  result <- tryCatch(
    LL3_fit_converged(model),
    error = function(e) {
      NA
    }
  )

  isTRUE(result)
}


# =========================================================================
# 5. Local model validation
#
# This validation does not depend on the names returned by check_LL3_fit().
#
# The artificial lower boundary is treated separately from the actual
# distribution support y > mu.
#
# A fitted mu equal to lower within floating-point tolerance is recorded as
# boundary contact, not as a material support violation.
# =========================================================================

validate_LL3_model <- function(
    model,
    y,
    lower,
    boundary_tolerance_multiplier = 100
) {
  output <- list(
    success = FALSE,
    parameters = NULL,
    error_message = NA_character_,

    checks = list(
      converged = FALSE,
      finite_parameters = FALSE,
      support_ok = FALSE,
      sigma_positive = FALSE,
      nu_positive = FALSE,
      lower_below_data = FALSE,

      mu_below_lower_exact = FALSE,
      mu_not_materially_above_lower = FALSE,
      mu_boundary_contact = FALSE,
      material_boundary_violation = TRUE,

      minimum_y_minus_mu = NA_real_,
      minimum_lower_minus_mu = NA_real_,
      minimum_y_minus_lower = NA_real_,
      boundary_tolerance = NA_real_
    )
  )

  if (
    inherits(model, "try-error") ||
      inherits(model, "error")
  ) {
    output$error_message <-
      if (inherits(model, "try-error")) {
        extract_try_error_message(model)
      } else {
        conditionMessage(model)
      }

    return(output)
  }


  parameters <- try(
    extract_LL3_parameters(model),
    silent = TRUE
  )

  if (inherits(parameters, "try-error")) {
    output$error_message <-
      extract_try_error_message(parameters)

    return(output)
  }


  required_parameter_names <- c(
    "mu",
    "sigma",
    "nu"
  )

  if (
    !all(
      required_parameter_names %in%
        names(parameters)
    )
  ) {
    output$error_message <-
      "The fitted parameter object does not contain mu, sigma, and nu."

    return(output)
  }


  y <- as.numeric(y)

  mu <- as.numeric(
    parameters$mu
  )

  sigma <- as.numeric(
    parameters$sigma
  )

  nu <- as.numeric(
    parameters$nu
  )


  if (
    length(mu) != length(y) ||
      length(sigma) != length(y) ||
      length(nu) != length(y)
  ) {
    output$error_message <-
      "Fitted parameter vectors do not match the response length."

    return(output)
  }


  finite_parameters <- all(
    is.finite(mu) &
      is.finite(sigma) &
      is.finite(nu)
  )


  converged <- get_model_convergence(
    model
  )


  support_gap <- y - mu

  support_ok <-
    finite_parameters &&
    all(
      support_gap > 0
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


  lower_below_data <-
    is.finite(lower) &&
    lower < min(y)


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


  mu_below_lower_exact <-
    finite_parameters &&
    all(
      lower_gap > 0
    )


  material_boundary_violation <-
    !finite_parameters ||
    any(
      mu >
        lower + boundary_tolerance
    )


  mu_not_materially_above_lower <-
    !material_boundary_violation


  mu_boundary_contact <-
    finite_parameters &&
    any(
      lower_gap <= boundary_tolerance
    )


  output$parameters <- parameters

  output$checks <- list(
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

    mu_below_lower_exact =
      mu_below_lower_exact,

    mu_not_materially_above_lower =
      mu_not_materially_above_lower,

    mu_boundary_contact =
      mu_boundary_contact,

    material_boundary_violation =
      material_boundary_violation,

    minimum_y_minus_mu =
      if (all(is.finite(support_gap))) {
        min(support_gap)
      } else {
        NA_real_
      },

    minimum_lower_minus_mu =
      if (all(is.finite(lower_gap))) {
        min(lower_gap)
      } else {
        NA_real_
      },

    minimum_y_minus_lower =
      if (
        is.finite(lower) &&
          all(is.finite(y))
      ) {
        min(y - lower)
      } else {
        NA_real_
      },

    boundary_tolerance =
      boundary_tolerance
  )


  output$success <- all(
    c(
      converged,
      finite_parameters,
      support_ok,
      sigma_positive,
      nu_positive,
      lower_below_data,
      mu_not_materially_above_lower
    )
  )


  if (!output$success) {
    failed_checks <- names(
      output$checks
    )[
      vapply(
        output$checks[
          c(
            "converged",
            "finite_parameters",
            "support_ok",
            "sigma_positive",
            "nu_positive",
            "lower_below_data",
            "mu_not_materially_above_lower"
          )
        ],
        function(value) {
          !isTRUE(value)
        },
        logical(1)
      )
    ]

    output$error_message <- paste0(
      "Failed checks: ",
      paste(
        failed_checks,
        collapse = ", "
      ),
      "."
    )
  }


  output
}


# =========================================================================
# 6. Direct-optimizer likelihood comparison
# =========================================================================

run_direct_optimizer_check <- function(
    joint_model,
    data,
    lower,
    seed,
    n_starts
) {
  output <- list(
    attempted = TRUE,
    success = FALSE,

    error_message = NA_character_,

    gamlss_logLik = NA_real_,
    direct_logLik = NA_real_,
    absolute_difference = NA_real_,
    agreement_1e3 = NA
  )


  comparison <- try(
    compare_LL3_gamlss_direct(
      gamlss_fit = joint_model,
      data = data,
      mu.formula = y ~ 1,
      sigma.formula = ~ time,
      nu.formula = ~ time,
      n_starts = n_starts,
      seed = seed
    ),
    silent = TRUE
  )


  if (inherits(comparison, "try-error")) {
    output$error_message <-
      extract_try_error_message(
        comparison
      )

    return(output)
  }


  direct_fit <-
    if (!is.null(comparison$direct)) {
      comparison$direct
    } else if (!is.null(comparison$direct_fit)) {
      comparison$direct_fit
    } else {
      NULL
    }


  if (
    is.null(direct_fit) ||
      !is.list(direct_fit)
  ) {
    output$error_message <-
      "The direct comparison did not return a direct-fit object."

    return(output)
  }


  joint_parameters <- try(
    extract_LL3_parameters(
      joint_model
    ),
    silent = TRUE
  )


  if (inherits(joint_parameters, "try-error")) {
    output$error_message <-
      "Could not extract parameters from the joint GAMLSS model."

    return(output)
  }


  gamlss_logLik <- sum(
    dLL3(
      data$y,
      mu = joint_parameters$mu,
      sigma = joint_parameters$sigma,
      nu = joint_parameters$nu,
      log = TRUE
    )
  )


  direct_logLik <-
    if (
      !is.null(direct_fit$logLik) &&
        length(direct_fit$logLik) == 1L
    ) {
      as.numeric(
        direct_fit$logLik
      )
    } else {
      NA_real_
    }


  direct_convergence <-
    if (
      !is.null(direct_fit$convergence) &&
        length(direct_fit$convergence) == 1L
    ) {
      as.integer(
        direct_fit$convergence
      )
    } else if (
      !is.null(direct_fit$optim$convergence) &&
        length(direct_fit$optim$convergence) == 1L
    ) {
      as.integer(
        direct_fit$optim$convergence
      )
    } else {
      NA_integer_
    }


  direct_mu <-
    if (
      !is.null(direct_fit$mu) &&
        length(direct_fit$mu) == 1L
    ) {
      as.numeric(
        direct_fit$mu
      )
    } else {
      NA_real_
    }


  direct_valid <-
    identical(
      direct_convergence,
      0L
    ) &&
    is.finite(direct_logLik) &&
    direct_logLik > -1e50 &&
    is.finite(direct_mu) &&
    direct_mu <=
      lower +
      100 *
      .Machine$double.eps *
      max(
        1,
        abs(lower),
        abs(direct_mu)
      )


  if (!direct_valid) {
    output$error_message <- paste(
      "The direct optimizer returned an invalid solution.",
      paste0(
        "convergence=",
        direct_convergence
      ),
      paste0(
        "logLik=",
        format(
          direct_logLik,
          digits = 12
        )
      ),
      paste0(
        "mu=",
        format(
          direct_mu,
          digits = 12
        )
      )
    )

    return(output)
  }


  if (!is.finite(gamlss_logLik)) {
    output$error_message <-
      "The GAMLSS log-likelihood was not finite."

    return(output)
  }


  absolute_difference <- abs(
    gamlss_logLik -
      direct_logLik
  )


  if (!is.finite(absolute_difference)) {
    output$error_message <-
      "The likelihood difference was not finite."

    return(output)
  }


  output$success <- TRUE
  output$gamlss_logLik <- gamlss_logLik
  output$direct_logLik <- direct_logLik
  output$absolute_difference <- absolute_difference
  output$agreement_1e3 <- absolute_difference < 1e-3

  output
}


# =========================================================================
# 7. Empty result row
# =========================================================================

new_joint_result_row <- function(
    n,
    seed
) {
  data.frame(
    n = as.integer(n),
    seed = as.integer(seed),

    stationary_fit_success = FALSE,
    sigma_fit_success = FALSE,
    nu_fit_success = FALSE,

    joint_converged = FALSE,
    joint_fit_success = FALSE,
    comparison_complete = FALSE,

    finite_parameters = NA,
    support_ok = NA,
    sigma_positive = NA,
    nu_positive = NA,
    lower_below_data = NA,

    mu_below_lower = NA,
    mu_below_lower_exact = NA,
    mu_not_materially_above_lower = NA,
    mu_boundary_contact = NA,
    material_boundary_violation = NA,

    minimum_y_minus_mu = NA_real_,
    minimum_lower_minus_mu = NA_real_,
    minimum_y_minus_lower = NA_real_,
    boundary_tolerance = NA_real_,

    error_stage = NA_character_,
    error_message = NA_character_,

    mu_hat = NA_real_,

    sigma_intercept_hat = NA_real_,
    sigma_slope_hat = NA_real_,
    sigma_curve_bias = NA_real_,
    sigma_curve_rmse = NA_real_,
    sigma_curve_correlation = NA_real_,

    nu_intercept_hat = NA_real_,
    nu_slope_hat = NA_real_,
    nu_curve_bias = NA_real_,
    nu_curve_rmse = NA_real_,
    nu_curve_correlation = NA_real_,

    AIC_selected_model = NA_character_,
    BIC_selected_model = NA_character_,
    AIC_selects_joint = NA,
    BIC_selects_joint = NA,

    residual_mean = NA_real_,
    residual_sd = NA_real_,
    residual_tail_rate_95 = NA_real_,

    direct_optimizer_attempted = FALSE,
    direct_optimizer_success = FALSE,
    direct_optimizer_error = NA_character_,
    direct_gamlss_logLik = NA_real_,
    direct_optimizer_logLik = NA_real_,
    direct_logLik_difference = NA_real_,
    direct_agreement_1e3 = NA,

    stringsAsFactors = FALSE
  )
}


# =========================================================================
# 8. Run simulations
# =========================================================================

number_of_rows <-
  repetitions *
  length(sample_sizes)

rows <- vector(
  "list",
  number_of_rows
)

row_index <- 0L


for (current_n in sample_sizes) {
  direct_checks_used <- 0L

  for (replication in seq_len(repetitions)) {
    row_index <- row_index + 1L

    current_seed <-
      simulation_seeds[
        row_index
      ]


    result <- new_joint_result_row(
      n = current_n,
      seed = current_seed
    )


    set.seed(
      current_seed
    )


    time <- seq(
      -1,
      1,
      length.out = current_n
    )


    sigma_true <- exp(
      truth$sigma_intercept +
        truth$sigma_slope *
        time
    )


    nu_true <- exp(
      truth$nu_intercept +
        truth$nu_slope *
        time
    )


    y <- rLL3(
      n = current_n,
      mu = truth$mu,
      sigma = sigma_true,
      nu = nu_true
    )


    simulation_data <- data.frame(
      y = y,
      time = time,
      sigma_true = sigma_true,
      nu_true = nu_true
    )


    lower <- LL3_support_lower(
      y
    )


    control <- LL3_default_control(
      n_cycles = 700,
      mu_step = 0.03,
      sigma_step = 0.04,
      nu_step = 0.04,
      trace = FALSE
    )


    # ---------------------------------------------------------------------
    # Stationary model
    # ---------------------------------------------------------------------

    stationary_model <- try(
      fit_LL3_gamlss(
        mu.formula = y ~ 1,
        sigma.formula = ~ 1,
        nu.formula = ~ 1,

        data = simulation_data,
        lower = lower,

        information = "opg",
        allow_nonstationary_mu = FALSE,

        control = control
      ),
      silent = TRUE
    )


    stationary_validation <- validate_LL3_model(
      model = stationary_model,
      y = y,
      lower = lower
    )


    result$stationary_fit_success <-
      stationary_validation$success


    if (!stationary_validation$success) {
      result$error_stage <-
        "stationary model"

      result$error_message <-
        stationary_validation$error_message

      rows[[row_index]] <- result

      next
    }


    stationary_parameters <-
      stationary_validation$parameters


    # ---------------------------------------------------------------------
    # Sigma-only model
    # ---------------------------------------------------------------------

    sigma_model <- try(
      fit_LL3_gamlss(
        mu.formula = y ~ 1,
        sigma.formula = ~ time,
        nu.formula = ~ 1,

        data = simulation_data,
        lower = lower,

        information = "opg",
        allow_nonstationary_mu = FALSE,

        mu.start =
          stationary_parameters$mu,

        sigma.start =
          stationary_parameters$sigma,

        nu.start =
          stationary_parameters$nu,

        control = control
      ),
      silent = TRUE
    )


    sigma_validation <- validate_LL3_model(
      model = sigma_model,
      y = y,
      lower = lower
    )


    result$sigma_fit_success <-
      sigma_validation$success


    # ---------------------------------------------------------------------
    # Nu-only model
    # ---------------------------------------------------------------------

    nu_model <- try(
      fit_LL3_gamlss(
        mu.formula = y ~ 1,
        sigma.formula = ~ 1,
        nu.formula = ~ time,

        data = simulation_data,
        lower = lower,

        information = "opg",
        allow_nonstationary_mu = FALSE,

        mu.start =
          stationary_parameters$mu,

        sigma.start =
          stationary_parameters$sigma,

        nu.start =
          stationary_parameters$nu,

        control = control
      ),
      silent = TRUE
    )


    nu_validation <- validate_LL3_model(
      model = nu_model,
      y = y,
      lower = lower
    )


    result$nu_fit_success <-
      nu_validation$success


    # ---------------------------------------------------------------------
    # Starting curves for joint model
    # ---------------------------------------------------------------------

    sigma_start <-
      stationary_parameters$sigma


    if (sigma_validation$success) {
      sigma_start <-
        sigma_validation$
          parameters$
          sigma
    }


    nu_start <-
      stationary_parameters$nu


    if (nu_validation$success) {
      nu_start <-
        nu_validation$
          parameters$
          nu
    }


    # ---------------------------------------------------------------------
    # Joint sigma-nu model
    # ---------------------------------------------------------------------

    joint_model <- try(
      fit_LL3_gamlss(
        mu.formula = y ~ 1,
        sigma.formula = ~ time,
        nu.formula = ~ time,

        data = simulation_data,
        lower = lower,

        information = "opg",
        allow_nonstationary_mu = FALSE,

        mu.start =
          stationary_parameters$mu,

        sigma.start =
          sigma_start,

        nu.start =
          nu_start,

        control = control
      ),
      silent = TRUE
    )


    joint_validation <- validate_LL3_model(
      model = joint_model,
      y = y,
      lower = lower
    )


    result$joint_converged <-
      isTRUE(
        joint_validation$
          checks$
          converged
      )


    result$joint_fit_success <-
      joint_validation$success


    result$finite_parameters <-
      joint_validation$
        checks$
        finite_parameters


    result$support_ok <-
      joint_validation$
        checks$
        support_ok


    result$sigma_positive <-
      joint_validation$
        checks$
        sigma_positive


    result$nu_positive <-
      joint_validation$
        checks$
        nu_positive


    result$lower_below_data <-
      joint_validation$
        checks$
        lower_below_data


    result$mu_below_lower <-
      joint_validation$
        checks$
        mu_below_lower_exact


    result$mu_below_lower_exact <-
      joint_validation$
        checks$
        mu_below_lower_exact


    result$mu_not_materially_above_lower <-
      joint_validation$
        checks$
        mu_not_materially_above_lower


    result$mu_boundary_contact <-
      joint_validation$
        checks$
        mu_boundary_contact


    result$material_boundary_violation <-
      joint_validation$
        checks$
        material_boundary_violation


    result$minimum_y_minus_mu <-
      joint_validation$
        checks$
        minimum_y_minus_mu


    result$minimum_lower_minus_mu <-
      joint_validation$
        checks$
        minimum_lower_minus_mu


    result$minimum_y_minus_lower <-
      joint_validation$
        checks$
        minimum_y_minus_lower


    result$boundary_tolerance <-
      joint_validation$
        checks$
        boundary_tolerance


    result$comparison_complete <-
      isTRUE(
        result$
          stationary_fit_success
      ) &&
      isTRUE(
        result$
          sigma_fit_success
      ) &&
      isTRUE(
        result$
          nu_fit_success
      ) &&
      isTRUE(
        result$
          joint_fit_success
      )


    if (!result$joint_fit_success) {
      result$error_stage <-
        "joint model"

      result$error_message <-
        joint_validation$error_message

      rows[[row_index]] <- result

      next
    }


    joint_parameters <-
      joint_validation$parameters


    # ---------------------------------------------------------------------
    # Extract parameter estimates
    # ---------------------------------------------------------------------

    sigma_coefficients <- try(
      stats::coef(
        joint_model,
        what = "sigma"
      ),
      silent = TRUE
    )


    nu_coefficients <- try(
      stats::coef(
        joint_model,
        what = "nu"
      ),
      silent = TRUE
    )


    if (
      inherits(
        sigma_coefficients,
        "try-error"
      ) ||
      inherits(
        nu_coefficients,
        "try-error"
      )
    ) {
      result$joint_fit_success <- FALSE
      result$comparison_complete <- FALSE

      result$error_stage <-
        "joint parameter extraction"

      result$error_message <-
        "Could not extract sigma or nu coefficients from the joint model."

      rows[[row_index]] <- result

      next
    }


    required_sigma_terms <- c(
      "(Intercept)",
      "time"
    )


    required_nu_terms <- c(
      "(Intercept)",
      "time"
    )


    if (
      !all(
        required_sigma_terms %in%
          names(sigma_coefficients)
      ) ||
      !all(
        required_nu_terms %in%
          names(nu_coefficients)
      )
    ) {
      result$joint_fit_success <- FALSE
      result$comparison_complete <- FALSE

      result$error_stage <-
        "joint coefficient names"

      result$error_message <-
        "Expected intercept and time coefficients were not found."

      rows[[row_index]] <- result

      next
    }


    result$mu_hat <-
      mean(
        joint_parameters$mu
      )


    result$sigma_intercept_hat <-
      unname(
        sigma_coefficients[
          "(Intercept)"
        ]
      )


    result$sigma_slope_hat <-
      unname(
        sigma_coefficients[
          "time"
        ]
      )


    sigma_error <-
      joint_parameters$sigma -
      sigma_true


    result$sigma_curve_bias <-
      mean(
        sigma_error
      )


    result$sigma_curve_rmse <-
      sqrt(
        mean(
          sigma_error^2
        )
      )


    result$sigma_curve_correlation <-
      safe_correlation(
        joint_parameters$sigma,
        sigma_true
      )


    result$nu_intercept_hat <-
      unname(
        nu_coefficients[
          "(Intercept)"
        ]
      )


    result$nu_slope_hat <-
      unname(
        nu_coefficients[
          "time"
        ]
      )


    nu_error <-
      joint_parameters$nu -
      nu_true


    result$nu_curve_bias <-
      mean(
        nu_error
      )


    result$nu_curve_rmse <-
      sqrt(
        mean(
          nu_error^2
        )
      )


    result$nu_curve_correlation <-
      safe_correlation(
        joint_parameters$nu,
        nu_true
      )


    # ---------------------------------------------------------------------
    # Residual calibration
    # ---------------------------------------------------------------------

    normalized_residuals <- LL3_index(
      y = y,
      mu = joint_parameters$mu,
      sigma = joint_parameters$sigma,
      nu = joint_parameters$nu
    )


    result$residual_mean <-
      mean(
        normalized_residuals
      )


    result$residual_sd <-
      stats::sd(
        normalized_residuals
      )


    result$residual_tail_rate_95 <-
      mean(
        abs(
          normalized_residuals
        ) >
          stats::qnorm(
            0.975
          )
      )


    # ---------------------------------------------------------------------
    # Complete-model AIC/BIC comparison
    # ---------------------------------------------------------------------

    if (result$comparison_complete) {
      candidate_models <- list(
        stationary = stationary_model,
        sigma = sigma_model,
        nu = nu_model,
        joint = joint_model
      )


      AIC_values <- vapply(
        candidate_models,
        stats::AIC,
        numeric(1)
      )


      BIC_values <- vapply(
        candidate_models,
        function(model) {
          stats::deviance(
            model
          ) +
            log(
              current_n
            ) *
            model$df.fit
        },
        numeric(1)
      )


      result$AIC_selected_model <-
        names(
          which.min(
            AIC_values
          )
        )


      result$BIC_selected_model <-
        names(
          which.min(
            BIC_values
          )
        )


      result$AIC_selects_joint <-
        identical(
          result$AIC_selected_model,
          "joint"
        )


      result$BIC_selects_joint <-
        identical(
          result$BIC_selected_model,
          "joint"
        )
    }


    # ---------------------------------------------------------------------
    # Independent direct optimizer
    # ---------------------------------------------------------------------

    if (
      direct_checks_used <
        direct_checks_per_sample_size
    ) {
      direct_checks_used <-
        direct_checks_used + 1L


      direct_check <- run_direct_optimizer_check(
        joint_model = joint_model,
        data = simulation_data,
        lower = lower,
        seed = current_seed,
        n_starts = direct_optimizer_starts
      )


      result$direct_optimizer_attempted <-
        direct_check$attempted


      result$direct_optimizer_success <-
        direct_check$success


      result$direct_optimizer_error <-
        direct_check$error_message


      result$direct_gamlss_logLik <-
        direct_check$gamlss_logLik


      result$direct_optimizer_logLik <-
        direct_check$direct_logLik


      result$direct_logLik_difference <-
        direct_check$absolute_difference


      result$direct_agreement_1e3 <-
        direct_check$agreement_1e3
    }


    # ---------------------------------------------------------------------
    # Final status classification
    # ---------------------------------------------------------------------

    if (!result$comparison_complete) {
      failed_auxiliary_models <- c(
        if (!result$sigma_fit_success) {
          "sigma-only"
        },

        if (!result$nu_fit_success) {
          "nu-only"
        }
      )


      result$error_stage <-
        "comparison model"


      result$error_message <-
        paste0(
          "The joint model succeeded, but the following auxiliary ",
          "comparison model(s) failed: ",
          paste(
            failed_auxiliary_models,
            collapse = ", "
          ),
          "."
        )
    } else {
      result$error_stage <-
        NA_character_

      result$error_message <-
        NA_character_
    }


    rows[[row_index]] <- result
  }


  message(
    "Completed n = ",
    current_n
  )
}


# =========================================================================
# 9. Combine raw simulation results
# =========================================================================

simulation_results <- do.call(
  rbind,
  rows
)

rownames(
  simulation_results
) <- NULL


# =========================================================================
# 10. Summarize by sample size
# =========================================================================

summary_rows <- lapply(
  sample_sizes,
  function(current_n) {
    all_runs <- simulation_results[
      simulation_results$n ==
        current_n,
      ,
      drop = FALSE
    ]


    joint_runs <- all_runs[
      all_runs$
        joint_fit_success %in% TRUE,
      ,
      drop = FALSE
    ]


    comparison_runs <- all_runs[
      all_runs$
        comparison_complete %in% TRUE,
      ,
      drop = FALSE
    ]


    direct_attempts <- all_runs[
      all_runs$
        direct_optimizer_attempted %in% TRUE,
      ,
      drop = FALSE
    ]


    direct_successes <- direct_attempts[
      direct_attempts$
        direct_optimizer_success %in% TRUE &
        is.finite(
          direct_attempts$
            direct_logLik_difference
        ),
      ,
      drop = FALSE
    ]


    attempted_count <-
      nrow(
        all_runs
      )


    joint_success_count <-
      nrow(
        joint_runs
      )


    comparison_count <-
      nrow(
        comparison_runs
      )


    direct_attempt_count <-
      nrow(
        direct_attempts
      )


    direct_success_count <-
      nrow(
        direct_successes
      )


    data.frame(
      n =
        current_n,

      attempted =
        attempted_count,

      stationary_success_rate =
        safe_rate(
          all_runs$
            stationary_fit_success
        ),

      sigma_only_success_rate =
        safe_rate(
          all_runs$
            sigma_fit_success
        ),

      nu_only_success_rate =
        safe_rate(
          all_runs$
            nu_fit_success
        ),

      joint_successful =
        joint_success_count,

      joint_success_rate =
        joint_success_count /
        attempted_count,

      comparisons_complete =
        comparison_count,

      comparison_complete_rate =
        comparison_count /
        attempted_count,

      exact_mu_below_lower_rate =
        safe_rate(
          joint_runs$
            mu_below_lower_exact
        ),

      boundary_contact_rate =
        safe_rate(
          joint_runs$
            mu_boundary_contact
        ),

      material_boundary_violation_rate =
        safe_rate(
          all_runs$
            material_boundary_violation
        ),

      minimum_recorded_lower_minus_mu =
        if (nrow(joint_runs) > 0L) {
          min(
            joint_runs$
              minimum_lower_minus_mu,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },

      maximum_boundary_tolerance =
        safe_max(
          joint_runs$
            boundary_tolerance
        ),

      mu_bias =
        safe_mean(
          joint_runs$mu_hat -
            truth$mu
        ),

      mu_rmse =
        safe_rmse(
          joint_runs$mu_hat,
          truth$mu
        ),

      sigma_intercept_bias =
        safe_mean(
          joint_runs$
            sigma_intercept_hat -
            truth$sigma_intercept
        ),

      sigma_slope_bias =
        safe_mean(
          joint_runs$
            sigma_slope_hat -
            truth$sigma_slope
        ),

      sigma_slope_rmse =
        safe_rmse(
          joint_runs$
            sigma_slope_hat,
          truth$sigma_slope
        ),

      mean_sigma_curve_bias =
        safe_mean(
          joint_runs$
            sigma_curve_bias
        ),

      mean_sigma_curve_rmse =
        safe_mean(
          joint_runs$
            sigma_curve_rmse
        ),

      mean_sigma_curve_correlation =
        safe_mean(
          joint_runs$
            sigma_curve_correlation
        ),

      nu_intercept_bias =
        safe_mean(
          joint_runs$
            nu_intercept_hat -
            truth$nu_intercept
        ),

      nu_slope_bias =
        safe_mean(
          joint_runs$
            nu_slope_hat -
            truth$nu_slope
        ),

      nu_slope_rmse =
        safe_rmse(
          joint_runs$
            nu_slope_hat,
          truth$nu_slope
        ),

      mean_nu_curve_bias =
        safe_mean(
          joint_runs$
            nu_curve_bias
        ),

      mean_nu_curve_rmse =
        safe_mean(
          joint_runs$
            nu_curve_rmse
        ),

      mean_nu_curve_correlation =
        safe_mean(
          joint_runs$
            nu_curve_correlation
        ),

      AIC_joint_selection_rate =
        safe_rate(
          comparison_runs$
            AIC_selects_joint
        ),

      BIC_joint_selection_rate =
        safe_rate(
          comparison_runs$
            BIC_selects_joint
        ),

      mean_residual_mean =
        safe_mean(
          joint_runs$
            residual_mean
        ),

      mean_residual_sd =
        safe_mean(
          joint_runs$
            residual_sd
        ),

      mean_residual_tail_rate_95 =
        safe_mean(
          joint_runs$
            residual_tail_rate_95
        ),

      direct_comparisons_attempted =
        direct_attempt_count,

      direct_comparisons_completed =
        direct_success_count,

      direct_comparison_success_rate =
        if (direct_attempt_count > 0L) {
          direct_success_count /
            direct_attempt_count
        } else {
          NA_real_
        },

      direct_agreement_rate_1e3 =
        safe_rate(
          direct_successes$
            direct_agreement_1e3
        ),

      mean_direct_logLik_difference =
        safe_mean(
          direct_successes$
            direct_logLik_difference
        ),

      maximum_direct_logLik_difference =
        safe_max(
          direct_successes$
            direct_logLik_difference
        ),

      stringsAsFactors = FALSE
    )
  }
)


simulation_summary <- do.call(
  rbind,
  summary_rows
)

rownames(
  simulation_summary
) <- NULL


# =========================================================================
# 11. Save outputs
# =========================================================================

results_directory <- file.path(
  LL3_PROJECT_ROOT,
  "validation",
  "results"
)


dir.create(
  results_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


utils::write.csv(
  simulation_results,
  file.path(
    results_directory,
    "joint_nonstationarity_raw.csv"
  ),
  row.names = FALSE
)


utils::write.csv(
  simulation_summary,
  file.path(
    results_directory,
    "joint_nonstationarity_summary.csv"
  ),
  row.names = FALSE
)


joint_failures <- simulation_results[
  !simulation_results$
    joint_fit_success,
  c(
    "n",
    "seed",
    "joint_converged",
    "support_ok",
    "sigma_positive",
    "nu_positive",
    "mu_below_lower_exact",
    "mu_not_materially_above_lower",
    "mu_boundary_contact",
    "material_boundary_violation",
    "minimum_y_minus_mu",
    "minimum_lower_minus_mu",
    "boundary_tolerance",
    "error_stage",
    "error_message"
  ),
  drop = FALSE
]


utils::write.csv(
  joint_failures,
  file.path(
    results_directory,
    "joint_nonstationarity_joint_failures.csv"
  ),
  row.names = FALSE
)


comparison_failures <- simulation_results[
  simulation_results$
    joint_fit_success %in% TRUE &
    !simulation_results$
      comparison_complete,
  c(
    "n",
    "seed",
    "stationary_fit_success",
    "sigma_fit_success",
    "nu_fit_success",
    "joint_fit_success",
    "error_stage",
    "error_message"
  ),
  drop = FALSE
]


utils::write.csv(
  comparison_failures,
  file.path(
    results_directory,
    "joint_nonstationarity_comparison_failures.csv"
  ),
  row.names = FALSE
)


direct_failures <- simulation_results[
  simulation_results$
    direct_optimizer_attempted %in% TRUE &
    !simulation_results$
      direct_optimizer_success,
  c(
    "n",
    "seed",
    "direct_optimizer_error"
  ),
  drop = FALSE
]


utils::write.csv(
  direct_failures,
  file.path(
    results_directory,
    "joint_nonstationarity_direct_failures.csv"
  ),
  row.names = FALSE
)


# =========================================================================
# 12. Print summary
# =========================================================================

print(
  simulation_summary,
  digits = 10,
  row.names = FALSE
)


if (nrow(joint_failures) > 0L) {
  message(
    "\nJoint-model failures: ",
    nrow(joint_failures)
  )

  print(
    joint_failures,
    digits = 12,
    row.names = FALSE
  )
} else {
  message(
    "\nAll joint-model fits completed successfully."
  )
}


if (nrow(comparison_failures) > 0L) {
  message(
    "\nSuccessful joint fits with incomplete auxiliary comparisons: ",
    nrow(comparison_failures)
  )

  print(
    comparison_failures,
    row.names = FALSE
  )
} else {
  message(
    "\nAll successful joint fits had complete auxiliary comparisons."
  )
}


if (nrow(direct_failures) > 0L) {
  message(
    "\nIndependent direct-optimizer failures: ",
    nrow(direct_failures),
    ". Failed comparisons were excluded from likelihood summaries."
  )

  print(
    direct_failures,
    row.names = FALSE
  )
} else {
  message(
    "\nAll attempted independent direct optimizations completed successfully."
  )
}


message(
  "\nJoint nonstationarity results written to:\n",
  results_directory
)