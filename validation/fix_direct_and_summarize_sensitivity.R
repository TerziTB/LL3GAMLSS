# validation/fix_direct_and_summarize_sensitivity.R
#
# Corrected targeted direct-optimizer validation for the two discrepant
# joint nonstationarity simulations.
#
# Also summarizes the existing lower-bound sensitivity results separately
# for boundary-contact and non-boundary cases.
#
# The independent optimizer below evaluates the LL3 log-likelihood directly
# from its mathematical formula. It does not call dLL3().


# =========================================================================
# 1. Project setup
# =========================================================================

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

source(file.path(project_dir, "R", "LL3_distribution.R"))
source(file.path(project_dir, "R", "LL3_gamlss_family.R"))
source(file.path(project_dir, "R", "LL3_fitting.R"))
source(file.path(project_dir, "R", "LL3_diagnostics.R"))

if (!requireNamespace("gamlss", quietly = TRUE)) {
  stop("The gamlss package must be installed.")
}


# =========================================================================
# 2. Joint simulation specification
# =========================================================================

truth_mu <- -50

truth_sigma_intercept <- log(70)
truth_sigma_slope <- 0.55

truth_nu_intercept <- log(1.8)
truth_nu_slope <- -0.20

discrepancy_threshold <- 1e-3

number_of_random_starts <- 50L

quantile_probabilities <- c(
  0.02,
  0.05,
  0.10,
  0.50,
  0.90,
  0.95,
  0.98
)


# =========================================================================
# 3. Utilities
# =========================================================================

stable_log1pexp <- function(x) {
  result <- numeric(length(x))

  positive <- x > 0

  result[positive] <-
    x[positive] + log1p(exp(-x[positive]))

  result[!positive] <-
    log1p(exp(x[!positive]))

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


safe_rate <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(as.logical(x))
}


safe_median <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  stats::median(x)
}


safe_max <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  max(x)
}


as_logical_safe <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  if (is.numeric(x)) {
    return(!is.na(x) & x != 0)
  }

  toupper(trimws(as.character(x))) %in%
    c("TRUE", "T", "YES", "Y", "1")
}


simulate_joint_dataset <- function(n, seed) {
  n <- as.integer(n)
  seed <- as.integer(seed)

  set.seed(seed)

  x <- seq(
    -1,
    1,
    length.out = n
  )

  sigma <- exp(
    truth_sigma_intercept +
      truth_sigma_slope * x
  )

  nu <- exp(
    truth_nu_intercept +
      truth_nu_slope * x
  )

  y <- rLL3(
    n = n,
    mu = truth_mu,
    sigma = sigma,
    nu = nu
  )

  list(
    data = data.frame(
      y = y,
      x = x
    ),
    y = y,
    x = x
  )
}


safe_mu_start <- function(mu, lower, y) {
  margin <- max(
    1e-8,
    sqrt(.Machine$double.eps) *
      max(
        1,
        abs(lower),
        abs(y)
      )
  )

  pmin(
    as.numeric(mu),
    lower - margin
  )
}


recover_linear_coefficients <- function(values, x) {
  design <- cbind(
    intercept = 1,
    slope = x
  )

  as.numeric(
    qr.coef(
      qr(design),
      values
    )
  )
}


quantile_matrix <- function(
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


# =========================================================================
# 4. GAMLSS fitting
# =========================================================================

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


fit_one_gamlss_model <- function(
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
      arguments$mu.start <- safe_mu_start(
        mu.start,
        lower = lower,
        y = data$y
      )
    }

    if (!is.null(sigma.start)) {
      arguments$sigma.start <- as.numeric(sigma.start)
    }

    if (!is.null(nu.start)) {
      arguments$nu.start <- as.numeric(nu.start)
    }

    arguments
  }

  run_attempt <- function(control) {
    warnings <- character()

    fit <- withCallingHandlers(
      tryCatch(
        do.call(
          fit_LL3_gamlss,
          construct_arguments(control)
        ),
        error = function(e) {
          e
        }
      ),
      warning = function(w) {
        warnings <<- c(
          warnings,
          conditionMessage(w)
        )

        invokeRestart("muffleWarning")
      }
    )

    list(
      fit = fit,
      warnings = unique(warnings)
    )
  }

  validate_attempt <- function(attempt) {
    if (inherits(attempt$fit, "error")) {
      return(
        list(
          success = FALSE,
          message = conditionMessage(attempt$fit)
        )
      )
    }

    parameters <- try(
      extract_LL3_parameters(attempt$fit),
      silent = TRUE
    )

    if (inherits(parameters, "try-error")) {
      return(
        list(
          success = FALSE,
          message = as.character(parameters)
        )
      )
    }

    mu <- as.numeric(parameters$mu)
    sigma <- as.numeric(parameters$sigma)
    nu <- as.numeric(parameters$nu)

    success <- isTRUE(attempt$fit$converged) &&
      all(is.finite(mu)) &&
      all(is.finite(sigma)) &&
      all(is.finite(nu)) &&
      all(data$y > mu) &&
      all(sigma > 0) &&
      all(nu > 0)

    list(
      success = success,
      parameters = parameters,
      message =
        if (success) {
          NA_character_
        } else {
          "The GAMLSS fit failed convergence or parameter validity checks."
        }
    )
  }

  first_attempt <- run_attempt(primary_control)
  first_validation <- validate_attempt(first_attempt)

  if (first_validation$success) {
    return(
      list(
        success = TRUE,
        fit = first_attempt$fit,
        parameters = first_validation$parameters,
        retried = FALSE
      )
    )
  }

  second_attempt <- run_attempt(retry_control)
  second_validation <- validate_attempt(second_attempt)

  list(
    success = second_validation$success,
    fit = second_attempt$fit,
    parameters = second_validation$parameters,
    retried = TRUE,
    message = second_validation$message
  )
}


fit_joint_gamlss <- function(data, lower) {
  stationary <- fit_one_gamlss_model(
    mu.formula = y ~ 1,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,
    data = data,
    lower = lower
  )

  if (!stationary$success) {
    stop(
      "Stationary GAMLSS fit failed: ",
      stationary$message
    )
  }

  stationary_parameters <- stationary$parameters

  sigma_only <- fit_one_gamlss_model(
    mu.formula = y ~ 1,
    sigma.formula = ~ x,
    nu.formula = ~ 1,
    data = data,
    lower = lower,
    mu.start = stationary_parameters$mu,
    sigma.start = stationary_parameters$sigma,
    nu.start = stationary_parameters$nu
  )

  nu_only <- fit_one_gamlss_model(
    mu.formula = y ~ 1,
    sigma.formula = ~ 1,
    nu.formula = ~ x,
    data = data,
    lower = lower,
    mu.start = stationary_parameters$mu,
    sigma.start = stationary_parameters$sigma,
    nu.start = stationary_parameters$nu
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

  joint <- fit_one_gamlss_model(
    mu.formula = y ~ 1,
    sigma.formula = ~ x,
    nu.formula = ~ x,
    data = data,
    lower = lower,
    mu.start = stationary_parameters$mu,
    sigma.start = sigma_start,
    nu.start = nu_start
  )

  if (!joint$success) {
    stop(
      "Joint GAMLSS fit failed: ",
      joint$message
    )
  }

  parameters <- extract_LL3_parameters(
    joint$fit
  )

  mu <- as.numeric(parameters$mu)
  sigma <- as.numeric(parameters$sigma)
  nu <- as.numeric(parameters$nu)

  sigma_coefficients <- recover_linear_coefficients(
    log(sigma),
    data$x
  )

  nu_coefficients <- recover_linear_coefficients(
    log(nu),
    data$x
  )

  list(
    fit = joint$fit,
    logLik = as.numeric(stats::logLik(joint$fit)),
    mu = mean(mu),
    sigma = sigma,
    nu = nu,
    sigma_intercept = sigma_coefficients[1L],
    sigma_slope = sigma_coefficients[2L],
    nu_intercept = nu_coefficients[1L],
    nu_slope = nu_coefficients[2L]
  )
}


# =========================================================================
# 5. Independent formula-only LL3 likelihood
#
# Important correction:
#
# The support value is named support_lower. This avoids collision with the
# optimizer's own lower parameter-bound argument.
# =========================================================================

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

  powered_log_z <-
    nu *
    log_z

  log_density <-
    eta_nu -
    eta_sigma +
    (nu - 1) *
    log_z -
    2 *
    stable_log1pexp(
      powered_log_z
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


direct_curves <- function(
    theta,
    x,
    support_lower
) {
  list(
    mu =
      support_lower -
      exp(theta[1L]),

    sigma =
      exp(
        theta[2L] +
          theta[3L] * x
      ),

    nu =
      exp(
        theta[4L] +
          theta[5L] * x
      ),

    sigma_intercept =
      theta[2L],

    sigma_slope =
      theta[3L],

    nu_intercept =
      theta[4L],

    nu_slope =
      theta[5L]
  )
}


# =========================================================================
# 6. Corrected multistart direct optimizer
# =========================================================================

run_direct_optimizer <- function(
    y,
    x,
    support_lower,
    gamlss_result,
    seed
) {
  scale_y <- data_scale(y)

  initial_gap <-
    support_lower -
    gamlss_result$mu

  if (
    !is.finite(initial_gap) ||
      initial_gap <= 0
  ) {
    initial_gap <- scale_y
  }

  base_start <- c(
    log(initial_gap),
    gamlss_result$sigma_intercept,
    gamlss_result$sigma_slope,
    gamlss_result$nu_intercept,
    gamlss_result$nu_slope
  )

  generic_start <- c(
    log(scale_y),
    log(max(stats::sd(y), 1)),
    0,
    log(2),
    0
  )

  parameter_lower_bounds <- c(
    log(max(1e-12 * scale_y, 1e-12)),
    -20,
    -10,
    -10,
    -10
  )

  parameter_upper_bounds <- c(
    log(max(1e4 * scale_y, 1e4)),
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

  deterministic_starts <- rbind(
    base_start,
    generic_start,
    base_start + c(0.25, 0, 0, 0, 0),
    base_start - c(0.25, 0, 0, 0, 0),
    base_start + c(0, 0.20, 0.10, 0, 0),
    base_start - c(0, 0.20, 0.10, 0, 0),
    base_start + c(0, 0, 0, 0.20, 0.10),
    base_start - c(0, 0, 0, 0.20, 0.10)
  )

  random_seed <- (
    as.double(seed) +
      91027
  ) %% .Machine$integer.max

  set.seed(as.integer(random_seed))

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
      0.75,
      0.50,
      0.40,
      0.50,
      0.40
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

  candidates <- list()

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
      candidates[[length(candidates) + 1L]] <<-
        list(
          theta = as.numeric(theta),
          objective = as.numeric(objective),
          converged = isTRUE(converged),
          method = method
        )
    }

    invisible(NULL)
  }

  for (start_index in seq_len(nrow(starts))) {
    current_start <- starts[
      start_index,
      ,
      drop = TRUE
    ]

    optim_result <- try(
      stats::optim(
        par = current_start,
        fn = direct_LL3_negative_loglikelihood,
        y = y,
        x = x,
        support_lower = support_lower,
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

    if (!inherits(optim_result, "try-error")) {
      add_candidate(
        theta = optim_result$par,
        objective = optim_result$value,
        converged = optim_result$convergence == 0,
        method = paste0(
          "L-BFGS-B start ",
          start_index
        )
      )
    }

    nlminb_result <- try(
      stats::nlminb(
        start = current_start,
        objective = direct_LL3_negative_loglikelihood,
        y = y,
        x = x,
        support_lower = support_lower,
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
      add_candidate(
        theta = nlminb_result$par,
        objective = nlminb_result$objective,
        converged = nlminb_result$convergence == 0,
        method = paste0(
          "nlminb start ",
          start_index
        )
      )
    }
  }

  if (length(candidates) == 0L) {
    stop(
      "Every corrected direct-optimizer attempt failed."
    )
  }

  converged_flags <- vapply(
    candidates,
    function(candidate) {
      candidate$converged
    },
    logical(1)
  )

  eligible_indices <-
    if (any(converged_flags)) {
      which(converged_flags)
    } else {
      seq_along(candidates)
    }

  objectives <- vapply(
    candidates[eligible_indices],
    function(candidate) {
      candidate$objective
    },
    numeric(1)
  )

  best_index <- eligible_indices[
    which.min(objectives)
  ]

  best <- candidates[[best_index]]

  polished <- try(
    stats::optim(
      par = best$theta,
      fn = direct_LL3_negative_loglikelihood,
      y = y,
      x = x,
      support_lower = support_lower,
      method = "L-BFGS-B",
      lower = parameter_lower_bounds,
      upper = parameter_upper_bounds,
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
    !inherits(polished, "try-error") &&
      is.finite(polished$value) &&
      polished$value < best$objective
  ) {
    best <- list(
      theta = polished$par,
      objective = polished$value,
      converged = polished$convergence == 0,
      method = paste0(
        best$method,
        " + polish"
      )
    )
  }

  list(
    theta = best$theta,
    logLik = -best$objective,
    converged = best$converged,
    method = best$method,
    attempts = length(candidates),
    converged_attempts = sum(converged_flags)
  )
}


# =========================================================================
# 7. Read and rerun discrepant cases
# =========================================================================

joint_raw_file <- file.path(
  results_dir,
  "joint_nonstationarity_raw.csv"
)

if (!file.exists(joint_raw_file)) {
  stop(
    "Missing file: ",
    joint_raw_file
  )
}

joint_raw <- utils::read.csv(
  joint_raw_file,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "n",
  "seed",
  "direct_optimizer_success",
  "direct_gamlss_logLik",
  "direct_optimizer_logLik",
  "direct_logLik_difference"
)

missing_columns <- setdiff(
  required_columns,
  names(joint_raw)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing columns in joint_nonstationarity_raw.csv: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

discrepant <- joint_raw[
  as_logical_safe(
    joint_raw$direct_optimizer_success
  ) &
    is.finite(
      joint_raw$direct_logLik_difference
    ) &
    abs(
      joint_raw$direct_logLik_difference
    ) >= discrepancy_threshold,
  ,
  drop = FALSE
]

if (nrow(discrepant) == 0L) {
  stop(
    "No discrepant direct-optimizer cases were found."
  )
}

direct_rows <- vector(
  "list",
  nrow(discrepant)
)

for (case_index in seq_len(nrow(discrepant))) {
  current_n <- as.integer(
    discrepant$n[case_index]
  )

  current_seed <- as.integer(
    discrepant$seed[case_index]
  )

  cat(
    "\nCorrected direct rerun ",
    case_index,
    "/",
    nrow(discrepant),
    ": n = ",
    current_n,
    ", seed = ",
    current_seed,
    "\n",
    sep = ""
  )

  simulated <- simulate_joint_dataset(
    n = current_n,
    seed = current_seed
  )

  support_lower <- LL3_support_lower(
    simulated$y
  )

  gamlss_result <- fit_joint_gamlss(
    data = simulated$data,
    lower = support_lower
  )

  direct_result <- run_direct_optimizer(
    y = simulated$y,
    x = simulated$x,
    support_lower = support_lower,
    gamlss_result = gamlss_result,
    seed = current_seed
  )

  direct_parameter_curves <- direct_curves(
    theta = direct_result$theta,
    x = simulated$x,
    support_lower = support_lower
  )

  gamlss_quantiles <- quantile_matrix(
    probabilities = quantile_probabilities,
    mu = gamlss_result$mu,
    sigma = gamlss_result$sigma,
    nu = gamlss_result$nu
  )

  direct_quantiles <- quantile_matrix(
    probabilities = quantile_probabilities,
    mu = direct_parameter_curves$mu,
    sigma = direct_parameter_curves$sigma,
    nu = direct_parameter_curves$nu
  )

  logLik_difference <-
    direct_result$logLik -
    gamlss_result$logLik

  scale_y <- data_scale(
    simulated$y
  )

  direct_rows[[case_index]] <- data.frame(
    n = current_n,
    seed = current_seed,

    original_gamlss_logLik =
      discrepant$direct_gamlss_logLik[
        case_index
      ],

    original_direct_logLik =
      discrepant$direct_optimizer_logLik[
        case_index
      ],

    original_absolute_difference =
      abs(
        discrepant$direct_logLik_difference[
          case_index
        ]
      ),

    refitted_gamlss_logLik =
      gamlss_result$logLik,

    corrected_direct_logLik =
      direct_result$logLik,

    direct_minus_gamlss =
      logLik_difference,

    absolute_logLik_difference =
      abs(logLik_difference),

    agreement_below_1e_3 =
      abs(logLik_difference) <
      1e-3,

    agreement_below_1e_6 =
      abs(logLik_difference) <
      1e-6,

    likelihood_winner =
      if (
        logLik_difference > 1e-7
      ) {
        "direct optimizer"
      } else if (
        logLik_difference < -1e-7
      ) {
        "GAMLSS"
      } else {
        "numerical tie"
      },

    gamlss_mu =
      gamlss_result$mu,

    direct_mu =
      direct_parameter_curves$mu,

    absolute_mu_difference =
      abs(
        gamlss_result$mu -
          direct_parameter_curves$mu
      ),

    gamlss_sigma_intercept =
      gamlss_result$sigma_intercept,

    direct_sigma_intercept =
      direct_parameter_curves$sigma_intercept,

    gamlss_sigma_slope =
      gamlss_result$sigma_slope,

    direct_sigma_slope =
      direct_parameter_curves$sigma_slope,

    gamlss_nu_intercept =
      gamlss_result$nu_intercept,

    direct_nu_intercept =
      direct_parameter_curves$nu_intercept,

    gamlss_nu_slope =
      gamlss_result$nu_slope,

    direct_nu_slope =
      direct_parameter_curves$nu_slope,

    maximum_absolute_log_sigma_difference =
      max(
        abs(
          log(gamlss_result$sigma) -
            log(direct_parameter_curves$sigma)
        )
      ),

    maximum_absolute_log_nu_difference =
      max(
        abs(
          log(gamlss_result$nu) -
            log(direct_parameter_curves$nu)
        )
      ),

    maximum_scaled_quantile_difference =
      max(
        abs(
          gamlss_quantiles -
            direct_quantiles
        )
      ) /
      scale_y,

    optimizer_converged =
      direct_result$converged,

    optimizer_method =
      direct_result$method,

    optimizer_candidates =
      direct_result$attempts,

    converged_candidates =
      direct_result$converged_attempts,

    stringsAsFactors = FALSE
  )
}

direct_results <- do.call(
  rbind,
  direct_rows
)

direct_output_file <- file.path(
  results_dir,
  "joint_direct_discrepancies_corrected.csv"
)

utils::write.csv(
  direct_results,
  direct_output_file,
  row.names = FALSE
)

cat(
  "\n============================================================\n",
  "Corrected direct-optimizer results\n",
  "============================================================\n",
  sep = ""
)

print(
  direct_results,
  digits = 15,
  row.names = FALSE
)


# =========================================================================
# 8. Stratify existing sensitivity analysis by default boundary status
# =========================================================================

sensitivity_file <- file.path(
  results_dir,
  "lower_bound_sensitivity_raw.csv"
)

if (!file.exists(sensitivity_file)) {
  stop(
    "Missing sensitivity file: ",
    sensitivity_file
  )
}

sensitivity <- utils::read.csv(
  sensitivity_file,
  stringsAsFactors = FALSE
)

sensitivity$complete <- as_logical_safe(
  sensitivity$complete
)

sensitivity$boundary_contact <- as_logical_safe(
  sensitivity$boundary_contact
)

sensitivity$AIC_selection_same_as_default <-
  as_logical_safe(
    sensitivity$AIC_selection_same_as_default
  )

sensitivity$BIC_selection_same_as_default <-
  as_logical_safe(
    sensitivity$BIC_selection_same_as_default
  )

default_rows <- sensitivity[
  sensitivity$lower_label ==
    "project_default",
  c(
    "n",
    "seed",
    "complete",
    "boundary_contact"
  ),
  drop = FALSE
]

names(default_rows)[
  names(default_rows) ==
    "complete"
] <- "default_complete"

names(default_rows)[
  names(default_rows) ==
    "boundary_contact"
] <- "default_boundary_contact"

sensitivity_stratified <- merge(
  sensitivity,
  default_rows,
  by = c(
    "n",
    "seed"
  ),
  all.x = TRUE,
  sort = FALSE
)

sensitivity_stratified$boundary_group <-
  ifelse(
    sensitivity_stratified$
      default_boundary_contact,
    "default boundary contact",
    "no default boundary contact"
  )

summary_groups <- split(
  sensitivity_stratified,
  interaction(
    sensitivity_stratified$lower_label,
    sensitivity_stratified$boundary_group,
    drop = TRUE
  )
)

stratified_rows <- lapply(
  summary_groups,
  function(current) {
    complete_current <- current[
      current$complete,
      ,
      drop = FALSE
    ]

    data.frame(
      lower_label =
        current$lower_label[1L],

      boundary_group =
        current$boundary_group[1L],

      cases =
        nrow(current),

      complete =
        nrow(complete_current),

      completion_rate =
        safe_rate(current$complete),

      AIC_selection_agreement =
        safe_rate(
          complete_current$
            AIC_selection_same_as_default
        ),

      BIC_selection_agreement =
        safe_rate(
          complete_current$
            BIC_selection_same_as_default
        ),

      median_absolute_logLik_change =
        safe_median(
          abs(
            complete_current$
              delta_joint_logLik_from_default
          )
        ),

      maximum_absolute_logLik_change =
        safe_max(
          abs(
            complete_current$
              delta_joint_logLik_from_default
          )
        ),

      median_scaled_quantile_change =
        safe_median(
          complete_current$
            maximum_scaled_quantile_difference
        ),

      maximum_scaled_quantile_change =
        safe_max(
          complete_current$
            maximum_scaled_quantile_difference
        ),

      stringsAsFactors = FALSE
    )
  }
)

stratified_summary <- do.call(
  rbind,
  stratified_rows
)

stratified_summary <- stratified_summary[
  order(
    stratified_summary$lower_label,
    stratified_summary$boundary_group
  ),
  ,
  drop = FALSE
]

rownames(stratified_summary) <- NULL

stratified_output_file <- file.path(
  results_dir,
  "lower_bound_sensitivity_stratified.csv"
)

utils::write.csv(
  stratified_summary,
  stratified_output_file,
  row.names = FALSE
)

cat(
  "\n============================================================\n",
  "Sensitivity stratified by default boundary contact\n",
  "============================================================\n",
  sep = ""
)

print(
  stratified_summary,
  digits = 12,
  row.names = FALSE
)


# =========================================================================
# 9. Final status
# =========================================================================

unresolved <- direct_results[
  !direct_results$agreement_below_1e_3,
  ,
  drop = FALSE
]

cat(
  "\nDirect cases remaining above 1e-3: ",
  nrow(unresolved),
  "\n",
  sep = ""
)

cat(
  "\nFiles written:\n",
  direct_output_file,
  "\n",
  stratified_output_file,
  "\n",
  sep = ""
)

cat(
  "\nCorrected validation completed.\n"
)
