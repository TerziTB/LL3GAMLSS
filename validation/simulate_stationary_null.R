# validation/simulate_stationary_null.R
#
# Monte Carlo validation under a stationary LL3 data-generating process.
#
# True model:
#
#   mu_t    = mu
#   sigma_t = sigma
#   nu_t    = nu
#
# Candidate models:
#
#   stationary: sigma ~ 1,    nu ~ 1
#   sigma:      sigma ~ time, nu ~ 1
#   nu:         sigma ~ 1,    nu ~ time
#   joint:      sigma ~ time, nu ~ time
#
# The script estimates false nonstationary model-selection rates for AIC
# and BIC.
#
# Important bookkeeping rule:
#
# AIC, BIC, and residual diagnostics are stored only when:
#
#   * every candidate model converges;
#   * every fitted model passes support and positivity checks;
#   * the stationary fitted model produces valid, nondegenerate PIT
#     residuals;
#   * every information-criterion value is finite.
#
# Incomplete simulations retain NA for selection and residual results.


# =========================================================================
# 1. Locate the LL3 project root
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


if (
  !requireNamespace(
    "gamlss",
    quietly = TRUE
  )
) {
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
    "LL3_REPS must be a positive integer."
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
  sigma = 70,
  nu = 1.8
)


# Change this only when deliberately starting a new simulation experiment.

master_seed <- 2027L

set.seed(master_seed)

simulation_seeds <- sample.int(
  .Machine$integer.max,
  size =
    repetitions *
    length(sample_sizes),
  replace = FALSE
)


# =========================================================================
# 4. GAMLSS controls
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
# =========================================================================
# 5. Utility functions
# =========================================================================

safe_mean <- function(x) {
  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(x)
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


capture_fit <- function(expression) {
  warning_messages <- character()

  fitted_object <- withCallingHandlers(
    tryCatch(
      force(expression),
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
    warnings = unique(
      warning_messages
    )
  )
}


# =========================================================================
# 6. Model validation
# =========================================================================

validate_LL3_model <- function(
    model,
    y
) {
  output <- list(
    success = FALSE,
    parameters = NULL,
    checks = NULL,
    error_message = NA_character_
  )

  if (
    inherits(model, "error") ||
      inherits(model, "try-error")
  ) {
    output$error_message <-
      if (inherits(model, "error")) {
        conditionMessage(model)
      } else {
        extract_try_error_message(model)
      }

    return(output)
  }


  parameters <- try(
    extract_LL3_parameters(
      model
    ),
    silent = TRUE
  )

  if (inherits(parameters, "try-error")) {
    output$error_message <-
      extract_try_error_message(
        parameters
      )

    return(output)
  }


  checks <- try(
    check_LL3_fit(
      model,
      y = y
    ),
    silent = TRUE
  )

  if (inherits(checks, "try-error")) {
    output$error_message <-
      extract_try_error_message(
        checks
      )

    return(output)
  }


  # Compatibility fallback for an older diagnostic function.

  if (
    !"mu_not_materially_above_lower" %in%
      names(checks) &&
      "mu_below_lower" %in%
      names(checks)
  ) {
    checks <- c(
      checks,
      mu_not_materially_above_lower =
        isTRUE(
          checks[
            "mu_below_lower"
          ]
        )
    )
  }


  required_checks <- c(
    "converged",
    "finite_parameters",
    "support_ok",
    "sigma_positive",
    "nu_positive",
    "lower_below_data",
    "mu_not_materially_above_lower"
  )


  missing_checks <- setdiff(
    required_checks,
    names(checks)
  )

  if (length(missing_checks) > 0L) {
    output$error_message <- paste0(
      "Missing fit checks: ",
      paste(
        missing_checks,
        collapse = ", "
      ),
      "."
    )

    output$parameters <- parameters
    output$checks <- checks

    return(output)
  }


  failed_checks <- required_checks[
    !vapply(
      checks[
        required_checks
      ],
      isTRUE,
      logical(1)
    )
  ]


  output$parameters <- parameters
  output$checks <- checks

  output$success <-
    length(failed_checks) == 0L


  if (!output$success) {
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
# 7. Candidate-model fitting with one controlled retry
# =========================================================================

fit_candidate_model <- function(
    mu.formula,
    sigma.formula,
    nu.formula,
    data,
    lower,
    y,
    mu.start = NULL,
    sigma.start = NULL,
    nu.start = NULL
) {
  make_arguments <- function(control) {
    arguments <- list(
      mu.formula =
        mu.formula,

      sigma.formula =
        sigma.formula,

      nu.formula =
        nu.formula,

      data =
        data,

      lower =
        lower,

      information =
        "opg",

      allow_nonstationary_mu =
        FALSE,

      control =
        control
    )

    if (!is.null(mu.start)) {
      arguments$mu.start <-
        mu.start
    }

    if (!is.null(sigma.start)) {
      arguments$sigma.start <-
        sigma.start
    }

    if (!is.null(nu.start)) {
      arguments$nu.start <-
        nu.start
    }

    arguments
  }


  first_attempt <- capture_fit(
    do.call(
      fit_LL3_gamlss,
      make_arguments(
        primary_control
      )
    )
  )


  first_validation <- validate_LL3_model(
    first_attempt$fit,
    y = y
  )


  if (first_validation$success) {
    return(
      list(
        fit =
          first_attempt$fit,

        validation =
          first_validation,

        retried =
          FALSE,

        warnings =
          collapse_messages(
            first_attempt$warnings
          )
      )
    )
  }


  second_attempt <- capture_fit(
    do.call(
      fit_LL3_gamlss,
      make_arguments(
        retry_control
      )
    )
  )


  second_validation <- validate_LL3_model(
    second_attempt$fit,
    y = y
  )


  list(
    fit =
      second_attempt$fit,

    validation =
      second_validation,

    retried =
      TRUE,

    warnings =
      collapse_messages(
        c(
          first_attempt$warnings,
          second_attempt$warnings
        )
      ),

    first_error =
      first_validation$error_message,

    retry_error =
      second_validation$error_message
  )
}


# =========================================================================
# 8. Empty simulation-result row
# =========================================================================

new_stationary_null_result <- function(
    n,
    seed
) {
  data.frame(
    n =
      as.integer(n),

    seed =
      as.integer(seed),


    complete =
      FALSE,


    stationary_fit_success =
      FALSE,

    sigma_fit_success =
      FALSE,

    nu_fit_success =
      FALSE,

    joint_fit_success =
      FALSE,


    stationary_retried =
      FALSE,

    sigma_retried =
      FALSE,

    nu_retried =
      FALSE,

    joint_retried =
      FALSE,


    stationary_warning =
      NA_character_,

    sigma_warning =
      NA_character_,

    nu_warning =
      NA_character_,

    joint_warning =
      NA_character_,


    AIC_stationary =
      NA_real_,

    AIC_sigma =
      NA_real_,

    AIC_nu =
      NA_real_,

    AIC_joint =
      NA_real_,


    BIC_stationary =
      NA_real_,

    BIC_sigma =
      NA_real_,

    BIC_nu =
      NA_real_,

    BIC_joint =
      NA_real_,


    AIC_selected =
      NA_character_,

    BIC_selected =
      NA_character_,


    AIC_false_nonstationary =
      NA,

    BIC_false_nonstationary =
      NA,


    residual_mean =
      NA_real_,

    residual_sd =
      NA_real_,

    residual_tail_rate_95 =
      NA_real_,

    pit_minimum =
      NA_real_,

    pit_maximum =
      NA_real_,

    pit_clipped_count =
      NA_integer_,


    error_stage =
      NA_character_,

    error_message =
      NA_character_,


    stringsAsFactors =
      FALSE
  )
}


# =========================================================================
# 9. Run simulations
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
  for (replication in seq_len(repetitions)) {
    row_index <- row_index + 1L

    current_seed <-
      simulation_seeds[
        row_index
      ]


    result <- new_stationary_null_result(
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


    y <- rLL3(
      n = current_n,
      mu = truth$mu,
      sigma = truth$sigma,
      nu = truth$nu
    )


    simulation_data <- data.frame(
      y = y,
      time = time
    )


    lower <- LL3_support_lower(
      y
    )


    # ---------------------------------------------------------------------
    # Stationary model
    # ---------------------------------------------------------------------

    stationary_result <- fit_candidate_model(
      mu.formula =
        y ~ 1,

      sigma.formula =
        ~ 1,

      nu.formula =
        ~ 1,

      data =
        simulation_data,

      lower =
        lower,

      y =
        y
    )


    result$stationary_fit_success <-
      isTRUE(
        stationary_result$
          validation$
          success
      )


    result$stationary_retried <-
      stationary_result$retried


    result$stationary_warning <-
      stationary_result$warnings


    if (!result$stationary_fit_success) {
      result$error_stage <-
        "stationary model"

      result$error_message <-
        stationary_result$
          validation$
          error_message

      rows[[row_index]] <- result

      next
    }


    stationary_parameters <-
      stationary_result$
        validation$
        parameters


    # ---------------------------------------------------------------------
    # Sigma-only nonstationary model
    # ---------------------------------------------------------------------

    sigma_result <- fit_candidate_model(
      mu.formula =
        y ~ 1,

      sigma.formula =
        ~ time,

      nu.formula =
        ~ 1,

      data =
        simulation_data,

      lower =
        lower,

      y =
        y,

      mu.start =
        stationary_parameters$mu,

      sigma.start =
        stationary_parameters$sigma,

      nu.start =
        stationary_parameters$nu
    )


    result$sigma_fit_success <-
      isTRUE(
        sigma_result$
          validation$
          success
      )


    result$sigma_retried <-
      sigma_result$retried


    result$sigma_warning <-
      sigma_result$warnings


    # ---------------------------------------------------------------------
    # Nu-only nonstationary model
    # ---------------------------------------------------------------------

    nu_result <- fit_candidate_model(
      mu.formula =
        y ~ 1,

      sigma.formula =
        ~ 1,

      nu.formula =
        ~ time,

      data =
        simulation_data,

      lower =
        lower,

      y =
        y,

      mu.start =
        stationary_parameters$mu,

      sigma.start =
        stationary_parameters$sigma,

      nu.start =
        stationary_parameters$nu
    )


    result$nu_fit_success <-
      isTRUE(
        nu_result$
          validation$
          success
      )


    result$nu_retried <-
      nu_result$retried


    result$nu_warning <-
      nu_result$warnings


    # ---------------------------------------------------------------------
    # Joint-model starting values
    # ---------------------------------------------------------------------

    sigma_start <-
      stationary_parameters$sigma


    if (result$sigma_fit_success) {
      sigma_start <-
        sigma_result$
          validation$
          parameters$
          sigma
    }


    nu_start <-
      stationary_parameters$nu


    if (result$nu_fit_success) {
      nu_start <-
        nu_result$
          validation$
          parameters$
          nu
    }


    # ---------------------------------------------------------------------
    # Joint sigma-nu nonstationary model
    # ---------------------------------------------------------------------

    joint_result <- fit_candidate_model(
      mu.formula =
        y ~ 1,

      sigma.formula =
        ~ time,

      nu.formula =
        ~ time,

      data =
        simulation_data,

      lower =
        lower,

      y =
        y,

      mu.start =
        stationary_parameters$mu,

      sigma.start =
        sigma_start,

      nu.start =
        nu_start
    )


    result$joint_fit_success <-
      isTRUE(
        joint_result$
          validation$
          success
      )


    result$joint_retried <-
      joint_result$retried


    result$joint_warning <-
      joint_result$warnings


    # ---------------------------------------------------------------------
    # Require every candidate model to pass validation
    # ---------------------------------------------------------------------

    candidate_models_valid <- all(
      c(
        result$stationary_fit_success,
        result$sigma_fit_success,
        result$nu_fit_success,
        result$joint_fit_success
      )
    )


    if (!candidate_models_valid) {
      failed_models <- c(
        if (!result$stationary_fit_success) {
          "stationary"
        },

        if (!result$sigma_fit_success) {
          "sigma-only"
        },

        if (!result$nu_fit_success) {
          "nu-only"
        },

        if (!result$joint_fit_success) {
          "joint"
        }
      )


      failure_messages <- c(
        if (!result$stationary_fit_success) {
          paste0(
            "stationary: ",
            stationary_result$
              validation$
              error_message
          )
        },

        if (!result$sigma_fit_success) {
          paste0(
            "sigma-only: ",
            sigma_result$
              validation$
              error_message
          )
        },

        if (!result$nu_fit_success) {
          paste0(
            "nu-only: ",
            nu_result$
              validation$
              error_message
          )
        },

        if (!result$joint_fit_success) {
          paste0(
            "joint: ",
            joint_result$
              validation$
              error_message
          )
        }
      )


      result$error_stage <-
        "candidate model fitting"


      result$error_message <- paste0(
        "Failed candidate model(s): ",
        paste(
          failed_models,
          collapse = ", "
        ),
        ". ",
        paste(
          failure_messages,
          collapse = " | "
        )
      )


      rows[[row_index]] <- result

      next
    }


    # ---------------------------------------------------------------------
    # Residual calibration under the correctly specified stationary model
    #
    # Calculate residuals before assigning AIC or BIC selections. Thus,
    # an invalid residual calculation cannot leave misleading model-
    # selection values in an incomplete row.
    # ---------------------------------------------------------------------

    stationary_parameters <-
      stationary_result$
        validation$
        parameters


    pit_raw <- try(
      pLL3(
        y,
        mu =
          stationary_parameters$mu,

        sigma =
          stationary_parameters$sigma,

        nu =
          stationary_parameters$nu
      ),
      silent = TRUE
    )


    if (
      inherits(pit_raw, "try-error") ||
        length(pit_raw) != length(y) ||
        any(!is.finite(pit_raw)) ||
        any(pit_raw < 0) ||
        any(pit_raw > 1)
    ) {
      result$error_stage <-
        "PIT residual calculation"

      result$error_message <-
        "The fitted stationary model produced invalid PIT probabilities."

      rows[[row_index]] <- result

      next
    }


    result$pit_minimum <-
      min(pit_raw)


    result$pit_maximum <-
      max(pit_raw)


    # Permit isolated floating-point probabilities at exactly zero or one,
    # but record every clipped value.

    pit_epsilon <- 1e-12


    pit <- pmin(
      pmax(
        pit_raw,
        pit_epsilon
      ),
      1 - pit_epsilon
    )


    result$pit_clipped_count <-
      sum(
        pit != pit_raw
      )


    normalized_residuals <- stats::qnorm(
      pit
    )


    residual_standard_deviation <- stats::sd(
      normalized_residuals
    )


    if (
      any(!is.finite(normalized_residuals)) ||
        !is.finite(residual_standard_deviation) ||
        residual_standard_deviation <=
          sqrt(.Machine$double.eps)
    ) {
      result$error_stage <-
        "residual validation"

      result$error_message <-
        paste(
          "The normalized residuals were non-finite",
          "or numerically degenerate."
        )

      rows[[row_index]] <- result

      next
    }


    result$residual_mean <-
      mean(
        normalized_residuals
      )


    result$residual_sd <-
      residual_standard_deviation


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
    # AIC and BIC
    # ---------------------------------------------------------------------

    candidate_models <- list(
      stationary =
        stationary_result$fit,

      sigma =
        sigma_result$fit,

      nu =
        nu_result$fit,

      joint =
        joint_result$fit
    )


    AIC_values <- try(
      vapply(
        candidate_models,
        stats::AIC,
        numeric(1)
      ),
      silent = TRUE
    )


    deviances <- try(
      vapply(
        candidate_models,
        stats::deviance,
        numeric(1)
      ),
      silent = TRUE
    )


    fitted_degrees_of_freedom <- try(
      vapply(
        candidate_models,
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
        inherits(
          fitted_degrees_of_freedom,
          "try-error"
        )
    ) {
      result$residual_mean <- NA_real_
      result$residual_sd <- NA_real_
      result$residual_tail_rate_95 <- NA_real_

      result$error_stage <-
        "information criteria"

      result$error_message <-
        "AIC or BIC calculation produced an error."

      rows[[row_index]] <- result

      next
    }


    BIC_values <-
      deviances +
      log(current_n) *
      fitted_degrees_of_freedom


    if (
      any(!is.finite(AIC_values)) ||
        any(!is.finite(BIC_values))
    ) {
      result$residual_mean <- NA_real_
      result$residual_sd <- NA_real_
      result$residual_tail_rate_95 <- NA_real_

      result$error_stage <-
        "information criteria"

      result$error_message <-
        "At least one AIC or BIC value was non-finite."

      rows[[row_index]] <- result

      next
    }


    result$AIC_stationary <-
      unname(
        AIC_values[
          "stationary"
        ]
      )

    result$AIC_sigma <-
      unname(
        AIC_values[
          "sigma"
        ]
      )

    result$AIC_nu <-
      unname(
        AIC_values[
          "nu"
        ]
      )

    result$AIC_joint <-
      unname(
        AIC_values[
          "joint"
        ]
      )


    result$BIC_stationary <-
      unname(
        BIC_values[
          "stationary"
        ]
      )

    result$BIC_sigma <-
      unname(
        BIC_values[
          "sigma"
        ]
      )

    result$BIC_nu <-
      unname(
        BIC_values[
          "nu"
        ]
      )

    result$BIC_joint <-
      unname(
        BIC_values[
          "joint"
        ]
      )


    result$AIC_selected <-
      names(
        which.min(
          AIC_values
        )
      )


    result$BIC_selected <-
      names(
        which.min(
          BIC_values
        )
      )


    result$AIC_false_nonstationary <-
      !identical(
        result$AIC_selected,
        "stationary"
      )


    result$BIC_false_nonstationary <-
      !identical(
        result$BIC_selected,
        "stationary"
      )


    # The replication is complete only after all fitting, residual, and
    # information-criterion checks have passed.

    result$complete <- TRUE


    result$error_stage <-
      NA_character_


    result$error_message <-
      NA_character_


    rows[[row_index]] <- result
  }


  message(
    "Completed n = ",
    current_n
  )
}


# =========================================================================
# 10. Combine raw results
# =========================================================================

simulation_results <- do.call(
  rbind,
  rows
)


rownames(
  simulation_results
) <- NULL


# Enforce bookkeeping invariants.

incomplete_rows <-
  !simulation_results$complete


stopifnot(
  all(
    is.na(
      simulation_results$
        AIC_selected[
          incomplete_rows
        ]
    )
  ),

  all(
    is.na(
      simulation_results$
        BIC_selected[
          incomplete_rows
        ]
    )
  ),

  all(
    is.na(
      simulation_results$
        AIC_false_nonstationary[
          incomplete_rows
        ]
    )
  ),

  all(
    is.na(
      simulation_results$
        BIC_false_nonstationary[
          incomplete_rows
        ]
    )
  ),

  all(
    is.na(
      simulation_results$
        residual_mean[
          incomplete_rows
        ]
    )
  ),

  all(
    is.na(
      simulation_results$
        residual_sd[
          incomplete_rows
        ]
    )
  ),

  all(
    is.na(
      simulation_results$
        residual_tail_rate_95[
          incomplete_rows
        ]
    )
  )
)


# =========================================================================
# 11. Summarize by sample size
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


    complete_runs <- all_runs[
      all_runs$complete %in% TRUE,
      ,
      drop = FALSE
    ]


    attempted_count <-
      nrow(all_runs)


    complete_count <-
      nrow(complete_runs)


    data.frame(
      n =
        current_n,

      attempted =
        attempted_count,

      complete =
        complete_count,

      completion_rate =
        complete_count /
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

      joint_success_rate =
        safe_rate(
          all_runs$
            joint_fit_success
        ),


      AIC_false_nonstationary_rate =
        safe_rate(
          complete_runs$
            AIC_false_nonstationary
        ),

      BIC_false_nonstationary_rate =
        safe_rate(
          complete_runs$
            BIC_false_nonstationary
        ),


      AIC_stationary_selection_rate =
        safe_rate(
          complete_runs$
            AIC_selected ==
            "stationary"
        ),

      BIC_stationary_selection_rate =
        safe_rate(
          complete_runs$
            BIC_selected ==
            "stationary"
        ),


      mean_residual_mean =
        safe_mean(
          complete_runs$
            residual_mean
        ),

      mean_residual_sd =
        safe_mean(
          complete_runs$
            residual_sd
        ),

      mean_residual_tail_rate_95 =
        safe_mean(
          complete_runs$
            residual_tail_rate_95
        ),


      mean_PIT_clipped_count =
        safe_mean(
          complete_runs$
            pit_clipped_count
        ),


      stationary_retry_rate =
        safe_rate(
          all_runs$
            stationary_retried
        ),

      sigma_retry_rate =
        safe_rate(
          all_runs$
            sigma_retried
        ),

      nu_retry_rate =
        safe_rate(
          all_runs$
            nu_retried
        ),

      joint_retry_rate =
        safe_rate(
          all_runs$
            joint_retried
        ),


      stringsAsFactors =
        FALSE
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
# 12. Failure and warning tables
# =========================================================================

failure_results <- simulation_results[
  !simulation_results$complete,
  c(
    "n",
    "seed",

    "stationary_fit_success",
    "sigma_fit_success",
    "nu_fit_success",
    "joint_fit_success",

    "stationary_retried",
    "sigma_retried",
    "nu_retried",
    "joint_retried",

    "error_stage",
    "error_message",

    "stationary_warning",
    "sigma_warning",
    "nu_warning",
    "joint_warning"
  ),
  drop = FALSE
]


warning_results <- simulation_results[
  (
    !is.na(
      simulation_results$
        stationary_warning
    ) |
      !is.na(
        simulation_results$
          sigma_warning
      ) |
      !is.na(
        simulation_results$
          nu_warning
      ) |
      !is.na(
        simulation_results$
          joint_warning
      )
  ),
  c(
    "n",
    "seed",
    "complete",

    "stationary_warning",
    "sigma_warning",
    "nu_warning",
    "joint_warning"
  ),
  drop = FALSE
]


# =========================================================================
# 13. Save results
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
    "stationary_null_raw.csv"
  ),
  row.names = FALSE
)


utils::write.csv(
  simulation_summary,
  file.path(
    results_directory,
    "stationary_null_summary.csv"
  ),
  row.names = FALSE
)


utils::write.csv(
  failure_results,
  file.path(
    results_directory,
    "stationary_null_failures.csv"
  ),
  row.names = FALSE
)


utils::write.csv(
  warning_results,
  file.path(
    results_directory,
    "stationary_null_warnings.csv"
  ),
  row.names = FALSE
)


# =========================================================================
# 14. Print results
# =========================================================================

print(
  simulation_summary,
  digits = 10,
  row.names = FALSE
)


if (nrow(failure_results) > 0L) {
  message(
    "\nIncomplete stationary-null simulations: ",
    nrow(failure_results)
  )

  print(
    failure_results,
    row.names = FALSE
  )
} else {
  message(
    "\nAll stationary-null simulations completed successfully."
  )
}


if (nrow(warning_results) > 0L) {
  message(
    "\nSimulations with recorded GAMLSS warnings: ",
    nrow(warning_results),
    ". See stationary_null_warnings.csv."
  )
} else {
  message(
    "\nNo GAMLSS warnings were recorded."
  )
}


message(
  "\nStationary-null results written to:\n",
  results_directory
)