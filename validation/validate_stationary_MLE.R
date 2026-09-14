# validation/validate_stationary_MLE.R
#
# Stationary LL3 maximum-likelihood validation.
#
# Compares three implementations:
#
#   1. Classic GAMLSS LL3 family
#   2. Direct numerical optimization using dLL3()
#   3. Independent numerical optimization using FAdist::dllog3()
#
# Validation focuses on:
#
#   * convergence;
#   * support validity;
#   * maximized log-likelihood agreement;
#   * parameter agreement on scale-normalized units;
#   * near-zero likelihood scores at the GAMLSS solution.


# =========================================================================
# 1. Load project code
# =========================================================================

source(
  file.path(
    "R",
    "LL3_distribution.R"
  )
)

source(
  file.path(
    "R",
    "LL3_gamlss_family.R"
  )
)

source(
  file.path(
    "R",
    "LL3_diagnostics.R"
  )
)

source(
  file.path(
    "R",
    "LL3_fitting.R"
  )
)


# =========================================================================
# 2. Required packages
# =========================================================================

if (!requireNamespace("gamlss", quietly = TRUE)) {
  stop(
    "Install the gamlss package before running this validation."
  )
}

if (!requireNamespace("FAdist", quietly = TRUE)) {
  stop(
    "Install the FAdist package before running this validation."
  )
}


# =========================================================================
# 3. Simulation settings
# =========================================================================

set.seed(42)

n <- 1000L

parameter_order <- c(
  "mu",
  "sigma",
  "nu"
)

truth <- c(
  mu = -50,
  sigma = 70,
  nu = 1.8
)

y <- rLL3(
  n = n,
  mu = unname(truth["mu"]),
  sigma = unname(truth["sigma"]),
  nu = unname(truth["nu"])
)

validation_data <- data.frame(
  y = y
)

lower <- LL3_support_lower(y)


# =========================================================================
# 4. Fit the stationary LL3 model using classic GAMLSS
# =========================================================================

strict_control <- gamlss::gamlss.control(
  c.crit = 1e-8,
  n.cyc = 2000,
  mu.step = 0.10,
  sigma.step = 0.10,
  nu.step = 0.10,
  gd.tol = Inf,
  autostep = TRUE,
  trace = FALSE
)

gamlss_fit <- fit_LL3_gamlss(
  mu.formula = y ~ 1,
  sigma.formula = ~1,
  nu.formula = ~1,
  data = validation_data,
  lower = lower,
  information = "opg",
  allow_nonstationary_mu = FALSE,
  control = strict_control
)

gamlss_parameter_rows <- extract_LL3_parameters(
  gamlss_fit
)

gamlss_parameters <- c(
  mu = unname(
    gamlss_parameter_rows$mu[1]
  ),

  sigma = unname(
    gamlss_parameter_rows$sigma[1]
  ),

  nu = unname(
    gamlss_parameter_rows$nu[1]
  )
)

stopifnot(
  identical(
    names(gamlss_parameters),
    parameter_order
  )
)

stopifnot(
  all(
    is.finite(gamlss_parameters)
  )
)

gamlss_logLik <- sum(
  dLL3(
    y,
    mu = gamlss_parameters["mu"],
    sigma = gamlss_parameters["sigma"],
    nu = gamlss_parameters["nu"],
    log = TRUE
  )
)

gamlss_checks <- check_LL3_fit(
  gamlss_fit,
  y = y
)


# =========================================================================
# 5. Direct numerical maximum likelihood using dLL3()
# =========================================================================

X_intercept <- matrix(
  1,
  nrow = n,
  ncol = 1
)

colnames(X_intercept) <- "(Intercept)"

gamlss_mu_coefficient <- stats::coef(
  gamlss_fit,
  what = "mu"
)

gamlss_sigma_coefficient <- stats::coef(
  gamlss_fit,
  what = "sigma"
)

gamlss_nu_coefficient <- stats::coef(
  gamlss_fit,
  what = "nu"
)

direct_start <- c(
  unname(
    gamlss_mu_coefficient[1]
  ),

  unname(
    gamlss_sigma_coefficient[1]
  ),

  unname(
    gamlss_nu_coefficient[1]
  )
)

stopifnot(
  length(direct_start) == 3L,
  all(is.finite(direct_start))
)

direct_result <- fit_LL3_direct_linear(
  y = y,
  X_sigma = X_intercept,
  X_nu = X_intercept,
  lower = lower,
  start = direct_start,
  n_starts = 20,
  seed = 42,
  maxit = 10000
)

direct_parameters <- c(
  mu = unname(
    direct_result$mu
  ),

  sigma = unname(
    exp(
      direct_result$beta_sigma[1]
    )
  ),

  nu = unname(
    exp(
      direct_result$beta_nu[1]
    )
  )
)

stopifnot(
  identical(
    names(direct_parameters),
    parameter_order
  )
)

stopifnot(
  all(
    is.finite(direct_parameters)
  )
)

stopifnot(
  direct_result$convergence == 0
)

direct_logLik <- unname(
  direct_result$logLik
)


# =========================================================================
# 6. Independent maximum likelihood using FAdist
# =========================================================================

FAdist_objective <- function(theta) {
  if (
    length(theta) != 3L ||
      any(!is.finite(theta))
  ) {
    return(1e100)
  }

  mu <- lower - exp(theta[1])
  sigma <- exp(theta[2])
  nu <- exp(theta[3])

  if (
    !is.finite(mu) ||
      mu >= lower ||
      !is.finite(sigma) ||
      sigma <= 0 ||
      !is.finite(nu) ||
      nu <= 0
  ) {
    return(1e100)
  }

  log_density <- suppressWarnings(
    FAdist::dllog3(
      y,
      shape = 1 / nu,
      scale = log(sigma),
      thres = mu,
      log = TRUE
    )
  )

  if (
    length(log_density) != length(y) ||
      any(!is.finite(log_density))
  ) {
    return(1e100)
  }

  -sum(log_density)
}


# Generic starting values independent of the direct optimizer

initial_values <- .LL3_start_values(y)

initial_mu <- min(
  unname(
    initial_values["mu"]
  ),
  lower - sqrt(.Machine$double.eps)
)

generic_start <- c(
  log(
    lower - initial_mu
  ),

  log(
    unname(
      initial_values["sigma"]
    )
  ),

  log(
    unname(
      initial_values["nu"]
    )
  )
)

stopifnot(
  length(generic_start) == 3L,
  all(is.finite(generic_start))
)


# Use several independent starting points.
#
# The first is generic, the second uses the direct optimizer only as a
# diagnostic start, and the remainder are random perturbations.

set.seed(43)

number_FAdist_starts <- 20L

FAdist_start_list <- vector(
  "list",
  number_FAdist_starts
)

FAdist_start_list[[1]] <- generic_start

FAdist_start_list[[2]] <- unname(
  direct_result$coefficients
)

if (number_FAdist_starts > 2L) {
  for (i in 3:number_FAdist_starts) {
    FAdist_start_list[[i]] <-
      generic_start +
      stats::rnorm(
        3,
        mean = 0,
        sd = 0.25
      )
  }
}


FAdist_fits <- lapply(
  FAdist_start_list,
  function(current_start) {
    try(
      stats::optim(
        par = current_start,
        fn = FAdist_objective,
        method = "BFGS",
        hessian = TRUE,
        control = list(
          maxit = 10000,
          reltol = 1e-12
        )
      ),
      silent = TRUE
    )
  }
)

FAdist_valid <- vapply(
  FAdist_fits,
  function(current_fit) {
    !inherits(
      current_fit,
      "try-error"
    ) &&
      current_fit$convergence == 0 &&
      is.finite(current_fit$value) &&
      all(is.finite(current_fit$par))
  },
  logical(1)
)

if (!any(FAdist_valid)) {
  stop(
    "No FAdist-based numerical optimization converged."
  )
}

FAdist_candidates <- FAdist_fits[
  FAdist_valid
]

FAdist_candidate_values <- vapply(
  FAdist_candidates,
  function(current_fit) {
    current_fit$value
  },
  numeric(1)
)

FAdist_best <- FAdist_candidates[[
  which.min(
    FAdist_candidate_values
  )
]]

FAdist_parameters <- c(
  mu = unname(
    lower -
      exp(
        FAdist_best$par[1]
      )
  ),

  sigma = unname(
    exp(
      FAdist_best$par[2]
    )
  ),

  nu = unname(
    exp(
      FAdist_best$par[3]
    )
  )
)

stopifnot(
  identical(
    names(FAdist_parameters),
    parameter_order
  )
)

stopifnot(
  all(
    is.finite(FAdist_parameters)
  )
)

stopifnot(
  FAdist_best$convergence == 0
)

FAdist_logLik <- unname(
  -FAdist_best$value
)


# =========================================================================
# 7. Verify all parameter names before constructing comparisons
# =========================================================================

stopifnot(
  all(
    parameter_order %in%
      names(truth)
  )
)

stopifnot(
  all(
    parameter_order %in%
      names(gamlss_parameters)
  )
)

stopifnot(
  all(
    parameter_order %in%
      names(direct_parameters)
  )
)

stopifnot(
  all(
    parameter_order %in%
      names(FAdist_parameters)
  )
)


# =========================================================================
# 8. Parameter comparison
# =========================================================================

parameter_comparison <- data.frame(
  parameter = parameter_order,

  truth = unname(
    truth[
      parameter_order
    ]
  ),

  GAMLSS = unname(
    gamlss_parameters[
      parameter_order
    ]
  ),

  direct_MLE = unname(
    direct_parameters[
      parameter_order
    ]
  ),

  FAdist_MLE = unname(
    FAdist_parameters[
      parameter_order
    ]
  ),

  stringsAsFactors = FALSE,
  row.names = NULL
)

parameter_comparison$GAMLSS_direct_difference <-
  abs(
    parameter_comparison$GAMLSS -
      parameter_comparison$direct_MLE
  )

parameter_comparison$GAMLSS_FAdist_difference <-
  abs(
    parameter_comparison$GAMLSS -
      parameter_comparison$FAdist_MLE
  )

parameter_comparison$direct_FAdist_difference <-
  abs(
    parameter_comparison$direct_MLE -
      parameter_comparison$FAdist_MLE
  )


# Normalize differences because mu, sigma, and nu use different units.

parameter_scale <- c(
  mu = max(
    stats::IQR(y),
    1
  ),

  sigma = max(
    abs(
      direct_parameters["sigma"]
    ),
    1
  ),

  nu = max(
    abs(
      direct_parameters["nu"]
    ),
    1
  )
)

parameter_comparison$scale_reference <- unname(
  parameter_scale[
    parameter_comparison$parameter
  ]
)

parameter_comparison$GAMLSS_direct_scaled_difference <-
  parameter_comparison$GAMLSS_direct_difference /
  parameter_comparison$scale_reference

parameter_comparison$GAMLSS_FAdist_scaled_difference <-
  parameter_comparison$GAMLSS_FAdist_difference /
  parameter_comparison$scale_reference

parameter_comparison$direct_FAdist_scaled_difference <-
  parameter_comparison$direct_FAdist_difference /
  parameter_comparison$scale_reference


comparison_numeric_columns <- c(
  "truth",
  "GAMLSS",
  "direct_MLE",
  "FAdist_MLE",
  "GAMLSS_direct_difference",
  "GAMLSS_FAdist_difference",
  "direct_FAdist_difference",
  "scale_reference",
  "GAMLSS_direct_scaled_difference",
  "GAMLSS_FAdist_scaled_difference",
  "direct_FAdist_scaled_difference"
)

stopifnot(
  all(
    is.finite(
      as.matrix(
        parameter_comparison[
          comparison_numeric_columns
        ]
      )
    )
  )
)


# =========================================================================
# 9. Log-likelihood comparison
# =========================================================================

likelihood_comparison <- data.frame(
  method = c(
    "Classic GAMLSS",
    "Direct dLL3 optimizer",
    "Independent FAdist optimizer"
  ),

  logLik = c(
    gamlss_logLik,
    direct_logLik,
    FAdist_logLik
  ),

  stringsAsFactors = FALSE,
  row.names = NULL
)

stopifnot(
  all(
    is.finite(
      likelihood_comparison$logLik
    )
  )
)

gamlss_direct_logLik_difference <- abs(
  gamlss_logLik -
    direct_logLik
)

gamlss_FAdist_logLik_difference <- abs(
  gamlss_logLik -
    FAdist_logLik
)

direct_FAdist_logLik_difference <- abs(
  direct_logLik -
    FAdist_logLik
)


# =========================================================================
# 10. Score check at the GAMLSS solution
# =========================================================================

natural_scores <- .LL3_score(
  y,
  mu = gamlss_parameters["mu"],
  sigma = gamlss_parameters["sigma"],
  nu = gamlss_parameters["nu"]
)

# Chain-rule derivatives:
#
#   mu = lower - exp(eta_mu)
#   dmu / deta_mu = -(lower - mu)
#
#   sigma = exp(eta_sigma)
#   dsigma / deta_sigma = sigma
#
#   nu = exp(eta_nu)
#   dnu / deta_nu = nu

link_scores <- c(
  eta_mu = sum(
    natural_scores[, "mu"] *
      (
        -(
          lower -
            gamlss_parameters["mu"]
        )
      )
  ),

  eta_sigma = sum(
    natural_scores[, "sigma"] *
      gamlss_parameters["sigma"]
  ),

  eta_nu = sum(
    natural_scores[, "nu"] *
      gamlss_parameters["nu"]
  )
)

score_comparison <- data.frame(
  parameter = names(link_scores),

  total_link_score = unname(
    link_scores
  ),

  mean_link_score = unname(
    link_scores / n
  ),

  stringsAsFactors = FALSE,
  row.names = NULL
)

stopifnot(
  all(
    is.finite(
      score_comparison$total_link_score
    )
  )
)


# =========================================================================
# 11. Print diagnostics before applying validation assertions
# =========================================================================

cat(
  "\nStationary LL3 parameter comparison\n"
)

print(
  parameter_comparison,
  digits = 14,
  row.names = FALSE
)

cat(
  "\nStationary LL3 likelihood comparison\n"
)

print(
  likelihood_comparison,
  digits = 16,
  row.names = FALSE
)

cat(
  "\nAbsolute log-likelihood differences\n"
)

print(
  c(
    GAMLSS_vs_direct =
      gamlss_direct_logLik_difference,

    GAMLSS_vs_FAdist =
      gamlss_FAdist_logLik_difference,

    direct_vs_FAdist =
      direct_FAdist_logLik_difference
  ),
  digits = 16
)

cat(
  "\nGAMLSS link-scale score check\n"
)

print(
  score_comparison,
  digits = 14,
  row.names = FALSE
)

cat(
  "\nGAMLSS validity checks\n"
)

print(
  gamlss_checks
)

cat(
  "\nSupport boundary diagnostics\n"
)

print(
  LL3_boundary_diagnostic(
    gamlss_fit,
    y = y
  ),
  digits = 14,
  row.names = FALSE
)


# =========================================================================
# 12. Validation assertions
# =========================================================================

# GAMLSS convergence and parameter validity

stopifnot(
  isTRUE(
    gamlss_checks["converged"]
  )
)

stopifnot(
  isTRUE(
    gamlss_checks["finite_parameters"]
  )
)

stopifnot(
  isTRUE(
    gamlss_checks["support_ok"]
  )
)

stopifnot(
  isTRUE(
    gamlss_checks["sigma_positive"]
  )
)

stopifnot(
  isTRUE(
    gamlss_checks["nu_positive"]
  )
)

stopifnot(
  isTRUE(
    gamlss_checks["lower_below_data"]
  )
)

stopifnot(
  isTRUE(
    gamlss_checks["mu_below_lower"]
  )
)


# The independently implemented direct dLL3 and FAdist likelihoods
# should converge to effectively the same maximum.

stopifnot(
  direct_FAdist_logLik_difference <
    1e-6
)

stopifnot(
  max(
    parameter_comparison$
      direct_FAdist_scaled_difference
  ) <
    1e-5
)


# GAMLSS should reach effectively the same likelihood maximum.

stopifnot(
  gamlss_direct_logLik_difference <
    1e-4
)

stopifnot(
  gamlss_FAdist_logLik_difference <
    1e-4
)


# Natural-parameter agreement is evaluated using scale-normalized
# differences rather than a single raw-unit tolerance.

stopifnot(
  max(
    parameter_comparison$
      GAMLSS_direct_scaled_difference
  ) <
    1e-3
)

stopifnot(
  max(
    parameter_comparison$
      GAMLSS_FAdist_scaled_difference
  ) <
    1e-3
)


# At a converged likelihood maximum, the average link-scale score should
# be close to zero.

stopifnot(
  max(
    abs(
      score_comparison$mean_link_score
    )
  ) <
    1e-5
)


message(
  "Stationary MLE validation passed."
)

output_dir <- file.path("validation", "results", "core_numerical")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  parameter_comparison,
  file.path(output_dir, "stationary_MLE_parameter_agreement.csv"),
  row.names = FALSE
)
utils::write.csv(
  likelihood_comparison,
  file.path(output_dir, "stationary_MLE_logLik_agreement.csv"),
  row.names = FALSE
)
utils::write.csv(
  score_comparison,
  file.path(output_dir, "stationary_MLE_score_checks.csv"),
  row.names = FALSE
)
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
