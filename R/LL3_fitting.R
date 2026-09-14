# Fitting utilities for the LL3 GAMLSS family.

LL3_support_lower <- function(y,
                              relative_margin = 1e-6,
                              absolute_margin = sqrt(.Machine$double.eps)) {
  y <- y[is.finite(y)]
  if (length(y) < 4L) stop("At least four finite observations are required.")

  spread <- stats::IQR(y)
  if (!is.finite(spread) || spread <= 0) spread <- stats::sd(y)
  if (!is.finite(spread) || spread <= 0) spread <- max(abs(y), 1)

  numerical_margin <- 100 * .Machine$double.eps * max(abs(y), 1)
  margin <- max(absolute_margin, relative_margin * spread, numerical_margin)
  min(y) - margin
}

.LL3_formula_has_terms <- function(formula) {
  length(attr(stats::terms(formula), "term.labels")) > 0L
}

.LL3_validate_model_data <- function(mu.formula, sigma.formula, nu.formula, data) {
  if (!is.data.frame(data)) stop("data must be a data.frame.")

  variables <- unique(c(
    all.vars(mu.formula), all.vars(sigma.formula), all.vars(nu.formula)
  ))
  missing_variables <- setdiff(variables, names(data))
  if (length(missing_variables)) {
    stop("Variables missing from data: ", paste(missing_variables, collapse = ", "))
  }
  if (any(!stats::complete.cases(data[, variables, drop = FALSE]))) {
    stop("Remove or impute missing values in all variables used by the model.")
  }

  mf <- stats::model.frame(mu.formula, data = data, na.action = stats::na.fail)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || any(!is.finite(y))) {
    stop("The response must be finite and numeric.")
  }
  if (length(y) != nrow(data)) {
    stop("The model response does not align with the rows of data.")
  }
  as.numeric(y)
}

LL3_default_control <- function(n_cycles = 700,
                                mu_step = 0.05,
                                sigma_step = 0.05,
                                nu_step = 0.05,
                                trace = FALSE) {
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Install gamlss before fitting LL3 models.")
  }
  gamlss::gamlss.control(
    n.cyc = n_cycles,
    mu.step = mu_step,
    sigma.step = sigma_step,
    nu.step = nu_step,
    gd.tol = Inf,
    autostep = TRUE,
    trace = trace
  )
}

fit_LL3_gamlss <- function(mu.formula,
                           sigma.formula = ~1,
                           nu.formula = ~1,
                           data,
                           lower = NULL,
                           information = "opg",
                           allow_nonstationary_mu = FALSE,
                           mu.start = NULL,
                           sigma.start = NULL,
                           nu.start = NULL,
                           control = NULL,
                           boundary_action = c("warning", "error", "ignore")) {
  boundary_action <- match.arg(boundary_action)
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Install gamlss before fitting LL3 models.")
  }
  if (!allow_nonstationary_mu && .LL3_formula_has_terms(mu.formula)) {
    stop(
      "The validated nonstationary scope keeps mu constant. ",
      "Set allow_nonstationary_mu = TRUE only for experimental work."
    )
  }

  y <- .LL3_validate_model_data(mu.formula, sigma.formula, nu.formula, data)
  n <- length(y)

  if (is.null(lower)) lower <- LL3_support_lower(y)
  if (!is.finite(lower) || lower >= min(y)) {
    stop("lower must be finite and strictly below every response value.")
  }

  start <- .LL3_start_values(y)
  margin <- min(y) - lower
  safe_mu <- min(unname(start["mu"]), lower - max(margin, sqrt(.Machine$double.eps)))

  if (is.null(mu.start)) mu.start <- rep(safe_mu, n)
  if (is.null(sigma.start)) sigma.start <- rep(unname(start["sigma"]), n)
  if (is.null(nu.start)) nu.start <- rep(unname(start["nu"]), n)

  mu.start <- rep_len(mu.start, n)
  sigma.start <- rep_len(sigma.start, n)
  nu.start <- rep_len(nu.start, n)

  if (any(!is.finite(mu.start) | mu.start >= lower)) {
    stop("Every mu.start value must be finite and strictly below lower.")
  }
  if (any(!is.finite(sigma.start) | sigma.start <= 0)) {
    stop("Every sigma.start value must be finite and positive.")
  }
  if (any(!is.finite(nu.start) | nu.start <= 0)) {
    stop("Every nu.start value must be finite and positive.")
  }
  if (is.null(control)) control <- LL3_default_control()

  family_object <- LL3(
    lower = lower,
    information = information
  )

  fit <- gamlss::gamlss(
    formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    family = family_object,
    data = data,
    mu.start = mu.start,
    sigma.start = sigma.start,
    nu.start = nu.start,
    control = control
  )

  attr(fit, "LL3_lower") <- lower
  attr(fit, "LL3_validated_scope") <- !allow_nonstationary_mu
  attr(fit, "LL3_family_object") <- family_object

  if (!identical(boundary_action, "ignore")) {
    LL3_assert_inference_ready(
      fit,
      y = y,
      action = if (identical(boundary_action, "error")) "error" else "warning"
    )
  }

  fit
}

fit_LL3_stationary <- function(y,
                               trace = FALSE,
                               n_cycles = 700,
                               information = "opg") {
  y <- y[is.finite(y)]
  data <- data.frame(y = y)
  lower <- LL3_support_lower(y)

  fit <- fit_LL3_gamlss(
    y ~ 1,
    sigma.formula = ~1,
    nu.formula = ~1,
    data = data,
    lower = lower,
    information = information,
    control = LL3_default_control(
      n_cycles = n_cycles,
      mu_step = 0.10,
      sigma_step = 0.10,
      nu_step = 0.10,
      trace = trace
    )
  )

  pars <- extract_LL3_parameters(fit)
  parameters <- c(mu = pars$mu[1], sigma = pars$sigma[1], nu = pars$nu[1])
  ll <- sum(dLL3(y, parameters["mu"], parameters["sigma"], parameters["nu"], log = TRUE))

  list(
    fit = fit,
    parameters = parameters,
    logLik = ll,
    support_ok = all(y > parameters["mu"]),
    lower = lower
  )
}

fit_LL3_gamlss_stationary <- fit_LL3_stationary
fit_LL3safe_gamlss_stationary <- fit_LL3_stationary

fit_LL3_candidate_models <- function(data,
                                     response,
                                     covariate,
                                     lower = NULL,
                                     information = "opg",
                                     control = NULL) {
  if (!all(c(response, covariate) %in% names(data))) {
    stop("response and covariate must name columns in data.")
  }
  mu_formula <- stats::as.formula(paste(response, "~ 1"))
  ns_formula <- stats::reformulate(covariate)
  y <- data[[response]]
  if (is.null(lower)) lower <- LL3_support_lower(y)
  if (is.null(control)) control <- LL3_default_control()

  M0 <- fit_LL3_gamlss(
    mu_formula, ~1, ~1, data = data, lower = lower,
    information = information, control = control
  )
  p0 <- extract_LL3_parameters(M0)

  M_sigma <- fit_LL3_gamlss(
    mu_formula, ns_formula, ~1, data = data, lower = lower,
    information = information,
    mu.start = p0$mu, sigma.start = p0$sigma, nu.start = p0$nu,
    control = control
  )

  M_nu <- fit_LL3_gamlss(
    mu_formula, ~1, ns_formula, data = data, lower = lower,
    information = information,
    mu.start = p0$mu, sigma.start = p0$sigma, nu.start = p0$nu,
    control = control
  )

  ps <- extract_LL3_parameters(M_sigma)
  pn <- extract_LL3_parameters(M_nu)

  M_joint <- fit_LL3_gamlss(
    mu_formula, ns_formula, ns_formula, data = data, lower = lower,
    information = information,
    mu.start = p0$mu, sigma.start = ps$sigma, nu.start = pn$nu,
    control = control
  )

  list(
    stationary = M0,
    sigma = M_sigma,
    nu = M_nu,
    joint = M_joint,
    lower = lower
  )
}

# Independent numerical maximum likelihood for constant mu and linear
# predictors in sigma and nu. No analytical gradient is supplied to optim().
fit_LL3_direct_linear <- function(y,
                                  X_sigma,
                                  X_nu,
                                  lower = NULL,
                                  start = NULL,
                                  n_starts = 5,
                                  seed = 1,
                                  maxit = 5000) {
  y <- as.numeric(y)
  X_sigma <- as.matrix(X_sigma)
  X_nu <- as.matrix(X_nu)
  n <- length(y)

  if (nrow(X_sigma) != n || nrow(X_nu) != n) {
    stop("Design matrices must have one row per response value.")
  }
  if (any(!is.finite(y)) || any(!is.finite(X_sigma)) || any(!is.finite(X_nu))) {
    stop("Response and design matrices must be finite.")
  }
  if (is.null(lower)) lower <- LL3_support_lower(y)
  if (lower >= min(y)) stop("lower must be strictly below min(y).")

  p_sigma <- ncol(X_sigma)
  p_nu <- ncol(X_nu)
  starts0 <- .LL3_start_values(y)

  if (is.null(start)) {
    beta_sigma0 <- rep(0, p_sigma)
    beta_nu0 <- rep(0, p_nu)
    beta_sigma0[1] <- log(starts0["sigma"])
    beta_nu0[1] <- log(starts0["nu"])
    start <- c(
      eta_mu = log(lower - min(starts0["mu"], lower - sqrt(.Machine$double.eps))),
      beta_sigma0,
      beta_nu0
    )
  }

  expected_length <- 1L + p_sigma + p_nu
  if (length(start) != expected_length || any(!is.finite(start))) {
    stop("start has the wrong length or contains non-finite values.")
  }

  unpack <- function(theta) {
    beta_sigma <- theta[1L + seq_len(p_sigma)]
    beta_nu <- theta[1L + p_sigma + seq_len(p_nu)]
    eta_sigma <- drop(X_sigma %*% beta_sigma)
    eta_nu <- drop(X_nu %*% beta_nu)
    list(
      mu = lower - exp(theta[1]),
      beta_sigma = beta_sigma,
      beta_nu = beta_nu,
      eta_sigma = eta_sigma,
      eta_nu = eta_nu,
      sigma = exp(eta_sigma),
      nu = exp(eta_nu)
    )
  }

  objective <- function(theta) {
    if (any(!is.finite(theta))) return(1e100)
    pars <- unpack(theta)
    if (!is.finite(pars$mu) || pars$mu >= lower ||
        any(!is.finite(pars$sigma)) || any(pars$sigma <= 0) ||
        any(!is.finite(pars$nu)) || any(pars$nu <= 0)) {
      return(1e100)
    }
    ld <- dLL3(y, pars$mu, pars$sigma, pars$nu, log = TRUE)
    if (any(!is.finite(ld))) return(1e100)
    -sum(ld)
  }

  set.seed(seed)
  start_list <- vector("list", max(1L, n_starts))
  start_list[[1]] <- start
  if (length(start_list) > 1L) {
    for (i in 2:length(start_list)) {
      start_list[[i]] <- start + stats::rnorm(length(start), sd = 0.10)
    }
  }

  fits <- lapply(start_list, function(s) {
    try(stats::optim(
      par = s, fn = objective, method = "BFGS", hessian = TRUE,
      control = list(maxit = maxit, reltol = 1e-10)
    ), silent = TRUE)
  })
  valid <- vapply(fits, function(z) {
    !inherits(z, "try-error") && z$convergence == 0 && is.finite(z$value)
  }, logical(1))
  if (!any(valid)) stop("No direct numerical optimization converged.")

  candidates <- fits[valid]
  values <- vapply(candidates, function(z) z$value, numeric(1))
  best <- candidates[[which.min(values)]]
  pars <- unpack(best$par)

  hessian_symmetric <- (best$hessian + t(best$hessian)) / 2
  hessian_eigenvalues <- try(
    eigen(hessian_symmetric, symmetric = TRUE, only.values = TRUE)$values,
    silent = TRUE
  )
  vcov <- NULL
  if (!inherits(hessian_eigenvalues, "try-error") &&
      all(is.finite(hessian_eigenvalues)) &&
      min(hessian_eigenvalues) > 0) {
    candidate_vcov <- try(solve(hessian_symmetric), silent = TRUE)
    if (!inherits(candidate_vcov, "try-error") && all(is.finite(candidate_vcov))) {
      vcov <- candidate_vcov
    }
  }

  list(
    coefficients = best$par,
    mu = pars$mu,
    beta_sigma = pars$beta_sigma,
    beta_nu = pars$beta_nu,
    fitted = data.frame(mu = rep(pars$mu, n), sigma = pars$sigma, nu = pars$nu),
    logLik = -best$value,
    convergence = best$convergence,
    hessian = hessian_symmetric,
    hessian_eigenvalues = if (inherits(hessian_eigenvalues, "try-error")) NULL else hessian_eigenvalues,
    vcov = vcov,
    lower = lower,
    objective = objective,
    optim = best
  )
}

compare_LL3_gamlss_direct <- function(gamlss_fit,
                                      data,
                                      mu.formula,
                                      sigma.formula = ~1,
                                      nu.formula = ~1,
                                      n_starts = 5,
                                      seed = 1) {
  if (.LL3_formula_has_terms(mu.formula)) {
    stop("This direct comparator currently validates constant-mu models only.")
  }
  y <- .LL3_validate_model_data(mu.formula, sigma.formula, nu.formula, data)
  X_sigma <- stats::model.matrix(sigma.formula, data = data)
  X_nu <- stats::model.matrix(nu.formula, data = data)
  lower <- attr(gamlss_fit, "LL3_lower")
  if (is.null(lower)) lower <- gamlss_fit$family[1] # deliberate fallback checked below
  if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower)) {
    stop("The GAMLSS fit does not contain a valid LL3 lower boundary.")
  }

  mu_coef <- stats::coef(gamlss_fit, what = "mu")
  sigma_coef <- stats::coef(gamlss_fit, what = "sigma")
  nu_coef <- stats::coef(gamlss_fit, what = "nu")
  start <- c(unname(mu_coef[1]), unname(sigma_coef), unname(nu_coef))

  direct <- fit_LL3_direct_linear(
    y, X_sigma, X_nu, lower = lower, start = start,
    n_starts = n_starts, seed = seed
  )
  gp <- extract_LL3_parameters(gamlss_fit)
  g_ll <- sum(dLL3(y, gp$mu, gp$sigma, gp$nu, log = TRUE))

  list(
    gamlss_logLik = g_ll,
    direct_logLik = direct$logLik,
    absolute_logLik_difference = abs(g_ll - direct$logLik),
    gamlss_mu = mean(gp$mu),
    direct_mu = direct$mu,
    sigma_coefficient_difference = unname(sigma_coef - direct$beta_sigma),
    nu_coefficient_difference = unname(nu_coef - direct$beta_nu),
    direct = direct
  )
}

LL3_parametric_bootstrap <- function(fit,
                                     data,
                                     response,
                                     mu.formula,
                                     sigma.formula = ~1,
                                     nu.formula = ~1,
                                     B = 499,
                                     seed = 1,
                                     control = NULL,
                                     show_progress = TRUE,
                                     boundary_action = c("error", "warning")) {
  boundary_action <- match.arg(boundary_action)
  if (!response %in% names(data)) stop("response must name a column in data.")
  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) stop("B must be a positive integer.")
  if (is.null(control)) control <- LL3_default_control(trace = FALSE)

  LL3_assert_inference_ready(
    fit,
    y = data[[response]],
    action = boundary_action
  )

  fitted_pars <- extract_LL3_parameters(fit)
  original_lower <- attr(fit, "LL3_lower")
  set.seed(seed)
  out <- vector("list", B)

  for (b in seq_len(B)) {
    if (show_progress && (b == 1L || b %% 25L == 0L || b == B)) {
      message("Bootstrap replicate ", b, " of ", B)
    }
    boot_data <- data
    boot_data[[response]] <- rLL3(
      nrow(data), fitted_pars$mu, fitted_pars$sigma, fitted_pars$nu
    )
    boot_lower <- LL3_support_lower(boot_data[[response]])

    fb <- try(fit_LL3_gamlss(
      mu.formula, sigma.formula, nu.formula,
      data = boot_data, lower = boot_lower,
      information = "opg", allow_nonstationary_mu = FALSE,
      boundary_action = "ignore",
      control = control
    ), silent = TRUE)

    if (inherits(fb, "try-error") || !isTRUE(fb$converged)) {
      out[[b]] <- NULL
      next
    }
    pb <- extract_LL3_parameters(fb)
    out[[b]] <- c(
      mu = mean(pb$mu),
      stats::setNames(stats::coef(fb, what = "sigma"),
               paste0("sigma:", names(stats::coef(fb, what = "sigma")))),
      stats::setNames(stats::coef(fb, what = "nu"),
               paste0("nu:", names(stats::coef(fb, what = "nu"))))
    )
  }

  valid <- !vapply(out, is.null, logical(1))
  if (!any(valid)) stop("No bootstrap refits converged.")
  estimates <- do.call(rbind, out[valid])
  list(
    estimates = estimates,
    convergence_rate = mean(valid),
    percentile_95 = apply(estimates, 2, stats::quantile,
                          probs = c(0.025, 0.975), na.rm = TRUE),
    B = B,
    original_lower = original_lower
  )
}


LL3_moving_block_bootstrap <- function(
    fit,
    data,
    response,
    mu.formula,
    sigma.formula = ~1,
    nu.formula = ~1,
    B = 499,
    block_length = NULL,
    seed = 1,
    control = NULL,
    show_progress = TRUE,
    boundary_action = c("error", "warning")
) {
  boundary_action <- match.arg(boundary_action)

  if (!response %in% names(data)) {
    stop("response must name a column in data.")
  }

  if (.LL3_formula_has_terms(mu.formula)) {
    stop("The validated moving-block bootstrap currently requires constant mu.")
  }

  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) {
    stop("B must be a positive integer.")
  }

  n <- nrow(data)
  if (is.null(block_length)) {
    block_length <- max(2L, as.integer(ceiling(n^(1 / 3))))
  }
  block_length <- as.integer(block_length)
  if (!is.finite(block_length) || block_length < 2L || block_length > n) {
    stop("block_length must be an integer between 2 and nrow(data).")
  }

  if (is.null(control)) {
    control <- LL3_default_control(trace = FALSE)
  }

  LL3_assert_inference_ready(
    fit,
    y = data[[response]],
    action = boundary_action
  )

  fitted_pars <- extract_LL3_parameters(fit)
  residuals <- LL3_residual_diagnostics(
    fit,
    y = data[[response]]
  )$residuals

  draw_circular_blocks <- function() {
    blocks_needed <- ceiling(n / block_length)
    starts <- sample.int(n, blocks_needed, replace = TRUE)
    indices <- unlist(
      lapply(
        starts,
        function(start) {
          ((start - 1L + seq_len(block_length) - 1L) %% n) + 1L
        }
      ),
      use.names = FALSE
    )
    indices[seq_len(n)]
  }

  set.seed(seed)
  out <- vector("list", B)
  boundary_failures <- logical(B)

  for (b in seq_len(B)) {
    if (show_progress && (b == 1L || b %% 25L == 0L || b == B)) {
      message("Moving-block bootstrap replicate ", b, " of ", B)
    }

    sampled_residuals <- residuals[draw_circular_blocks()]
    probabilities <- stats::pnorm(sampled_residuals)
    probabilities <- pmin(pmax(probabilities, 1e-12), 1 - 1e-12)

    boot_data <- data
    boot_data[[response]] <- qLL3(
      probabilities,
      mu = fitted_pars$mu,
      sigma = fitted_pars$sigma,
      nu = fitted_pars$nu
    )
    boot_lower <- LL3_support_lower(boot_data[[response]])

    fb <- try(
      fit_LL3_gamlss(
        mu.formula,
        sigma.formula,
        nu.formula,
        data = boot_data,
        lower = boot_lower,
        information = "opg",
        allow_nonstationary_mu = FALSE,
        boundary_action = "ignore",
        control = control
      ),
      silent = TRUE
    )

    if (inherits(fb, "try-error") || !isTRUE(fb$converged)) {
      next
    }

    diagnostic <- LL3_boundary_diagnostic(
      fb,
      y = boot_data[[response]]
    )
    if (!isTRUE(diagnostic$inference_ready[1])) {
      boundary_failures[b] <- TRUE
      next
    }

    pb <- extract_LL3_parameters(fb)
    out[[b]] <- c(
      mu = mean(pb$mu),
      stats::setNames(
        stats::coef(fb, what = "sigma"),
        paste0("sigma:", names(stats::coef(fb, what = "sigma")))
      ),
      stats::setNames(
        stats::coef(fb, what = "nu"),
        paste0("nu:", names(stats::coef(fb, what = "nu")))
      )
    )
  }

  valid <- !vapply(out, is.null, logical(1))
  if (!any(valid)) {
    stop("No moving-block bootstrap refits produced interior solutions.")
  }

  estimates <- do.call(rbind, out[valid])

  list(
    estimates = estimates,
    convergence_rate = mean(valid),
    boundary_failure_rate = mean(boundary_failures),
    percentile_95 = apply(
      estimates,
      2,
      stats::quantile,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    ),
    B = B,
    block_length = block_length,
    method = "circular moving-block residual bootstrap"
  )
}
