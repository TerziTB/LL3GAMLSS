# Publication-scale short-record and robustness stress tests for the validated
# sigma-nonstationary LL3 model. Configure with LL3_STRESS_REPS,
# LL3_STRESS_CORES, and LL3_STRESS_WRITE environment variables.

source(file.path("R", "LL3_distribution.R"))
source(file.path("R", "LL3_gamlss_family.R"))
source(file.path("R", "LL3_diagnostics.R"))
source(file.path("R", "LL3_fitting.R"))

reps <- as.integer(Sys.getenv("LL3_STRESS_REPS", "100"))
cores <- as.integer(Sys.getenv("LL3_STRESS_CORES", "1"))
write_results <- identical(toupper(Sys.getenv("LL3_STRESS_WRITE", "FALSE")), "TRUE")
if (!is.finite(reps) || reps < 1L) stop("LL3_STRESS_REPS must be a positive integer.")
if (!is.finite(cores) || cores < 1L) stop("LL3_STRESS_CORES must be a positive integer.")
cores <- min(cores, parallel::detectCores(logical = TRUE))

sample_sizes <- c(30L, 60L, 100L)
scenarios <- c("iid", "ar1", "contaminated")
truth <- list(mu = -50, sigma_intercept = log(70), sigma_slope = 0.35, nu = 1.8)
master_seed <- 20260817L

simulate_probabilities <- function(n, scenario) {
  if (scenario %in% c("iid", "contaminated")) return(stats::runif(n))
  rho <- 0.55
  innovations <- stats::rnorm(n, sd = sqrt(1 - rho^2))
  z <- numeric(n)
  z[1] <- stats::rnorm(1)
  for (i in 2:n) z[i] <- rho * z[i - 1L] + innovations[i]
  stats::pnorm(z)
}

run_one <- function(n, scenario, seed) {
  set.seed(seed)
  x <- seq(-1, 1, length.out = n)
  sigma <- exp(truth$sigma_intercept + truth$sigma_slope * x)
  y <- qLL3(simulate_probabilities(n, scenario), truth$mu, sigma, truth$nu)
  if (identical(scenario, "contaminated")) {
    contaminated <- sample.int(n, max(1L, round(0.05 * n)))
    y[contaminated] <- y[contaminated] + 4 * sigma[contaminated]
  }

  dat <- data.frame(y = y, x = x)
  result <- data.frame(
    n = n, scenario = scenario, seed = seed,
    fits_returned = FALSE, stationary_converged = FALSE,
    sigma_converged = FALSE, inference_ready = FALSE, complete = FALSE,
    boundary_contact = NA, sigma_slope = NA_real_, sigma_slope_error = NA_real_,
    AIC_selects_sigma = NA, BIC_selects_sigma = NA,
    max_residual_acf = NA_real_, failure_reason = "fit error",
    stringsAsFactors = FALSE
  )

  fitted <- try(suppressWarnings({
    lower <- LL3_support_lower(y)
    control <- LL3_default_control()
    stationary <- fit_LL3_gamlss(
      y ~ 1, ~1, ~1, data = dat, lower = lower, control = control
    )
    starts <- extract_LL3_parameters(stationary)
    sigma_model <- fit_LL3_gamlss(
      y ~ 1, ~x, ~1, data = dat, lower = lower,
      mu.start = starts$mu, sigma.start = starts$sigma,
      nu.start = starts$nu, control = control
    )
    list(stationary = stationary, sigma = sigma_model)
  }), silent = TRUE)
  if (inherits(fitted, "try-error")) return(result)

  result$fits_returned <- TRUE
  result$stationary_converged <- LL3_fit_converged(fitted$stationary)
  result$sigma_converged <- LL3_fit_converged(fitted$sigma)

  diagnostic <- try(LL3_boundary_diagnostic(fitted$sigma, y), silent = TRUE)
  if (!inherits(diagnostic, "try-error")) {
    result$boundary_contact <- diagnostic$exact_boundary_contact[1]
    result$inference_ready <- isTRUE(diagnostic$inference_ready[1])
  }

  if (!result$stationary_converged || !result$sigma_converged) {
    result$failure_reason <- "non-convergence"
    return(result)
  }
  if (!result$inference_ready) {
    result$failure_reason <- "boundary or invalid fit"
    return(result)
  }

  selected <- try(
    LL3_model_table(fitted[c("stationary", "sigma")], n = n),
    silent = TRUE
  )
  residual <- try(
    LL3_residual_diagnostics(fitted$sigma, y, max_lag = min(12, n - 1L)),
    silent = TRUE
  )
  slope <- try(unname(stats::coef(fitted$sigma, what = "sigma")["x"]), silent = TRUE)
  if (inherits(selected, "try-error") || inherits(residual, "try-error") ||
      inherits(slope, "try-error") || length(slope) != 1L || !is.finite(slope)) {
    result$failure_reason <- "post-fit diagnostic error"
    return(result)
  }

  result$sigma_slope <- slope
  result$sigma_slope_error <- slope - truth$sigma_slope
  result$AIC_selects_sigma <- selected$model[which.min(selected$AIC)] == "sigma"
  result$BIC_selects_sigma <- selected$model[which.min(selected$BIC)] == "sigma"
  result$max_residual_acf <- residual$max_absolute_acf
  result$complete <- all(is.finite(c(
    result$sigma_slope_error, result$max_residual_acf
  )))
  result$failure_reason <- if (result$complete) "" else "non-finite diagnostic"
  result
}

safe_rate <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

selection_se <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  p <- mean(x)
  sqrt(p * (1 - p) / length(x))
}

design <- expand.grid(
  n = sample_sizes, scenario = scenarios, replicate = seq_len(reps),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
set.seed(master_seed)
design$seed <- sample.int(.Machine$integer.max, nrow(design))

started <- Sys.time()
indices <- seq_len(nrow(design))
runner <- function(i) run_one(design$n[i], design$scenario[i], design$seed[i])

if (cores > 1L) {
  cluster <- parallel::makeCluster(cores)
  parallel::clusterEvalQ(cluster, {
    source(file.path("R", "LL3_distribution.R"))
    source(file.path("R", "LL3_gamlss_family.R"))
    source(file.path("R", "LL3_diagnostics.R"))
    source(file.path("R", "LL3_fitting.R"))
    NULL
  })
  parallel::clusterExport(
    cluster,
    c("design", "truth", "simulate_probabilities", "run_one"),
    envir = environment()
  )
  rows <- parallel::parLapply(cluster, indices, runner)
  parallel::stopCluster(cluster)
  cluster <- NULL
} else {
  rows <- lapply(indices, runner)
}

raw <- do.call(rbind, rows)
elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

summary <- do.call(rbind, lapply(split(raw, interaction(raw$n, raw$scenario)), function(z) {
  complete <- z[z$complete, , drop = FALSE]
  errors <- complete$sigma_slope_error
  data.frame(
    n = z$n[1], scenario = z$scenario[1], attempted = nrow(z),
    fits_returned = sum(z$fits_returned),
    stationary_convergence_rate = mean(z$stationary_converged),
    sigma_convergence_rate = mean(z$sigma_converged),
    inference_ready_rate = mean(z$inference_ready),
    complete = nrow(complete), completion_rate = mean(z$complete),
    boundary_contact_rate = safe_rate(z$boundary_contact),
    sigma_slope_bias = if (length(errors)) mean(errors) else NA_real_,
    sigma_slope_rmse = if (length(errors)) sqrt(mean(errors^2)) else NA_real_,
    AIC_sigma_selection_rate = safe_rate(complete$AIC_selects_sigma),
    AIC_selection_mcse = selection_se(complete$AIC_selects_sigma),
    BIC_sigma_selection_rate = safe_rate(complete$BIC_selects_sigma),
    BIC_selection_mcse = selection_se(complete$BIC_selects_sigma),
    mean_max_residual_acf = if (nrow(complete)) mean(complete$max_residual_acf) else NA_real_
  )
}))
rownames(summary) <- NULL
summary <- summary[order(summary$scenario, summary$n), ]
print(summary, digits = 4, row.names = FALSE)
cat(sprintf("Completed %d simulations in %.1f seconds using %d core(s).\n",
            nrow(raw), elapsed_seconds, cores))

if (write_results) {
  result_dir <- file.path("validation", "results")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(raw, file.path(result_dir, "publication_stress_raw.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(result_dir, "publication_stress_summary.csv"), row.names = FALSE)
  manifest <- data.frame(
    run_date_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    master_seed = master_seed, replicates_per_cell = reps,
    sample_sizes = paste(sample_sizes, collapse = ";"),
    scenarios = paste(scenarios, collapse = ";"),
    attempted = nrow(raw), cores = cores,
    elapsed_seconds = elapsed_seconds,
    R_version = R.version.string,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    manifest, file.path(result_dir, "publication_stress_manifest.csv"),
    row.names = FALSE
  )
}
