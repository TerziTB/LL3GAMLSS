source(file.path("R", "LL3_distribution.R"))
source(file.path("R", "LL3_gamlss_family.R"))
source(file.path("R", "LL3_diagnostics.R"))
source(file.path("R", "LL3_fitting.R"))

repetitions <- as.integer(Sys.getenv("LL3_REPS", "200"))
sample_sizes <- c(60, 100, 250, 1000)
master_seed <- 2026
truth <- list(mu = -50, sigma_intercept = log(70), sigma_slope = 0.55, nu = 1.8)
set.seed(master_seed)
seeds <- sample.int(.Machine$integer.max, repetitions * length(sample_sizes))
empty_row <- function(n, seed) data.frame(
  n = n, seed = seed, success = FALSE,
  mu_hat = NA_real_, sigma_intercept_hat = NA_real_, sigma_slope_hat = NA_real_,
  sigma_curve_rmse = NA_real_, sigma_curve_correlation = NA_real_,
  delta_AIC = NA_real_, delta_BIC = NA_real_,
  residual_mean = NA_real_, residual_sd = NA_real_, residual_tail_rate_95 = NA_real_,
  direct_logLik_difference = NA_real_
)
rows <- list(); index <- 0L

for (n in sample_sizes) {
  for (b in seq_len(repetitions)) {
    index <- index + 1L
    set.seed(seeds[index])
    time <- seq(-1, 1, length.out = n)
    sigma_true <- exp(truth$sigma_intercept + truth$sigma_slope * time)
    y <- rLL3(n, truth$mu, sigma_true, truth$nu)
    dat <- data.frame(y = y, time = time)
    lower <- LL3_support_lower(y)
    control <- LL3_default_control(trace = FALSE)

    m0 <- try(fit_LL3_gamlss(y ~ 1, ~1, ~1, dat, lower = lower, control = control), silent = TRUE)
    ms <- if (!inherits(m0, "try-error")) {
      p0 <- extract_LL3_parameters(m0)
      try(fit_LL3_gamlss(y ~ 1, ~time, ~1, dat, lower = lower,
                         mu.start = p0$mu, sigma.start = p0$sigma, nu.start = p0$nu,
                         control = control), silent = TRUE)
    } else m0

    success <- !inherits(ms, "try-error") && LL3_fit_converged(ms) &&
      all(check_LL3_fit(ms, y)[c("support_ok", "sigma_positive", "nu_positive")])
    row <- empty_row(n, seeds[index])
    row$success <- success

    if (success) {
      p <- extract_LL3_parameters(ms)
      bs <- stats::coef(ms, what = "sigma")
      z <- LL3_index(y, p$mu, p$sigma, p$nu)
      direct_diff <- NA_real_
      if (b <= min(20L, repetitions)) {
        cmp <- try(compare_LL3_gamlss_direct(
          ms, dat, y ~ 1, ~time, ~1, n_starts = 3, seed = seeds[index]
        ), silent = TRUE)
        if (!inherits(cmp, "try-error")) direct_diff <- cmp$absolute_logLik_difference
      }
      row$mu_hat <- mean(p$mu)
      row$sigma_intercept_hat <- unname(bs["(Intercept)"])
      row$sigma_slope_hat <- unname(bs["time"])
      row$sigma_curve_rmse <- sqrt(mean((p$sigma - sigma_true)^2))
      row$sigma_curve_correlation <- stats::cor(p$sigma, sigma_true)
      row$delta_AIC <- stats::AIC(m0) - stats::AIC(ms)
      row$delta_BIC <- (stats::deviance(m0) + log(n) * m0$df.fit) -
        (stats::deviance(ms) + log(n) * ms$df.fit)
      row$residual_mean <- mean(z)
      row$residual_sd <- stats::sd(z)
      row$residual_tail_rate_95 <- mean(abs(z) > stats::qnorm(0.975))
      row$direct_logLik_difference <- direct_diff
      # Values are assigned explicitly so failed and successful rows share a schema.
      invisible(NULL)
      #
    }
    rows[[index]] <- row
  }
  message("Completed n = ", n)
}

results <- do.call(rbind, rows)
summary_rows <- lapply(sample_sizes, function(n) {
  d <- results[results$n == n, , drop = FALSE]
  s <- d[d$success, , drop = FALSE]
  data.frame(
    n = n, attempted = nrow(d), successful = nrow(s),
    success_rate = mean(d$success),
    mu_bias = mean(s$mu_hat - truth$mu),
    sigma_intercept_bias = mean(s$sigma_intercept_hat - truth$sigma_intercept),
    sigma_slope_bias = mean(s$sigma_slope_hat - truth$sigma_slope),
    sigma_slope_rmse = sqrt(mean((s$sigma_slope_hat - truth$sigma_slope)^2)),
    mean_sigma_curve_rmse = mean(s$sigma_curve_rmse),
    mean_sigma_curve_correlation = mean(s$sigma_curve_correlation),
    AIC_selection_rate = mean(s$delta_AIC > 0),
    BIC_selection_rate = mean(s$delta_BIC > 0),
    mean_residual_mean = mean(s$residual_mean),
    mean_residual_sd = mean(s$residual_sd),
    mean_residual_tail_rate_95 = mean(s$residual_tail_rate_95),
    maximum_direct_logLik_difference = if (all(is.na(s$direct_logLik_difference))) {
      NA_real_
    } else max(s$direct_logLik_difference, na.rm = TRUE)
  )
})
summary <- do.call(rbind, summary_rows)
dir.create(file.path("validation", "results"), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(results, file.path("validation", "results", "sigma_nonstationarity_raw.csv"), row.names = FALSE)
utils::write.csv(summary, file.path("validation", "results", "sigma_nonstationarity_summary.csv"), row.names = FALSE)
print(summary)
