source(file.path("R", "LL3_distribution.R"))
source(file.path("R", "LL3_gamlss_family.R"))
source(file.path("R", "LL3_diagnostics.R"))
source(file.path("R", "LL3_fitting.R"))

set.seed(10)
y <- rLL3(500, mu = -50, sigma = 70, nu = 1.8)
result <- fit_LL3_stationary(y, trace = TRUE)
print(result$parameters)
print(check_LL3_fit(result$fit, y))
print(LL3_boundary_diagnostic(result$fit, y))
print(LL3_residual_diagnostics(result$fit, y)$summary)

# Independent numerical-likelihood comparison.
dat <- data.frame(y = y)
comparison <- compare_LL3_gamlss_direct(
  result$fit, dat, y ~ 1, ~1, ~1, n_starts = 5, seed = 10
)
print(comparison[c(
  "gamlss_logLik", "direct_logLik", "absolute_logLik_difference",
  "gamlss_mu", "direct_mu"
)])
