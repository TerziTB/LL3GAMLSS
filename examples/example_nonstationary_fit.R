source(file.path("R", "LL3_distribution.R"))
source(file.path("R", "LL3_gamlss_family.R"))
source(file.path("R", "LL3_diagnostics.R"))
source(file.path("R", "LL3_fitting.R"))

set.seed(20)
n <- 500
time <- seq(-1, 1, length.out = n)
mu <- rep(-50, n)
sigma <- exp(log(70) + 0.55 * time)
nu <- exp(log(1.8) - 0.20 * time)
y <- rLL3(n, mu, sigma, nu)
dat <- data.frame(y = y, time = time)

models <- fit_LL3_candidate_models(
  dat, response = "y", covariate = "time",
  control = LL3_default_control(trace = TRUE)
)
print(LL3_model_table(models[c("stationary", "sigma", "nu", "joint")], n = n))

selected <- models$joint
print(check_LL3_fit(selected, y))
pars <- extract_LL3_parameters(selected)
dat$LL3_probability <- LL3_pit(y, pars$mu, pars$sigma, pars$nu)
dat$LL3_index <- LL3_index(y, pars$mu, pars$sigma, pars$nu)
print(LL3_residual_diagnostics(selected, y)$summary)

# Independent optimizer check for the same linear predictors.
comparison <- compare_LL3_gamlss_direct(
  selected, dat, y ~ 1, ~time, ~time, n_starts = 5, seed = 20
)
print(comparison[c(
  "gamlss_logLik", "direct_logLik", "absolute_logLik_difference",
  "sigma_coefficient_difference", "nu_coefficient_difference"
)])

# Small demonstration bootstrap. Use B >= 499 for reported intervals.
bootstrap_demo <- LL3_parametric_bootstrap(
  selected, dat, response = "y",
  mu.formula = y ~ 1,
  sigma.formula = ~time,
  nu.formula = ~time,
  B = 25, seed = 20, show_progress = TRUE
)
print(bootstrap_demo$convergence_rate)
print(bootstrap_demo$percentile_95)
