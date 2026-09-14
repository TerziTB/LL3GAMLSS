library(LL3GAMLSS)

set.seed(2026)
n <- 66L * 12L
climate <- data.frame(
  date = seq(as.Date("1960-01-01"), by = "month", length.out = n),
  time_scaled = as.numeric(scale(seq_len(n)))
)
climate$water_balance <- rLL3(
  n,
  mu = -50,
  sigma = exp(log(70) + 0.20 * climate$time_scaled),
  nu = 2.5
)

fit <- fit_LL3_spei(
  climate,
  response = "water_balance",
  date = "date",
  scale = 3,
  sigma.formula = ~time_scaled,
  nu.formula = ~1
)

stopifnot(
  inherits(fit, "LL3_spei"),
  sum(is.finite(fit$index)) == n - 2L,
  all(vapply(fit$fits, LL3_fit_converged, logical(1)))
)
