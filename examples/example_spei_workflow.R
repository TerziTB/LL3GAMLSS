library(LL3GAMLSS)

set.seed(2026)
n <- 60L * 12L
date <- seq(as.Date("1961-01-01"), by = "month", length.out = n)
time_scaled <- as.numeric(scale(seq_len(n)))

climate <- data.frame(date = date, time_scaled = time_scaled)
climate$water_balance <- rLL3(
  n,
  mu = -50,
  sigma = exp(log(70) + 0.15 * time_scaled),
  nu = 1.8
)

fit <- fit_LL3_spei(
  climate,
  response = "water_balance",
  date = "date",
  scale = 3,
  sigma.formula = ~ time_scaled,
  min_observations = 40
)

print(fit)
head(fit$results)
plot(fit$results$date, fit$index, type = "l", ylab = "LL3-SPEI", xlab = "")
