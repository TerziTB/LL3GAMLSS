test_that("accumulation uses trailing calendar observations", {
  expect_equal(
    LL3_accumulate(1:5, scale = 3),
    c(NA, NA, 6, 9, 12)
  )
  expect_equal(
    LL3_accumulate(1:4, scale = 2, weights = c(0.25, 0.75)),
    c(NA, 1.75, 2.75, 3.75)
  )
})

test_that("monthly LL3-SPEI workflow produces calibrated finite indices", {
  skip_if_not_installed("gamlss")
  set.seed(2026)
  n <- 60L * 12L
  dates <- seq(as.Date("1961-01-01"), by = "month", length.out = n)
  time_scaled <- as.numeric(scale(seq_len(n)))
  sigma <- exp(log(70) + 0.15 * time_scaled)
  dat <- data.frame(
    date = dates,
    time_scaled = time_scaled,
    balance = rLL3(n, mu = -50, sigma = sigma, nu = 1.8)
  )

  fit <- suppressWarnings(
    fit_LL3_spei(
      dat,
      response = "balance",
      date = "date",
      scale = 3,
      sigma.formula = ~time_scaled,
      min_observations = 40,
      boundary_action = "warning"
    )
  )

  expect_s3_class(fit, "LL3_spei")
  expect_length(fit$fits, 12)
  expect_equal(sum(is.na(fit$index)), 2)
  expect_true(all(is.finite(fit$index[-c(1, 2)])))
  expect_lt(abs(mean(fit$index, na.rm = TRUE)), 0.15)
  expect_lt(abs(sd(fit$index, na.rm = TRUE) - 1), 0.15)
})

test_that("moving-block bootstrap returns coefficient intervals", {
  skip_if_not_installed("gamlss")
  set.seed(88)
  n <- 180L
  x <- seq(-1, 1, length.out = n)
  dat <- data.frame(
    y = rLL3(n, -50, exp(log(70) + 0.25 * x), 1.8),
    x = x
  )
  fit <- fit_LL3_gamlss(
    y ~ 1,
    sigma.formula = ~x,
    data = dat,
    boundary_action = "error"
  )

  result <- LL3_moving_block_bootstrap(
    fit,
    data = dat,
    response = "y",
    mu.formula = y ~ 1,
    sigma.formula = ~x,
    B = 3,
    block_length = 8,
    seed = 99,
    show_progress = FALSE
  )

  expect_equal(result$B, 3)
  expect_equal(result$block_length, 8)
  expect_true(nrow(result$estimates) >= 1)
  expect_true(all(is.finite(result$percentile_95)))
})
