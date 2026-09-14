test_that("mu link rejects material support violations", {
  family <- LL3(lower = 0)

  expect_error(
    family$mu.linkfun(1),
    "exceeds the LL3 upper boundary"
  )

  tolerance <- family$mu_link_tolerance
  expect_equal(
    family$mu.linkfun(tolerance / 2),
    log(tolerance)
  )
  expect_lt(family$mu.linkinv(0), 0)
})

test_that("stationary GAMLSS agrees with direct optimization", {
  skip_if_not_installed("gamlss")
  set.seed(42)
  y <- rLL3(300, mu = -50, sigma = 70, nu = 1.8)
  fit <- fit_LL3_gamlss(
    y ~ 1,
    data = data.frame(y = y),
    boundary_action = "error"
  )
  comparison <- compare_LL3_gamlss_direct(
    fit,
    data = data.frame(y = y),
    mu.formula = y ~ 1,
    n_starts = 8,
    seed = 42
  )

  expect_true(LL3_fit_converged(fit))
  expect_lt(comparison$absolute_logLik_difference, 1e-3)
  expect_silent(LL3_assert_inference_ready(fit, y))
})

test_that("boundary contact blocks ordinary inference", {
  fake <- list(
    y = c(1, 2, 3, 4),
    mu.fv = rep(1 - 10 * .Machine$double.eps, 4),
    sigma.fv = rep(1, 4),
    nu.fv = rep(2, 4),
    converged = TRUE
  )
  attr(fake, "LL3_lower") <- 1

  expect_error(
    LL3_assert_inference_ready(fake, fake$y),
    "not ready for ordinary likelihood inference"
  )
})

test_that("model tables withhold information criteria at the boundary", {
  skip_if_not_installed("gamlss")
  set.seed(2026)
  y <- rLL3(120, mu = -40, sigma = 55, nu = 2.1)
  fit <- fit_LL3_gamlss(
    y ~ 1,
    data = data.frame(y = y),
    boundary_action = "error"
  )
  boundary_fit <- fit
  lower <- attr(boundary_fit, "LL3_lower")
  boundary_fit$mu.fv <- rep(
    lower - 10 * .Machine$double.eps * max(abs(lower), 1),
    length(boundary_fit$y)
  )

  model_table <- LL3_model_table(
    list(interior = fit, boundary = boundary_fit)
  )
  boundary_row <- model_table[model_table$model == "boundary", ]
  interior_row <- model_table[model_table$model == "interior", ]

  expect_true(interior_row$inference_ready)
  expect_false(boundary_row$inference_ready)
  expect_true(all(is.na(unlist(boundary_row[c(
    "deviance", "AIC", "BIC", "delta_AIC", "delta_BIC"
  )], use.names = FALSE))))
  expect_true(all(is.finite(unlist(interior_row[c(
    "deviance", "AIC", "BIC", "delta_AIC", "delta_BIC"
  )], use.names = FALSE))))
})
