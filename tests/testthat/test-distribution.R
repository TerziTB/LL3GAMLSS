test_that("probability functions invert and integrate", {
  pars <- c(mu = -50, sigma = 70, nu = 1.8)
  probabilities <- seq(0.001, 0.999, length.out = 199)
  quantiles <- qLL3(
    probabilities,
    pars["mu"],
    pars["sigma"],
    pars["nu"]
  )

  expect_equal(
    pLL3(quantiles, pars["mu"], pars["sigma"], pars["nu"]),
    probabilities,
    tolerance = 1e-12
  )

  area <- integrate(
    function(x) dLL3(x, pars["mu"], pars["sigma"], pars["nu"]),
    lower = pars["mu"],
    upper = Inf,
    rel.tol = 1e-9
  )$value
  expect_equal(area, 1, tolerance = 1e-8)
})

test_that("FAdist parameter conversion round-trips", {
  converted <- LL3_to_FAdist(-50, 70, 1.8)
  restored <- FAdist_to_LL3(
    converted$shape,
    converted$scale,
    converted$thres
  )
  expect_equal(restored$mu, -50)
  expect_equal(restored$sigma, 70)
  expect_equal(restored$nu, 1.8)
})

test_that("moments and support endpoints are correct", {
  expect_equal(LL3_median(-50, 70, 1.8), 20)
  expect_true(is.infinite(LL3_variance(-50, 70, 1.8)))
  expect_equal(pLL3(-50, -50, 70, 1.8), 0)
  expect_equal(qLL3(0, -50, 70, 1.8), -50)
  expect_true(is.infinite(qLL3(1, -50, 70, 1.8)))
})
