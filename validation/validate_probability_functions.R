source(file.path("R", "LL3_distribution.R"))

tol_inverse <- 1e-9
tol_integral <- 1e-7
parameter_sets <- list(
  c(mu = -100, sigma = 80, nu = 1.2),
  c(mu = -50, sigma = 70, nu = 1.8),
  c(mu = -20, sigma = 35, nu = 3.0),
  c(mu = 0, sigma = 1, nu = 5.0)
)
p <- seq(0.001, 0.999, length.out = 999)
rows <- list()

for (i in seq_along(parameter_sets)) {
  pars <- parameter_sets[[i]]
  x <- qLL3(p, pars["mu"], pars["sigma"], pars["nu"])
  inversion_error <- max(abs(pLL3(x, pars["mu"], pars["sigma"], pars["nu"]) - p))
  area <- stats::integrate(
    function(z) dLL3(z, pars["mu"], pars["sigma"], pars["nu"]),
    lower = pars["mu"], upper = Inf,
    subdivisions = 3000L, rel.tol = 1e-9
  )$value
  rows[[i]] <- data.frame(
    set = i, inversion_error = inversion_error,
    density_integral = area, integral_error = abs(area - 1)
  )
}
result <- do.call(rbind, rows)
stopifnot(max(result$inversion_error) < tol_inverse)
stopifnot(max(result$integral_error) < tol_integral)

if (requireNamespace("FAdist", quietly = TRUE)) {
  pars <- c(mu = -50, sigma = 70, nu = 1.8)
  fad <- LL3_to_FAdist(pars["mu"], pars["sigma"], pars["nu"])
  x <- qLL3(p, pars["mu"], pars["sigma"], pars["nu"])
  result_FAdist <- data.frame(
    density = max(abs(
      dLL3(x, pars["mu"], pars["sigma"], pars["nu"]) -
        FAdist::dllog3(x, fad$shape, fad$scale, fad$thres)
    )),
    CDF = max(abs(
      pLL3(x, pars["mu"], pars["sigma"], pars["nu"]) -
        FAdist::pllog3(x, fad$shape, fad$scale, fad$thres)
    )),
    quantile = max(abs(
      qLL3(p, pars["mu"], pars["sigma"], pars["nu"]) -
        FAdist::qllog3(p, fad$shape, fad$scale, fad$thres)
    ))
  )
  stopifnot(result_FAdist$density < 1e-10)
  stopifnot(result_FAdist$CDF < 1e-10)
  stopifnot(result_FAdist$quantile < 1e-8)
  print(result_FAdist)
} else {
  message("FAdist is not installed; independent function comparison skipped.")
}

print(result)
output_dir <- file.path("validation", "results", "core_numerical")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, file.path(output_dir, "probability_function_checks.csv"), row.names = FALSE)
if (exists("result_FAdist")) {
  utils::write.csv(result_FAdist, file.path(output_dir, "FAdist_probability_agreement.csv"), row.names = FALSE)
}
message("Probability-function validation passed.")
