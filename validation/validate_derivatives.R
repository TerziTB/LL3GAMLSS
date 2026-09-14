source(file.path("R", "LL3_distribution.R"))
source(file.path("R", "LL3_gamlss_family.R"))
if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Install numDeriv.")

cases <- data.frame(
  y = c(-20, 25, 200),
  mu = c(-50, -50, -50),
  sigma = c(70, 70, 70),
  nu = c(1.2, 1.8, 3.0)
)
rows <- vector("list", nrow(cases))

for (i in seq_len(nrow(cases))) {
  z <- cases[i, ]
  ll <- function(par) dLL3(z$y, par[1], par[2], par[3], log = TRUE)
  numerical_gradient <- numDeriv::grad(ll, c(z$mu, z$sigma, z$nu))
  analytical_gradient <- as.numeric(.LL3_score(z$y, z$mu, z$sigma, z$nu)[1, ])
  numerical_hessian <- numDeriv::hessian(ll, c(z$mu, z$sigma, z$nu))
  h <- .LL3_observed_hessian(z$y, z$mu, z$sigma, z$nu)
  analytical_hessian <- matrix(c(
    h$mm, h$ms, h$mn,
    h$ms, h$ss, h$sn,
    h$mn, h$sn, h$nn
  ), nrow = 3, byrow = TRUE)
  rows[[i]] <- data.frame(
    case = i,
    maximum_gradient_error = max(abs(numerical_gradient - analytical_gradient)),
    maximum_hessian_error = max(abs(numerical_hessian - analytical_hessian))
  )
}
result <- do.call(rbind, rows)
print(result)
stopifnot(max(result$maximum_gradient_error) < 1e-6)
stopifnot(max(result$maximum_hessian_error) < 1e-5)
output_dir <- file.path("validation", "results", "core_numerical")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, file.path(output_dir, "derivative_checks.csv"), row.names = FALSE)
message("Derivative validation passed.")
