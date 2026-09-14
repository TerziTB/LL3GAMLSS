# Three-parameter shifted log-logistic distribution
#
# Parameterization:
#   mu    = threshold/shift, real
#   sigma = scale, strictly positive
#   nu    = shape, strictly positive
#
# For x > mu:
#   F(x) = 1 / {1 + [sigma / (x - mu)]^nu}
#   Q(p) = mu + sigma * [p / (1 - p)]^(1/nu)

.LL3_recycle <- function(...) {
  values <- list(...)
  n <- max(lengths(values))
  lapply(values, rep_len, length.out = n)
}

.LL3_log1pexp <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

.LL3_validate_parameters <- function(sigma, nu) {
  if (any(!is.finite(sigma) | sigma <= 0, na.rm = TRUE)) {
    stop("sigma must contain only finite, strictly positive values.")
  }
  if (any(!is.finite(nu) | nu <= 0, na.rm = TRUE)) {
    stop("nu must contain only finite, strictly positive values.")
  }
  invisible(TRUE)
}

dLL3 <- function(x, mu = 0, sigma = 1, nu = 2, log = FALSE) {
  a <- .LL3_recycle(x, mu, sigma, nu)
  x <- a[[1]]
  mu <- a[[2]]
  sigma <- a[[3]]
  nu <- a[[4]]

  .LL3_validate_parameters(sigma, nu)

  log_density <- rep(NA_real_, length(x))
  valid <- is.finite(mu) & is.finite(sigma) & is.finite(nu)
  outside <- !is.na(x) & valid & x <= mu
  pos_inf <- is.infinite(x) & x > 0 & valid
  inside <- is.finite(x) & valid & x > mu

  log_density[outside | pos_inf] <- -Inf

  if (any(inside)) {
    log_z <- log(x[inside] - mu[inside]) - log(sigma[inside])
    eta <- nu[inside] * log_z
    log_density[inside] <-
      log(nu[inside]) - log(sigma[inside]) +
      (nu[inside] - 1) * log_z -
      2 * .LL3_log1pexp(eta)
  }

  if (log) log_density else exp(log_density)
}

pLL3 <- function(q, mu = 0, sigma = 1, nu = 2,
                 lower.tail = TRUE, log.p = FALSE) {
  a <- .LL3_recycle(q, mu, sigma, nu)
  q <- a[[1]]
  mu <- a[[2]]
  sigma <- a[[3]]
  nu <- a[[4]]

  .LL3_validate_parameters(sigma, nu)

  ans <- rep(NA_real_, length(q))
  valid <- is.finite(mu) & is.finite(sigma) & is.finite(nu)
  below <- !is.na(q) & valid & q <= mu
  pos_inf <- is.infinite(q) & q > 0 & valid
  inside <- is.finite(q) & valid & q > mu

  if (log.p) {
    ans[below] <- if (lower.tail) -Inf else 0
    ans[pos_inf] <- if (lower.tail) 0 else -Inf
  } else {
    ans[below] <- if (lower.tail) 0 else 1
    ans[pos_inf] <- if (lower.tail) 1 else 0
  }

  if (any(inside)) {
    eta <- nu[inside] *
      (log(q[inside] - mu[inside]) - log(sigma[inside]))

    if (log.p) {
      ans[inside] <- if (lower.tail) {
        -.LL3_log1pexp(-eta)
      } else {
        -.LL3_log1pexp(eta)
      }
    } else {
      ans[inside] <- if (lower.tail) {
        stats::plogis(eta)
      } else {
        stats::plogis(-eta)
      }
    }
  }

  ans
}

qLL3 <- function(p, mu = 0, sigma = 1, nu = 2,
                 lower.tail = TRUE, log.p = FALSE) {
  a <- .LL3_recycle(p, mu, sigma, nu)
  p <- a[[1]]
  mu <- a[[2]]
  sigma <- a[[3]]
  nu <- a[[4]]

  .LL3_validate_parameters(sigma, nu)

  if (log.p) {
    if (any(p > 0, na.rm = TRUE)) {
      stop("Log-probabilities must be less than or equal to zero.")
    }
    p <- if (lower.tail) exp(p) else -expm1(p)
  } else {
    if (any(p < 0 | p > 1, na.rm = TRUE)) {
      stop("p must lie in [0, 1].")
    }
    if (!lower.tail) p <- 1 - p
  }

  ans <- mu + sigma * exp(stats::qlogis(p) / nu)
  ans[p == 0] <- mu[p == 0]
  ans[p == 1] <- Inf
  ans
}

rLL3 <- function(n, mu = 0, sigma = 1, nu = 2) {
  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1])
  if (!is.finite(n) || n < 1L) stop("n must be a positive integer.")

  qLL3(
    stats::runif(n),
    mu = rep_len(mu, n),
    sigma = rep_len(sigma, n),
    nu = rep_len(nu, n)
  )
}

LL3_mean <- function(mu, sigma, nu) {
  .LL3_validate_parameters(sigma, nu)
  theta <- pi / nu
  ifelse(nu > 1, mu + sigma * theta / sin(theta), Inf)
}

LL3_variance <- function(mu, sigma, nu) {
  .LL3_validate_parameters(sigma, nu)
  theta <- pi / nu
  ifelse(
    nu > 2,
    sigma^2 * theta * (2 / sin(2 * theta) - theta / sin(theta)^2),
    Inf
  )
}

LL3_median <- function(mu, sigma, nu) {
  .LL3_validate_parameters(sigma, nu)
  mu + sigma
}

LL3_to_FAdist <- function(mu, sigma, nu) {
  .LL3_validate_parameters(sigma, nu)
  data.frame(shape = 1 / nu, scale = log(sigma), thres = mu)
}

FAdist_to_LL3 <- function(shape, scale, thres) {
  if (any(!is.finite(shape) | shape <= 0, na.rm = TRUE)) {
    stop("FAdist shape must be finite and positive.")
  }
  data.frame(mu = thres, sigma = exp(scale), nu = 1 / shape)
}
