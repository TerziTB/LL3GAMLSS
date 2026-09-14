\name{LL3-fitting}
\alias{LL3_support_lower}
\alias{LL3_default_control}
\alias{fit_LL3_gamlss}
\alias{fit_LL3_stationary}
\alias{fit_LL3_gamlss_stationary}
\alias{fit_LL3safe_gamlss_stationary}
\alias{fit_LL3_candidate_models}
\alias{fit_LL3_direct_linear}
\alias{compare_LL3_gamlss_direct}
\alias{LL3_parametric_bootstrap}
\alias{LL3_moving_block_bootstrap}
\title{Fit and Validate LL3 Distributional Regression Models}
\description{
Fitting helpers, an independent direct optimizer, candidate-model comparison,
and iid or moving-block uncertainty procedures for LL3 models.
}
\usage{
LL3_support_lower(y, relative_margin = 1e-6,
                  absolute_margin = sqrt(.Machine$double.eps))
LL3_default_control(n_cycles = 700, mu_step = 0.05,
                    sigma_step = 0.05, nu_step = 0.05,
                    trace = FALSE)
fit_LL3_gamlss(mu.formula, sigma.formula = ~1, nu.formula = ~1,
               data, lower = NULL, information = "opg",
               allow_nonstationary_mu = FALSE, mu.start = NULL,
               sigma.start = NULL, nu.start = NULL, control = NULL,
               boundary_action = c("warning", "error", "ignore"))
fit_LL3_stationary(y, trace = FALSE, n_cycles = 700,
                   information = "opg")
fit_LL3_gamlss_stationary(y, trace = FALSE, n_cycles = 700,
                          information = "opg")
fit_LL3safe_gamlss_stationary(y, trace = FALSE, n_cycles = 700,
                              information = "opg")
fit_LL3_candidate_models(data, response, covariate, lower = NULL,
                         information = "opg", control = NULL)
fit_LL3_direct_linear(y, X_sigma, X_nu, lower = NULL, start = NULL,
                      n_starts = 5, seed = 1, maxit = 5000)
compare_LL3_gamlss_direct(gamlss_fit, data, mu.formula,
                          sigma.formula = ~1, nu.formula = ~1,
                          n_starts = 5, seed = 1)
LL3_parametric_bootstrap(fit, data, response, mu.formula,
                         sigma.formula = ~1, nu.formula = ~1,
                         B = 499, seed = 1, control = NULL,
                         show_progress = TRUE,
                         boundary_action = c("error", "warning"))
LL3_moving_block_bootstrap(fit, data, response, mu.formula,
                           sigma.formula = ~1, nu.formula = ~1,
                           B = 499, block_length = NULL, seed = 1,
                           control = NULL, show_progress = TRUE,
                           boundary_action = c("error", "warning"))
}
\arguments{
  \item{y}{Finite numeric response vector.}
  \item{relative_margin, absolute_margin}{Margins used to place the numerical
  threshold boundary below the sample minimum.}
  \item{n_cycles, mu_step, sigma_step, nu_step, trace}{GAMLSS control values.}
  \item{mu.formula, sigma.formula, nu.formula}{Parameter predictor formulas.}
  \item{data}{Model data frame.}
  \item{lower}{Numerical threshold boundary.}
  \item{information}{Working curvature passed to \code{LL3}.}
  \item{allow_nonstationary_mu}{Whether to permit the experimental
  covariate-dependent threshold formulation.}
  \item{mu.start, sigma.start, nu.start, start}{Optional starting values.}
  \item{control}{A \code{gamlss.control} object.}
  \item{boundary_action}{Action when the fitted threshold contacts the
  numerical boundary.}
  \item{response, covariate}{Column names.}
  \item{X_sigma, X_nu}{Direct-optimizer design matrices.}
  \item{n_starts, seed, maxit}{Direct optimization controls.}
  \item{gamlss_fit, fit}{A fitted LL3 GAMLSS model.}
  \item{B}{Number of bootstrap replications.}
  \item{show_progress}{Whether to report bootstrap progress.}
  \item{block_length}{Circular moving-block length; defaults to a cube-root
  rule.}
}
\value{
Fit functions return GAMLSS objects or structured fit summaries. Comparison
and bootstrap functions return lists containing estimates and diagnostics.
}
