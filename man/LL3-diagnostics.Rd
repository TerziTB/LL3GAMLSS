\name{LL3-diagnostics}
\alias{extract_LL3_parameters}
\alias{LL3_fit_converged}
\alias{check_LL3_fit}
\alias{LL3_pit}
\alias{LL3_index}
\alias{LL3_residual_diagnostics}
\alias{LL3_boundary_diagnostic}
\alias{LL3_assert_inference_ready}
\alias{LL3_model_table}
\alias{plot_LL3_diagnostics}
\title{Diagnostics for LL3 GAMLSS Fits}
\description{
Parameter extraction, validity and support checks, probability-integral
transform diagnostics, model tables, and an explicit inference-readiness gate.
}
\usage{
extract_LL3_parameters(object)
LL3_fit_converged(fit)
check_LL3_fit(object, y = NULL, boundary_tolerance_multiplier = 100)
LL3_pit(y, mu, sigma, nu, epsilon = 1e-10)
LL3_index(y, mu, sigma, nu, epsilon = 1e-10)
LL3_residual_diagnostics(object, y = NULL, max_lag = 24)
LL3_boundary_diagnostic(object, y = NULL,
                        boundary_tolerance_multiplier = 100)
LL3_assert_inference_ready(object, y = NULL,
                           action = c("error", "warning"))
LL3_model_table(models, n = NULL)
plot_LL3_diagnostics(object, y = NULL)
}
\arguments{
  \item{object, fit}{A fitted LL3 GAMLSS model.}
  \item{y}{Optional response vector.}
  \item{boundary_tolerance_multiplier}{Floating-point boundary tolerance
  multiplier.}
  \item{mu, sigma, nu}{LL3 parameters.}
  \item{epsilon}{Probability clipping constant.}
  \item{max_lag}{Maximum residual ACF lag.}
  \item{action}{Whether a non-interior fit raises an error or warning.}
  \item{models}{Named list of fitted models. Boundary-contact and otherwise
  non-regular fits are retained in the table but do not receive ordinary
  deviance, AIC, BIC, or delta values.}
  \item{n}{Optional common sample size for BIC.}
}
\value{
The functions return parameter data frames, logical diagnostics, normalized
indices, diagnostic lists or tables. \code{LL3_model_table()} includes an
\code{inference_ready} column and reports information criteria only for fits
that pass \code{LL3_assert_inference_ready()}. The plotting helper returns its
diagnostic list invisibly.
}
