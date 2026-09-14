\name{LL3-spei}
\alias{LL3_accumulate}
\alias{fit_LL3_spei}
\alias{print.LL3_spei}
\title{Calendar-Month LL3-Based SPEI Workflow}
\description{
Accumulates climatic water balance and fits independent conditional LL3
distributions for each calendar month before Gaussian standardization.
}
\usage{
LL3_accumulate(x, scale = 1L, weights = NULL)
fit_LL3_spei(data, response, date, scale = 1L, weights = NULL,
             sigma.formula = ~1, nu.formula = ~1,
             reference_start = NULL, reference_end = NULL,
             min_observations = 30L, lower = NULL,
             epsilon = 1e-10, information = "opg", control = NULL,
             boundary_action = c("error", "warning"))
\method{print}{LL3_spei}(x, ...)
}
\arguments{
  \item{x}{For accumulation, a numeric balance vector; for printing, an
  \code{LL3_spei} object.}
  \item{scale}{Positive integer accumulation scale in observations/months.}
  \item{weights}{Optional trailing accumulation weights.}
  \item{data}{Chronologically ordered data frame.}
  \item{response, date}{Names of the balance and date columns.}
  \item{sigma.formula, nu.formula}{Parametric scale and shape formulas.}
  \item{reference_start, reference_end}{Optional reference-period dates.}
  \item{min_observations}{Minimum reference observations required per month.}
  \item{lower}{Optional common or month-specific numerical threshold boundary.}
  \item{epsilon}{Probability clipping constant before Gaussian transformation.}
  \item{information, control}{LL3/GAMLSS fitting controls.}
  \item{boundary_action}{Action for numerical threshold contact.}
  \item{...}{Unused print arguments.}
}
\details{
Dates must be unique and strictly increasing. Accumulation is trailing, and
the calendar month of the final observation identifies the monthly model.
Only parametric predictor terms are currently supported for prediction. This
function implements an LL3-based standardized index; it does not claim that
the log-logistic distribution is the only admissible SPEI distribution.
}
\value{
\code{LL3_accumulate} returns a numeric vector. \code{fit_LL3_spei} returns an
object of class \code{LL3_spei} containing the index, fitted parameter curves,
monthly fits and settings.
}
\examples{
LL3_accumulate(1:6, scale = 3)
}
