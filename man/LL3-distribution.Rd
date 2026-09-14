\name{LL3-distribution}
\alias{dLL3}
\alias{pLL3}
\alias{qLL3}
\alias{rLL3}
\alias{LL3_mean}
\alias{LL3_variance}
\alias{LL3_median}
\alias{LL3_to_FAdist}
\alias{FAdist_to_LL3}
\title{Three-Parameter Shifted Log-Logistic Distribution}
\description{
Density, distribution, quantile, random-generation, moment, and
parameter-conversion functions for the shifted log-logistic distribution.
}
\usage{
dLL3(x, mu = 0, sigma = 1, nu = 2, log = FALSE)
pLL3(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
qLL3(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
rLL3(n, mu = 0, sigma = 1, nu = 2)
LL3_mean(mu, sigma, nu)
LL3_variance(mu, sigma, nu)
LL3_median(mu, sigma, nu)
LL3_to_FAdist(mu, sigma, nu)
FAdist_to_LL3(shape, scale, thres)
}
\arguments{
  \item{x, q}{Numeric quantiles.}
  \item{p}{Numeric probabilities.}
  \item{n}{Number of random values.}
  \item{mu}{Threshold parameter; observations have support above \code{mu}.}
  \item{sigma}{Strictly positive scale parameter.}
  \item{nu}{Strictly positive shape parameter.}
  \item{log, log.p}{Return logarithmic density or accept/return log probabilities.}
  \item{lower.tail}{Whether probabilities refer to the lower tail.}
  \item{shape, scale, thres}{The corresponding \pkg{FAdist} parameters.}
}
\details{
For \eqn{y>\mu}, the distribution function is
\deqn{F(y) = [1 + \{\sigma/(y-\mu)\}^{\nu}]^{-1}.}
The mapping to \pkg{FAdist} is \code{shape = 1 / nu},
\code{scale = log(sigma)}, and \code{thres = mu}.
}
\value{
The probability functions return numeric vectors. Moment functions return
numeric vectors, using infinity when a moment does not exist. Conversion
functions return data frames.
}
\examples{
p <- c(0.1, 0.5, 0.9)
x <- qLL3(p, mu = -50, sigma = 70, nu = 1.8)
pLL3(x, mu = -50, sigma = 70, nu = 1.8)
}
