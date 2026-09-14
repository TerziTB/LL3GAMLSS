\name{LL3-family}
\alias{LL3}
\alias{LL3safe}
\title{Classic GAMLSS Family for the Shifted Log-Logistic Distribution}
\description{
Constructs a three-parameter \code{gamlss.family} object with a support-safe
threshold link.
}
\usage{
LL3(lower, sigma.link = "log", nu.link = "log",
    information = c("opg", "observed"))
LL3safe(lower, sigma.link = "log", nu.link = "log",
        information = c("opg", "observed"))
}
\arguments{
  \item{lower}{Finite upper boundary imposed on the threshold parameter. It
  must be strictly below every fitting response.}
  \item{sigma.link, nu.link}{GAMLSS links for positive scale and shape.}
  \item{information}{Working curvature: stable outer-product-of-gradients
  curvature or the observed Hessian.}
}
\details{
The threshold link is \eqn{\eta_\mu=\log(lower-\mu)}. Values exceeding the
boundary by no more than floating-point tolerance are projected for numerical
safety; material violations raise an error. The OPG curvature is an
optimization device and is not a reported covariance estimator.
}
\value{A \code{gamlss.family} object.}
