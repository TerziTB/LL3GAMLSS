# R/LL3_gamlss_family.R
#
# Classic-GAMLSS family for the three-parameter shifted log-logistic
# distribution.
#
# Distribution:
#
#   F(y) = 1 / {1 + [sigma / (y - mu)]^nu},   y > mu
#
# Parameters:
#
#   mu    threshold / shift
#   sigma positive scale
#   nu    positive shape
#
# The distribution functions dLL3(), pLL3(), qLL3(), and rLL3() are defined
# in R/LL3_distribution.R. The starting-value helper is defined below.
#
# This file defines:
#
#   .LL3_terms()
#   .LL3_score()
#   .LL3_observed_hessian()
#   .LL3_opg_hessian()
#   LL3()
#
# Support-safe mu link:
#
#   eta_mu = log(lower - mu)
#   mu     = lower - exp(eta_mu)
#
# The inverse link maps finite eta_mu values below lower. Floating-point
# contact with lower is tolerated because lower itself must be below all
# observations.

# Conservative LL3 starting values
#
# These values are intended only to start the GAMLSS optimizer.
# They are not treated as parameter estimates.
#
# The threshold is placed safely below the sample minimum, the scale is
# based on the sample median and robust spread, and the shape starts from
# the moderate generic value nu = 2.

.LL3_start_values <- function(y) {
  y <- as.numeric(y)

  if (
    length(y) < 3L ||
    any(!is.finite(y))
  ) {
    stop(
      "y must contain at least three finite observations."
    )
  }

  if (length(unique(y)) < 2L) {
    stop(
      "LL3 starting values cannot be constructed from a constant response."
    )
  }

  sample_quantiles <- stats::quantile(
    y,
    probabilities = c(0.25, 0.50, 0.75),
    names = FALSE,
    type = 8
  )

  q25 <- unname(sample_quantiles[1])
  q50 <- unname(sample_quantiles[2])
  q75 <- unname(sample_quantiles[3])

  robust_spread <- q75 - q25

  if (
    !is.finite(robust_spread) ||
    robust_spread <= 0
  ) {
    robust_spread <- stats::mad(
      y,
      center = q50,
      constant = 1,
      na.rm = TRUE
    )
  }

  if (
    !is.finite(robust_spread) ||
    robust_spread <= 0
  ) {
    robust_spread <- diff(range(y)) / 4
  }

  if (
    !is.finite(robust_spread) ||
    robust_spread <= 0
  ) {
    robust_spread <- max(
      0.1 * abs(q50),
      1
    )
  }

  numerical_scale <- max(
    1,
    abs(min(y)),
    abs(q50),
    robust_spread
  )

  threshold_gap <- max(
    0.50 * robust_spread,
    sqrt(.Machine$double.eps) * numerical_scale,
    1e-6
  )

  mu_start <- min(y) - threshold_gap

  sigma_start <- max(
    q50 - mu_start,
    0.50 * robust_spread,
    sqrt(.Machine$double.eps) * numerical_scale,
    1e-6
  )

  nu_start <- 2

  starting_values <- c(
    mu = unname(mu_start),
    sigma = unname(sigma_start),
    nu = unname(nu_start)
  )

  stopifnot(
    identical(
      names(starting_values),
      c("mu", "sigma", "nu")
    ),
    all(is.finite(starting_values)),
    starting_values["mu"] < min(y),
    starting_values["sigma"] > 0,
    starting_values["nu"] > 0
  )

  starting_values
}
# Dependency check

.LL3_check_family_dependencies <- function() {
  required_functions <- c(
    "dLL3",
    "pLL3",
    "qLL3",
    "rLL3"
  )

  lookup_environment <- environment()

  missing_functions <- required_functions[
    !vapply(
      required_functions,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE,
      envir = lookup_environment
    )
  ]

  if (length(missing_functions) > 0L) {
    stop(
      "The following LL3 distribution functions are missing: ",
      paste(
        missing_functions,
        collapse = ", "
      ),
      ". Source R/LL3_distribution.R before ",
      "R/LL3_gamlss_family.R."
    )
  }

  invisible(TRUE)
}


# Internal derivative helpers

.LL3_family_recycle <- function(...) {
  arguments <- list(...)

  output_length <- max(
    lengths(arguments)
  )

  lapply(
    arguments,
    rep_len,
    length.out = output_length
  )
}


.LL3_terms <- function(
    y,
    mu,
    sigma,
    nu
) {
  recycled <- .LL3_family_recycle(
    y,
    mu,
    sigma,
    nu
  )

  y <- recycled[[1]]
  mu <- recycled[[2]]
  sigma <- recycled[[3]]
  nu <- recycled[[4]]

  numerical_scale <- pmax(
    abs(y),
    abs(mu),
    abs(sigma),
    1
  )

  distance <- y - mu

  # This barrier is used only by the score and Hessian calculations.
  # The probability density itself remains zero outside y > mu.
  distance_safe <- pmax(
    distance,
    sqrt(.Machine$double.eps) *
      numerical_scale
  )

  log_standardized_distance <- log(
    distance_safe / sigma
  )

  linear_term <-
    nu *
    log_standardized_distance

  probability <- stats::plogis(
    linear_term
  )

  signed_probability <-
    2 * probability - 1

  probability_curvature <-
    2 *
    probability *
    (1 - probability)

  list(
    y = y,
    mu = mu,
    sigma = sigma,
    nu = nu,

    distance = distance,
    distance_safe = distance_safe,

    log_standardized_distance =
      log_standardized_distance,

    linear_term =
      linear_term,

    probability =
      probability,

    signed_probability =
      signed_probability,

    probability_curvature =
      probability_curvature
  )
}


# Per-observation score

.LL3_score <- function(
    y,
    mu,
    sigma,
    nu
) {
  terms <- .LL3_terms(
    y = y,
    mu = mu,
    sigma = sigma,
    nu = nu
  )

  score <- cbind(
    mu =
      (
        1 +
          terms$nu *
          terms$signed_probability
      ) /
      terms$distance_safe,

    sigma =
      terms$nu *
      terms$signed_probability /
      terms$sigma,

    nu =
      1 / terms$nu -
      terms$log_standardized_distance *
      terms$signed_probability
  )

  colnames(score) <- c(
    "mu",
    "sigma",
    "nu"
  )

  score
}


# Exact observed Hessian
#
# These are per-observation second derivatives of the log likelihood with
# respect to the natural parameters mu, sigma, and nu.

.LL3_observed_hessian <- function(
    y,
    mu,
    sigma,
    nu
) {
  terms <- .LL3_terms(
    y = y,
    mu = mu,
    sigma = sigma,
    nu = nu
  )

  hessian_mu_mu <-
    (
      1 +
        terms$nu *
        terms$signed_probability -
        terms$nu^2 *
        terms$probability_curvature
    ) /
    terms$distance_safe^2

  hessian_sigma_sigma <-
    -terms$nu *
    (
      terms$signed_probability +
        terms$nu *
        terms$probability_curvature
    ) /
    terms$sigma^2

  hessian_nu_nu <-
    -1 / terms$nu^2 -
    terms$probability_curvature *
    terms$log_standardized_distance^2

  hessian_mu_sigma <-
    -terms$nu^2 *
    terms$probability_curvature /
    (
      terms$sigma *
        terms$distance_safe
    )

  hessian_mu_nu <-
    (
      terms$signed_probability +
        terms$nu *
        terms$probability_curvature *
        terms$log_standardized_distance
    ) /
    terms$distance_safe

  hessian_sigma_nu <-
    (
      terms$signed_probability +
        terms$nu *
        terms$probability_curvature *
        terms$log_standardized_distance
    ) /
    terms$sigma

  list(
    mm = hessian_mu_mu,
    ss = hessian_sigma_sigma,
    nn = hessian_nu_nu,

    ms = hessian_mu_sigma,
    mn = hessian_mu_nu,
    sn = hessian_sigma_nu
  )
}


# Outer-product-of-gradients working curvature
#
# This is used as stable working curvature by the GAMLSS fitting algorithm.
# The returned values are negative per-observation score products because
# GAMLSS expects second derivatives of the log likelihood.

.LL3_opg_hessian <- function(
    y,
    mu,
    sigma,
    nu
) {
  score <- .LL3_score(
    y = y,
    mu = mu,
    sigma = sigma,
    nu = nu
  )

  if (
    !is.matrix(score) ||
      !all(
        c(
          "mu",
          "sigma",
          "nu"
        ) %in% colnames(score)
      )
  ) {
    stop(
      ".LL3_score() must return a matrix containing ",
      "mu, sigma, and nu columns."
    )
  }

  score_mu <- score[, "mu"]
  score_sigma <- score[, "sigma"]
  score_nu <- score[, "nu"]

  if (
    any(!is.finite(score_mu)) ||
      any(!is.finite(score_sigma)) ||
      any(!is.finite(score_nu))
  ) {
    stop(
      "Non-finite LL3 scores were encountered while ",
      "constructing the OPG working curvature."
    )
  }

  curvature_floor <- 1e-12

  list(
    mm =
      -pmax(
        score_mu^2,
        curvature_floor
      ),

    ss =
      -pmax(
        score_sigma^2,
        curvature_floor
      ),

    nn =
      -pmax(
        score_nu^2,
        curvature_floor
      ),

    ms =
      -score_mu *
      score_sigma,

    mn =
      -score_mu *
      score_nu,

    sn =
      -score_sigma *
      score_nu
  )
}


# Classic-GAMLSS LL3 family

LL3 <- function(
    lower,
    sigma.link = "log",
    nu.link = "log",
    information = c(
      "opg",
      "observed"
    )
) {
  information <- match.arg(
    information
  )

  .LL3_check_family_dependencies()

  if (
    !is.numeric(lower) ||
      length(lower) != 1L ||
      !is.finite(lower)
  ) {
    stop(
      "lower must be one finite numeric value."
    )
  }

  lower <- unname(
    as.numeric(lower)
  )

  if (
    !requireNamespace(
      "gamlss.dist",
      quietly = TRUE
    )
  ) {
    stop(
      "Install the gamlss.dist package before using LL3()."
    )
  }


  # -----------------------------------------------------------------------
  # Standard sigma and nu links
  # -----------------------------------------------------------------------

  checklink <- gamlss.dist::checklink

  sigma_link <- checklink(
    "sigma.link",
    "Three-parameter shifted log-logistic",
    substitute(sigma.link),
    c(
      "log",
      "identity",
      "inverse",
      "own"
    )
  )

  nu_link <- checklink(
    "nu.link",
    "Three-parameter shifted log-logistic",
    substitute(nu.link),
    c(
      "log",
      "identity",
      "inverse",
      "own"
    )
  )


  # -----------------------------------------------------------------------
  # Support-safe mu link
  #
  #   eta_mu = log(lower - mu)
  #   mu     = lower - exp(eta_mu)
  #
  # The forward link does not stop when mu touches or very slightly exceeds
  # lower during internal GAMLSS calculations. Instead, the gap is clipped
  # to a small positive tolerance. The inverse link remains support-safe.
  # -----------------------------------------------------------------------

  mu_link_tolerance <- max(
    100 *
      .Machine$double.eps *
      max(
        1,
        abs(lower)
      ),
    .Machine$double.xmin
  )

  mu_initial_gap <- max(
    sqrt(.Machine$double.eps) *
      max(
        1,
        abs(lower)
      ),
    1000 * mu_link_tolerance,
    1e-10
  )


  mu_linkfun <- eval(
    substitute(
      function(mu) {
        if (any(!is.finite(mu))) {
          stop(
            "Every mu value supplied to the LL3 link must be finite."
          )
        }

        gap <- L - mu

        if (any(mu > L + T)) {
          stop(
            "mu exceeds the LL3 upper boundary by more than numerical ",
            "tolerance. Every mu value must be less than or equal to lower."
          )
        }

        log(
          pmax(
            gap,
            T
          )
        )
      },
      list(
        L = lower,
        T = mu_link_tolerance
      )
    )
  )


  mu_linkinv <- eval(
    substitute(
      function(eta) {
        L - exp(eta)
      },
      list(
        L = lower
      )
    )
  )


  mu_mueta <- function(eta) {
    -exp(eta)
  }


  mu_valid <- eval(
    substitute(
      function(mu) {
        all(
          is.finite(mu) &
            mu <= L + T
        )
      },
      list(
        L = lower,
        T = mu_link_tolerance
      )
    )
  )


  y_valid <- eval(
    substitute(
      function(y) {
        is.numeric(y) &&
          all(is.finite(y)) &&
          all(y > L)
      },
      list(
        L = lower
      )
    )
  )


  # -----------------------------------------------------------------------
  # Initial values
  # -----------------------------------------------------------------------

  mu_initial <- as.expression(
    substitute(
      {
        start <- base::get(
          ".LL3_start_values",
          envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
            base::asNamespace("LL3GAMLSS")
          } else {
            base::globalenv()
          },
          inherits = TRUE
        )(y)

        candidate <- unname(
          start["mu"]
        )

        if (!is.finite(candidate)) {
          candidate <- L - G
        }

        candidate <- min(
          candidate,
          L - G
        )

        mu <- rep(
          candidate,
          length(y)
        )
      },
      list(
        L = lower,
        G = mu_initial_gap
      )
    )
  )


  sigma_initial <- expression({
    start <- base::get(
      ".LL3_start_values",
      envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
        base::asNamespace("LL3GAMLSS")
      } else {
        base::globalenv()
      },
      inherits = TRUE
    )(y)

    sigma <- rep(
      unname(
        start["sigma"]
      ),
      length(y)
    )
  })


  nu_initial <- expression({
    start <- base::get(
      ".LL3_start_values",
      envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
        base::asNamespace("LL3GAMLSS")
      } else {
        base::globalenv()
      },
      inherits = TRUE
    )(y)

    nu <- rep(
      unname(
        start["nu"]
      ),
      length(y)
    )
  })


  # -----------------------------------------------------------------------
  # Curvature functions
  #
  # Functions are defined explicitly so they remain available after the
  # family object is processed by GAMLSS.
  # -----------------------------------------------------------------------

  if (identical(
    information,
    "opg"
  )) {
    d2_mu <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_opg_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$mm
    }


    d2_sigma <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_opg_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$ss
    }


    d2_nu <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_opg_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$nn
    }


    d2_mu_sigma <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_opg_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$ms
    }


    d2_mu_nu <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_opg_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$mn
    }


    d2_sigma_nu <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_opg_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$sn
    }
  } else {
    d2_mu <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_observed_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$mm
    }


    d2_sigma <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_observed_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$ss
    }


    d2_nu <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_observed_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$nn
    }


    d2_mu_sigma <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_observed_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$ms
    }


    d2_mu_nu <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_observed_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$mn
    }


    d2_sigma_nu <- function(
        y,
        mu,
        sigma,
        nu
    ) {
      base::get(
        ".LL3_observed_hessian",
        envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
          base::asNamespace("LL3GAMLSS")
        } else {
          base::globalenv()
        },
        inherits = TRUE
      )(
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      )$sn
    }
  }


  # -----------------------------------------------------------------------
  # Construct family object
  # -----------------------------------------------------------------------

  rqres_pfun <- if (base::isNamespaceLoaded("LL3GAMLSS")) {
    "LL3GAMLSS::pLL3"
  } else {
    "pLL3"
  }

  rqres_expression <- as.expression(
    substitute(
      rqres(
        pfun = PFUN,
        type = "Continuous",
        y = y,
        mu = mu,
        sigma = sigma,
        nu = nu
      ),
      list(PFUN = rqres_pfun)
    )
  )

  structure(
    list(
      family = c(
        "LL3",
        "Support-safe three-parameter shifted log-logistic"
      ),

      parameters = list(
        mu = TRUE,
        sigma = TRUE,
        nu = TRUE
      ),

      nopar = 3,

      type = "Continuous",


      # -------------------------------------------------------------------
      # Link descriptions
      # -------------------------------------------------------------------

      mu.link = paste0(
        "below(",
        format(
          lower,
          digits = 16
        ),
        ")"
      ),

      sigma.link = as.character(
        substitute(
          sigma.link
        )
      ),

      nu.link = as.character(
        substitute(
          nu.link
        )
      ),


      # -------------------------------------------------------------------
      # Link functions
      # -------------------------------------------------------------------

      mu.linkfun =
        mu_linkfun,

      sigma.linkfun =
        sigma_link$linkfun,

      nu.linkfun =
        nu_link$linkfun,


      mu.linkinv =
        mu_linkinv,

      sigma.linkinv =
        sigma_link$linkinv,

      nu.linkinv =
        nu_link$linkinv,


      mu.dr =
        mu_mueta,

      sigma.dr =
        sigma_link$mu.eta,

      nu.dr =
        nu_link$mu.eta,


      # -------------------------------------------------------------------
      # First derivatives of log likelihood
      # -------------------------------------------------------------------

      dldm = function(
          y,
          mu,
          sigma,
          nu
      ) {
        base::get(
          ".LL3_score",
          envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
            base::asNamespace("LL3GAMLSS")
          } else {
            base::globalenv()
          },
          inherits = TRUE
        )(
          y = y,
          mu = mu,
          sigma = sigma,
          nu = nu
        )[, "mu"]
      },


      dldd = function(
          y,
          mu,
          sigma,
          nu
      ) {
        base::get(
          ".LL3_score",
          envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
            base::asNamespace("LL3GAMLSS")
          } else {
            base::globalenv()
          },
          inherits = TRUE
        )(
          y = y,
          mu = mu,
          sigma = sigma,
          nu = nu
        )[, "sigma"]
      },


      dldv = function(
          y,
          mu,
          sigma,
          nu
      ) {
        base::get(
          ".LL3_score",
          envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
            base::asNamespace("LL3GAMLSS")
          } else {
            base::globalenv()
          },
          inherits = TRUE
        )(
          y = y,
          mu = mu,
          sigma = sigma,
          nu = nu
        )[, "nu"]
      },


      # -------------------------------------------------------------------
      # Second and cross derivatives
      # -------------------------------------------------------------------

      d2ldm2 =
        d2_mu,

      d2ldd2 =
        d2_sigma,

      d2ldv2 =
        d2_nu,

      d2ldmdd =
        d2_mu_sigma,

      d2ldmdv =
        d2_mu_nu,

      d2ldddv =
        d2_sigma_nu,


      # -------------------------------------------------------------------
      # Deviance contribution
      # -------------------------------------------------------------------

      G.dev.incr = function(
          y,
          mu,
          sigma,
          nu,
          ...
      ) {
        deviance_increment <- -2 *
          base::get(
            "dLL3",
            envir = if (base::isNamespaceLoaded("LL3GAMLSS")) {
              base::asNamespace("LL3GAMLSS")
            } else {
              base::globalenv()
            },
            inherits = TRUE
          )(
            y,
            mu = mu,
            sigma = sigma,
            nu = nu,
            log = TRUE
          )

        # A finite penalty lets GAMLSS shorten or reject an invalid step
        # rather than terminating because of an infinite global deviance.
        deviance_increment[
          !is.finite(
            deviance_increment
          )
        ] <- 1e12

        deviance_increment
      },


      # -------------------------------------------------------------------
      # Randomized quantile residuals
      # -------------------------------------------------------------------

      rqres = rqres_expression,


      # -------------------------------------------------------------------
      # Starting values
      # -------------------------------------------------------------------

      mu.initial =
        mu_initial,

      sigma.initial =
        sigma_initial,

      nu.initial =
        nu_initial,


      # -------------------------------------------------------------------
      # Parameter and response validity
      # -------------------------------------------------------------------

      mu.valid =
        mu_valid,

      sigma.valid = function(sigma) {
        all(
          is.finite(sigma) &
            sigma > 0
        )
      },

      nu.valid = function(nu) {
        all(
          is.finite(nu) &
            nu > 0
        )
      },

      y.valid =
        y_valid,


      # -------------------------------------------------------------------
      # Distribution moments
      # -------------------------------------------------------------------

      mean = function(
          mu,
          sigma,
          nu
      ) {
        angle <- pi / nu

        ifelse(
          nu > 1,
          mu +
            sigma *
            angle /
            sin(angle),
          Inf
        )
      },


      variance = function(
          mu,
          sigma,
          nu
      ) {
        angle <- pi / nu

        ifelse(
          nu > 2,
          sigma^2 *
            angle *
            (
              2 /
                sin(
                  2 * angle
                ) -
                angle /
                sin(angle)^2
            ),
          Inf
        )
      },


      median = function(
          mu,
          sigma,
          nu
      ) {
        mu + sigma
      },


      # -------------------------------------------------------------------
      # LL3-specific metadata
      # -------------------------------------------------------------------

      lower =
        lower,

      mu_upper_boundary =
        lower,

      mu_link_tolerance =
        mu_link_tolerance,

      mu_initial_gap =
        mu_initial_gap,

      information =
        information
    ),

    class = c(
      "gamlss.family",
      "family"
    )
  )
}


# Backward-compatible alias.
LL3safe <- LL3
