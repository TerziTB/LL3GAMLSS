# Operational helpers for stationary and nonstationary LL3-based SPEI.

LL3_accumulate <- function(x, scale = 1L, weights = NULL) {
  x <- as.numeric(x)
  scale <- as.integer(scale)

  if (!is.finite(scale) || scale < 1L) {
    stop("scale must be a positive integer.")
  }

  if (is.null(weights)) {
    weights <- rep(1, scale)
  }
  weights <- as.numeric(weights)

  if (length(weights) != scale || any(!is.finite(weights))) {
    stop("weights must contain one finite value for every accumulation lag.")
  }
  if (all(weights == 0)) {
    stop("At least one accumulation weight must be non-zero.")
  }

  result <- rep(NA_real_, length(x))
  if (length(x) < scale) {
    return(result)
  }

  for (i in seq.int(scale, length(x))) {
    window <- x[seq.int(i - scale + 1L, i)]
    if (all(is.finite(window))) {
      result[i] <- sum(window * weights)
    }
  }

  result
}


.LL3_as_date <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }

  converted <- as.Date(x)
  if (anyNA(converted)) {
    stop("date values must be Date objects or unambiguous ISO date strings.")
  }
  converted
}


.LL3_predict_parameter <- function(object, what, newdata) {
  smooth_matrix <- object[[paste0(what, ".s")]]
  if (!is.null(smooth_matrix)) {
    stop(
      "fit_LL3_spei() currently supports parametric predictor terms only; ",
      "smooth-term prediction requires an explicit prediction method."
    )
  }

  parameter_formula <- stats::formula(object, what)
  if (length(parameter_formula) == 3L) {
    parameter_formula[[2L]] <- NULL
  }

  parameter_terms <- stats::terms(parameter_formula)
  model_frame <- stats::model.frame(
    parameter_terms,
    newdata,
    na.action = stats::na.fail,
    xlev = object[[paste0(what, ".xlevels")]]
  )
  design <- stats::model.matrix(
    parameter_terms,
    model_frame,
    contrasts.arg = object$contrasts
  )
  coefficients <- stats::coef(object, what = what)

  missing_columns <- setdiff(names(coefficients), colnames(design))
  if (length(missing_columns)) {
    stop(
      "Prediction design matrix is missing coefficient columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  eta <- drop(
    design[, names(coefficients), drop = FALSE] %*%
      unname(coefficients)
  )

  family_object <- attr(object, "LL3_family_object")
  if (is.null(family_object)) {
    family_object <- LL3(
      lower = attr(object, "LL3_lower")
    )
  }

  inverse_link <- family_object[[paste0(what, ".linkinv")]]
  as.numeric(inverse_link(eta))
}


fit_LL3_spei <- function(
    data,
    response,
    date,
    scale = 1L,
    weights = NULL,
    sigma.formula = ~1,
    nu.formula = ~1,
    reference_start = NULL,
    reference_end = NULL,
    min_observations = 30L,
    lower = NULL,
    epsilon = 1e-10,
    information = "opg",
    control = NULL,
    boundary_action = c("error", "warning")
) {
  boundary_action <- match.arg(boundary_action)

  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }
  if (!is.character(response) || length(response) != 1L ||
      !response %in% names(data)) {
    stop("response must name one column in data.")
  }
  if (!is.character(date) || length(date) != 1L || !date %in% names(data)) {
    stop("date must name one column in data.")
  }
  if (!is.numeric(data[[response]]) || any(!is.finite(data[[response]]))) {
    stop("The response column must contain only finite numeric values.")
  }

  dates <- .LL3_as_date(data[[date]])
  if (is.unsorted(dates, strictly = TRUE)) {
    stop("date values must be unique and strictly increasing.")
  }

  min_observations <- as.integer(min_observations)
  if (!is.finite(min_observations) || min_observations < 4L) {
    stop("min_observations must be an integer of at least four.")
  }
  if (!is.finite(epsilon) || epsilon <= 0 || epsilon >= 0.5) {
    stop("epsilon must lie strictly between zero and 0.5.")
  }

  accumulated <- LL3_accumulate(
    data[[response]],
    scale = scale,
    weights = weights
  )
  month <- as.POSIXlt(dates, tz = "UTC")$mon + 1L

  reference <- rep(TRUE, nrow(data))
  if (!is.null(reference_start)) {
    reference <- reference & dates >= .LL3_as_date(reference_start)
  }
  if (!is.null(reference_end)) {
    reference <- reference & dates <= .LL3_as_date(reference_end)
  }

  if (!any(reference & is.finite(accumulated))) {
    stop("The reference period contains no complete accumulated observations.")
  }

  model_data <- data
  model_data$.LL3_response <- accumulated
  model_data$.LL3_month <- month
  model_data$.LL3_reference <- reference

  fitted_values <- data.frame(
    mu = rep(NA_real_, nrow(data)),
    sigma = rep(NA_real_, nrow(data)),
    nu = rep(NA_real_, nrow(data))
  )
  fits <- list()

  represented_months <- sort(unique(month[is.finite(accumulated)]))

  resolve_lower <- function(current_month, y) {
    if (is.null(lower)) {
      return(LL3_support_lower(y))
    }
    if (length(lower) == 1L) {
      return(as.numeric(lower))
    }
    if (!is.null(names(lower)) && as.character(current_month) %in% names(lower)) {
      return(as.numeric(lower[as.character(current_month)]))
    }
    if (length(lower) == 12L) {
      return(as.numeric(lower[current_month]))
    }
    stop("lower must be NULL, one value, or twelve month-specific values.")
  }

  for (current_month in represented_months) {
    fit_rows <-
      month == current_month &
      reference &
      is.finite(accumulated)
    prediction_rows <-
      month == current_month &
      is.finite(accumulated)

    if (sum(fit_rows) < min_observations) {
      stop(
        "Month ", current_month, " has only ", sum(fit_rows),
        " usable reference observations; at least ", min_observations,
        " are required."
      )
    }

    fit_data <- model_data[fit_rows, , drop = FALSE]
    prediction_data <- model_data[prediction_rows, , drop = FALSE]
    current_lower <- resolve_lower(
      current_month,
      fit_data$.LL3_response
    )

    current_fit <- fit_LL3_gamlss(
      .LL3_response ~ 1,
      sigma.formula = sigma.formula,
      nu.formula = nu.formula,
      data = fit_data,
      lower = current_lower,
      information = information,
      allow_nonstationary_mu = FALSE,
      control = control,
      boundary_action = boundary_action
    )

    fits[[as.character(current_month)]] <- current_fit
    fitted_values$mu[prediction_rows] <- .LL3_predict_parameter(
      current_fit,
      "mu",
      prediction_data
    )
    fitted_values$sigma[prediction_rows] <- .LL3_predict_parameter(
      current_fit,
      "sigma",
      prediction_data
    )
    fitted_values$nu[prediction_rows] <- .LL3_predict_parameter(
      current_fit,
      "nu",
      prediction_data
    )
  }

  complete <-
    is.finite(accumulated) &
    is.finite(fitted_values$mu) &
    is.finite(fitted_values$sigma) &
    is.finite(fitted_values$nu)

  probability <- rep(NA_real_, nrow(data))
  probability[complete] <- pLL3(
    accumulated[complete],
    mu = fitted_values$mu[complete],
    sigma = fitted_values$sigma[complete],
    nu = fitted_values$nu[complete]
  )

  outside_support <- complete & accumulated <= fitted_values$mu
  if (any(outside_support)) {
    warning(
      sum(outside_support),
      " accumulated observations lie outside the fitted reference-period support; ",
      "their probabilities are clipped before Gaussian transformation.",
      call. = FALSE
    )
  }

  index <- rep(NA_real_, nrow(data))
  index[complete] <- stats::qnorm(
    pmin(pmax(probability[complete], epsilon), 1 - epsilon)
  )

  result <- list(
    index = index,
    results = data.frame(
      date = dates,
      month = month,
      accumulated_balance = accumulated,
      mu = fitted_values$mu,
      sigma = fitted_values$sigma,
      nu = fitted_values$nu,
      probability = probability,
      LL3_SPEI = index,
      reference = reference
    ),
    fits = fits,
    call = match.call(),
    settings = list(
      response = response,
      date = date,
      scale = as.integer(scale),
      weights = if (is.null(weights)) rep(1, as.integer(scale)) else weights,
      reference_start = reference_start,
      reference_end = reference_end,
      sigma.formula = sigma.formula,
      nu.formula = nu.formula,
      epsilon = epsilon
    )
  )

  class(result) <- "LL3_spei"
  result
}


print.LL3_spei <- function(x, ...) {
  cat("LL3-based monthly SPEI fit\n")
  cat("  accumulation scale:", x$settings$scale, "month(s)\n")
  cat("  fitted calendar months:", paste(names(x$fits), collapse = ", "), "\n")
  cat("  finite index values:", sum(is.finite(x$index)), "of", length(x$index), "\n")
  invisible(x)
}
