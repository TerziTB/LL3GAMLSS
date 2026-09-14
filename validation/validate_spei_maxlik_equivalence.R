project_dir <- local({
  source(file.path("validation", "validation_helpers.R"), local = TRUE)
  resolve_project_dir()
})
setwd(project_dir)
source(file.path("validation", "validation_helpers.R"))
load_LL3_source(project_dir)

if (!requireNamespace("SPEI", quietly = TRUE) || !requireNamespace("lmom", quietly = TRUE)) {
  stop("SPEI and lmom are required for this validation.")
}

output_dir <- file.path(project_dir, "validation", "results", "spei_maxlik_equivalence")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
control <- empirical_LL3_control()

glo_to_LL3 <- function(parameters) {
  xi <- unname(parameters[[1L]])
  alpha <- unname(parameters[[2L]])
  kappa <- unname(parameters[[3L]])
  if (!is.finite(kappa) || kappa >= 0) {
    return(c(mu = NA_real_, sigma = NA_real_, nu = NA_real_))
  }
  c(mu = xi + alpha / kappa, sigma = -alpha / kappa, nu = -1 / kappa)
}

spei_internal_maxlik <- function(y) {
  namespace <- asNamespace("SPEI")
  pwm <- get("PWM", namespace)(y, order = 0:2)
  lmoments <- get("pwm2lmom", namespace)(pwm)
  fortran <- c(lmoments$lambdas[1:2], lmoments$ratios[3])
  initial <- tryCatch(
    get("pelglo", namespace)(fortran),
    error = function(e) get("parglo", namespace)(lmoments)$para
  )
  get("parglo.maxlik", namespace)(y, initial)
}

LL3_to_glo <- function(mu, sigma, nu) {
  c(xi = mu + sigma, alpha = sigma / nu, kappa = -1 / nu)
}

fit_and_compare_month <- function(
    y,
    spei_parameters,
    spei_index,
    source,
    case_id,
    month,
    scale_months,
    direct_call_status
) {
  fit_data <- data.frame(response = y)
  fit <- suppressWarnings(fit_LL3_gamlss(
    response ~ 1,
    sigma.formula = ~1,
    nu.formula = ~1,
    data = fit_data,
    control = control,
    boundary_action = "error"
  ))
  estimated <- extract_LL3_parameters(fit)[1L, ]
  mapped <- glo_to_LL3(spei_parameters)
  polished_fit <- get("parglo.maxlik", asNamespace("SPEI"))(
    y,
    LL3_to_glo(estimated$mu, estimated$sigma, estimated$nu)
  )
  polished_parameters <- polished_fit$para
  probability_gamlss <- pLL3(y, estimated$mu, estimated$sigma, estimated$nu)
  probability_spei <- lmom::cdfglo(y, spei_parameters)
  probability_polished <- lmom::cdfglo(y, polished_parameters)
  index_gamlss <- stats::qnorm(pmin(pmax(probability_gamlss, 1e-10), 1 - 1e-10))
  index_polished <- stats::qnorm(pmin(pmax(probability_polished, 1e-10), 1 - 1e-10))
  keep <- is.finite(index_gamlss) & is.finite(spei_index)
  data.frame(
    source = source,
    case_id = case_id,
    month = month,
    scale_months = scale_months,
    n = length(y),
    SPEI_xi = unname(spei_parameters[[1L]]),
    SPEI_alpha = unname(spei_parameters[[2L]]),
    SPEI_kappa = unname(spei_parameters[[3L]]),
    LL3_compatible_kappa = is.finite(spei_parameters[[3L]]) && spei_parameters[[3L]] < 0,
    GAMLSS_mu = estimated$mu,
    SPEI_mapped_mu = mapped[["mu"]],
    GAMLSS_sigma = estimated$sigma,
    SPEI_mapped_sigma = mapped[["sigma"]],
    GAMLSS_nu = estimated$nu,
    SPEI_mapped_nu = mapped[["nu"]],
    max_absolute_probability_difference = max(abs(probability_gamlss - probability_spei)),
    SPEI_MAE = mean(abs(index_gamlss[keep] - spei_index[keep])),
    SPEI_max_absolute_difference = max(abs(index_gamlss[keep] - spei_index[keep])),
    SPEI_correlation = stats::cor(index_gamlss[keep], spei_index[keep]),
    polished_max_absolute_probability_difference = max(abs(probability_gamlss - probability_polished)),
    polished_SPEI_MAE = mean(abs(index_gamlss - index_polished)),
    polished_SPEI_max_absolute_difference = max(abs(index_gamlss - index_polished)),
    polished_SPEI_correlation = stats::cor(index_gamlss, index_polished),
    polished_SPEI_convergence = polished_fit$conv,
    SPEI_direct_maxlik_call_status = direct_call_status,
    stringsAsFactors = FALSE
  )
}

compare_series <- function(balance, start_year, case_id, source, scale_months = 1L, months = seq_len(12L)) {
  series <- stats::ts(balance, start = c(start_year, 1), frequency = 12)
  reference <- suppressWarnings(SPEI::spei(
    series,
    scale = scale_months,
    distribution = "log-Logistic",
    fit = "max-lik",
    na.rm = TRUE,
    verbose = FALSE
  ))
  accumulated <- LL3_accumulate(balance, scale_months)
  calendar_month <- rep(seq_len(12L), length.out = length(balance))
  direct_parameters_valid <- all(is.finite(reference$coefficients))
  do.call(rbind, lapply(months, function(month) {
    rows <- calendar_month == month & is.finite(accumulated)
    parameters <- if (direct_parameters_valid) {
      reference$coefficients[, 1L, month]
    } else {
      spei_internal_maxlik(accumulated[rows])$para
    }
    index <- if (direct_parameters_valid) {
      as.numeric(reference$fitted)[rows]
    } else {
      stats::qnorm(lmom::cdfglo(accumulated[rows], parameters))
    }
    fit_and_compare_month(
      accumulated[rows],
      parameters,
      index,
      source,
      case_id,
      month,
      scale_months,
      if (direct_parameters_valid) "returned_finite_parameters" else "returned_NA_used_internal_engine"
    )
  }))
}

simulated_results <- list()
for (case in seq_len(3L)) {
  set.seed(951000L + case)
  years <- 66L
  month <- rep(seq_len(12L), years)
  angle <- 2 * pi * (month - 1) / 12
  mu <- -130 + 20 * cos(angle) + 5 * case
  sigma <- (60 + 10 * case) * exp(0.18 * sin(angle))
  nu <- (1.7 + 0.15 * case) * exp(0.06 * cos(angle))
  balance <- rLL3(length(month), mu, sigma, nu)
  simulated_results[[case]] <- compare_series(
    balance,
    1960,
    sprintf("simulated_%d", case),
    "simulated",
    scale_months = 1L
  )
}

empirical_results <- list()
empirical_path <- file.path(project_dir, "analysis", "seyhan_results", "spei_values.csv")
if (file.exists(empirical_path)) {
  empirical <- utils::read.csv(empirical_path, stringsAsFactors = FALSE)
  empirical$date <- as.Date(empirical$date)
  cases <- data.frame(
    station = c("17351", "17802", "17840"),
    scale_months = c(1L, 6L, 12L),
    month = c(1L, 7L, 12L),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(cases))) {
    subset <- empirical[
      as.character(empirical$station) == cases$station[i] &
        empirical$scale_months == cases$scale_months[i],
      , drop = FALSE
    ]
    subset <- subset[order(subset$date), ]
    if (!nrow(subset)) next
    empirical_results[[i]] <- compare_series(
      subset$water_balance_mm,
      as.integer(format(min(subset$date), "%Y")),
      sprintf("station_%s_scale_%d", cases$station[i], cases$scale_months[i]),
      "empirical",
      scale_months = cases$scale_months[i],
      months = cases$month[i]
    )
  }
}

raw <- do.call(rbind, c(simulated_results, empirical_results))
summary_groups <- split(raw, interaction(raw$source, raw$LL3_compatible_kappa, drop = TRUE))
summary <- do.call(rbind, lapply(summary_groups, function(d) {
  data.frame(
    source = unique(d$source),
    LL3_compatible_kappa = unique(d$LL3_compatible_kappa),
    comparisons = nrow(d),
    maximum_probability_difference = max(d$max_absolute_probability_difference),
    median_SPEI_MAE = stats::median(d$SPEI_MAE),
    maximum_SPEI_difference = max(d$SPEI_max_absolute_difference),
    minimum_SPEI_correlation = min(d$SPEI_correlation),
    polished_maximum_probability_difference = max(d$polished_max_absolute_probability_difference),
    polished_median_SPEI_MAE = stats::median(d$polished_SPEI_MAE),
    polished_maximum_SPEI_difference = max(d$polished_SPEI_max_absolute_difference),
    polished_minimum_SPEI_correlation = min(d$polished_SPEI_correlation),
    polished_convergence_rate = mean(d$polished_SPEI_convergence == 0),
    direct_call_finite_rate = mean(d$SPEI_direct_maxlik_call_status == "returned_finite_parameters"),
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(raw, file.path(output_dir, "spei_maxlik_equivalence_raw.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(output_dir, "spei_maxlik_equivalence_summary.csv"), row.names = FALSE)
writeLines(c(
  "SPEI generalized-logistic to LL3 mapping:",
  "nu=-1/kappa; sigma=-alpha/kappa; mu=xi+alpha/kappa, requiring kappa<0.",
  "Installed SPEI 1.8.1 direct fit=max-lik calls returned NA coefficients because the fitting-method initialization switch omitted max-lik. The benchmark therefore called SPEI's own parglo.maxlik engine using its unbiased-PWM generalized-logistic start and cdfglo standardization.",
  "A second SPEI-engine run started at the parameter-equivalent GAMLSS solution to distinguish objective/parameterization equivalence from the default SPEI optimizer's initialization and local convergence.",
  "The SPEI generalized-logistic optimizer does not constrain kappa below zero. Rows with kappa>=0 have an upper rather than lower finite support endpoint and are flagged as not parameter-equivalent to the lower-threshold LL3 family.",
  "The benchmark compares fitted probabilities and standardized values. M0 remains the stationary comparator in model selection; this script is implementation validation only.",
  if (file.exists(empirical_path)) {
    "Three representative empirical station-scale-month series were included."
  } else {
    "Empirical input was unavailable; simulated validation completed and the empirical benchmark was skipped."
  }
), file.path(output_dir, "spei_maxlik_equivalence_settings.txt"))
write_session_info(file.path(output_dir, "sessionInfo.txt"))
print(summary)
