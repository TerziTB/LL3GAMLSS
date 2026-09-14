#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args

input_json <- file.path("artifact_work", "inspection", "workbook_values.json")
output_dir <- file.path("analysis", "seyhan_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_json)) {
  stop(
    "The access-controlled station input is not included in the public repository. ",
    "Create artifact_work/inspection/workbook_values.json from the provider workbook ",
    "using the documented Date, Precipitation, Average_Temp, Min_Temp, and Max_Temp columns."
  )
}

.libPaths(c(file.path(getwd(), ".rlib-realdata"), .libPaths()))
library(LL3GAMLSS)

station_metadata <- data.frame(
  station = c("17351", "17802", "17837", "17840", "17906", "17934", "17936", "17981"),
  station_name = c("Adana", "Pinarbasi", "Tomarza", "Sariz", "Ulukisla", "Pozanti", "Karaisali", "Karatas"),
  latitude_deg_n = c(37.0041, 38.7251, 38.4522, 38.4781, 37.5480, 37.4758, 37.2505, 36.5683),
  longitude_deg_e = c(35.3443, 36.3904, 35.7912, 36.5035, 34.4867, 34.9022, 35.0628, 35.3894),
  elevation_m = c(23, 1542, 1402, 1599, 1453, 1080, 240, 22),
  stringsAsFactors = FALSE
)

raw <- jsonlite::fromJSON(input_json, simplifyVector = FALSE)

excel_date <- function(x) {
  as.Date(as.numeric(x), origin = "1899-12-30")
}

spei_hargreaves_pet_monthly <- function(date, tmin, tmax, latitude_deg) {
  if (!requireNamespace("SPEI", quietly = TRUE)) {
    stop("Install the SPEI package before running the real-data analysis.")
  }
  start <- c(
    as.integer(format(min(date), "%Y")),
    as.integer(format(min(date), "%m"))
  )
  tmin_ts <- stats::ts(tmin, start = start, frequency = 12)
  tmax_ts <- stats::ts(tmax, start = start, frequency = 12)
  as.numeric(SPEI::hargreaves(
    Tmin = tmin_ts,
    Tmax = tmax_ts,
    lat = latitude_deg,
    na.rm = FALSE,
    verbose = FALSE
  ))
}

parse_station <- function(station) {
  rows <- raw[[station]]
  header <- unlist(rows[[1]], use.names = FALSE)
  values <- do.call(rbind, lapply(rows[-1], function(x) unlist(x, use.names = FALSE)))
  values <- as.data.frame(values, stringsAsFactors = FALSE)
  names(values) <- header
  for (nm in setdiff(names(values), "Date")) {
    values[[nm]] <- as.numeric(values[[nm]])
  }
  values$Date <- excel_date(values$Date)
  values$station <- station
  values
}

station_data <- lapply(station_metadata$station, parse_station)
names(station_data) <- station_metadata$station

screen_annual_repetition <- function(d, exclusion_threshold_months = 12L) {
  years <- as.integer(format(d$Date, "%Y"))
  months <- as.integer(format(d$Date, "%m"))
  prior_dates <- as.Date(sprintf("%04d-%02d-01", years - 1L, months))
  prior_index <- match(prior_dates, d$Date)
  comparable <- !is.na(prior_index)
  same_as_prior_year <- rep(FALSE, nrow(d))
  same_as_prior_year[comparable] <-
    d$Precipitation[comparable] == d$Precipitation[prior_index[comparable]] &
    d$Min_Temp[comparable] == d$Min_Temp[prior_index[comparable]] &
    d$Max_Temp[comparable] == d$Max_Temp[prior_index[comparable]]

  runs <- rle(same_as_prior_year)
  true_runs <- which(runs$values)
  if (!length(true_runs)) {
    return(data.frame(
      max_consecutive_prior_year_matches = 0L,
      repeated_run_start = as.Date(NA),
      repeated_run_end = as.Date(NA),
      inferred_fill_start = as.Date(NA),
      excluded = FALSE
    ))
  }

  longest_run <- true_runs[which.max(runs$lengths[true_runs])]
  run_end <- cumsum(runs$lengths)[longest_run]
  run_start <- run_end - runs$lengths[longest_run] + 1L
  inferred_start <- as.Date(sprintf(
    "%04d-%02d-01",
    as.integer(format(d$Date[run_start], "%Y")) - 1L,
    as.integer(format(d$Date[run_start], "%m"))
  ))
  data.frame(
    max_consecutive_prior_year_matches = runs$lengths[longest_run],
    repeated_run_start = d$Date[run_start],
    repeated_run_end = d$Date[run_end],
    inferred_fill_start = inferred_start,
    excluded = runs$lengths[longest_run] >= exclusion_threshold_months
  )
}

repetition_screen <- do.call(rbind, lapply(names(station_data), function(station) {
  result <- screen_annual_repetition(station_data[[station]])
  cbind(station = station, result)
}))
repetition_screen$reason <- ifelse(
  repetition_screen$excluded,
  "At least 12 consecutive precipitation-Tmin-Tmax vectors exactly matched the preceding year",
  "Retained"
)
write.csv(
  repetition_screen,
  file.path(output_dir, "annual_repetition_screen.csv"),
  row.names = FALSE
)

excluded_stations <- repetition_screen$station[repetition_screen$excluded]
station_metadata <- station_metadata[!station_metadata$station %in% excluded_stations, , drop = FALSE]
station_data <- station_data[station_metadata$station]

quality_rows <- lapply(station_metadata$station, function(station) {
  d <- station_data[[station]]
  expected <- seq(min(d$Date), max(d$Date), by = "month")
  data.frame(
    station = station,
    first_date = min(d$Date),
    last_date = max(d$Date),
    observations = nrow(d),
    expected_months = length(expected),
    missing_months = sum(!expected %in% d$Date),
    duplicate_dates = sum(duplicated(d$Date)),
    missing_measurements = sum(!stats::complete.cases(d[c("Precipitation", "Average_Temp", "Min_Temp", "Max_Temp")])),
    negative_precipitation = sum(d$Precipitation < 0, na.rm = TRUE),
    min_above_max = sum(d$Min_Temp > d$Max_Temp, na.rm = TRUE),
    mean_outside_min_max = sum(
      d$Average_Temp < d$Min_Temp | d$Average_Temp > d$Max_Temp,
      na.rm = TRUE
    )
  )
})
quality <- do.call(rbind, quality_rows)

common_start <- max(vapply(station_data, function(x) as.numeric(min(x$Date)), numeric(1)))
common_end <- min(vapply(station_data, function(x) as.numeric(max(x$Date)), numeric(1)))
common_start <- as.Date(common_start, origin = "1970-01-01")
common_end <- as.Date(common_end, origin = "1970-01-01")

for (station in station_metadata$station) {
  d <- station_data[[station]]
  meta <- station_metadata[station_metadata$station == station, ]
  d$PET_Hargreaves_mm <- spei_hargreaves_pet_monthly(
    d$Date,
    d$Min_Temp,
    d$Max_Temp,
    meta$latitude_deg_n
  )
  d$water_balance_mm <- d$Precipitation - d$PET_Hargreaves_mm
  d <- d[d$Date >= common_start & d$Date <= common_end, , drop = FALSE]
  d$time_scaled <- as.numeric(scale(seq_len(nrow(d))))
  station_data[[station]] <- d
}

quality$analysis_common_start <- common_start
quality$analysis_common_end <- common_end
quality$analysis_common_months <- as.integer(
  (as.integer(format(common_end, "%Y")) - as.integer(format(common_start, "%Y"))) * 12 +
    as.integer(format(common_end, "%m")) - as.integer(format(common_start, "%m")) + 1
)

write.csv(station_metadata, file.path(output_dir, "station_metadata.csv"), row.names = FALSE)
write.csv(quality, file.path(output_dir, "data_quality.csv"), row.names = FALSE)

if (quick_mode) {
  stations_to_run <- "17351"
  scales_to_run <- 12L
} else {
  stations_to_run <- station_metadata$station
  scales_to_run <- c(1L, 3L, 6L, 12L)
}

model_names <- c("stationary", "sigma_time", "nu_time", "sigma_nu_time")
model_metrics <- list()
all_values <- list()
progress_file <- file.path(output_dir, "progress.log")
writeLines(
  paste("Started", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "quick_mode=", quick_mode),
  progress_file
)

fit_metrics <- function(fits, model_name, station, scale) {
  ok <- vapply(fits, function(x) inherits(x, "gamlss"), logical(1))
  valid_fits <- fits[ok]
  if (!length(valid_fits)) {
    return(data.frame(
      station = station, scale_months = scale, model = model_name,
      fitted_months = 0L, converged_months = 0L, inference_ready_months = 0L,
      boundary_contact_months = NA_integer_, n = NA_integer_, df = NA_real_,
      deviance = NA_real_, AIC = NA_real_, AICc = NA_real_, BIC = NA_real_,
      min_lower_minus_mu = NA_real_
    ))
  }
  diagnostics <- lapply(valid_fits, function(fit) {
    tryCatch(LL3_boundary_diagnostic(fit), error = function(e) NULL)
  })
  inference_ready <- vapply(diagnostics, function(x) {
    !is.null(x) && isTRUE(x$inference_ready[1])
  }, logical(1))
  boundary_contact <- vapply(diagnostics, function(x) {
    !is.null(x) && isTRUE(x$exact_boundary_contact[1])
  }, logical(1))
  lower_gap <- vapply(diagnostics, function(x) {
    if (is.null(x)) NA_real_ else x$minimum_lower_minus_mu[1]
  }, numeric(1))
  dev <- sum(vapply(valid_fits, stats::deviance, numeric(1)))
  df <- sum(vapply(valid_fits, function(x) x$df.fit, numeric(1)))
  n <- sum(vapply(valid_fits, function(x) length(x$y), integer(1)))
  aic <- dev + 2 * df
  aicc <- if (n > df + 1) aic + 2 * df * (df + 1) / (n - df - 1) else NA_real_
  data.frame(
    station = station,
    scale_months = scale,
    model = model_name,
    fitted_months = length(valid_fits),
    converged_months = sum(vapply(valid_fits, LL3_fit_converged, logical(1))),
    inference_ready_months = sum(inference_ready),
    boundary_contact_months = sum(boundary_contact),
    n = n,
    df = df,
    deviance = dev,
    AIC = aic,
    AICc = aicc,
    BIC = dev + log(n) * df,
    min_lower_minus_mu = suppressWarnings(min(lower_gap, na.rm = TRUE))
  )
}

index_from_parameters <- function(y, parameter_frame, epsilon = 1e-10) {
  complete <- is.finite(y) & stats::complete.cases(parameter_frame)
  probability <- rep(NA_real_, length(y))
  probability[complete] <- pLL3(
    y[complete],
    parameter_frame$mu[complete],
    parameter_frame$sigma[complete],
    parameter_frame$nu[complete]
  )
  result <- rep(NA_real_, length(y))
  result[complete] <- stats::qnorm(
    pmin(pmax(probability[complete], epsilon), 1 - epsilon)
  )
  result
}

for (station in stations_to_run) {
  d <- station_data[[station]]
  for (scale_months in scales_to_run) {
    message("Fitting station ", station, ", scale ", scale_months)
    cat(
      paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), station, scale_months),
      file = progress_file, append = TRUE, sep = "\n"
    )
    accumulated <- LL3_accumulate(d$water_balance_mm, scale = scale_months)
    month <- as.integer(format(d$Date, "%m"))
    parameter_store <- setNames(lapply(model_names, function(x) {
      data.frame(
        mu = rep(NA_real_, nrow(d)),
        sigma = rep(NA_real_, nrow(d)),
        nu = rep(NA_real_, nrow(d))
      )
    }), model_names)
    monthly_fits <- setNames(lapply(model_names, function(x) vector("list", 12L)), model_names)

    for (current_month in seq_len(12L)) {
      rows <- month == current_month & is.finite(accumulated)
      fit_data <- data.frame(
        response = accumulated[rows],
        time_scaled = d$time_scaled[rows]
      )
      candidates <- tryCatch(
        suppressWarnings(
          fit_LL3_candidate_models(
            fit_data,
            response = "response",
            covariate = "time_scaled",
            control = LL3_default_control(
              n_cycles = 2000,
              mu_step = 0.03,
              sigma_step = 0.03,
              nu_step = 0.03,
              trace = FALSE
            )
          )
        ),
        error = function(e) {
          structure(list(message = conditionMessage(e)), class = "fit_failure")
        }
      )
      if (inherits(candidates, "fit_failure")) {
        cat(
          paste("FAILED", station, scale_months, current_month, candidates$message),
          file = progress_file, append = TRUE, sep = "\n"
        )
        next
      }
      mapped <- list(
        stationary = candidates$stationary,
        sigma_time = candidates$sigma,
        nu_time = candidates$nu,
        sigma_nu_time = candidates$joint
      )
      for (model_name in model_names) {
        fit <- mapped[[model_name]]
        monthly_fits[[model_name]][[current_month]] <- fit
        pars <- extract_LL3_parameters(fit)
        parameter_store[[model_name]][rows, ] <- data.frame(
          mu = pars$mu,
          sigma = pars$sigma,
          nu = pars$nu
        )
      }
    }

    indices <- lapply(model_names, function(model_name) {
      index_from_parameters(accumulated, parameter_store[[model_name]])
    })
    names(indices) <- model_names

    for (model_name in model_names) {
      model_metrics[[length(model_metrics) + 1L]] <- fit_metrics(
        monthly_fits[[model_name]], model_name, station, scale_months
      )
    }

    traditional <- tryCatch(
      as.numeric(
        SPEI::spei(
          d$water_balance_mm,
          scale = scale_months,
          distribution = "log-Logistic",
          fit = "ub-pwm",
          na.rm = TRUE,
          verbose = FALSE
        )$fitted
      ),
      error = function(e) rep(NA_real_, nrow(d))
    )

    all_values[[length(all_values) + 1L]] <- data.frame(
      station = station,
      date = d$Date,
      scale_months = scale_months,
      precipitation_mm = d$Precipitation,
      mean_temp_c = d$Average_Temp,
      min_temp_c = d$Min_Temp,
      max_temp_c = d$Max_Temp,
      PET_Hargreaves_mm = d$PET_Hargreaves_mm,
      water_balance_mm = d$water_balance_mm,
      accumulated_balance_mm = accumulated,
      SPEI_PWM_stationary = traditional,
      LL3_stationary = indices$stationary,
      LL3_sigma_time = indices$sigma_time,
      LL3_nu_time = indices$nu_time,
      LL3_sigma_nu_time = indices$sigma_nu_time
    )
  }
}

model_comparison <- do.call(rbind, model_metrics)
spei_values <- do.call(rbind, all_values)

summary_rows <- list()
diagnostic_rows <- list()

safe_acf1 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  as.numeric(stats::acf(x, plot = FALSE, lag.max = 1)$acf[2])
}

safe_trend <- function(x, date) {
  keep <- is.finite(x)
  if (sum(keep) < 3L) return(c(slope_decade = NA_real_, kendall_tau = NA_real_, kendall_p = NA_real_))
  decimal_year <- as.integer(format(date[keep], "%Y")) +
    (as.integer(format(date[keep], "%m")) - 0.5) / 12
  slope <- unname(stats::coef(stats::lm(x[keep] ~ decimal_year))[2]) * 10
  kt <- suppressWarnings(stats::cor.test(x[keep], decimal_year, method = "kendall", exact = FALSE))
  c(slope_decade = slope, kendall_tau = unname(kt$estimate), kendall_p = kt$p.value)
}

for (station in stations_to_run) {
  for (scale_months in scales_to_run) {
    mc <- model_comparison[
      model_comparison$station == station & model_comparison$scale_months == scale_months,
      , drop = FALSE
    ]
    stat <- mc[mc$model == "stationary", , drop = FALSE]
    ns <- mc[mc$model != "stationary" & mc$fitted_months == 12L, , drop = FALSE]
    best_ns <- if (nrow(ns)) ns[which.min(ns$AICc), , drop = FALSE] else ns
    best_name <- if (nrow(best_ns)) best_ns$model[1] else NA_character_
    delta <- if (nrow(best_ns) && nrow(stat)) stat$AICc[1] - best_ns$AICc[1] else NA_real_
    valid_inference <- nrow(best_ns) && nrow(stat) &&
      best_ns$inference_ready_months[1] == 12L && stat$inference_ready_months[1] == 12L
    evidence <- if (!valid_inference) {
      "exploratory_boundary_or_fit_issue"
    } else if (!is.finite(delta) || delta < 2) {
      "no_material_support"
    } else if (delta < 4) {
      "weak_support"
    } else if (delta < 10) {
      "moderate_support"
    } else {
      "strong_support"
    }

    vals <- spei_values[
      spei_values$station == station & spei_values$scale_months == scale_months,
      , drop = FALSE
    ]
    ns_col <- if (is.na(best_name)) NA_character_ else switch(
      best_name,
      sigma_time = "LL3_sigma_time",
      nu_time = "LL3_nu_time",
      sigma_nu_time = "LL3_sigma_nu_time"
    )
    selected_ns <- if (is.na(ns_col)) rep(NA_real_, nrow(vals)) else vals[[ns_col]]
    vals$LL3_selected_nonstationary <- selected_ns
    vals$selected_nonstationary_model <- best_name
    spei_values[
      spei_values$station == station & spei_values$scale_months == scale_months,
      "LL3_selected_nonstationary"
    ] <- selected_ns
    spei_values[
      spei_values$station == station & spei_values$scale_months == scale_months,
      "selected_nonstationary_model"
    ] <- best_name

    keep <- is.finite(vals$LL3_stationary) & is.finite(selected_ns)
    stat_trend <- safe_trend(vals$LL3_stationary, vals$date)
    ns_trend <- safe_trend(selected_ns, vals$date)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      station = station,
      scale_months = scale_months,
      best_nonstationary_model = best_name,
      stationary_AICc = if (nrow(stat)) stat$AICc[1] else NA_real_,
      best_nonstationary_AICc = if (nrow(best_ns)) best_ns$AICc[1] else NA_real_,
      delta_AICc_favoring_nonstationary = delta,
      inference_status = evidence,
      stationary_boundary_months = if (nrow(stat)) stat$boundary_contact_months[1] else NA_integer_,
      nonstationary_boundary_months = if (nrow(best_ns)) best_ns$boundary_contact_months[1] else NA_integer_,
      correlation_stationary_nonstationary = if (sum(keep) > 2) stats::cor(vals$LL3_stationary[keep], selected_ns[keep]) else NA_real_,
      mean_absolute_difference = if (any(keep)) mean(abs(vals$LL3_stationary[keep] - selected_ns[keep])) else NA_real_,
      max_absolute_difference = if (any(keep)) max(abs(vals$LL3_stationary[keep] - selected_ns[keep])) else NA_real_,
      observations_abs_difference_ge_0_5 = sum(abs(vals$LL3_stationary[keep] - selected_ns[keep]) >= 0.5),
      stationary_extreme_dry_count = sum(vals$LL3_stationary <= -2, na.rm = TRUE),
      nonstationary_extreme_dry_count = sum(selected_ns <= -2, na.rm = TRUE),
      stationary_slope_per_decade = stat_trend["slope_decade"],
      nonstationary_slope_per_decade = ns_trend["slope_decade"],
      stationary_kendall_p = stat_trend["kendall_p"],
      nonstationary_kendall_p = ns_trend["kendall_p"]
    )

    index_series <- list(
      SPEI_PWM_stationary = vals$SPEI_PWM_stationary,
      LL3_stationary = vals$LL3_stationary,
      LL3_selected_nonstationary = selected_ns
    )
    for (index_name in names(index_series)) {
      z <- index_series[[index_name]]
      z <- z[is.finite(z)]
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
        station = station,
        scale_months = scale_months,
        index = index_name,
        n = length(z),
        mean = mean(z),
        sd = stats::sd(z),
        skewness = mean((z - mean(z))^3) / stats::sd(z)^3,
        excess_kurtosis = mean((z - mean(z))^4) / stats::sd(z)^4 - 3,
        acf_lag1 = safe_acf1(z),
        minimum = min(z),
        maximum = max(z)
      )
    }
  }
}

selected_summary <- do.call(rbind, summary_rows)
index_diagnostics <- do.call(rbind, diagnostic_rows)

write.csv(model_comparison, file.path(output_dir, "model_comparison.csv"), row.names = FALSE)
write.csv(selected_summary, file.path(output_dir, "selected_model_summary.csv"), row.names = FALSE)
write.csv(index_diagnostics, file.path(output_dir, "index_diagnostics.csv"), row.names = FALSE)
write.csv(spei_values, file.path(output_dir, "spei_values.csv"), row.names = FALSE)

methodology <- c(
  "Real-data LL3-GAMLSS SPEI comparison",
  paste("Common analysis period:", common_start, "to", common_end),
  paste0(
    "PET: monthly Hargreaves reference evapotranspiration from SPEI::hargreaves ",
    "(SPEI version ", as.character(utils::packageVersion("SPEI")), "), using monthly ",
    "mean daily minimum/maximum temperature and station latitude; the function ",
    "internally estimates extraterrestrial radiation and returns monthly millimetres."
  ),
  "Water balance: monthly precipitation minus Hargreaves PET.",
  paste0(
    "Input-quality exclusion: stations with at least 12 consecutive monthly precipitation-Tmin-Tmax vectors ",
    "exactly matching the preceding year were excluded before fitting. Excluded station(s): ",
    paste(excluded_stations, collapse = ", "), "."
  ),
  "Accumulation scales: 1, 3, 6, and 12 months (quick mode uses one station and scale).",
  "LL3 fitting: independent model for each calendar month with a constant estimated threshold.",
  "Candidate structures: stationary; scale~time; shape~time; scale+shape~time.",
  "Comparison criterion: joint monthly AICc, with boundary diagnostics. AICc is treated as exploratory whenever any fitted threshold contacts its numerical boundary.",
  "External stationary reference: SPEI package, three-parameter log-logistic distribution, unbiased probability-weighted moments.",
  "Station coordinates are recorded in station_metadata.csv and should be verified against the final authoritative metadata supplied with the meteorological dataset before publication."
)
writeLines(methodology, file.path(output_dir, "methodology.txt"))

png(file.path(output_dir, "scale12_stationary_vs_nonstationary.png"), width = 1800, height = 1400, res = 160)
op <- par(mfrow = c(ceiling(length(stations_to_run) / 2), 2), mar = c(3, 4, 2, 1), oma = c(2, 1, 3, 0))
for (station in stations_to_run) {
  vals <- spei_values[spei_values$station == station & spei_values$scale_months == max(scales_to_run), ]
  plot(vals$date, vals$LL3_stationary, type = "l", col = "#64748B", lwd = 1,
       xlab = "", ylab = "SPEI", main = station, ylim = range(c(vals$LL3_stationary, vals$LL3_selected_nonstationary), na.rm = TRUE))
  lines(vals$date, vals$LL3_selected_nonstationary, col = "#0F766E", lwd = 1)
  abline(h = c(-2, -1, 0), col = c("#B91C1C", "#F59E0B", "#CBD5E1"), lty = c(2, 2, 3))
}
mtext(paste0("LL3 stationary vs selected nonstationary, ", max(scales_to_run), "-month scale"), outer = TRUE, cex = 1.2)
legend("bottom", legend = c("Stationary", "Selected nonstationary"), col = c("#64748B", "#0F766E"), lty = 1, horiz = TRUE, bty = "n", xpd = NA)
par(op)
dev.off()

writeLines(
  c(readLines(progress_file), paste("Completed", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
  progress_file
)

message("Completed. Results written to ", normalizePath(output_dir, winslash = "/"))
