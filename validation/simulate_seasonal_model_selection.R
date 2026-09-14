args <- commandArgs(trailingOnly = TRUE)
argument_value <- function(prefix, default) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit)) sub(prefix, "", hit[[1L]]) else default
}

project_dir <- local({
  source(file.path("validation", "validation_helpers.R"), local = TRUE)
  resolve_project_dir()
})
setwd(project_dir)
source(file.path("validation", "validation_helpers.R"))
load_LL3_source(project_dir)

replicates <- as.integer(argument_value("--replicates=", "30"))
workers <- max(1L, as.integer(argument_value("--workers=", "4")))
n_years <- as.integer(argument_value("--years=", "66"))
if (!is.finite(replicates) || replicates < 1L) stop("replicates must be positive.")
if (!is.finite(n_years) || n_years < 60L || n_years > 75L) {
  stop("years must be between 60 and 75.")
}

output_dir <- file.path(project_dir, "validation", "results", "seasonal_selection")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
truths <- c("M0", "M_sigma")
dependence <- c(iid = 0, ar1 = 0.4)
scales <- c(1L, 3L, 6L, 12L)
sigma_slope_truth <- 0.25
base_seed <- 940000L
control <- empirical_LL3_control()

simulate_balance <- function(n_years, truth, phi, seed) {
  set.seed(seed)
  n <- 12L * n_years
  month <- rep(seq_len(12L), n_years)
  time_scaled <- as.numeric(scale(seq_len(n)))
  angle <- 2 * pi * (month - 1) / 12
  mu <- -150 + 25 * cos(angle)
  sigma <- 75 * exp(0.22 * sin(angle - pi / 6))
  nu <- 4.0 * exp(0.05 * cos(angle + pi / 4))
  if (truth == "M_sigma") sigma <- sigma * exp(sigma_slope_truth * time_scaled)
  innovations_sd <- sqrt(1 - phi^2)
  latent <- if (phi == 0) {
    stats::rnorm(n)
  } else {
    as.numeric(stats::arima.sim(model = list(ar = phi), n = n, sd = innovations_sd))
  }
  probability <- pmin(pmax(stats::pnorm(latent), 1e-10), 1 - 1e-10)
  data.frame(
    date = seq(as.Date("1960-01-01"), by = "month", length.out = n),
    month = month,
    time_scaled = time_scaled,
    balance = qLL3(probability, mu, sigma, nu),
    stringsAsFactors = FALSE
  )
}

fit_scale <- function(data, scale_months) {
  accumulated <- LL3_accumulate(data$balance, scale_months)
  models <- c("M0", "M_sigma", "M_nu", "M_sigma_nu")
  fits <- setNames(lapply(models, function(unused) vector("list", 12L)), models)
  errors <- character()
  for (current_month in seq_len(12L)) {
    rows <- data$month == current_month & is.finite(accumulated)
    fit_data <- data.frame(
      response = accumulated[rows],
      time_scaled = data$time_scaled[rows]
    )
    candidates <- tryCatch(
      suppressWarnings(fit_LL3_candidate_models(
        fit_data,
        response = "response",
        covariate = "time_scaled",
        control = control
      )),
      error = function(e) e
    )
    if (inherits(candidates, "error")) {
      errors <- c(errors, sprintf("month %d: %s", current_month, conditionMessage(candidates)))
      next
    }
    fits$M0[[current_month]] <- candidates$stationary
    fits$M_sigma[[current_month]] <- candidates$sigma
    fits$M_nu[[current_month]] <- candidates$nu
    fits$M_sigma_nu[[current_month]] <- candidates$joint
  }
  metrics <- do.call(rbind, lapply(models, function(model) {
    seasonal_candidate_metrics(fits[[model]], model)
  }))
  eligible <- with(metrics,
    fitted_months == 12L & converged_months == 12L &
      inference_ready_months == 12L & is.finite(AICc)
  )
  stationary_row <- which(metrics$model == "M0")
  ns_rows <- which(metrics$model != "M0" & eligible)
  best_ns <- if (length(ns_rows)) ns_rows[[which.min(metrics$AICc[ns_rows])]] else integer()
  delta <- if (eligible[[stationary_row]] && length(best_ns)) {
    metrics$AICc[[stationary_row]] - metrics$AICc[[best_ns]]
  } else {
    NA_real_
  }
  best_ns_model <- if (length(best_ns)) metrics$model[[best_ns]] else NA_character_
  selected_2 <- if (is.finite(delta) && delta >= 2) best_ns_model else if (eligible[[stationary_row]]) "M0" else NA_character_
  selected_10 <- if (is.finite(delta) && delta >= 10) best_ns_model else if (eligible[[stationary_row]]) "M0" else NA_character_
  sigma_row <- which(metrics$model == "M_sigma")
  data.frame(
    scale_months = scale_months,
    total_monthly_observations = nrow(data),
    minimum_calendar_month_n = min(vapply(fits$M0, function(x) if (inherits(x, "gamlss")) length(x$y) else 0L, integer(1))),
    all_candidates_complete = all(eligible),
    any_boundary_contact = any(metrics$boundary_contact_months > 0, na.rm = TRUE),
    delta_AICc = delta,
    best_nonstationary_model = best_ns_model,
    selected_delta2 = selected_2,
    selected_delta10 = selected_10,
    M_sigma_slope_median = metrics$sigma_slope_median[[sigma_row]],
    fitting_errors = paste(errors, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

run_replicate <- function(task) {
  started <- proc.time()[["elapsed"]]
  simulated <- simulate_balance(task$n_years, task$truth, task$phi, task$seed)
  result <- do.call(rbind, lapply(scales, function(scale_months) {
    fit_scale(simulated, scale_months)
  }))
  result$truth <- task$truth
  result$dependence <- task$dependence
  result$phi <- task$phi
  result$replicate <- task$replicate
  result$seed <- task$seed
  result$runtime_seconds <- proc.time()[["elapsed"]] - started
  result
}

tasks <- expand.grid(
  truth = truths,
  dependence = names(dependence),
  replicate = seq_len(replicates),
  stringsAsFactors = FALSE
)
tasks$phi <- unname(dependence[tasks$dependence])
tasks$n_years <- n_years
tasks$seed <- base_seed + seq_len(nrow(tasks))
task_list <- split(tasks, seq_len(nrow(tasks)))

run_started <- Sys.time()
if (workers == 1L) {
  raw_list <- lapply(task_list, run_replicate)
} else {
  cluster <- parallel::makePSOCKcluster(min(workers, length(task_list)))
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterExport(cluster, c("project_dir"), envir = environment())
  parallel::clusterEvalQ(cluster, {
    setwd(project_dir)
    source(file.path("validation", "validation_helpers.R"))
    load_LL3_source(project_dir)
    NULL
  })
  parallel::clusterExport(
    cluster,
    c(
      "simulate_balance", "fit_scale", "run_replicate", "scales", "control",
      "sigma_slope_truth", "seasonal_candidate_metrics", "rate_with_mcse"
    ),
    envir = environment()
  )
  raw_list <- parallel::parLapplyLB(cluster, task_list, run_replicate)
  parallel::stopCluster(cluster)
  on.exit(NULL, add = FALSE)
}
raw <- do.call(rbind, raw_list)
raw <- raw[order(raw$truth, raw$dependence, raw$replicate, raw$scale_months), ]
rownames(raw) <- NULL

summarize_group <- function(d) {
  false_2 <- rate_with_mcse(d$truth == "M0" & d$selected_delta2 != "M0")
  false_10 <- rate_with_mcse(d$truth == "M0" & d$selected_delta10 != "M0")
  detection_2 <- rate_with_mcse(d$truth == "M_sigma" & d$selected_delta2 != "M0")
  detection_10 <- rate_with_mcse(d$truth == "M_sigma" & d$selected_delta10 != "M0")
  correct_2 <- rate_with_mcse(d$truth == "M_sigma" & d$selected_delta2 == "M_sigma")
  completion <- rate_with_mcse(d$all_candidates_complete)
  boundary <- rate_with_mcse(d$any_boundary_contact)
  slope_target <- if (unique(d$truth) == "M_sigma" && unique(d$scale_months) == 1L) sigma_slope_truth else NA_real_
  slope_bias <- if (is.finite(slope_target)) mean(d$M_sigma_slope_median - slope_target, na.rm = TRUE) else NA_real_
  slope_rmse <- if (is.finite(slope_target)) sqrt(mean((d$M_sigma_slope_median - slope_target)^2, na.rm = TRUE)) else NA_real_
  data.frame(
    truth = unique(d$truth),
    dependence = unique(d$dependence),
    phi = unique(d$phi),
    scale_months = unique(d$scale_months),
    replicates = nrow(d),
    total_monthly_observations = unique(d$total_monthly_observations),
    calendar_month_n = paste(range(d$minimum_calendar_month_n), collapse = "-"),
    false_selection_delta2 = if (unique(d$truth) == "M0") false_2[["estimate"]] else NA_real_,
    false_selection_delta2_mcse = if (unique(d$truth) == "M0") false_2[["mcse"]] else NA_real_,
    false_selection_delta10 = if (unique(d$truth) == "M0") false_10[["estimate"]] else NA_real_,
    false_selection_delta10_mcse = if (unique(d$truth) == "M0") false_10[["mcse"]] else NA_real_,
    detection_power_delta2 = if (unique(d$truth) == "M_sigma") detection_2[["estimate"]] else NA_real_,
    detection_power_delta2_mcse = if (unique(d$truth) == "M_sigma") detection_2[["mcse"]] else NA_real_,
    detection_power_delta10 = if (unique(d$truth) == "M_sigma") detection_10[["estimate"]] else NA_real_,
    detection_power_delta10_mcse = if (unique(d$truth) == "M_sigma") detection_10[["mcse"]] else NA_real_,
    correct_Msigma_delta2 = if (unique(d$truth) == "M_sigma") correct_2[["estimate"]] else NA_real_,
    correct_Msigma_delta2_mcse = if (unique(d$truth) == "M_sigma") correct_2[["mcse"]] else NA_real_,
    slope_target = slope_target,
    slope_bias = slope_bias,
    slope_rmse = slope_rmse,
    boundary_contact_rate = boundary[["estimate"]],
    boundary_contact_mcse = boundary[["mcse"]],
    fit_completion_rate = completion[["estimate"]],
    fit_completion_mcse = completion[["mcse"]],
    stringsAsFactors = FALSE
  )
}

groups <- split(raw, interaction(raw$truth, raw$dependence, raw$scale_months, drop = TRUE))
summary <- do.call(rbind, lapply(groups, summarize_group))
rownames(summary) <- NULL
run_finished <- Sys.time()

utils::write.csv(raw, file.path(output_dir, "seasonal_selection_raw.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(output_dir, "seasonal_selection_summary.csv"), row.names = FALSE)
utils::write.csv(tasks, file.path(output_dir, "seasonal_selection_seeds.csv"), row.names = FALSE)
settings <- c(
  sprintf("replicates_per_truth_dependence=%d", replicates),
  sprintf("years=%d", n_years),
  sprintf("monthly_observations=%d", 12L * n_years),
  sprintf("calendar_month_observations=%d", n_years),
  sprintf("scales=%s", paste(scales, collapse = ",")),
  sprintf("sigma_slope_truth_at_scale1=%.3f", sigma_slope_truth),
  "ar1_phi=0.4",
  "AICc aggregation=sum monthly deviance and df across 12 independent calendar-month fits",
  "selection=M0 unless best non-stationary model improves AICc by the stated threshold",
  "fit_completion_rate=all four candidates complete and inference-ready for all 12 months",
  "selection_completion_rate=M0 and at least one non-stationary candidate permit a finite selection",
  "slope bias and RMSE are reported only at scale 1 because sums of LL3 variables are not generally LL3",
  "control=n.cyc 2000; c.crit 0.001; autostep TRUE; gd.tol Inf; steps 0.03",
  sprintf("workers=%d", min(workers, length(task_list))),
  sprintf("started=%s", format(run_started, "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("finished=%s", format(run_finished, "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("wall_clock_seconds=%.3f", as.numeric(difftime(run_finished, run_started, units = "secs")))
)
writeLines(settings, file.path(output_dir, "seasonal_selection_settings.txt"))
write_session_info(file.path(output_dir, "sessionInfo.txt"))
source(file.path(project_dir, "validation", "summarize_seasonal_model_selection.R"))
print(summary)
