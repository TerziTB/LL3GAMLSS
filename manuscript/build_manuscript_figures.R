options(stringsAsFactors = FALSE)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "manuscript", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

navy <- "#17365D"
blue <- "#2F6B9A"
teal <- "#168A88"
orange <- "#D97706"
red <- "#B42318"
gray <- "#667085"
light <- "#EEF3F7"

open_png <- function(name, width = 2400, height = 1600) {
  grDevices::png(
    file.path(out_dir, name), width = width, height = height,
    res = 300, type = "cairo-png", bg = "white"
  )
}

panel_label <- function(label) {
  graphics::mtext(label, side = 3, line = 0.2, adj = 0, font = 2, cex = 1.05)
}

# Fig. 3: simulation validation -------------------------------------------------
stress <- read.csv(file.path(root, "validation", "results", "publication_stress_summary.csv"))
null <- read.csv(file.path(root, "validation", "results", "stationary_null_summary.csv"))
coverage <- read.csv(file.path(root, "validation", "results", "interval_coverage_recovered_summary.csv"))

scenario_labels <- c(ar1 = "AR(1), correlation 0.55", contaminated = "5% contamination", iid = "IID")
scenario_cols <- c(ar1 = red, contaminated = orange, iid = teal)

open_png("Fig3_simulation_validation.png", 2400, 1900)
par(mfrow = c(2, 2), mar = c(4.6, 4.8, 2.0, 1.0), oma = c(0, 0, 0.4, 0),
    las = 1, cex.axis = 0.85, cex.lab = 0.95)

plot(NA, xlim = range(stress$n), ylim = c(0, max(stress$sigma_slope_rmse) * 1.08),
     xlab = "Monthly sample size", ylab = "Scale-slope RMSE")
for (s in names(scenario_labels)) {
  z <- stress[stress$scenario == s, ]
  lines(z$n, z$sigma_slope_rmse, type = "b", pch = 16, lwd = 2,
        col = scenario_cols[s])
}
legend("topright", legend = scenario_labels, col = scenario_cols,
       lwd = 2, pch = 16, bty = "n", cex = 0.78)
panel_label("a")

plot(NA, xlim = range(stress$n), ylim = c(0, 0.075),
     xlab = "Monthly sample size", ylab = "Boundary-contact rate")
for (s in names(scenario_labels)) {
  z <- stress[stress$scenario == s, ]
  lines(z$n, z$boundary_contact_rate, type = "b", pch = 16, lwd = 2,
        col = scenario_cols[s])
}
abline(h = 0, col = "#98A2B3", lty = 3)
panel_label("b")

matplot(null$n,
        cbind(null$AIC_false_nonstationary_rate, null$BIC_false_nonstationary_rate),
        type = "b", pch = c(16, 17), lty = 1, lwd = 2,
        col = c(blue, orange), ylim = c(0, 0.40),
        xlab = "Sample size", ylab = "False non-stationary selection rate")
legend("topright", legend = c("Minimum AIC", "Minimum BIC"),
       col = c(blue, orange), pch = c(16, 17), lwd = 2, bty = "n", cex = 0.80)
panel_label("c")

plot(coverage$n, coverage$sigma_slope_coverage, type = "b", pch = 16,
     lwd = 2, col = teal, ylim = c(0.90, 0.98), xlim = c(90, 260),
     xlab = "Sample size", ylab = "Empirical 95% interval coverage")
lines(coverage$n, coverage$nu_slope_coverage, type = "b", pch = 17,
      lwd = 2, col = blue)
abline(h = 0.95, lty = 2, col = red, lwd = 1.5)
legend("bottomleft", legend = c("Scale slope", "Shape slope", "Nominal 0.95"),
       pch = c(16, 17, NA), lty = c(1, 1, 2), lwd = 2,
       col = c(teal, blue, red), bty = "n", cex = 0.78)
panel_label("d")

dev.off()

# Shared station metadata for Figs. 4 and 5 -------------------------------------
stations <- read.csv(file.path(root, "analysis", "manuscript_inputs", "station_metadata.csv"))
stations <- stations[as.character(stations$station) != "17934", , drop = FALSE]

# Fig. 4: AICc support heat map --------------------------------------------------
sel <- read.csv(file.path(root, "analysis", "manuscript_inputs", "station_selection.csv"))
sel$station <- as.character(sel$station)
sel <- sel[sel$station != "17934", , drop = FALSE]
station_order <- as.character(stations$station)
scale_order <- c(1, 3, 6, 12)
zmat <- matrix(NA_real_, nrow = length(station_order), ncol = length(scale_order),
               dimnames = list(station_order, paste0(scale_order, " month")))
for (i in seq_len(nrow(sel))) {
  zmat[sel$station[i], paste0(sel$scale_months[i], " month")] <-
    sel$delta_AICc_favoring_nonstationary[i]
}
row_labels <- paste0(stations$station_name, " (", stations$station, ")")

palette <- colorRampPalette(c("#355C7D", "#D9E4EC", "#FFFFFF", "#F7D7A8", "#B42318"))(101)
limit <- max(abs(zmat), na.rm = TRUE)
open_png("Fig4_model_support_heatmap.png", 2200, 1700)
par(mar = c(5.5, 10.5, 2.2, 4.5), las = 1)
image(seq_along(scale_order), seq_along(station_order), t(zmat),
      col = palette, zlim = c(-limit, limit), axes = FALSE,
      xlab = "Accumulation scale", ylab = "")
axis(1, at = seq_along(scale_order), labels = paste0(scale_order, " months"))
axis(2, at = seq_along(station_order), labels = row_labels)
box()
for (i in seq_along(station_order)) {
  for (j in seq_along(scale_order)) {
    v <- zmat[i, j]
    text(j, i, sprintf("%.1f%s", v, ifelse(v >= 10, "*", "")),
         cex = 0.83, font = ifelse(v >= 10, 2, 1),
         col = ifelse(abs(v) > 0.55 * limit, "white", "#101828"))
  }
}
mtext("Asterisks denote strong support at a Delta AICc threshold of 10",
      side = 1, line = 4.0, cex = 0.78, col = gray)
dev.off()

# Fig. 5: selected 12-month time series -----------------------------------------
values <- read.csv(file.path(root, "analysis", "manuscript_inputs", "figure5_stationary_nonstationary.csv"))
values$date <- as.Date(values$date)
show_stations <- c("17351", "17802", "17837", "17840")
show_names <- setNames(stations$station_name, as.character(stations$station))

open_png("Fig5_spei12_comparison.png", 2400, 1800)
par(mfrow = c(2, 2), mar = c(3.8, 4.2, 2.4, 1.2), oma = c(4.2, 0.5, 0.4, 0.5),
    las = 1, cex.axis = 0.82, cex.lab = 0.90)
for (s in show_stations) {
  z <- values[as.character(values$station) == s & values$scale_months == 12, ]
  sm <- sel[sel$station == s & sel$scale_months == 12, ]
  yr <- range(c(z$LL3_stationary, z$LL3_selected_nonstationary), finite = TRUE)
  plot(z$date, z$LL3_stationary, type = "l", lwd = 1.15,
       col = "#1F2937", lty = 2,
       xlab = "", ylab = "SPEI-12", ylim = yr,
       main = sprintf("%s (%s), Delta AICc: %.1f",
                      show_names[s], s, sm$delta_AICc_favoring_nonstationary))
  lines(z$date, z$LL3_selected_nonstationary, col = "#0072B2",
        lwd = 1.65, lty = 1)
  abline(h = -1, col = orange, lty = 3, lwd = 1)
  abline(h = -2, col = red, lty = 4, lwd = 1)
  abline(h = 0, col = "#D0D5DD", lty = 3)
}
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()
legend("bottom", inset = 0.012, horiz = TRUE, bty = "n",
       legend = c("Stationary LL3 (M0)", "Best non-stationary candidate", "SPEI threshold -1", "SPEI threshold -2"),
       col = c("#1F2937", "#0072B2", orange, red),
       lty = c(2, 1, 3, 4), lwd = c(1.4, 2.0, 1.1, 1.1),
       cex = 0.72)
dev.off()

cat("Created manuscript figures in", out_dir, "\n")
