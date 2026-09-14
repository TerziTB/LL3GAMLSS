#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(terra))

args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args)) args[[1L]] else file.path(
  "analysis", "era5_results", "final_full_195001_202512"
)

cell_summary <- read.csv(file.path(result_dir, "cell_model_selection.csv"))
area_fractions <- read.csv(file.path(result_dir, "area_weighted_model_fractions.csv"))
basin_nspei <- read.csv(file.path(result_dir, "basin_nspei.csv"))
basin_nspei$date <- as.Date(basin_nspei$date)
basin <- vect(file.path("analysis", "era5_input", "seyhan_hydrobasins_lev08.gpkg"))
stations <- read.csv(file.path("analysis", "era5_input", "station_basin_membership.csv"))
station_points <- vect(stations, geom = c("longitude", "latitude"), crs = "EPSG:4326")

grid <- rast(
  ncols = round((max(cell_summary$longitude) - min(cell_summary$longitude)) / 0.1) + 1L,
  nrows = round((max(cell_summary$latitude) - min(cell_summary$latitude)) / 0.1) + 1L,
  xmin = min(cell_summary$longitude) - 0.05,
  xmax = max(cell_summary$longitude) + 0.05,
  ymin = min(cell_summary$latitude) - 0.05,
  ymax = max(cell_summary$latitude) + 0.05,
  crs = "EPSG:4326"
)

model_levels <- c("stationary", "sigma_time", "nu_time", "sigma_nu_time")
model_colors <- c("#D1D5DB", "#0F766E", "#D97706", "#7C3AED")
model_labels <- c(
  "M0: stationary",
  expression(M[sigma] * ": scale ~ time"),
  expression(M[nu] * ": shape ~ time"),
  expression(M[sigma*nu] * ": scale + shape ~ time")
)

make_raster <- function(d, value) {
  output <- grid
  values(output) <- NA_real_
  cells <- cellFromXY(output, d[, c("longitude", "latitude")])
  values(output)[cells] <- value
  output
}

png(
  file.path(result_dir, "era5_all_scales_selected_models.png"),
  width = 3000,
  height = 2800,
  res = 240,
  pointsize = 15
)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.20))
par(oma = c(0.4, 0.8, 0.5, 0.6))
for (scale_months in c(1L, 3L, 6L, 12L)) {
  d <- cell_summary[cell_summary$scale_months == scale_months, ]
  raster <- make_raster(d, match(d$selected_model, model_levels))
  par(mar = c(4.0, 5.8, 3.4, 1.2), mgp = c(3.5, 1.15, 0),
      cex.axis = 1.35, cex.lab = 1.65, font.lab = 2)
  plot(
    raster,
    col = model_colors,
    breaks = seq(0.5, 4.5, by = 1),
    legend = FALSE,
    axes = TRUE,
    xlab = "",
    ylab = "Latitude (degrees N)",
    main = paste0(scale_months, "-month scale"),
    cex.main = 1.48,
    font.main = 2
  )
  lines(basin, lwd = 1.5, col = "#111827")
  points(
    station_points,
    pch = ifelse(stations$inside_hydrobasins_boundary, 21, 1),
    bg = ifelse(stations$inside_hydrobasins_boundary, "#FFFFFF", NA),
    col = "#111827",
    cex = 0.65
  )
}
par(mar = rep(0, 4))
plot.new()
text(0.5, 0.86, "Longitude (degrees E)", cex = 1.25, font = 2)
legend(
  "bottom",
  inset = 0.02,
  legend = c(model_labels, "In-basin station", "Outside station"),
  col = c(model_colors, "#111827", "#111827"),
  pch = c(rep(15, 4), 21, 1),
  pt.bg = c(rep(NA, 4), "#FFFFFF", NA),
  bty = "n",
  ncol = 3,
  cex = 1.08
)
dev.off()

fraction_matrix <- xtabs(
  basin_area_percent ~ selected_model + scale_months,
  data = area_fractions[area_fractions$selected_model != "fit_failure", ]
)
fraction_matrix <- fraction_matrix[model_levels, , drop = FALSE]
png(
  file.path(result_dir, "era5_all_scales_model_fractions.png"),
  width = 2100,
  height = 1400,
  res = 210
)
par(mar = c(4.8, 4.8, 3.4, 1.0))
barplot(
  fraction_matrix,
  col = model_colors,
  border = "white",
  ylim = c(0, 100),
  xlab = "Accumulation scale (months)",
  ylab = "Basin area (%)",
  main = "Area-weighted LL3 model selection",
  legend.text = model_labels,
  args.legend = list(x = "top", inset = c(0, -0.01), bty = "n", ncol = 2, cex = 0.85)
)
abline(h = seq(0, 100, 20), col = "#FFFFFF80", lwd = 0.8)
box()
dev.off()

png(
  file.path(result_dir, "era5_basin_nspei_timeseries.png"),
  width = 2400,
  height = 1650,
  res = 210
)
par(mfrow = c(2, 2), mar = c(3.7, 4.2, 3.0, 1.0), oma = c(1.2, 1.0, 2.4, 0))
for (scale_months in c(1L, 3L, 6L, 12L)) {
  d <- basin_nspei[basin_nspei$scale_months == scale_months, ]
  plot(
    d$date,
    d$nSPEI,
    type = "n",
    xlab = "",
    ylab = "Selected LL3 index",
    main = paste0(scale_months, "-month scale"),
    ylim = c(-3.5, 3.5)
  )
  rect(
    par("usr")[1], -3.5, par("usr")[2], -2,
    col = "#B91C1C18", border = NA
  )
  rect(
    par("usr")[1], -2, par("usr")[2], -1,
    col = "#D9770618", border = NA
  )
  lines(d$date, d$nSPEI, col = "#0F766E", lwd = 0.75)
  abline(h = c(-2, -1, 0), col = c("#B91C1C", "#D97706", "#9CA3AF"), lty = c(2, 2, 3))
  mtext(
    paste("Selected model:", unique(d$selected_model)),
    side = 3,
    line = 0.25,
    cex = 0.72
  )
}
mtext(
  "Area-weighted Seyhan nSPEI (ERA5-Land, 1950-2025)",
  outer = TRUE,
  side = 3,
  cex = 1.2,
  line = 0.5
)
dev.off()

cat("All-scale figures written to", normalizePath(result_dir, winslash = "/"), "\n")
