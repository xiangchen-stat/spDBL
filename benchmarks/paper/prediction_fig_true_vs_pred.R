# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Predicted (median) against true y + u for the local prediction runs.
#   Rscript benchmarks/paper/prediction_fig_true_vs_pred.R benchmarks/results/e3/sites benchmarks/results/e3/offgrid true-vs-predicted.png
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(TRUE)
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"; surface <- "#fcfcfb"
inside <- "#2a78d6"; outside <- "#eb6834"; plain <- "#6e6d68"
theme_set(theme_minimal(base_size = 10) + theme(
  text = element_text(colour = ink), axis.text = element_text(colour = ink2),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
  plot.background = element_rect(fill = surface, colour = NA),
  strip.text = element_text(colour = ink, face = "bold", hjust = 0, size = 9),
  legend.position = "top", legend.title = element_blank(),
  plot.title = element_text(face = "bold", size = 10.5)))

read_details <- function(dir, prefix) lapply(list.files(file.path(dir, "raw"),
  paste0("^", prefix, "_r.*rds$"), full.names = TRUE), function(f) readRDS(f)$detail)
coverage_table <- function(dir, file) {
  tab <- utils::read.csv(file.path(dir, file))
  tab <- tab[tab$status == "OK" & tab$target == "y+u", ]
  tapply(tab$coverage95, tab$predictor, mean)
}
# One row per prediction: truth, median and (when known) whether the 95%
# interval covers the truth.
rows <- function(truth, med, lo = NULL, hi = NULL, predictor) data.frame(
  predictor = predictor, truth = as.vector(truth), pred = as.vector(med),
  status = if (is.null(lo)) "No interval shown" else
    ifelse(as.vector(truth) >= as.vector(lo) & as.vector(truth) <= as.vector(hi),
           "Inside spDBL's 95% interval", "Outside spDBL's 95% interval"))

sites <- read_details(args[1], "sites")
cov_s <- coverage_table(args[1], "e3_sites_results.csv")
on <- do.call(rbind, lapply(sites, function(dt) rbind(
  rows(dt$truth_f, dt$f[2, , ], dt$f[1, , ], dt$f[3, , ], "spDBL"),
  rows(dt$truth_f, dt$z_sim[2, , ], predictor = "Emulator, no discrepancy"),
  rows(dt$truth_f, dt$idw, predictor = "Inverse-distance interpolation"))))
on$predictor <- factor(on$predictor, c("spDBL", "Emulator, no discrepancy", "Inverse-distance interpolation"))
on_cov <- c(spDBL = cov_s[["spdbl"]], "Emulator, no discrepancy" = cov_s[["simulator_only"]],
            "Inverse-distance interpolation" = NA)

off <- read_details(args[2], "offgrid")
cov_o <- coverage_table(args[2], "e3_offgrid_results.csv")
og <- do.call(rbind, lapply(off, function(dt) rbind(
  rows(dt$truth_f, dt$bilinear_f[2, , ], dt$bilinear_f[1, , ], dt$bilinear_f[3, , ], "spDBL, bilinear"),
  rows(dt$truth_f, dt$cell_f[2, , ], dt$cell_f[1, , ], dt$cell_f[3, , ], "spDBL, cell"),
  rows(dt$truth_f, dt$sim_bilinear[2, , ], predictor = "Emulator, no discrepancy"),
  rows(dt$truth_f, dt$idw, predictor = "Inverse-distance interpolation"))))
og$predictor <- factor(og$predictor, c("spDBL, bilinear", "spDBL, cell",
  "Emulator, no discrepancy", "Inverse-distance interpolation"))
og_cov <- c("spDBL, bilinear" = cov_o[["spdbl_bilinear"]], "spDBL, cell" = cov_o[["spdbl_cell"]],
            "Emulator, no discrepancy" = cov_o[["sim_bilinear"]], "Inverse-distance interpolation" = NA)

panel <- function(d, cov, title) {
  lim <- range(c(d$truth, d$pred))
  lab <- do.call(rbind, lapply(split(d, d$predictor), function(x) data.frame(
    predictor = x$predictor[1],
    text = sprintf("RMSE %.2f\n95%% coverage %s", sqrt(mean((x$pred - x$truth)^2)),
      if (is.na(cov[[as.character(x$predictor[1])]])) "-" else
        sprintf("%.0f%%", 100 * cov[[as.character(x$predictor[1])]])))))
  ggplot(d, aes(truth, pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = ink2, linewidth = 0.4) +
    geom_point(aes(colour = status), size = 0.7, alpha = 0.45, stroke = 0) +
    geom_text(data = lab, aes(x = lim[1], y = lim[2], label = text), hjust = 0, vjust = 1,
              size = 3, colour = ink, lineheight = 0.95) +
    scale_colour_manual(values = c("Inside spDBL's 95% interval" = inside,
      "Outside spDBL's 95% interval" = outside, "No interval shown" = plain),
      breaks = c("Inside spDBL's 95% interval", "Outside spDBL's 95% interval")) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    facet_wrap(~predictor, nrow = 1) + coord_equal(xlim = lim, ylim = lim) +
    labs(x = "True y + u", y = "Predicted (median)", title = title)
}
pa <- panel(on, on_cov, sprintf("(a) On-grid: held-out grid sites (%d predictions, %d replicates)",
  sum(on$predictor == "spDBL"), length(sites)))
pb <- panel(og, og_cov, sprintf("(b) Off-grid: points between cell centres (%d predictions, %d replicates)",
  sum(og$predictor == "spDBL, bilinear"), length(off)))

png(args[3], width = 11, height = 7.6, units = "in", res = 150, bg = surface)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 8)))
print(pa, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2:7))
print(pb, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:8))
invisible(dev.off())
cat("wrote", args[3], "\n")
