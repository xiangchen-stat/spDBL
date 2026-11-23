# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Figure 3 of the MLOSS paper (Appendix B.6): spDBL's predicted (posterior
# predictive median) against true y + u at locations without field data, on
# the grid (results/e3/sites) and off it (results/e3/offgrid, bilinear field
# map), in the 30 E3 replicates at all 26 times.
# Run from the paper folder (it writes figures/spdbl-true-vs-predicted.png
# there):
#   Rscript <spDBL>/benchmarks/paper/prediction_figure.R [results directory]
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(TRUE)
# The results directory defaults to benchmarks/results next to this script.
script_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else "."
})
res <- if (length(args)) args[1] else file.path(script_dir, "..", "results")
out <- file.path("figures", "spdbl-true-vs-predicted.png")

ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"; surface <- "#ffffff"
inside <- "#2a78d6"; outside <- "#eb6834"
theme_set(theme_minimal(base_size = 9.5) + theme(
  text = element_text(colour = ink), axis.text = element_text(colour = ink2),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
  plot.background = element_rect(fill = surface, colour = NA),
  strip.text = element_text(colour = ink, face = "bold", hjust = 0, size = 9),
  legend.position = "bottom", legend.title = element_blank(),
  plot.title = element_text(face = "bold", size = 11.5)))

# Per-replicate detail saved by e3_sites_worker.R and e3_offgrid_worker.R:
# truth_f is y + u at the targets (times x targets), and each band holds the
# 2.5%, 50% and 97.5% predictive quantiles (3 x times x targets).
details <- function(dir, prefix) lapply(list.files(file.path(dir, "raw"),
  paste0("^", prefix, "_r.*rds$"), full.names = TRUE), function(f) readRDS(f)$detail)
rows <- function(truth, band, panel) data.frame(panel = panel,
  truth = as.vector(truth), pred = as.vector(band[2, , ]),
  status = ifelse(as.vector(truth) >= as.vector(band[1, , ]) & as.vector(truth) <= as.vector(band[3, , ]),
                  "Truth inside the 95% interval", "Truth outside the 95% interval"))
on <- details(file.path(res, "e3", "sites"), "sites")
off <- details(file.path(res, "e3", "offgrid"), "offgrid")
stopifnot(length(on) > 0, length(off) > 0)
panels <- c(sprintf("On grid: held-out sites (n = %s)", format(6L * 26L * length(on), big.mark = ",")),
            sprintf("Off grid: bilinear (n = %s)", format(24L * 26L * length(off), big.mark = ",")))
d <- rbind(do.call(rbind, lapply(on, function(dt) rows(dt$truth_f, dt$f, panels[1]))),
           do.call(rbind, lapply(off, function(dt) rows(dt$truth_f, dt$bilinear_f, panels[2]))))
d$panel <- factor(d$panel, panels)
# Random drawing order, so the colours' proportions match the coverage.
set.seed(1)
d <- d[sample.int(nrow(d)), ]
lim <- range(c(d$truth, d$pred))
lab <- do.call(rbind, lapply(split(d, d$panel), function(x) data.frame(panel = x$panel[1],
  text = sprintf("RMSE %.2f\n95%% coverage %.0f%%", sqrt(mean((x$pred - x$truth)^2)),
                 100 * mean(x$status == "Truth inside the 95% interval")))))
p <- ggplot(d, aes(truth, pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = ink2, linewidth = 0.45) +
  geom_point(aes(colour = status), size = 0.55, alpha = 0.4, stroke = 0) +
  geom_text(data = lab, aes(x = lim[1], y = lim[2], label = text), hjust = 0, vjust = 1,
            size = 3.1, colour = ink, lineheight = 0.95) +
  scale_colour_manual(values = c("Truth inside the 95% interval" = inside,
                                 "Truth outside the 95% interval" = outside)) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  facet_wrap(~panel, nrow = 1) + coord_equal(xlim = lim, ylim = lim) +
  labs(x = "True field value", y = "spDBL prediction (posterior predictive median)")
dir.create("figures", showWarnings = FALSE)
ggsave(out, p, width = 7, height = 4.1, dpi = 300, bg = surface)
cat("wrote", out, "\n")
