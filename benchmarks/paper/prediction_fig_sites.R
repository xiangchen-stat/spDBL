# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Figure for the on-grid site-prediction study (results/e3/<profile>/).
#   Rscript benchmarks/paper/prediction_fig_sites.R benchmarks/results/e3/sites 1 sites.png
# (a) Held-out sites of one replicate: truth y + u, held-out observations z,
#     spDBL median and 95% band, and the two baselines' point predictions.
# (b) Per-replicate RMSE for the held-out observations, by predictor.
# (c) Per-replicate 95% coverage of the held-out observations.
suppressPackageStartupMessages({library(ggplot2)})
args <- commandArgs(TRUE)
dir <- args[1]; rep_id <- as.integer(args[2]); out <- args[3]

ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"; surface <- "#fcfcfb"
pal <- c(spDBL = "#2a78d6", "Emulator, no discrepancy" = "#eb6834",
         "Inverse-distance interpolation" = "#1baf7a")
theme_set(theme_minimal(base_size = 10) + theme(
  text = element_text(colour = ink), axis.text = element_text(colour = ink2),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
  plot.background = element_rect(fill = surface, colour = NA),
  strip.text = element_text(colour = ink, face = "bold", hjust = 0),
  legend.position = "top", legend.title = element_blank(),
  plot.title = element_text(face = "bold", size = 10.5)))

res <- readRDS(file.path(dir, "raw", sprintf("sites_r%03d.rds", rep_id)))
dt <- res$detail
K <- length(dt$held)
grid_rc <- function(s) sprintf("site %d (row %d, col %d)", s, (s - 1) %% 6 + 1, (s - 1) %/% 6 + 1)
long <- do.call(rbind, lapply(seq_len(K), function(k) data.frame(
  site = grid_rc(dt$held[k]), time = dt$time,
  truth = dt$truth_f[, k], obs = dt$truth_z[, k],
  lo = dt$f[1, , k], med = dt$f[2, , k], hi = dt$f[3, , k],
  sim = dt$z_sim[2, , k], idw = dt$idw[, k])))
long$site <- factor(long$site, unique(long$site))
lines <- rbind(
  data.frame(long[, c("site", "time")], value = long$med, what = "spDBL"),
  data.frame(long[, c("site", "time")], value = long$sim, what = "Emulator, no discrepancy"),
  data.frame(long[, c("site", "time")], value = long$idw, what = "Inverse-distance interpolation"))
lines$what <- factor(lines$what, names(pal))
pa <- ggplot(long, aes(time)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = pal[["spDBL"]], alpha = 0.18) +
  geom_line(data = lines, aes(y = value, colour = what), linewidth = 0.7) +
  geom_line(aes(y = truth, linetype = "Truth y + u"), colour = ink, linewidth = 0.6) +
  geom_point(aes(y = obs, shape = "Held-out observation z"), colour = ink2, size = 1.3) +
  scale_colour_manual(values = pal) +
  scale_linetype_manual(values = c("Truth y + u" = "solid")) +
  scale_shape_manual(values = c("Held-out observation z" = 1)) +
  facet_wrap(~site, ncol = 3) +
  labs(x = "Time", y = "log(I + 1)",
       title = sprintf("(a) Replicate %d: the %d held-out sites (band: spDBL 95%% interval for y + u)", rep_id, K)) +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2),
         shape = guide_legend(order = 3))

tab <- utils::read.csv(file.path(dir, "e3_sites_results.csv"))
tab <- tab[tab$status == "OK" & tab$target == "z", ]
tab$predictor <- factor(c(spdbl = "spDBL", simulator_only = "Emulator, no discrepancy",
  idw = "Inverse-distance interpolation")[tab$predictor], names(pal))
pb <- ggplot(tab, aes(predictor, rmse, colour = predictor)) +
  geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.8, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, colour = ink, linewidth = 0.3) +
  scale_colour_manual(values = pal) + scale_x_discrete(labels = function(x) sub(", ", ",\n", sub(" interpolation", "\ninterpolation", x))) +
  labs(x = NULL, y = "RMSE of held-out z", title = sprintf("(b) RMSE per replicate (n = %d)", length(unique(tab$rep))))
cov <- tab[!is.na(tab$coverage95), ]
pc <- ggplot(cov, aes(predictor, coverage95, colour = predictor)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = ink2, linewidth = 0.4) +
  geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.8, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, colour = ink, linewidth = 0.3) +
  scale_colour_manual(values = pal) + scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = function(x) sub(", ", ",\n", x)) +
  labs(x = NULL, y = "95% interval coverage of held-out z", title = "(c) Coverage per replicate (dashed: nominal)")

png(out, width = 9, height = 8.2, units = "in", res = 150, bg = surface)
grid::grid.newpage()
lay <- grid::grid.layout(2, 2, heights = grid::unit(c(0.58, 0.42), "npc"))
grid::pushViewport(grid::viewport(layout = lay))
print(pa, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(pb, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pc, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
invisible(dev.off())
cat("wrote", out, "\n")
