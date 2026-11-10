# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Figure 2 of the MLOSS paper: benchmarks E1-E3 (E2 from the fixed-domain,
# 10-design profile, with means and 95% t intervals over designs).
# Run from the paper folder (it writes figures/benchmarks.pdf there):
#   Rscript <spDBL>/benchmarks/paper/benchmarks.R [results directory]
# cairo_pdf embeds the fonts.
args <- commandArgs(TRUE)
# The results directory defaults to benchmarks/results next to this script.
script_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else "."
})
res <- if (length(args)) args[1] else file.path(script_dir, "..", "results")
read <- function(e, f) {
  x <- utils::read.csv(file.path(res, e, "paper", f))
  x[x$status == "OK", ]
}
e1 <- read("e1", "e1_results.csv")
# Coefficient draws regenerated on demand (emulator_learn(keep_draws = FALSE)):
# faster fits and less memory, identical results.
e1_lazy <- utils::read.csv(file.path(res, "e1", "paper_lazy", "e1_results.csv"))
e1_lazy <- e1_lazy[e1_lazy$status == "OK", ]
e2 <- utils::read.csv(file.path(res, "e2", "fixed", "e2_per_design.csv"))
e2$S <- e2$grid^2
e3 <- utils::read.csv(file.path(res, "e3", "mode", "e3_results.csv")); e3 <- e3[e3$status == "OK", ]

# Colour follows the package: spDBL blue, Gu's RobustGaSP/RobustCalibration
# orange, mlegp aqua; shapes repeat the encoding for greyscale print.
col <- c(spdbl = "#2a78d6", gu = "#eb6834", mlegp = "#1baf7a")
win_col <- c("1" = "#86b6ef", "2" = "#3987e5", "5" = "#184f95")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"

dir.create("figures", showWarnings = FALSE)
out <- file.path("figures", "benchmarks.pdf")
cairo_pdf(out, width = 6.0, height = 1.68, family = "Helvetica", pointsize = 9)
par(mfrow = c(1, 4), mar = c(2.5, 2.7, 1.5, 0.5), mgp = c(1.45, 0.35, 0),
    tcl = -0.2, las = 1, cex.axis = 0.8, cex.lab = 0.85, col.axis = ink2,
    col.lab = ink, fg = ink2, bty = "n", lwd = 0.8)
panel_title <- function(s) mtext(s, side = 3, line = 0.35, adj = 0, cex = 0.8, font = 2, col = ink)
hgrid <- function(at) abline(h = at, col = grid_col, lwd = 0.6)

# (a) E1: emulator fit time against the number of sites, by episode window.
xs <- c(100, 10000)
plot(NA, xlim = c(100, 10000 * 9), ylim = c(0.1, 1500), log = "xy", axes = FALSE,
     xlab = expression("sites " * italic(S)), ylab = "fit time (s)")
yt <- 10^(-1:3); hgrid(yt)
axis(1, at = c(1e2, 1e3, 1e4), labels = c("10\u00b2", "10\u00b3", "10\u2074"))
axis(2, at = yt, labels = c("0.1", "1", "10", "100", "1000"))
for (w in names(win_col)) {
  stored <- e1[e1$episode_window == as.numeric(w), ]; stored <- stored[order(stored$S), ]
  lines(stored$S, stored$fit_s, col = adjustcolor(win_col[w], 0.45), lwd = 1, lty = 2)
  d <- e1_lazy[e1_lazy$episode_window == as.numeric(w), ]; d <- d[order(d$S), ]
  lines(d$S, d$fit_s, col = win_col[w], lwd = 1.6)
  points(d$S, d$fit_s, col = win_col[w], pch = c("1" = 16, "2" = 17, "5" = 15)[w], cex = 0.75)
  last <- d[nrow(d), ]
  # Nudge the w = 1 and w = 2 labels apart; their end points are close on a log scale.
  nudge <- c("1" = 0.55, "2" = 1.6, "5" = 1.15)[[w]]
  text(last$S * 1.35, last$fit_s * nudge, sprintf("w=%s\n%.0f GiB", w, last$peak_rss_mb / 1024),
       adj = c(0, 0.5), cex = 0.75, col = ink2, xpd = NA, lineheight = 0.85)
}
panel_title("(a) E1: fit time")

# (b), (c) E2: held-out RMSE and 95% coverage by grid size, mean over the 10
# designs with 95% t intervals; methods are offset slightly to keep the bars apart.
e2m <- list(spdbl = "spdbl_s2r_designgp", gu = "ppgasp_tuned", mlegp = "mlegp_pca")
pch2 <- c(spdbl = 16, gu = 17, mlegp = 15)
offset <- c(spdbl = 1, gu = 0.9, mlegp = 1.1)
e2_panel <- function(metric, ylim, yat, ref = NULL, title, ylab) {
  plot(NA, xlim = c(80, 3200), ylim = ylim, log = "x", axes = FALSE,
       xlab = expression("sites " * italic(S)), ylab = ylab)
  hgrid(yat)
  if (!is.null(ref)) abline(h = ref, lty = 2, col = ink2, lwd = 0.7)
  axis(1, at = c(100, 400, 2500), labels = c("10\u00b2", "20\u00b2", "50\u00b2"))
  axis(2, at = yat)
  for (m in names(e2m)) {
    d <- e2[e2$method == e2m[[m]], ]
    s <- do.call(rbind, lapply(split(d[[metric]], d$S), function(v) data.frame(
      mean = mean(v), half = stats::qt(0.975, length(v) - 1) * stats::sd(v) / sqrt(length(v)))))
    S <- as.numeric(rownames(s)) * offset[[m]]
    segments(S, s$mean - s$half, S, s$mean + s$half, col = col[m], lwd = 1.1)
    lines(S, s$mean, col = col[m], lwd = 1.6)
    points(S, s$mean, col = col[m], pch = pch2[m], cex = 0.8)
  }
  panel_title(title)
}
e2_panel("rmse", c(0.3, 1.1), seq(0.3, 1.1, 0.2), title = "(b) E2: RMSE", ylab = "held-out RMSE")
e2_panel("coverage95", c(0.7, 1), seq(0.7, 1, 0.1), ref = 0.95,
         title = "(c) E2: 95% coverage", ylab = "coverage")
legend("bottomright", legend = c("spDBL", "RobustGaSP (tuned)", "PCA-GP (mlegp)"),
       col = col, pch = pch2, lwd = 1.6, bty = "n", cex = 0.8, text.col = ink,
       seg.len = 1.4, inset = c(0, 0.02))

# (d) E3: relative error of calibrated parameters over 30 replicates x 3
# parameters, with 95% interval coverage above each method.
e3m <- c(spdbl = "spDBL", rc_sgasp = "S-GaSP", rc_gasp = "GaSP")
e3c <- c(spdbl = col[["spdbl"]], rc_sgasp = col[["gu"]], rc_gasp = col[["gu"]])
plot(NA, xlim = c(0.5, 3.5), ylim = c(-0.34, 0.64), axes = FALSE, xlab = "",
     ylab = "relative error")
hgrid(seq(-0.2, 0.4, 0.2)); abline(h = 0, col = ink2, lwd = 0.7)
axis(2, at = seq(-0.2, 0.4, 0.2))
axis(1, at = 1:3, labels = e3m, tick = FALSE, line = -0.4)
mtext("RobustCalibration", side = 1, at = 2.5, line = 1.15, cex = 0.72, col = ink2)
set.seed(1)
for (i in seq_along(e3m)) {
  d <- e3[e3$method == names(e3m)[i], ]
  points(i + runif(nrow(d), -0.22, 0.22), d$rel_error, pch = 16, cex = 0.45,
         col = adjustcolor(e3c[[i]], 0.55))
  q <- stats::quantile(d$rel_error, c(0.25, 0.5, 0.75))
  segments(i - 0.3, q[2], i + 0.3, q[2], col = ink, lwd = 1.4)
  segments(i, q[1], i, q[3], col = ink, lwd = 1.4)
  text(i, 0.6, sprintf("%.0f%%", 100 * mean(d$covered)), cex = 0.78, col = ink)
}
panel_title("(d) E3: calibration")

invisible(dev.off())
message("Wrote ", out)
