# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Figure 1 of the MLOSS paper: workflow functions, the objects they return,
# and the exported building blocks underneath.
# Run from the paper folder (it writes figures/architecture.pdf there):
#   Rscript <spDBL>/benchmarks/paper/architecture.R
# Drawn at the JMLR text width (6 in) so font sizes are the printed sizes;
# cairo_pdf embeds the fonts. R rounds text sizes to whole points, so sizes
# are 7 pt (titles, headings) and 6 pt (details). Label widths are checked on
# a pdf() device, whose font metrics are exact; Courier's advance (0.600 em)
# stands in for Menlo's (0.602 em).
dir.create("figures", showWarnings = FALSE)
out <- file.path("figures", "architecture.pdf")
W <- 6.0; H <- 2.75
sans <- "Helvetica"; mono <- "Menlo"
big <- 1; small <- 6 / 7   # cex for 7 pt and 6 pt at pointsize 7

ink <- "#0b0b0b"; ink2 <- "#3d3c39"; ink3 <- "#6b6a65"
blue <- "#2a78d6"; blue_dark <- "#1c5cab"; band_fill <- "#eef4fc"; grey <- "#9a9993"

pdf(NULL, width = W, height = H, pointsize = 7)
plot.new(); measure_dev <- dev.cur()
cairo_pdf(out, width = W, height = H, family = sans, pointsize = 7)
draw_dev <- dev.cur()
par(mar = c(0, 0, 0, 0), xpd = NA)
plot.new()
plot.window(xlim = c(0, W), ylim = c(0, H), xaxs = "i", yaxs = "i")

# One unit is one inch: four columns across the full width.
x0 <- 0.04; x1 <- W - 0.04; gap <- 0.09; pad <- 0.04
box_w <- (x1 - x0 - 2 * 0.06 - 3 * gap) / 4
cx <- x0 + 0.06 + box_w / 2 + (0:3) * (box_w + gap)
bh <- 0.56; line_step <- 0.112

rrect <- function(xa, ya, xb, yb, r = 0.045, ...) {
  a <- seq(0, pi / 2, length.out = 8)
  xs <- c(xb - r + r * cos(a), xa + r - r * sin(a), xa + r - r * cos(a), xb - r + r * sin(a))
  ys <- c(yb - r + r * sin(a), yb - r + r * cos(a), ya + r - r * sin(a), ya + r - r * cos(a))
  polygon(xs, ys, ...)
}
fits <- function(s, cex, family, font = 1, width = box_w - 2 * pad) {
  dev.set(measure_dev)
  w <- strwidth(s, units = "inches", cex = cex, font = font,
                family = if (family == mono) "Courier" else "Helvetica") * 1.004
  dev.set(draw_dev)
  if (w > width) stop(sprintf("label too wide (%.3f > %.3f in): %s", w, width, s))
}
# A box: title, description lines and, in blue, the class of the returned object.
fbox <- function(x, y, title, desc, ret = NULL, border = blue, lty = 1,
                 title_family = mono, fill = "white") {
  rrect(x - box_w / 2, y - bh / 2, x + box_w / 2, y + bh / 2, col = fill,
        border = border, lwd = 0.9, lty = lty)
  ls <- list(list(title, title_family, 2, big, ink))
  for (d in desc) ls <- c(ls, list(list(d, sans, 1, small, ink2)))
  if (!is.null(ret)) ls <- c(ls, list(list(ret, mono, 1, small, blue_dark)))
  n <- length(ls)
  ys <- y + ((n - 1) / 2 - seq_len(n) + 1) * line_step
  for (i in seq_len(n)) {
    l <- ls[[i]]
    fits(l[[1]], l[[4]], l[[2]], l[[3]])
    text(x, ys[i], l[[1]], family = l[[2]], font = l[[3]], cex = l[[4]], col = l[[5]])
  }
}
arrow <- function(xa, ya, xb, yb, lty = 1)
  arrows(xa, ya, xb, yb, length = 0.045, angle = 22, lwd = 0.9, col = ink2, lty = lty)
band <- function(lo, hi, label, note, dashed = FALSE) {
  if (dashed) rrect(x0, lo, x1, hi, r = 0.03, col = NA, border = grey, lwd = 0.7, lty = 2)
  else rect(x0, lo, x1, hi, col = band_fill, border = NA)
  text(x0 + 0.06, hi - 0.09, label, adj = c(0, 0.5), family = sans, font = 2, cex = big, col = ink)
  dev.set(measure_dev)
  lw <- strwidth(label, units = "inches", family = "Helvetica", font = 2, cex = big)
  dev.set(draw_dev)
  text(x0 + 0.06 + lw + 0.08, hi - 0.09, note, adj = c(0, 0.5), family = sans, font = 3,
       cex = small, col = ink3)
}

# Bands, top to bottom; boxes sit below each band's heading.
b1 <- c(1.88, 2.71); b2 <- c(0.97, 1.80); b3 <- c(0.04, 0.89)
band(b1[1], b1[2], "Emulation", "closed form, no MCMC")
band(b2[1], b2[2], "Calibration", "no new simulator runs")
band(b3[1], b3[2], "Building blocks", "exported, for new emulators", dashed = TRUE)
y1 <- b1[1] + 0.05 + bh / 2; y2 <- b2[1] + 0.05 + bh / 2; y3 <- b3[1] + 0.06 + 0.31

# Emulation row
fbox(cx[1], y1, "prepare_data()", c("reads simulator runs,", "training and test split"))
fbox(cx[2], y1, "emulator_learn()", c("closed-form filter, smoother"), "spDBL_emulator")
fbox(cx[3], y1, "emulator_predict()", c("fields at new inputs,", "exact or Monte Carlo"))
fbox(cx[4], y1, "S3 methods", c("print(), summary(),", "predict(), plot()"),
     border = grey, title_family = sans)
arrow(cx[1] + box_w / 2, y1, cx[2] - box_w / 2, y1)
arrow(cx[2] + box_w / 2, y1, cx[3] - box_w / 2, y1)

# Calibration row
fbox(cx[1], y2, "calibration_validate()", c("optional held-out check", "of predicted fields"), lty = 2)
fbox(cx[2], y2, "calibration_prepare()", c("design GP, Matérn 5/2"), "spDBL_calibration_input")
fbox(cx[3], y2, "calibration_fit()", c("field data, NUTS per draw"), "spDBL_calibration")
fbox(cx[4], y2, "calibration_predict()", c("draws of fields y and z,", "discrepancy u"))
arrow(cx[2] - box_w / 2, y2, cx[1] + box_w / 2, y2)
arrow(cx[2] + box_w / 2, y2, cx[3] - box_w / 2, y2)
arrow(cx[3] + box_w / 2, y2, cx[4] - box_w / 2, y2)
# Posterior emulator draws feed the calibration.
arrow(cx[2], y1 - bh / 2, cx[2], y2 + bh / 2)
text(cx[2] + 0.05, (y1 - bh / 2 + y2 + bh / 2) / 2 + 0.01, "posterior draws",
     adj = c(0, 0.5), family = sans, font = 3, cex = small, col = ink3)

# Building blocks
bb <- list(
  list("C++ kernels", c("FF(), BS(), FFBS()", "via RcppEigen")),
  list("Out of core", c("FF/BS_bigdata_R()", "stream CSV blocks")),
  list("Stan programs", c("marginal, explicit u;", "R reference densities")),
  list("Samplers, kernels", c("MNIW_sampler()", "gen_gp_kernel()")))
for (k in 1:4) {
  b <- bb[[k]]
  rrect(cx[k] - box_w / 2, y3 - 0.31, cx[k] + box_w / 2, y3 + 0.31, col = "white",
        border = grey, lwd = 0.7)
  fits(b[[1]], big, sans, 2)
  text(cx[k], y3 + line_step, b[[1]], family = sans, font = 2, cex = big, col = ink)
  for (i in 1:2) {
    fam <- if (grepl("\\(", b[[2]][i])) mono else sans
    fits(b[[2]][i], small, fam)
    text(cx[k], y3 - (i - 1) * line_step, b[[2]][i], family = fam, cex = small, col = ink2)
  }
}

invisible(dev.off(draw_dev)); invisible(dev.off(measure_dev))
message("Wrote ", out)
