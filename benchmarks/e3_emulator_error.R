# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E3 ablation, no MCMC: how well does the calibration emulator predict the
# simulator at each replicate's true parameters, and are its errors
# correlated across time? calibration_fit() treats the design-row (design
# GP) error as independent across times; if the errors persist over time,
# the field likelihood counts the same error many times and the parameter
# posterior is too narrow.
#   Rscript benchmarks/e3_emulator_error.R [--profile=mode]
# Uses the emulator, design GP and replicates of the E3 profile (as
# e3_worker.R does). For each replicate, predict(cal, eta, interval =
# "prediction") at the true eta gives pointwise predictive means and
# standard deviations at every time and site (moment-matched over the
# emulator draws, with the design-GP variance), which are compared with the
# simulator output y used to generate the field data.
# Output in results/e3_ablation/emulator_error/: errors.csv (one row per
# replicate, time and site) and summary.csv.

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))
args <- bench_parse_args()
if (is.null(args$profile)) args$profile <- "mode"
bench_setup_lib(Sys.getenv("SPDBL_BENCH_LIB"))
cfg <- bench_config("e3", args)
suppressPackageStartupMessages(library(spDBL))
out <- file.path(bench_dir(), "results", "e3_ablation", "emulator_error")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(sir_calibration_replicates(cfg$n_rep, cfg$seed, cfg$box, cfg$cache,
  inner = cfg$inner, n0 = cfg$field_n0, d0 = cfg$field_d0))
em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
  covariance = "sigma2R", episode_window = cfg$episode_window,
  nsam = cfg$n_draws, seed = cfg$seed)
cal <- calibration_prepare(em, d$lower, d$upper)

rows <- list()
for (r in seq_along(d$replicates)) {
  rep <- d$replicates[[r]]
  p <- predict(cal, matrix(rep$eta, nrow = 1), interval = "prediction")
  n_t <- length(p$mean)
  # Predictions cover the last n_t times of the simulator output.
  times <- seq(nrow(rep$y) - n_t + 1, nrow(rep$y))
  for (k in seq_len(n_t)) {
    y <- rep$y[times[k], ]
    rows[[length(rows) + 1L]] <- data.frame(rep = r, time = times[k] - 1L,
      site = seq_along(y), y = y, mean = p$mean[[k]][1, ], sd = p$sd[[k]][1, ])
  }
}
e <- do.call(rbind, rows)
e$z <- (e$y - e$mean) / e$sd
e <- e[is.finite(e$z) & e$sd > 0, ]
e <- e[order(e$rep, e$site, e$time), ]
utils::write.csv(e, file.path(out, "errors.csv"), row.names = FALSE)

# Lag-k correlation of the standardised errors within each (replicate, site)
# series, pooled over series.
lag_cor <- function(k) {
  a <- b <- numeric()
  for (g in split(e, list(e$rep, e$site), drop = TRUE)) {
    g <- g[order(g$time), ]
    if (nrow(g) > k) { a <- c(a, g$z[seq_len(nrow(g) - k)]); b <- c(b, g$z[-seq_len(k)]) }
  }
  stats::cor(a, b)
}
rho <- vapply(1:5, lag_cor, numeric(1))
# Variance of a sum of T errors with these correlations, relative to the
# independent case (AR(1) approximation with the lag-1 correlation).
n_time <- length(unique(e$time))
r1 <- rho[1]
inflation <- 1 + 2 * sum((1 - seq_len(n_time - 1) / n_time) * r1^seq_len(n_time - 1))
summ <- data.frame(n_rep = length(unique(e$rep)), n_time = n_time,
  n_errors = nrow(e), coverage95 = mean(abs(e$z) <= stats::qnorm(0.975)),
  sd_z = stats::sd(e$z), mean_z = mean(e$z),
  lag1 = rho[1], lag2 = rho[2], lag3 = rho[3], lag4 = rho[4], lag5 = rho[5],
  variance_inflation_ar1 = inflation)

# Per emulator draw, as each conditional fit sees it: y_t given draw m is
# Gaussian with mean mu_t and covariance Sigma + design-GP term
# (.cal_time_covariance), independently across t. Whitened errors
# w_t = U_t^{-T} (y_t - mu_t), with Cov_t = U_t'U_t, should be iid N(0, 1).
grid <- d$grid
adjacent <- which(outer(seq_len(nrow(grid)), seq_len(nrow(grid)), function(i, j)
  abs(grid$row[i] - grid$row[j]) + abs(grid$col[i] - grid$col[j]) == 1 & i < j), arr.ind = TRUE)
per <- list()
for (r in seq_along(d$replicates)) {
  rep <- d$replicates[[r]]
  for (m in seq_len(cal$M)) {
    pred <- spDBL:::predict_calibration(rep$eta, m, cal)
    n_t <- nrow(pred$mu); times <- seq(nrow(rep$y) - n_t + 1, nrow(rep$y))
    zs <- matrix(NA_real_, n_t, ncol(pred$mu))
    for (k in seq_len(n_t)) {
      C <- spDBL:::.cal_time_covariance(pred, k)
      err <- rep$y[times[k], ] - pred$mu[k, ]
      if (all(abs(err) < 1e-10)) next   # the fixed initial field
      U <- chol(C)
      w <- backsolve(U, err, transpose = TRUE)
      sdv <- sqrt(diag(C)); zs[k, ] <- err / sdv
      R <- C / outer(sdv, sdv)
      per[[length(per) + 1L]] <- data.frame(rep = r, draw = m, time = times[k] - 1L,
        mean_w2 = mean(w^2), coverage95 = mean(abs(err / sdv) <= stats::qnorm(0.975)),
        adj_cor_model = mean(R[adjacent]),
        adj_cor_observed = stats::cor(zs[k, adjacent[, 1]], zs[k, adjacent[, 2]]))
    }
    ok <- stats::complete.cases(zs)
    per[[length(per)]]$lag1_draw <- if (sum(ok) > 2)
      stats::cor(as.vector(zs[ok, ][-1, ]), as.vector(zs[ok, ][-sum(ok), ])) else NA_real_
  }
}
pd <- do.call(rbind, lapply(per, function(x) { if (is.null(x$lag1_draw)) x$lag1_draw <- NA_real_; x }))
utils::write.csv(pd, file.path(out, "per_draw.csv"), row.names = FALSE)
summ$draw_mean_w2 <- mean(pd$mean_w2)
summ$draw_coverage95 <- mean(pd$coverage95)
summ$draw_lag1 <- mean(pd$lag1_draw, na.rm = TRUE)
summ$adjacent_cor_model <- mean(pd$adj_cor_model)
summ$adjacent_cor_observed <- mean(pd$adj_cor_observed)
utils::write.csv(summ, file.path(out, "summary.csv"), row.names = FALSE)
cfg$out <- out
bench_write_meta(cfg, "meta")
print(t(summ), digits = 3)
