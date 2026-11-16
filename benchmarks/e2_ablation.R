# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E2 ablation: what does the dynamic emulator add to simulator-free
# prediction? On the paper-profile E2 data (same training and held-out runs),
# compares at times 3, ..., T - 1:
#   spdbl_designgp   spDBL's simulator-free prediction (from the E2 raw results)
#   designgp_only    the design GP alone: the lag-1 block of its predicted
#                    design row at time t + 1 is its prediction of the field at
#                    time t (no prediction exists for time T). Design rows are
#                    cbind(lag 2, lag 1), so lag 1 is columns C + 1, ..., 2C.
#   ppgasp_joint     RobustGaSP ppgasp, one kernel over all S x T outputs (E2)
#   ppgasp_per_time  ppgasp fitted separately at each time step
#   mlegp_pca        principal-component GP (E2)
# The emulator and design GP are refitted exactly as in E2 (seed 20261101).
#   Rscript benchmarks/e2_ablation.R [--grids=10,20,50] [--no-ppgasp] [--out=dir]
# Environment: SPDBL_BENCH_LIB, SPDBL_BENCH_CACHE (the E2 simulated data).
# Output: results/e2_ablation/paper/{e2_ablation.csv, e2_ablation_per_time.csv}

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
bench <- getOption("spdbl_bench_dir")
source(file.path(bench, "R", "common.R"))
args <- commandArgs(TRUE)
grids <- if (any(grepl("^--grids=", args)))
  as.numeric(strsplit(sub("^--grids=", "", grep("^--grids=", args, value = TRUE)), ",")[[1]]) else
  c(10, 20, 50)
use_ppgasp <- !("--no-ppgasp" %in% args)
bench_setup_lib(Sys.getenv("SPDBL_BENCH_LIB"))
suppressPackageStartupMessages(library(spDBL))
if (use_ppgasp) suppressPackageStartupMessages(library(RobustGaSP))
cache <- Sys.getenv("SPDBL_BENCH_CACHE", file.path(bench, "cache"))
e2_raw <- file.path(bench, "results", "e2", "paper", "raw")
out_dir <- if (any(grepl("^--out=", args))) sub("^--out=", "", grep("^--out=", args, value = TRUE)) else
  file.path(bench, "results", "e2_ablation", "paper")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seed <- 20261101

score <- function(truth, mean, sd, times) {
  y <- truth[, , times]; m <- mean[, , times]; s <- sd[, , times]
  ok <- is.finite(s) & s > 0
  data.frame(rmse = sqrt(mean((y - m)^2)),
    coverage95 = if (all(is.na(s))) NA else mean(abs(y - m) <= stats::qnorm(0.975) * s),
    crps = if (all(is.na(s))) NA else mean(ifelse(ok, crps_gaussian(y, m, pmax(s, 1e-12)), abs(y - m))))
}

rows <- list(); per_time <- list()
for (g in grids) {
  d <- readRDS(file.path(cache, sprintf("sir_g%03d_n40-50_T26_s%d.rds", g, seed)))
  T <- length(d$train); S <- g^2; n <- nrow(d$input_test)
  times <- 3:(T - 1)
  truth <- list_to_array(d$test)
  t0 <- proc.time()[["elapsed"]]
  em <- emulator_learn(d$input_train, d$train, Nx = g, Ny = g, nsam = 20, AR_choice = 2,
    episode_window = 2, covariance = "sigma2R", seed = seed)
  cal <- calibration_prepare(em, d$lower, d$upper)
  s <- em$setup; pr <- cal$predictor
  B <- s$nT / s$nT_ori; C <- s$bncol; p <- pr$p
  stopifnot(p == 2 * C, length(pr$center) == s$nT * p)
  # Design GP alone: episode (t + 1, b) has lag-1 columns C + 1..2C = field at t.
  M <- V <- array(NA_real_, c(n, S, T))
  ks <- seq_len((T - 1) * B)
  for (i in seq_len(n)) {
    des <- spDBL:::predict_offline_ar2(pr, d$input_test[i, ])
    ep_m <- lapply(ks, function(k) matrix(des$F[[B + k]][C + seq_len(C)], 1))
    ep_v <- lapply(ks, function(k) matrix(des$variance[(B + k - 1) * p + C + seq_len(C)], 1))
    fm <- recover_from_EP_exact(ep_m, T - 1, (T - 1) * B, block_cols = s$block_cols)
    fv <- recover_from_EP_exact(ep_v, T - 1, (T - 1) * B, block_cols = s$block_cols)
    for (t in seq_len(T - 1)) { M[i, , t] <- fm[[t]]; V[i, , t] <- fv[[t]] }
  }
  dg_s <- proc.time()[["elapsed"]] - t0
  sc <- score(truth, M, sqrt(V), times)
  rows[[length(rows) + 1]] <- cbind(grid = g, method = "designgp_only", sc, seconds = dg_s)
  per_time[[length(per_time) + 1]] <- data.frame(grid = g, method = "designgp_only",
    time = seq_len(T), rmse = sqrt(apply((truth - M)^2, 3, mean)))
  rm(em, cal); invisible(gc())

  if (use_ppgasp) {
    t0 <- proc.time()[["elapsed"]]
    active <- d$upper > d$lower
    X <- d$input_train[, active, drop = FALSE]; Xn <- d$input_test[, active, drop = FALSE]
    PM <- PS <- array(NA_real_, c(n, S, T))
    for (t in seq_len(T)) {
      Y <- as.matrix(d$train[[t]])
      const <- apply(Y, 2, function(v) max(v) - min(v) < 1e-12)
      PM[, , t] <- matrix(rep(Y[1, ], each = n), n, S); PS[, , t] <- 0
      if (any(!const)) {
        utils::capture.output(fit <- ppgasp(design = X, response = Y[, !const, drop = FALSE]))
        pp <- predict(fit, Xn)
        PM[, !const, t] <- pp$mean; PS[, !const, t] <- pp$sd
      }
    }
    sc <- score(truth, PM, PS, times)
    rows[[length(rows) + 1]] <- cbind(grid = g, method = "ppgasp_per_time", sc,
      seconds = proc.time()[["elapsed"]] - t0)
    per_time[[length(per_time) + 1]] <- data.frame(grid = g, method = "ppgasp_per_time",
      time = seq_len(T), rmse = sqrt(apply((truth - PM)^2, 3, mean)))
  }

  # E2's own per-time RMSE for the other methods (same data, held-out runs, cells).
  for (m in c("spdbl_s2r_designgp", "ppgasp_default", "mlegp_pca")) {
    f <- file.path(e2_raw, sprintf("e2_g%03d_%s.rds", g, m))
    if (!file.exists(f)) next
    pt <- readRDS(f)$per_time_rmse
    lab <- c(spdbl_s2r_designgp = "spdbl_designgp", ppgasp_default = "ppgasp_joint",
             mlegp_pca = "mlegp_pca")[[m]]
    rows[[length(rows) + 1]] <- data.frame(grid = g, method = lab,
      rmse = sqrt(mean(pt[times]^2)), coverage95 = NA, crps = NA, seconds = NA)
    per_time[[length(per_time) + 1]] <- data.frame(grid = g, method = lab,
      time = seq_len(T), rmse = pt)
  }
  message(sprintf("grid %d done", g))
}
tab <- do.call(rbind, rows)
utils::write.csv(tab, file.path(out_dir, "e2_ablation.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, per_time), file.path(out_dir, "e2_ablation_per_time.csv"),
  row.names = FALSE)
print(tab[order(tab$grid, tab$rmse), ], row.names = FALSE, digits = 4)
