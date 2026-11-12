# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E2 worker: one (grid, method) emulator fit + held-out prediction.
# Called by e2_accuracy.R as: Rscript e2_worker.R <task.rds>
#
# Methods
#   spdbl_iw_known          spDBL IW emulator, held-out AR rows built from the
#                           true held-out trajectories (conditional/oracle).
#   spdbl_s2r_known         as above with covariance = "sigma2R".
#   spdbl_s2r_designgp_lazy the same with emulator_learn(keep_draws = FALSE)
#   spdbl_s2r_designgp      sigma2R emulator + calibration_prepare(); AR rows
#                           predicted from eta by the offline design GP
#                           (PDE-free, comparable to ppgasp). Package defaults.
#   spdbl_s2r_designgp_unit tuning pass: same, inputs rescaled to [0, 1].
#   designgp_only           ablation: the design GP's own prediction of each
#                           field (no dynamic model); no prediction at time T.
#   ppgasp_per_time         ablation: ppgasp fitted separately at each time.
#   ppgasp_default          RobustGaSP::ppgasp, package defaults.
#   ppgasp_tuned            tuning pass: nugget.est = TRUE and
#                           optimization = "nelder-mead" (ppgasp itself warns
#                           that the derivative-based optimiser can be
#                           inaccurate for many outputs).
#   mlegp_pca               mlegp GPs on principal-component weights of the
#                           centred outputs (Higdon-style PCA-GP).

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "spdbl_predictive.R"))

# Runs x (S * T) matrix, column (t - 1) * S + s, from a time-indexed list.
flatten_runs <- function(x) do.call(cbind, lapply(x, as.matrix))

fit_spdbl <- function(d, cfg, task, covariance, design, unit = FALSE, keep_draws = TRUE) {
  suppressPackageStartupMessages(library(spDBL))
  X_train <- d$input_train
  X_test <- d$input_test
  lower <- d$lower
  upper <- d$upper
  if (unit) {
    X_train <- unit_inputs(X_train, d$lower, d$upper)
    X_test <- unit_inputs(X_test, d$lower, d$upper)
    lower <- stats::setNames(rep(0, length(d$lower)), names(d$lower))
    upper <- stats::setNames(as.numeric(d$upper > d$lower), names(d$lower))
  }
  fit <- timed({
    em <- emulator_learn(X_train, d$train, Nx = d$grid, Ny = d$grid,
      nsam = cfg$nsam, AR_choice = cfg$ar_choice,
      episode_window = cfg$episode_window, covariance = covariance,
      seed = task$seed, keep_draws = keep_draws)
    if (design == "designgp") list(em = em, cal = calibration_prepare(em, lower, upper)) else
      list(em = em)
  })
  heap_fit <- gc_peak_mb()
  invisible(gc(reset = TRUE))
  pred <- timed(if (design == "designgp")
    spdbl_designgp_predictive(fit$value$cal, X_test) else
    spdbl_known_predictive(fit$value$em, X_test,
      spdbl_known_designs(fit$value$em, d$test, task$seed)))
  extra <- list()
  if (design == "designgp") {
    cal <- fit$value$cal
    extra$design_gp_loo_rmse <- cal$predictor$loo_rmse
    # One row of lengthscales per kernel group (time step); medians recorded.
    ell <- matrix(cal$predictor$lengthscale, ncol = length(cal$predictor$active))
    extra$design_gp_lengthscale <- paste(signif(apply(ell, 2, stats::median), 3),
                                         collapse = ";")
    extra$design_gp_kernel <- if (is.null(cal$design_kernel)) "sqexp" else cal$design_kernel
    extra$design_gp_groups <- if (is.null(cal$design_groups)) "shared" else cal$design_groups
    extra$design_uncertainty <- isTRUE(cal$design_uncertainty)
  }
  extra$fit_object_mb <- as.numeric(utils::object.size(fit$value)) / 2^20
  list(pred = pred$value, fit_s = fit$elapsed, fit_cpu_s = fit$cpu,
       predict_s = pred$elapsed, predict_cpu_s = pred$cpu,
       heap_fit = heap_fit, extra = extra)
}

# Ablation: the design GP alone. Its predicted design row for episode
# (t + 1, b) is cbind(lag 2, lag 1); the lag-1 half is its prediction of the
# field at time t (none exists for the last time, which stays NA).
fit_designgp_only <- function(d, cfg, task) {
  suppressPackageStartupMessages(library(spDBL))
  fit <- timed({
    em <- emulator_learn(d$input_train, d$train, Nx = d$grid, Ny = d$grid,
      nsam = cfg$nsam, AR_choice = cfg$ar_choice,
      episode_window = cfg$episode_window, covariance = "sigma2R", seed = task$seed)
    list(em = em, cal = calibration_prepare(em, d$lower, d$upper))
  })
  s <- fit$value$em$setup
  pr <- fit$value$cal$predictor
  n <- nrow(d$input_test); S <- d$grid^2; T <- length(d$train)
  B <- s$nT / s$nT_ori; C <- s$bncol
  stopifnot(pr$p == 2 * C)
  block_cols <- s$block_cols
  pred <- timed({
    M <- V <- array(NA_real_, c(n, S, T))
    ks <- seq_len((T - 1) * B)
    for (i in seq_len(n)) {
      des <- spDBL:::predict_offline_ar2(pr, d$input_test[i, ])
      ep_m <- lapply(ks, function(k) matrix(des$F[[B + k]][C + seq_len(C)], 1))
      ep_v <- lapply(ks, function(k) matrix(des$variance[(B + k - 1) * pr$p + C + seq_len(C)], 1))
      fm <- recover_from_EP_exact(ep_m, T - 1, (T - 1) * B, block_cols = block_cols)
      fv <- recover_from_EP_exact(ep_v, T - 1, (T - 1) * B, block_cols = block_cols)
      for (t in seq_len(T - 1)) { M[i, , t] <- fm[[t]]; V[i, , t] <- fv[[t]] }
    }
    sd <- sqrt(V)
    list(mean = M, sd = sd, lower = M - stats::qnorm(0.975) * sd,
         upper = M + stats::qnorm(0.975) * sd)
  })
  list(pred = pred$value, fit_s = fit$elapsed, fit_cpu_s = fit$cpu,
       predict_s = pred$elapsed, predict_cpu_s = pred$cpu, heap_fit = NA_real_,
       extra = list(last_time_predicted = FALSE))
}

# Ablation: ppgasp fitted separately at each time step (same defaults).
fit_ppgasp_per_time <- function(d, cfg, task) {
  suppressPackageStartupMessages(library(RobustGaSP))
  active <- d$upper > d$lower
  X <- d$input_train[, active, drop = FALSE]
  X_new <- d$input_test[, active, drop = FALSE]
  n <- nrow(X_new); S <- d$grid^2; T <- length(d$train)
  out <- list(mean = array(NA_real_, c(n, S, T)))
  out$sd <- out$lower <- out$upper <- out$mean
  tm <- timed(for (t in seq_len(T)) {
    Y <- as.matrix(d$train[[t]])
    const <- apply(Y, 2, function(v) max(v) - min(v) < 1e-12)
    out$mean[, , t] <- out$lower[, , t] <- out$upper[, , t] <- matrix(rep(Y[1, ], each = n), n, S)
    out$sd[, , t] <- 0
    if (any(!const)) {
      utils::capture.output(fit <- ppgasp(design = X, response = Y[, !const, drop = FALSE]))
      pp <- predict(fit, X_new)
      out$mean[, !const, t] <- pp$mean; out$sd[, !const, t] <- pp$sd
      out$lower[, !const, t] <- pp$lower95; out$upper[, !const, t] <- pp$upper95
    }
  })
  list(pred = out, fit_s = tm$elapsed, fit_cpu_s = tm$cpu, predict_s = NA_real_,
       predict_cpu_s = NA_real_, heap_fit = NA_real_, extra = list())
}

fit_ppgasp <- function(d, cfg, task, tuned) {
  suppressPackageStartupMessages(library(RobustGaSP))
  active <- d$upper > d$lower
  Y <- flatten_runs(d$train)
  const <- apply(Y, 2, function(v) max(v) - min(v) < 1e-12)
  X <- d$input_train[, active, drop = FALSE]
  X_new <- d$input_test[, active, drop = FALSE]
  fit <- timed(if (tuned)
    ppgasp(design = X, response = Y[, !const, drop = FALSE], nugget.est = TRUE,
           optimization = "nelder-mead") else
    ppgasp(design = X, response = Y[, !const, drop = FALSE]))
  heap_fit <- gc_peak_mb()
  invisible(gc(reset = TRUE))
  pred <- timed(predict(fit$value, X_new))
  # ppgasp refuses constant output columns (identical initial condition at
  # t = 0), so they are removed before fitting and predicted exactly here.
  n <- nrow(X_new)
  fill <- function(v, constant_value) {
    out <- matrix(rep(Y[1, ], each = n), n, ncol(Y))
    out[, !const] <- v
    if (!is.null(constant_value)) out[, const] <- constant_value
    out
  }
  S <- d$grid^2
  T <- length(d$train)
  shape <- function(m) array(m, c(n, S, T))
  p <- pred$value
  list(pred = list(mean = shape(fill(p$mean, NULL)), sd = shape(fill(p$sd, 0)),
                   lower = shape(fill(p$lower95, NULL)),
                   upper = shape(fill(p$upper95, NULL))),
       fit_s = fit$elapsed, fit_cpu_s = fit$cpu, predict_s = pred$elapsed,
       predict_cpu_s = pred$cpu, heap_fit = heap_fit,
       extra = list(n_constant_columns = sum(const),
         ppgasp_beta = paste(signif(fit$value@beta_hat, 3), collapse = ";"),
         ppgasp_nugget = fit$value@nugget,
         fit_object_mb = as.numeric(utils::object.size(fit$value)) / 2^20))
}

fit_mlegp <- function(d, cfg, task) {
  suppressPackageStartupMessages(library(mlegp))
  Y <- flatten_runs(d$train)
  const <- apply(Y, 2, function(v) max(v) - min(v) < 1e-12)
  active <- d$upper > d$lower
  X <- unit_inputs(d$input_train, d$lower, d$upper)[, active, drop = FALSE]
  X_new <- unit_inputs(d$input_test, d$lower, d$upper)[, active, drop = FALSE]
  centre <- colMeans(Y[, !const, drop = FALSE])
  fit <- timed({
    Z <- t(sweep(Y[, !const, drop = FALSE], 2, centre))
    utils::capture.output(m <- mlegp(X, Z, PC.percent = cfg$mlegp_pc_percent,
                                     verbose = 0, seed = task$seed))
    m
  })
  heap_fit <- gc_peak_mb()
  invisible(gc(reset = TRUE))
  pred <- timed({
    m <- fit$value
    w <- lapply(seq_len(m$numGPs), function(j)
      predict(m[[j]], newData = X_new, se.fit = TRUE))
    W <- sapply(w, function(o) as.numeric(o$fit))
    W_se <- sapply(w, function(o) as.numeric(o$se.fit))
    W <- matrix(W, nrow(X_new)); W_se <- matrix(W_se, nrow(X_new))
    list(mean = sweep(W %*% t(m$UD), 2, centre, "+"),
         sd = sqrt(W_se^2 %*% t(m$UD^2)))
  })
  n <- nrow(X_new)
  S <- d$grid^2
  T <- length(d$train)
  fill <- function(v, constant_value) {
    out <- matrix(rep(Y[1, ], each = n), n, ncol(Y))
    out[, !const] <- v
    if (!is.null(constant_value)) out[, const] <- constant_value
    array(out, c(n, S, T))
  }
  mean <- fill(pred$value$mean, NULL)
  sd <- fill(pred$value$sd, 0)
  list(pred = list(mean = mean, sd = sd, lower = mean - 1.96 * sd,
                   upper = mean + 1.96 * sd),
       fit_s = fit$elapsed, fit_cpu_s = fit$cpu, predict_s = pred$elapsed,
       predict_cpu_s = pred$cpu, heap_fit = heap_fit,
       extra = list(n_pc = fit$value$numGPs,
         fit_object_mb = as.numeric(utils::object.size(fit$value)) / 2^20))
}

bench_worker_main(function(task) {
  cfg <- task$cfg
  bench_setup_lib()
  d <- readRDS(task$data_file)
  invisible(gc(reset = TRUE))
  res <- switch(task$method,
    spdbl_iw_known = fit_spdbl(d, cfg, task, "iw", "known"),
    spdbl_s2r_known = fit_spdbl(d, cfg, task, "sigma2R", "known"),
    spdbl_s2r_designgp = fit_spdbl(d, cfg, task, "sigma2R", "designgp"),
    spdbl_s2r_designgp_unit = fit_spdbl(d, cfg, task, "sigma2R", "designgp", unit = TRUE),
    spdbl_s2r_designgp_lazy = fit_spdbl(d, cfg, task, "sigma2R", "designgp", keep_draws = FALSE),
    ppgasp_default = fit_ppgasp(d, cfg, task, tuned = FALSE),
    ppgasp_tuned = fit_ppgasp(d, cfg, task, tuned = TRUE),
    mlegp_pca = fit_mlegp(d, cfg, task),
    designgp_only = fit_designgp_only(d, cfg, task),
    ppgasp_per_time = fit_ppgasp_per_time(d, cfg, task),
    stop("Unknown method ", task$method))
  truth <- list_to_array(d$test)
  mask <- array(FALSE, dim(truth))
  mask[, , seq_len(dim(truth)[3]) >= cfg$eval_from] <- TRUE
  p <- res$pred
  if (!identical(dim(p$mean), dim(truth))) stop("Prediction has the wrong shape")
  # Methods without a prediction at some times (designgp_only at time T) are
  # scored where they predict; per-time scores allow a common time range.
  mask <- mask & !is.na(p$mean)
  scores <- score_predictions(truth, p$mean, p$sd, p$lower, p$upper, mask)
  row <- cbind(data.frame(grid = d$grid, S = d$grid^2,
    design = if (is.null(task$design)) 1L else task$design, method = task$method,
    n_train = nrow(d$input_train), n_test = nrow(d$input_test),
    n_time = length(d$train), eval_from = cfg$eval_from), scores,
    data.frame(fit_s = res$fit_s, fit_cpu_s = res$fit_cpu_s,
      predict_s = res$predict_s, predict_cpu_s = res$predict_cpu_s,
      heap_peak_fit_mb = res$heap_fit))
  for (k in names(res$extra)) row[[k]] <- res$extra[[k]]
  # Per-time scores over held-out runs and sites, kept in the raw result.
  per_time <- sqrt(apply((truth - p$mean)^2, 3, mean))
  per_time_scores <- do.call(rbind, lapply(seq_len(dim(truth)[3]), function(t)
    if (all(is.na(p$mean[, , t]))) data.frame(time = t, rmse = NA, coverage95 = NA,
      width95 = NA, crps = NA) else cbind(time = t, score_predictions(
        truth[, , t, drop = FALSE], p$mean[, , t, drop = FALSE], p$sd[, , t, drop = FALSE],
        p$lower[, , t, drop = FALSE], p$upper[, , t, drop = FALSE])[,
        c("rmse", "coverage95", "width95", "crps")])))
  list(row = row, per_time_rmse = per_time, per_time = per_time_scores)
})
