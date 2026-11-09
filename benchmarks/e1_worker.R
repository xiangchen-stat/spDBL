# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E1 worker: one (grid, episode_window) configuration in a fresh process.
# Called by e1_scaling.R as: Rscript e1_worker.R <task.rds>

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "spdbl_predictive.R"))

bench_worker_main(function(task) {
  cfg <- task$cfg
  bench_setup_lib()
  suppressPackageStartupMessages(library(spDBL))
  d <- readRDS(task$data_file)
  g <- task$grid
  w <- task$window
  invisible(gc(reset = TRUE))
  rss_start <- rss_now_mb()

  fit <- timed(emulator_learn(d$input_train, d$train, Nx = g, Ny = g,
    nsam = cfg$nsam, AR_choice = cfg$ar_choice, episode_window = w,
    covariance = cfg$covariance, seed = task$seed,
    keep_draws = if (is.null(cfg$keep_draws)) TRUE else isTRUE(cfg$keep_draws)))
  heap_fit <- gc_peak_mb()
  rss_after_fit <- rss_now_mb()
  em <- fit$value
  s <- em$setup
  size_total <- as.numeric(utils::object.size(em)) / 2^20
  size_rfield <- as.numeric(utils::object.size(s$R_field)) / 2^20

  invisible(gc(reset = TRUE))
  F_new <- spdbl_known_designs(em, d$test, task$seed)
  pred <- timed(emulator_predict(em, d$input_test, F_new_ls = F_new))
  heap_pred <- gc_peak_mb()
  truth <- list_to_array(d$test)
  mean <- list_to_array(pred$value)
  keep <- seq_len(dim(truth)[3]) >= cfg$eval_from
  rmse <- sqrt(mean((truth[, , keep] - mean[, , keep])^2))

  pred_mc_s <- NA_real_
  if (isTRUE(cfg$predict_mc)) {
    set.seed(task$seed)
    pred_mc_s <- timed(emulator_predict(em, d$input_test, F_new_ls = F_new,
                                        MC = TRUE))$elapsed
  }

  # Cost of the dense S x S field covariance that spDBL 1.1.1 built and
  # factorised in every emulator_learn() call. Since 1.2.0 it is built only
  # for calibration, at the observed sites; timed here to show the saving.
  rfield_s <- NA_real_
  if (isTRUE(cfg$measure_rfield)) {
    fit$value <- pred$value <- em$setup$R_field <- NULL
    invisible(gc())
    rfield_s <- timed({
      R <- s$gp_sigma2 * exp(-s$phi_sp * as.matrix(stats::dist(s$ind_sp))) +
        diag(s$gp_tau2, nrow(s$ind_sp))
      chol(R)
      NULL
    })$elapsed
  }

  list(row = data.frame(grid = g, S = g * g, episode_window = w,
    block_cols = s$bncol, n_blocks = s$n_b, episodes = s$nT,
    covariance = cfg$covariance, n_train = nrow(d$input_train),
    n_test = nrow(d$input_test), n_time = length(d$train), nsam = cfg$nsam,
    fit_s = fit$elapsed, fit_cpu_s = fit$cpu, predict_s = pred$elapsed,
    predict_cpu_s = pred$cpu, predict_mc_s = pred_mc_s,
    rfield_build_chol_s = rfield_s, rfield_share = rfield_s / fit$elapsed,
    object_mb = size_total, rfield_mb = size_rfield,
    heap_peak_fit_mb = heap_fit, heap_peak_predict_mb = heap_pred,
    rss_start_mb = rss_start, rss_after_fit_mb = rss_after_fit,
    rmse_known_ar = rmse))
})
