# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E3 site-prediction worker: one replicate. Called by e3_sites.R as:
#   Rscript e3_sites_worker.R <task.rds>
#
# Leaves cfg$site_holdout grid sites out of the calibration, predicts them
# with calibration_predict(sites = ...) at every observation time, and scores
# the predictions of the held-out observations z and of the noise-free field
# y + u (the simulator at the true parameters plus the true discrepancy)
# against two baselines:
#   spdbl           composition draws of y, u and z at the held-out sites
#   simulator_only  the calibrated emulator at the held-out sites without the
#                   discrepancy: y from predict_calibration() at each
#                   posterior (eta, emulator draw), z = y + N(0, tau2_t)
#   idw             inverse-distance-weighted (power 2) interpolation of the
#                   observed field at the same time; a point prediction

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))

simulate_runs <- function(X, d, n_people) {
  runs <- lapply(seq_len(nrow(X)), function(i)
    spDBL::gen_pde(X[i, ], d$Nx, d$Ny, n_people, length(d$time)))
  lapply(seq_along(d$time), function(t)
    do.call(rbind, lapply(runs, function(y) y[t, , drop = FALSE])))
}

run_sites <- function(d, rep, cfg, task) {
  suppressPackageStartupMessages(library(spDBL))
  S <- d$Nx * d$Ny
  T <- length(d$time)
  set.seed(task$seed + 100000L)
  held <- sort(sample.int(S, cfg$site_holdout))
  observed <- setdiff(seq_len(S), held)
  training <- if (is.null(cfg$train_n)) d$training else {
    X <- sir_lhs(as.integer(cfg$train_n), d$lower, d$upper, cfg$train_seed)
    list(input = X, responses = simulate_runs(X, d, cfg$box$n_people))
  }
  em <- timed(emulator_learn(training$input, training$responses, d$Nx, d$Ny,
    covariance = "sigma2R", episode_window = cfg$episode_window,
    nsam = cfg$n_draws, seed = cfg$seed))
  prep <- timed(calibration_prepare(em$value, d$lower, d$upper))
  fit <- timed(calibration_fit(prep$value, rep$z[, observed, drop = FALSE], d$U,
    locations = observed, chains = cfg$chains, iter = cfg$iter,
    warmup = cfg$warmup, draws_per_emulator = cfg$draws_per_emulator,
    n0 = d$prior$n0, d0 = d$prior$d0, seed = task$seed, cores = cfg$cores,
    refresh = 0L, init = cfg$init, keep_discrepancy = TRUE,
    cache_dir = cfg$stan_cache))
  f <- fit$value
  base <- timed(calibration_predict(f, seed = task$seed))
  pred <- timed(calibration_predict(f, seed = task$seed, sites = held))
  stopifnot(identical(pred$value[c("y", "u", "z_rep")], base$value))
  p <- pred$value$sites
  n <- nrow(f$mcmc_eta)

  # Calibrated emulator without the discrepancy, at the same posterior draws.
  y_sim <- array(NA_real_, c(n, T, length(held)))
  set.seed(task$seed + 200000L)
  for (i in seq_len(n)) {
    g <- spDBL:::predict_calibration(f$mcmc_eta[i, ], f$emulator_draw[i],
      f$emulator, held, covariance = "marginal")
    y_sim[i, , ] <- g$mu + sqrt(g$variance) * stats::rnorm(length(g$mu))
  }
  # Noise sd tau_t of each draw, recycled over the sites (element [i, t, k]).
  z_sim <- y_sim + array(stats::rnorm(length(y_sim)), dim(y_sim)) *
    array(sqrt(f$mcmc_tau2), dim(y_sim))

  # Inverse-distance weighting of the observed field at each time.
  xy <- as.matrix(d$grid[, c("row", "col")])
  dist2 <- outer(rowSums(xy[held, , drop = FALSE]^2), rowSums(xy[observed, , drop = FALSE]^2), "+") -
    2 * xy[held, , drop = FALSE] %*% t(xy[observed, , drop = FALSE])
  w <- 1 / dist2
  w <- w / rowSums(w)
  idw <- rep$z[, observed, drop = FALSE] %*% t(w)

  truth_z <- rep$z[, held, drop = FALSE]
  truth_f <- rep$y[, held, drop = FALSE] + rep$u[, held, drop = FALSE]
  scores <- rbind(
    cbind(target = "z", predictor = "spdbl", score_draws(p$z, truth_z)),
    cbind(target = "z", predictor = "simulator_only", score_draws(z_sim, truth_z)),
    cbind(target = "z", predictor = "idw", score_draws(idw, truth_z)),
    cbind(target = "y+u", predictor = "spdbl", score_draws(p$y + p$u, truth_f)),
    cbind(target = "y+u", predictor = "simulator_only", score_draws(y_sim, truth_f)),
    cbind(target = "y+u", predictor = "idw", score_draws(idw, truth_f)),
    cbind(target = "u", predictor = "spdbl",
          score_draws(p$u, rep$u[, held, drop = FALSE])))
  scores$held_sites <- paste(held, collapse = ";")
  scores$n_posterior <- n
  scores$fit_s <- em$elapsed + prep$elapsed + fit$elapsed
  scores$predict_s <- base$elapsed
  scores$predict_sites_s <- pred$elapsed
  scores$settings <- sprintf("train=%s M=%d chains=%d iter=%d warmup=%d per_draw=%d init=%s holdout=%d",
    if (is.null(cfg$train_n)) "bundled40" else paste0("lhs", cfg$train_n),
    cfg$n_draws, cfg$chains, cfg$iter, cfg$warmup, cfg$draws_per_emulator,
    paste(cfg$init, collapse = ","), cfg$site_holdout)
  # Pointwise predictive summaries for plots (held-out sites by time).
  band <- function(x) apply(x, c(2, 3), stats::quantile, c(0.025, 0.5, 0.975), names = FALSE)
  detail <- list(held = held, time = d$time, truth_z = truth_z, truth_f = truth_f,
    z = band(p$z), f = band(p$y + p$u), z_sim = band(z_sim), idw = idw)
  list(row = scores, detail = detail)
}

bench_worker_main(function(task) {
  cfg <- task$cfg
  bench_setup_lib()
  d <- readRDS(task$data_file)
  rep <- d$replicates[[task$rep]]
  out <- run_sites(d, rep, cfg, task)
  out$row <- cbind(data.frame(rep = task$rep, method = task$method), out$row)
  out
})
