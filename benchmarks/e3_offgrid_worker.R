# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Off-grid prediction worker: one replicate of sir_offgrid_replicates().
# Called by e3_offgrid.R as:  Rscript e3_offgrid_worker.R <task.rds>
#
# Calibrates spDBL on the field data at the 36 coarse cell centres, predicts
# the target points between them with calibration_predict(coords = ...), and
# scores the predictions of the held-out observations z and of the true field
# y + u (fine-grid simulator at the true parameters plus the discrepancy):
#   spdbl_cell, spdbl_bilinear  calibration_predict(field = "cell"/"bilinear")
#   sim_cell, sim_bilinear      the calibrated emulator at the same field map,
#                               without the discrepancy (z adds N(0, tau2_t))
#   idw                         inverse-distance-weighted (power 2)
#                               interpolation of the observed field

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

# Weight matrix (points x cells) of a field map from .cal_point_cells().
field_weights <- function(map, cells) {
  W <- matrix(0, nrow(map$cells), length(cells))
  for (m in seq_len(ncol(map$cells))) {
    use <- which(!is.na(map$cells[, m]))
    W[cbind(use, match(map$cells[use, m], cells))] <- map$weights[use, m]
  }
  W
}

run_offgrid <- function(d, rep, cfg, task) {
  suppressPackageStartupMessages(library(spDBL))
  S <- d$Nx * d$Ny
  T <- length(d$time)
  grid <- unname(as.matrix(d$grid))
  coords <- rep$target_coords
  K <- nrow(coords)
  kernel <- sir_offgrid_kernel(d$kernel_phi)
  training <- if (is.null(cfg$train_n)) d$training else {
    X <- sir_lhs(as.integer(cfg$train_n), d$lower, d$upper, cfg$train_seed)
    list(input = X, responses = simulate_runs(X, d, cfg$box$n_people))
  }
  em <- timed(emulator_learn(training$input, training$responses, d$Nx, d$Ny,
    covariance = "sigma2R", episode_window = cfg$episode_window,
    nsam = cfg$n_draws, seed = cfg$seed))
  prep <- timed(calibration_prepare(em$value, d$lower, d$upper))
  fit <- timed(calibration_fit(prep$value, rep$z, d$U, locations = seq_len(S),
    chains = cfg$chains, iter = cfg$iter, warmup = cfg$warmup,
    draws_per_emulator = cfg$draws_per_emulator, n0 = d$prior$n0,
    d0 = d$prior$d0, seed = task$seed, cores = cfg$cores, refresh = 0L,
    init = cfg$init, keep_discrepancy = TRUE, cache_dir = cfg$stan_cache))
  f <- fit$value
  n <- nrow(f$mcmc_eta)
  base <- timed(calibration_predict(f, seed = task$seed))
  pc <- timed(calibration_predict(f, seed = task$seed, coords = coords,
    discrepancy_kernel = kernel, field = "cell"))
  pb <- timed(calibration_predict(f, seed = task$seed, coords = coords,
    discrepancy_kernel = kernel, field = "bilinear"))
  stopifnot(identical(pc$value[c("y", "u", "z_rep")], base$value),
            identical(pb$value[c("y", "u", "z_rep")], base$value))

  # Calibrated emulator without the discrepancy, through the same field maps.
  maps <- list(cell = spDBL:::.cal_point_cells(grid, coords, "cell"),
               bilinear = spDBL:::.cal_point_cells(grid, coords, "bilinear"))
  cells <- sort(unique(stats::na.omit(c(maps$cell$cells, maps$bilinear$cells))))
  W <- lapply(maps, field_weights, cells = cells)
  y_sim <- lapply(W, function(w) array(NA_real_, c(n, T, K)))
  set.seed(task$seed + 200000L)
  for (i in seq_len(n)) {
    g <- spDBL:::predict_calibration(f$mcmc_eta[i, ], f$emulator_draw[i],
      f$emulator, cells, covariance = "full")
    for (t in seq_len(T)) {
      S_t <- spDBL:::.cal_time_covariance(g, t)
      for (m in names(W)) {
        mean <- as.numeric(W[[m]] %*% g$mu[t, ])
        sd <- sqrt(pmax(rowSums((W[[m]] %*% S_t) * W[[m]]), 0))
        y_sim[[m]][i, t, ] <- mean + sd * stats::rnorm(K)
      }
    }
  }
  noise_sd <- array(sqrt(f$mcmc_tau2), c(n, T, K))
  z_sim <- lapply(y_sim, function(y) y + array(stats::rnorm(length(y)), dim(y)) * noise_sd)

  # Inverse-distance weighting of the observed field at each time.
  dist2 <- outer(rowSums(coords^2), rowSums(grid^2), "+") - 2 * coords %*% t(grid)
  w <- 1 / dist2
  w <- w / rowSums(w)
  idw <- rep$z %*% t(w)

  truth_z <- rep$z_target
  truth_f <- rep$y_target + rep$u_target
  p <- list(spdbl_cell = pc$value$points, spdbl_bilinear = pb$value$points)
  scores <- rbind(
    cbind(target = "z", predictor = "spdbl_cell", score_draws(p$spdbl_cell$z, truth_z)),
    cbind(target = "z", predictor = "spdbl_bilinear", score_draws(p$spdbl_bilinear$z, truth_z)),
    cbind(target = "z", predictor = "sim_cell", score_draws(z_sim$cell, truth_z)),
    cbind(target = "z", predictor = "sim_bilinear", score_draws(z_sim$bilinear, truth_z)),
    cbind(target = "z", predictor = "idw", score_draws(idw, truth_z)),
    cbind(target = "y+u", predictor = "spdbl_cell",
          score_draws(p$spdbl_cell$y + p$spdbl_cell$u, truth_f)),
    cbind(target = "y+u", predictor = "spdbl_bilinear",
          score_draws(p$spdbl_bilinear$y + p$spdbl_bilinear$u, truth_f)),
    cbind(target = "y+u", predictor = "sim_cell", score_draws(y_sim$cell, truth_f)),
    cbind(target = "y+u", predictor = "sim_bilinear", score_draws(y_sim$bilinear, truth_f)),
    cbind(target = "y+u", predictor = "idw", score_draws(idw, truth_f)),
    cbind(target = "u", predictor = "spdbl_cell", score_draws(p$spdbl_cell$u, rep$u_target)),
    cbind(target = "u", predictor = "spdbl_bilinear", score_draws(p$spdbl_bilinear$u, rep$u_target)))
  scores$n_points <- K
  scores$n_posterior <- n
  scores$fit_s <- em$elapsed + prep$elapsed + fit$elapsed
  scores$predict_s <- base$elapsed
  scores$predict_cell_s <- pc$elapsed
  scores$predict_bilinear_s <- pb$elapsed
  scores$settings <- sprintf("train=%s M=%d chains=%d iter=%d warmup=%d per_draw=%d init=%s refine=%d points=%d",
    if (is.null(cfg$train_n)) "bundled40" else paste0("lhs", cfg$train_n),
    cfg$n_draws, cfg$chains, cfg$iter, cfg$warmup, cfg$draws_per_emulator,
    paste(cfg$init, collapse = ","), d$refine, K)
  band <- function(x) apply(x, c(2, 3), stats::quantile, c(0.025, 0.5, 0.975), names = FALSE)
  detail <- list(coords = coords, time = d$time, truth_z = truth_z, truth_f = truth_f,
    cell_f = band(p$spdbl_cell$y + p$spdbl_cell$u),
    bilinear_f = band(p$spdbl_bilinear$y + p$spdbl_bilinear$u),
    sim_bilinear = band(z_sim$bilinear), idw = idw)
  list(row = scores, detail = detail)
}

bench_worker_main(function(task) {
  cfg <- task$cfg
  bench_setup_lib()
  d <- readRDS(task$data_file)
  rep <- d$replicates[[task$rep]]
  out <- run_offgrid(d, rep, cfg, task)
  out$row <- cbind(data.frame(rep = task$rep, method = task$method), out$row)
  out
})
