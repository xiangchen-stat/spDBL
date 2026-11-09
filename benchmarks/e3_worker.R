# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E3 worker: calibrate one replicate field dataset with one method.
# Called by e3_calibration.R as: Rscript e3_worker.R <task.rds>
#
# Methods
#   spdbl     emulator_learn(covariance = "sigma2R") -> calibration_prepare()
#             -> calibration_fit() (RStan NUTS; explicit chains/iter/warmup).
#   rc_sgasp  RobustCalibration::rcalibration, simul_type = 0 (ppgasp
#             emulator of the T x S output), S-GaSP discrepancy (its default).
#   rc_gasp   as rc_sgasp with the Kennedy-O'Hagan GaSP discrepancy.

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))

# Rank-normalised split R-hat and bulk/tail ESS (rstan implementations, as in
# Vehtari et al. 2021), applied identically to both packages' draws.
chain_diagnostics <- function(mat) {
  c(rhat = rstan::Rhat(mat), ess_bulk = rstan::ess_bulk(mat),
    ess_tail = rstan::ess_tail(mat))
}

summarise_parameter <- function(draws, truth, lower, upper) {
  q <- stats::quantile(draws, c(0.025, 0.5, 0.975), names = FALSE)
  range <- upper - lower
  data.frame(truth = truth, post_mean = mean(draws), post_median = q[2],
    q025 = q[1], q975 = q[3], covered = truth >= q[1] && truth <= q[3],
    width = q[3] - q[1], rel_width = (q[3] - q[1]) / range,
    error = mean(draws) - truth, rel_error = (mean(draws) - truth) / range)
}

# Field data for the abl_self ablation: the replicate's discrepancy and noise
# with the simulator output y replaced by a draw from the emulator's own model
# at the true parameters (one emulator draw chosen at random; at each time
# Gaussian with the conditional fit's mean and covariance, independently
# across times). The fixed initial field is kept.
self_consistent_field <- function(rep, cal, seed) {
  set.seed(seed + 500000L)
  m <- sample.int(cal$M, 1L)
  pred <- spDBL:::predict_calibration(rep$eta, m, cal)
  n_t <- nrow(pred$mu)
  times <- seq(nrow(rep$y) - n_t + 1L, nrow(rep$y))
  y <- rep$y
  for (k in seq_len(n_t)) {
    if (all(abs(rep$y[times[k], ] - pred$mu[k, ]) < 1e-10)) next
    C <- spDBL:::.cal_time_covariance(pred, k)
    y[times[k], ] <- pred$mu[k, ] + drop(crossprod(chol(C), stats::rnorm(ncol(C))))
  }
  rep$z - rep$y + y
}

# Training runs for the training-size ablations: the bundled 40 runs by
# default; cfg$train_n new Latin-hypercube runs over the whole box
# (abl_n*); or, with cfg$train_add, the bundled runs plus that many runs
# placed in a widened first-stage posterior (abl_two).
simulate_runs <- function(X, d, n_people) {
  runs <- lapply(seq_len(nrow(X)), function(i)
    spDBL::gen_pde(X[i, ], d$Nx, d$Ny, n_people, length(d$time)))
  lapply(seq_along(d$time), function(t)
    do.call(rbind, lapply(runs, function(y) y[t, , drop = FALSE])))
}

spdbl_training <- function(d, rep, cfg, task) {
  if (!is.null(cfg$train_n)) {
    X <- sir_lhs(as.integer(cfg$train_n), d$lower, d$upper, cfg$train_seed)
    return(list(input = X, responses = simulate_runs(X, d, cfg$box$n_people)))
  }
  if (is.null(cfg$train_add)) return(c(d$training, list(stage1 = NULL)))
  # Stage 1: a short calibration on the bundled runs.
  active <- which(d$upper > d$lower)
  range <- (d$upper - d$lower)[active]
  t1 <- proc.time()
  em1 <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
    covariance = "sigma2R", episode_window = cfg$episode_window,
    nsam = cfg$stage1_draws, seed = cfg$seed)
  cal1 <- calibration_prepare(em1, d$lower, d$upper)
  fit1 <- calibration_fit(cal1, rep$z, d$U, chains = cfg$stage1_chains,
    iter = cfg$stage1_iter, warmup = cfg$stage1_iter %/% 2L, n0 = d$prior$n0,
    d0 = d$prior$d0, seed = task$seed + 9000L, cores = cfg$cores, refresh = 0L,
    init = cfg$init, keep_discrepancy = FALSE, cache_dir = cfg$stan_cache)
  stage1_cpu <- sum((proc.time() - t1)[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)
  # Candidates from the stage-1 posterior with its spread tripled, inside the
  # box; then greedy maximin on the unit scale, away from the bundled runs.
  draws <- sweep(sweep(fit1$mcmc_eta[, active, drop = FALSE], 2, d$lower[active]), 2, range, "/")
  centre <- colMeans(draws)
  L <- t(chol(stats::cov(draws) + diag(1e-10, length(active))))
  set.seed(task$seed + 700000L)
  cand <- t(centre + 3 * L %*% matrix(stats::rnorm(length(active) * 20000L), length(active)))
  cand <- cand[apply(cand >= 0 & cand <= 1, 1, all), , drop = FALSE]
  if (nrow(cand) < 4L * cfg$train_add) stop("Too few stage-2 candidates inside the box")
  old <- sweep(sweep(d$training$input[, active, drop = FALSE], 2, d$lower[active]), 2, range, "/")
  dmin <- apply(cand, 1, function(x) min(colSums((t(old) - x)^2)))
  chosen <- integer()
  for (k in seq_len(cfg$train_add)) {
    j <- which.max(dmin); chosen <- c(chosen, j)
    dmin <- pmin(dmin, colSums((t(cand) - cand[j, ])^2))
  }
  X <- matrix(rep(d$lower, each = cfg$train_add), cfg$train_add, length(d$lower),
              dimnames = list(NULL, names(d$lower)))
  X[, active] <- sweep(sweep(cand[chosen, , drop = FALSE], 2, range, "*"), 2, d$lower[active], "+")
  added <- simulate_runs(X, d, cfg$box$n_people)
  list(input = rbind(d$training$input, X),
       responses = Map(rbind, d$training$responses, added),
       stage1 = list(cpu_s = stage1_cpu, mean = d$lower[active] + centre * range,
                     sd = apply(fit1$mcmc_eta[, active, drop = FALSE], 2, stats::sd)))
}

run_spdbl <- function(d, rep, cfg, task) {
  suppressPackageStartupMessages(library(spDBL))
  active <- which(d$upper > d$lower)
  tr <- timed(spdbl_training(d, rep, cfg, task))
  training <- tr$value
  em <- timed(emulator_learn(training$input, training$responses, d$Nx, d$Ny,
    covariance = "sigma2R", episode_window = cfg$episode_window,
    nsam = cfg$n_draws, seed = cfg$seed))
  prep <- timed(calibration_prepare(em$value, d$lower, d$upper,
    sampling = cfg$sampling,
    design_uncertainty = if (is.null(cfg$design_uncertainty)) TRUE else
      isTRUE(cfg$design_uncertainty)))
  z <- if (identical(cfg$field_source, "emulator"))
    self_consistent_field(rep, prep$value, task$seed) else rep$z
  keep <- cfg$chains * (cfg$iter - cfg$warmup)
  fit <- timed(calibration_fit(prep$value, z, d$U,
    chains = cfg$chains, iter = cfg$iter, warmup = cfg$warmup,
    draws_per_emulator = keep, n0 = d$prior$n0, d0 = d$prior$d0,
    seed = task$seed, cores = cfg$cores, refresh = 0L, init = cfg$init,
    keep_discrepancy = FALSE, cache_dir = cfg$stan_cache))
  f <- fit$value
  out <- lapply(active, function(j) {
    name <- colnames(f$mcmc_eta)[j]
    diag <- t(sapply(unique(f$emulator_draw), function(m) {
      sel <- f$emulator_draw == m
      mat <- matrix(NA_real_, cfg$iter - cfg$warmup, cfg$chains)
      mat[cbind(f$iteration[sel], f$chain[sel])] <- f$mcmc_eta[sel, j]
      chain_diagnostics(mat)
    }))
    cbind(data.frame(parameter = name),
      summarise_parameter(f$mcmc_eta[, j], rep$eta[j], d$lower[j], d$upper[j]),
      data.frame(rhat_max = max(diag[, "rhat"]), ess_bulk_min = min(diag[, "ess_bulk"]),
        ess_tail_min = min(diag[, "ess_tail"]), ess_bulk_total = sum(diag[, "ess_bulk"])))
  })
  out <- do.call(rbind, out)
  pkg_diag <- do.call(rbind, lapply(f$diagnostics, function(x)
    data.frame(div = x$divergences, td = x$max_treedepth_hits,
               rhat = x$max_rhat, neff = x$min_n_eff)))
  out$n_samples <- nrow(f$mcmc_eta)
  out$divergences <- sum(pkg_diag$div)
  out$treedepth_hits <- sum(pkg_diag$td)
  out$pkg_max_rhat_all <- max(pkg_diag$rhat)
  # init = "mode": the fewest starts that reached the best optimum (within one
  # log-density unit) over the conditional fits.
  near <- vapply(f$diagnostics, function(x)
    if (is.null(x$mode_search)) NA_real_ else as.numeric(x$mode_search$near_best), numeric(1))
  out$mode_near_best_min <- if (all(is.na(near))) NA_real_ else min(near, na.rm = TRUE)
  out$pkg_min_neff_all <- min(pkg_diag$neff)
  out$emulator_s <- em$elapsed + prep$elapsed
  out$n_train <- nrow(training$input)
  out$training_s <- tr$elapsed
  out$stage1_cpu_s <- if (is.null(training$stage1)) NA_real_ else training$stage1$cpu_s
  out$sampling_s <- fit$elapsed
  out$sampling_cpu_s <- fit$cpu
  out$wall_total_s <- em$elapsed + prep$elapsed + fit$elapsed
  # ESS/s: the slowest parameter's ESS, summed over independent conditionals,
  # per wall second and per CPU second (rstan's chain processes included).
  out$ess_per_s <- min(out$ess_bulk_total) / fit$elapsed
  out$ess_per_cpu_s <- min(out$ess_bulk_total) / fit$cpu
  out$settings <- sprintf("M=%d chains=%d iter=%d warmup=%d cores=%d init=%s sampling=%s design_uncertainty=%s field=%s",
    cfg$n_draws, cfg$chains, cfg$iter, cfg$warmup, cfg$cores,
    paste(cfg$init, collapse = ","), prep$value$sampling_method,
    prep$value$design_uncertainty,
    if (is.null(cfg$field_source)) "simulator" else cfg$field_source)
  out
}

run_rc <- function(d, rep, cfg, task, discrepancy) {
  suppressPackageStartupMessages({library(RobustCalibration); library(RobustGaSP)})
  active <- which(d$upper > d$lower)
  S <- d$Nx * d$Ny
  T <- length(d$time)
  Y <- do.call(cbind, lapply(d$training$responses, as.matrix))  # column (t-1)S+s
  # ppgasp rejects constant columns: all runs share the t = 0 field, so the
  # t = 0 observations are dropped for RobustCalibration (spDBL uses them).
  keep <- apply(Y, 2, function(v) max(v) - min(v) > 1e-10)
  obs <- as.vector(t(rep$z))[keep]
  x <- cbind(time = rep(d$time, each = S), row = rep(d$grid$row, T),
             col = rep(d$grid$col, T))[keep, , drop = FALSE]
  X_sim <- d$training$input[, active, drop = FALSE]
  theta_range <- cbind(d$lower[active], d$upper[active])
  emu <- timed(utils::capture.output(ppgasp(design = X_sim, response = Y[, keep])))
  # init = "mode": RobustCalibration's own maximum-likelihood estimate from
  # rc_mle_starts initial values, the counterpart of spDBL's mode search; its
  # time counts towards the calibration time.
  mode_s <- 0
  theta_mode <- NULL
  if (identical(cfg$init, "mode")) {
    set.seed(task$seed)
    mle <- timed(utils::capture.output(fit_mle <- rcalibration(design = x,
      observations = obs, p_theta = length(active), simul_type = 0,
      input_simul = X_sim, output_simul = Y[, keep],
      loc_index_emulator = seq_len(sum(keep)), theta_range = theta_range,
      discrepancy_type = discrepancy, method = "mle",
      num_initial_values = cfg$rc_mle_starts)))
    theta_mode <- fit_mle@param_est[seq_along(active)]
    mode_s <- mle$elapsed
  }
  chains <- vector("list", cfg$rc_chains)
  times <- accept <- numeric(cfg$rc_chains)
  for (c in seq_len(cfg$rc_chains)) {
    set.seed(task$seed + 100L * c)
    # Same starting rules as spDBL: the default (centre of the box) for
    # init = "midpoint", a uniform draw in the box for init = "random", and
    # the mode jittered on the logit scale (sd 0.05) for init = "mode".
    start <- if (!is.null(theta_mode)) {
      width <- theta_range[, 2] - theta_range[, 1]
      u <- pmin(pmax((theta_mode - theta_range[, 1]) / width, 1e-6), 1 - 1e-6)
      theta_range[, 1] + width * stats::plogis(stats::qlogis(u) +
        stats::rnorm(length(active), sd = 0.05))
    } else if (identical(cfg$init, "random"))
      stats::runif(length(active), theta_range[, 1], theta_range[, 2]) else NULL
    run <- timed(rcalibration(design = x, observations = obs,
      p_theta = length(active), simul_type = 0, input_simul = X_sim,
      output_simul = Y[, keep], loc_index_emulator = seq_len(sum(keep)),
      theta_range = theta_range, S = cfg$rc_S, S_0 = cfg$rc_S0,
      discrepancy_type = discrepancy, initial_values = start))
    chains[[c]] <- run$value@post_sample[, seq_along(active), drop = FALSE]
    times[c] <- run$elapsed
    accept[c] <- run$value@accept_S[1] / cfg$rc_S
  }
  n_keep <- min(vapply(chains, nrow, integer(1)))
  out <- lapply(seq_along(active), function(k) {
    j <- active[k]
    mat <- sapply(chains, function(ch) ch[seq_len(n_keep), k])
    dg <- chain_diagnostics(mat)
    cbind(data.frame(parameter = names(d$lower)[j]),
      summarise_parameter(as.vector(mat), rep$eta[j], d$lower[j], d$upper[j]),
      data.frame(rhat_max = dg[["rhat"]], ess_bulk_min = dg[["ess_bulk"]],
        ess_tail_min = dg[["ess_tail"]], ess_bulk_total = dg[["ess_bulk"]]))
  })
  out <- do.call(rbind, out)
  out$n_samples <- n_keep * cfg$rc_chains
  out$accept_rate_theta <- mean(accept)
  out$n_obs_used <- length(obs)
  out$emulator_s <- emu$elapsed
  out$mode_search_s <- mode_s
  out$sampling_s <- sum(times) + mode_s
  out$sampling_cpu_s <- sum(times) + mode_s
  out$wall_total_s <- sum(times) + mode_s
  out$ess_per_s <- min(out$ess_bulk_total) / (sum(times) + mode_s)
  out$ess_per_cpu_s <- out$ess_per_s  # chains run one after another, one core
  out$settings <- sprintf("%s S=%d S_0=%d chains=%d init=%s", discrepancy,
    cfg$rc_S, cfg$rc_S0, cfg$rc_chains, paste(cfg$init, collapse = ","))
  out
}

bench_worker_main(function(task) {
  cfg <- task$cfg
  bench_setup_lib()
  d <- readRDS(task$data_file)
  rep <- d$replicates[[task$rep]]
  out <- switch(task$method,
    spdbl = run_spdbl(d, rep, cfg, task),
    rc_sgasp = run_rc(d, rep, cfg, task, "S-GaSP"),
    rc_gasp = run_rc(d, rep, cfg, task, "GaSP"),
    stop("Unknown method ", task$method))
  list(row = cbind(data.frame(rep = task$rep, method = task$method), out))
})
