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

run_spdbl <- function(d, rep, cfg, task) {
  suppressPackageStartupMessages(library(spDBL))
  active <- which(d$upper > d$lower)
  em <- timed(emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
    covariance = "sigma2R", episode_window = cfg$episode_window,
    nsam = cfg$n_draws, seed = cfg$seed))
  prep <- timed(calibration_prepare(em$value, d$lower, d$upper))
  keep <- cfg$chains * (cfg$iter - cfg$warmup)
  fit <- timed(calibration_fit(prep$value, rep$z, d$U,
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
  out$pkg_min_neff_all <- min(pkg_diag$neff)
  out$emulator_s <- em$elapsed + prep$elapsed
  out$sampling_s <- fit$elapsed
  out$sampling_cpu_s <- fit$cpu
  out$wall_total_s <- em$elapsed + prep$elapsed + fit$elapsed
  # ESS/s: the slowest parameter's ESS, summed over independent conditionals,
  # per wall second and per CPU second (rstan's chain processes included).
  out$ess_per_s <- min(out$ess_bulk_total) / fit$elapsed
  out$ess_per_cpu_s <- min(out$ess_bulk_total) / fit$cpu
  out$settings <- sprintf("M=%d chains=%d iter=%d warmup=%d cores=%d init=%s",
    cfg$n_draws, cfg$chains, cfg$iter, cfg$warmup, cfg$cores,
    paste(cfg$init, collapse = ","))
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
  chains <- vector("list", cfg$rc_chains)
  times <- accept <- numeric(cfg$rc_chains)
  for (c in seq_len(cfg$rc_chains)) {
    set.seed(task$seed + 100L * c)
    run <- timed(rcalibration(design = x, observations = obs,
      p_theta = length(active), simul_type = 0, input_simul = X_sim,
      output_simul = Y[, keep], loc_index_emulator = seq_len(sum(keep)),
      theta_range = theta_range, S = cfg$rc_S, S_0 = cfg$rc_S0,
      discrepancy_type = discrepancy))
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
  out$sampling_s <- sum(times)
  out$sampling_cpu_s <- sum(times)
  out$wall_total_s <- sum(times)
  out$ess_per_s <- min(out$ess_bulk_total) / sum(times)
  out$ess_per_cpu_s <- out$ess_per_s  # chains run one after another, one core
  out$settings <- sprintf("%s S=%d S_0=%d chains=%d", discrepancy, cfg$rc_S,
                          cfg$rc_S0, cfg$rc_chains)
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
