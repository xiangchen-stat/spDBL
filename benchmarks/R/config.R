# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Benchmark profiles. "validation" is the small configuration used to check
# the harness on a laptop; "paper" is the configuration intended for the
# manuscript on dedicated hardware. Any field can be overridden on the
# command line, e.g. --grids=10,20 --jobs=2. Unknown keys are rejected.
# Every grid size must be divisible by every episode window (1, 2 and 5).

# Parameter box shared by all SIR experiments (taken from dt_calibration):
# eta1 (transmission) in [2, 4], eta2 (recovery) in [0.2, 0.4],
# alpha2 (diffusion of I) in [0, 0.2]; alpha1 = alpha3 = 0 are fixed.
bench_sir_box <- list(
  lower = c(eta1 = 2, eta2 = 0.2, alpha1 = 0, alpha2 = 0, alpha3 = 0),
  upper = c(eta1 = 4, eta2 = 0.4, alpha1 = 0, alpha2 = 0.2, alpha3 = 0),
  n_people = 10000)

bench_profiles <- list(
  e1 = list(
    validation = list(
      label = "validation only", grids = c(10, 20, 30, 40, 50),
      windows = c(1, 2, 5), covariance = "iw", n_train = 40, n_test = 10,
      n_time = 26, nsam = 10, ar_choice = 2, eval_from = 3,
      predict_mc = FALSE, measure_rfield = TRUE, seed = 20261101,
      timeout_sec = 1800, jobs = 1, threads = 1),
    paper = list(
      # 100^2 with window 5 may exceed 64 GB (peak RSS about 68 GB projected
      # from 50^2); it runs last, and its peak memory is recorded either way.
      label = "paper", grids = c(10, 20, 30, 40, 50, 70, 100),
      windows = c(1, 2, 5), covariance = "iw", n_train = 40, n_test = 10,
      n_time = 26, nsam = 10, ar_choice = 2, eval_from = 3, predict_mc = TRUE, measure_rfield = TRUE,
      seed = 20261101, timeout_sec = 6 * 3600, jobs = 1, threads = 1)),
  e2 = list(
    validation = list(
      label = "validation only", grids = c(10, 20), episode_window = 2,
      methods = c("spdbl_iw_known", "spdbl_s2r_known", "spdbl_s2r_designgp",
                  "spdbl_s2r_designgp_unit", "ppgasp_default", "ppgasp_tuned",
                  "mlegp_pca"),
      n_train = 40, n_test = 10, n_time = 26, nsam = 20, ar_choice = 2,
      eval_from = 3, mlegp_pc_percent = 99.9, seed = 20261101,
      timeout_sec = 1800, jobs = 2, threads = 1),
    paper = list(
      label = "paper", grids = c(10, 20, 50), episode_window = 2,
      methods = c("spdbl_iw_known", "spdbl_s2r_known", "spdbl_s2r_designgp",
                  "spdbl_s2r_designgp_unit", "ppgasp_default", "ppgasp_tuned",
                  "mlegp_pca"),
      # nsam = 20: spDBL stores every coefficient draw (nsam x episodes x
      # 2C x C doubles). Peak RSS grows about as grid^3 at a fixed window:
      # about 15 GB at 50^2 but about 120 GB at 100^2, so E2 stops at 50^2.
      # 50 held-out runs: 10 sampled the input box too thinly for coverage.
      n_train = 40, n_test = 50, n_time = 26, nsam = 20, ar_choice = 2,
      eval_from = 3, mlegp_pc_percent = 99.9, seed = 20261101,
      timeout_sec = 6 * 3600, jobs = 1, threads = 1)),
  e3 = list(
    validation = list(
      label = "validation only", n_rep = 3,
      methods = c("spdbl", "rc_sgasp", "rc_gasp"), episode_window = 2,
      n_draws = 4, chains = 2, iter = 500, warmup = 250, cores = 1,
      init = "midpoint",
      rc_S = 2000, rc_S0 = 500, rc_chains = 2, field_n0 = 30, field_d0 = 2,
      inner = 0.8, seed = 20261101, timeout_sec = 3600, jobs = 3, threads = 1),
    paper = list(
      # Every task is single-core (cores = 1; RobustCalibration has no
      # parallel chains), so wall and CPU time agree and both packages are
      # timed alike. Throughput comes from running tasks side by side: run
      # with --jobs=14 on a 16-core machine (about 13 h for 30 replicates).
      label = "paper", n_rep = 30,
      methods = c("spdbl", "rc_sgasp", "rc_gasp"), episode_window = 2,
      n_draws = 20, chains = 4, iter = 2000, warmup = 1000, cores = 1,
      init = "midpoint",
      rc_S = 10000, rc_S0 = 2000, rc_chains = 4, field_n0 = 30, field_d0 = 2,
      inner = 0.8, seed = 20261101, timeout_sec = 24 * 3600, jobs = 1,
      threads = 1)))
