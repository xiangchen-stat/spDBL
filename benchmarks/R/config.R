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
      timeout_sec = 6 * 3600, jobs = 1, threads = 1),
    fixed = list(
      # Fixed physical domain: the 10^2, 20^2 and 50^2 grids all discretise
      # the 10 x 10 domain of the 10^2 grid (sir_pde_domain()), so they
      # simulate the same epidemic. 10 independent training/test designs per
      # grid; design d uses seed + 1000 (d - 1), so design 1 at 10^2 equals
      # the paper profile's data. designgp_only and ppgasp_per_time are the
      # ablations; the [0, 1] input-rescaling pass of spDBL is dropped (it
      # changed RMSE by under 0.001 in the paper profile). Fixed before any
      # fixed-domain result was seen.
      label = "fixed domain, 10 designs", grids = c(10, 20, 50), episode_window = 2,
      domain = 10, designs = 10,
      methods = c("spdbl_iw_known", "spdbl_s2r_known", "spdbl_s2r_designgp",
                  "designgp_only", "ppgasp_default", "ppgasp_tuned",
                  "ppgasp_per_time", "mlegp_pca"),
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
      # init = "random": every chain of both packages starts at a uniform
      # draw in the box. Some replicates have separated posterior modes;
      # dispersed starts let R-hat detect them, while chains that all start
      # at the centre tend to stay in the same mode.
      n_draws = 20, chains = 4, iter = 2000, warmup = 1000, cores = 1,
      init = "random",
      rc_S = 10000, rc_S0 = 2000, rc_chains = 4, field_n0 = 30, field_d0 = 2,
      inner = 0.8, seed = 20261101, timeout_sec = 24 * 3600, jobs = 1,
      threads = 1),
    mode = list(
      # The paper profile with every chain started at its package's posterior
      # mode. With random starts, chains settled in local modes 35 to 586
      # log-density units below the main mode (results/e3_diagnostic), so the
      # pooled draws were unreliable. spDBL: calibration_fit(init = "mode"),
      # L-BFGS from 20 random starts in each conditional fit. RobustCalibration:
      # its maximum-likelihood estimate (method = "mle", rc_mle_starts initial
      # values), each chain started there with a small jitter, and twice the
      # iterations of the paper profile for its random-walk Metropolis chains.
      label = "paper, mode starts", n_rep = 30,
      methods = c("spdbl", "rc_sgasp", "rc_gasp"), episode_window = 2,
      n_draws = 20, chains = 4, iter = 2000, warmup = 1000, cores = 1,
      init = "mode", rc_mle_starts = 10,
      rc_S = 20000, rc_S0 = 5000, rc_chains = 4, field_n0 = 30, field_d0 = 2,
      inner = 0.8, seed = 20261101, timeout_sec = 24 * 3600, jobs = 1,
      threads = 1)))

# E3 ablations of the under-coverage (results/e3/abl_*): the mode profile with
# one change each, spDBL only. Run them on replicates 1-10 (--reps=1,...,10)
# of the same 30-replicate data, so they pair with the mode run's replicates
# and task seeds.
#   abl_joint  emulator coefficient paths sampled jointly (sampling = "joint")
#              rather than per-episode marginals
#   abl_nodu   design-GP uncertainty not propagated (design_uncertainty = FALSE)
#   abl_m60    60 emulator draws instead of 20 (4 cores per task)
#   abl_n40, abl_n80, abl_n160, abl_two  training runs (see below)
#   abl_self   field data with y drawn from the emulator's own model at the
#              true parameters (one random emulator draw, Gaussian across
#              sites with the design-GP term, independent across times), the
#              same u and noise: spDBL is then exactly specified, so this
#              checks the inference rather than the emulator
local({
  base <- bench_profiles$e3$mode
  base$methods <- "spdbl"
  abl <- function(label, ...) utils::modifyList(base, c(list(label = label), list(...)))
  bench_profiles$e3$abl_joint <<- abl("ablation: joint coefficient paths", sampling = "joint")
  bench_profiles$e3$abl_nodu <<- abl("ablation: no design-GP uncertainty",
                                      design_uncertainty = FALSE)
  bench_profiles$e3$abl_m60 <<- abl("ablation: 60 emulator draws", n_draws = 60, cores = 4)
  bench_profiles$e3$abl_self <<- abl("ablation: field data from the emulator's model",
                                      field_source = "emulator")
  # Training-size ablations: new Latin-hypercube designs over the whole box
  # (the same design for every replicate), and a two-stage design that adds
  # 40 runs where a short first-stage calibration puts the posterior
  # (widened threefold), chosen by greedy maximin.
  lhs_seed <- base$seed + 3000L
  bench_profiles$e3$abl_n40 <<- abl("ablation: 40 new Latin-hypercube training runs",
                                     train_n = 40, train_seed = lhs_seed)
  bench_profiles$e3$abl_n80 <<- abl("ablation: 80 Latin-hypercube training runs",
                                     train_n = 80, train_seed = lhs_seed, cores = 4)
  bench_profiles$e3$abl_n160 <<- abl("ablation: 160 Latin-hypercube training runs",
                                      train_n = 160, train_seed = lhs_seed, cores = 4)
  bench_profiles$e3$abl_two <<- abl("ablation: two-stage design, 40 bundled + 40 targeted runs",
                                     train_add = 40, stage1_draws = 10, stage1_chains = 2,
                                     stage1_iter = 1000, cores = 4)
})

# Timing reruns after the 1.2.0 speed-ups (batched simulator-free prediction,
# keep_draws = FALSE): E1 with regenerated draws, and E2 design 1 with every
# method, plus spDBL with regenerated draws, on the same node.
local({
  bench_profiles$e1$paper_lazy <<- utils::modifyList(bench_profiles$e1$paper,
    list(label = "paper, keep_draws = FALSE", keep_draws = FALSE))
  bench_profiles$e2$fixed_timing <<- utils::modifyList(bench_profiles$e2$fixed,
    list(label = "fixed domain, design 1, timing after the 1.2.0 speed-ups",
         methods = c(bench_profiles$e2$fixed$methods, "spdbl_s2r_designgp_lazy")))
})

# E3 with a larger training design for every method: 100 Latin-hypercube runs
# (the same design for all methods and replicates) instead of the bundled 40,
# all 30 replicates, otherwise the mode profile. RobustCalibration's four
# chains run at the same time (rc_cores = 4) to stay within a 24-hour job.
local({
  bench_profiles$e3$n100 <<- utils::modifyList(bench_profiles$e3$mode,
    list(label = "100 Latin-hypercube training runs, all methods",
         train_n = 100, train_seed = bench_profiles$e3$mode$seed + 3000L,
         rc_cores = 4))
})
