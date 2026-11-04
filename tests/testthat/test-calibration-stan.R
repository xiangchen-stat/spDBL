# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

test_that("installed Stan programs match R and pool conditionals equally", {
  skip_if(Sys.getenv("SPDBL_RUN_STAN_TESTS") != "true",
          "Set SPDBL_RUN_STAN_TESTS=true for compiler-dependent integration checks")
  skip_if_not_installed("rstan", minimum_version = "2.32.7")
  f <- calibration_fixture()
  f$U <- 1.1001 * f$U
  e <- f$emulator
  # The default design GP of calibration_prepare(): Matern 5/2, one kernel
  # per time step (here two groups).
  per_time <- calibration_fixture(design_kernel = "matern52", design_groups = "time")
  per_time$U <- 1.1001 * per_time$U
  for (marginal in c(TRUE, FALSE)) {
    model <- spDBL:::.cal_stan_model(marginal, cache = TRUE,
      cache_dir = Sys.getenv("SPDBL_STAN_CACHE_DIR", tempdir()))
    for (fx in list(f, per_time)) {
      em <- fx$emulator
      initial <- list(eta_unit = fx$eta, precision = c(5, 6, 7))
      if (!marginal) initial$u_std <- matrix(seq(-0.2, 0.2, length.out = 12), 3)
      data <- spDBL:::make_stan_calibration_data(em, 1L, fx$z, fx$U)
      fixed <- rstan::sampling(model, data = data, chains = 1L, iter = 1L,
        warmup = 0L, algorithm = "Fixed_param", init = list(initial), refresh = 0L)
      values <- rstan::extract(fixed)
      prediction <- spDBL:::predict_calibration(fx$eta, 1L, em)
      expect_equal(matrix(values$mu, em$T, em$S), prediction$mu, tolerance = 1e-9)
      expect_equal(matrix(values$Sigma, em$S, em$S), unname(prediction$Sigma),
                   tolerance = 1e-9)
      expect_equal(as.numeric(values$design_variance), prediction$design_variance,
                   tolerance = 1e-9)
      for (t in seq_len(em$T))
        expect_equal(as.numeric(values$design_variance)[t] * data$design_covariance[t, , ],
          design_cov(prediction, t), tolerance = 1e-9, ignore_attr = TRUE)
      alternative <- initial
      alternative$eta_unit <- c(0.7, 0.35)
      alternative$precision <- c(4, 8, 6)
      stan_density <- function(v) rstan::log_prob(fixed,
        rstan::unconstrain_pars(fixed, v), adjust_transform = FALSE)
      r_density <- if (marginal) function(v)
        spDBL:::calibration_marginal_log_density(v$eta_unit, v$precision, em, 1L,
          fx$z, fx$U) else
        function(v) spDBL:::calibration_log_density(v$eta_unit, v$precision,
          v$u_std, em, 1L, fx$z, fx$U)
      expect_equal(stan_density(alternative) - stan_density(initial),
                   r_density(alternative) - r_density(initial), tolerance = 1e-8)
    }
    # The two time steps have their own design-GP variance factors.
    expect_false(isTRUE(all.equal(prediction$design_variance[1],
                                  prediction$design_variance[2])))
    fit <- suppressWarnings(calibration_fit(e, f$z, f$U,
      chains = 1L, iter = 40L, draws_per_emulator = 4L,
      compiled_model = model, marginal_discrepancy = marginal,
      refresh = 0L, keep_discrepancy = FALSE))
    expect_s3_class(fit, "spDBL_calibration")
    expect_equal(as.integer(table(fit$emulator_draw)), c(4L, 4L))
    expect_null(fit$mcmc_u)
    path <- tempfile(fileext = ".rds")
    saveRDS(fit, path)
    roundtrip <- readRDS(path)
    predictions <- calibration_predict(roundtrip, seed = 12L)
    expect_identical(dim(predictions$y), c(8L, 2L, 4L))
    expect_true(all(is.finite(unlist(predictions))))
    expect_equal(summary(roundtrip)$n_samples, 8L)
    expect_equal(fit$settings$warmup, 20L)
    expect_identical(fit$settings$control$adapt_delta, 0.95)
    if (marginal) {
      # This zero-warmup transition checks bookkeeping, not convergence.
      one <- suppressWarnings(calibration_fit(e, f$z, f$U, iter = 1L, compiled_model = model,
        keep_discrepancy = FALSE, refresh = 0L))
      expect_true(one$execution_only)
      # rstan calls the start function once per chain with its chain_id.
      random <- suppressWarnings(calibration_fit(e, f$z, f$U, compiled_model = model,
        draws = 1L, chains = 2L, iter = 2L, warmup = 0L, draws_per_emulator = 2L,
        init = "random", keep_discrepancy = FALSE, refresh = 0L))
      expect_identical(random$settings$init, "random")
      expect_true(all(is.finite(random$mcmc_eta)))
      expect_false(isTRUE(all.equal(random$mcmc_eta[1, ], random$mcmc_eta[2, ])))
      expect_equal(one$settings$iter, 1L)
      expect_equal(one$settings$chains, 1L)
      expect_equal(one$settings$warmup, 0L)
      expect_equal(one$settings$draws_per_emulator, 1L)
      expect_false(one$settings$control$adapt_engaged)
      expect_equal(one$settings$control$stepsize, 0.001)
      expect_equal(as.integer(table(one$emulator_draw)), c(1L, 1L))
      expect_true(all(is.na(summary(one)$diagnostics$max_rhat)))
      expect_true(all(is.na(summary(one)$diagnostics$min_n_eff)))
    }
  }
})
