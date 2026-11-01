test_that("installed Stan programs match R and pool conditionals equally", {
  skip_if(Sys.getenv("SPDBL_RUN_STAN_TESTS") != "true",
          "Set SPDBL_RUN_STAN_TESTS=true for compiler-dependent integration checks")
  skip_if_not_installed("rstan", minimum_version = "2.32.7")
  f <- calibration_fixture()
  f$U <- 1.1001 * f$U
  e <- f$emulator
  for (marginal in c(TRUE, FALSE)) {
    model <- spDBL:::.cal_stan_model(marginal, cache = TRUE,
      cache_dir = Sys.getenv("SPDBL_STAN_CACHE_DIR", tempdir()))
    initial <- list(eta_unit = f$eta, precision = c(5, 6, 7))
    if (!marginal) initial$u_std <- matrix(seq(-0.2, 0.2, length.out = 12), 3)
    data <- spDBL:::make_stan_calibration_data(e, 1L, f$z, f$U)
    fixed <- rstan::sampling(model, data = data, chains = 1L, iter = 1L,
      warmup = 0L, algorithm = "Fixed_param", init = list(initial), refresh = 0L)
    values <- rstan::extract(fixed)
    prediction <- spDBL:::predict_calibration(f$eta, 1L, e)
    expect_equal(matrix(values$mu, e$T, e$S), prediction$mu, tolerance = 1e-9)
    expect_equal(matrix(values$Sigma, e$S, e$S), unname(prediction$Sigma), tolerance = 1e-9)
    alternative <- initial
    alternative$eta_unit <- c(0.7, 0.35)
    alternative$precision <- c(4, 8, 6)
    stan_density <- function(v) rstan::log_prob(fixed,
      rstan::unconstrain_pars(fixed, v), adjust_transform = FALSE)
    r_density <- if (marginal) function(v)
      spDBL:::calibration_marginal_log_density(v$eta_unit, v$precision, e, 1L, f$z, f$U) else
      function(v) spDBL:::calibration_log_density(v$eta_unit, v$precision,
        v$u_std, e, 1L, f$z, f$U)
    expect_equal(stan_density(alternative) - stan_density(initial),
                 r_density(alternative) - r_density(initial), tolerance = 1e-8)
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
      one <- suppressWarnings(calibration_fit(e, f$z, f$U, compiled_model = model,
        keep_discrepancy = FALSE, refresh = 0L))
      expect_true(one$execution_only)
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
