# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

method_data <- function() {
  input <- cbind(a = seq(0.05, 0.95, length.out = 8),
                 b = c(0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4, 0.8))
  responses <- lapply(1:3, function(t)
    outer(input[, 1] + t * input[, 2], 1:4, function(x, y) sin(x) + y / 10))
  list(input = input, responses = responses)
}

test_that("print and summary describe the fitted emulator", {
  d <- method_data()
  for (covariance in c("iw", "sigma2R")) {
    em <- emulator_learn(d$input, d$responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
      covariance = covariance, nsam = 4000L, seed = 5L)
    out <- capture_output(expect_identical(print(em), em))
    expect_match(out, "Training runs: 8 with 2 inputs")
    expect_match(out, "Grid: 2 x 2 \\(4 sites\\) over 3 times")
    expect_match(out, "Stored posterior draws: 4000")
    s <- summary(em)
    expect_s3_class(s, "summary.spDBL_emulator")
    expect_identical(s$settings$covariance, covariance)
    expect_identical(s$settings$sampling, "marginal")
    # In-sample RMSE of the smoothed fit, computed time by time.
    st <- em$para_ffbs$bs$st
    fitted_error <- vapply(seq_len(em$setup$nT), function(e)
      sum((em$setup$Y[[e]] - em$setup$F_ls[[e]] %*% st[, , e])^2), numeric(1))
    by_time <- tapply(fitted_error, (seq_len(em$setup$nT) - 1) %/% em$setup$n_b, sum)
    expect_equal(s$rmse_by_time, sqrt(as.numeric(by_time) / (8 * 4)))
    expect_equal(s$rmse, sqrt(sum(fitted_error) / (8 * 4 * 3)))
    expect_equal(s$response_sd, sd(unlist(d$responses)))
    # The residual variance is the posterior mean behind the stored draws.
    draws <- em$res_ffbs$Sigma
    mc <- if (covariance == "iw") mean(apply(draws, 3, function(S) mean(diag(S)))) else
      mean(draws) * mean(diag(em$setup$R_episode))
    expect_equal(s$residual_variance, mc, tolerance = 0.05)
    expect_match(capture_output(print(s)), "In-sample RMSE of the smoothed fit")
  }
})

test_that("predict methods match emulator_predict", {
  d <- method_data()
  em <- emulator_learn(d$input, d$responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
    covariance = "sigma2R", nsam = 3L, seed = 6L)
  new <- d$input[c(2, 5), ] + 0.01
  held_out <- lapply(d$responses, function(y) y[c(2, 5), ] + 0.01)
  expect_identical(predict(em, new, dt_pde_test = held_out),
                   emulator_predict(em, new, dt_pde_test = held_out))
  expect_identical(predict(em, new, dt_pde_test = held_out, MC = TRUE),
                   emulator_predict(em, new, dt_pde_test = held_out, MC = TRUE))
  cal <- calibration_prepare(em, c(0, 0), c(1, 1))
  expect_identical(predict(cal, new), emulator_predict(cal, new))
  expect_identical(predict(cal, new, MC = TRUE), emulator_predict(cal, new, MC = TRUE))
  out <- capture_output(expect_identical(print(cal), cal))
  expect_match(out, "3 emulator draws \\(marginal sampling\\), 4 sites, 3 times")
  expect_match(out, "Design-row GP variance is propagated")
  expect_match(out, "Matern 5/2 kernel, hyperparameters per time step \\(3\\)")
  shared <- calibration_prepare(em, c(0, 0), c(1, 1), design_kernel = "sqexp",
                                design_groups = "shared")
  expect_match(capture_output(print(shared)),
               "squared-exponential kernel, one set of hyperparameters")
  # A parameter held fixed in every training run is reported as fixed.
  d$input[, "b"] <- 0.5
  fixed_em <- emulator_learn(d$input, d$responses, Nx = 2L, Ny = 2L,
    AR_choice = 1L, covariance = "sigma2R", nsam = 3L, seed = 6L)
  fixed <- calibration_prepare(fixed_em, c(0, 0.5), c(1, 0.5))
  expect_match(capture_output(print(fixed)), "b +0.5 +0.5 +fixed")
})

test_that("methods accept objects saved without the newer metadata", {
  d <- method_data()
  em <- emulator_learn(d$input, d$responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
    nsam = 3L, seed = 7L)
  old <- em
  old$setup$covariance <- NULL
  old$setup$sampling_method <- NULL
  expect_identical(summary(old)$settings$covariance, "iw")
  expect_identical(summary(old)$settings$sampling, "marginal")
  attr(old$res_ffbs, "sampling_method") <- NULL
  expect_identical(summary(old)$settings$sampling, "unknown")
  expect_match(capture_output(print(old)), "inverse-Wishart covariance, unknown")
  # A degenerate filter with too few degrees of freedom has no finite mean.
  last <- length(old$para_ffbs$ff) - 1L
  old$para_ffbs$ff[[last]]$nt <- 2
  expect_true(is.na(summary(old)$residual_variance))
  s2r <- emulator_learn(d$input, d$responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
    covariance = "sigma2R", nsam = 3L, seed = 7L)
  s2r$para_ffbs$ff[[last]]$nt <- 1
  expect_true(is.na(summary(s2r)$residual_variance))
  cal <- calibration_prepare(s2r, c(0, 0), c(1, 1))
  cal$sampling_method <- NULL
  cal$parameter_names <- NULL
  out <- capture_output(print(cal))
  expect_match(out, "unknown sampling")
  expect_match(out, "eta1")
})
