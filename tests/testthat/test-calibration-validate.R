# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

test_that("sampling metadata is read, upgraded, or rejected", {
  method <- spDBL:::.cal_sampling_method
  expect_identical(method("independent_smoothed_marginals"), "marginal")
  expect_identical(method(NULL, TRUE), "joint")
  expect_identical(method(NULL, FALSE), "marginal")
  expect_identical(method("joint", TRUE), "joint")
  expect_error(method(NULL, NULL), "sampling-method metadata")
  expect_error(method("joint", FALSE), "metadata disagree")
  expect_error(method("other"), "should be one of")
})

test_that("covariance and seed helpers validate their inputs", {
  expect_error(spDBL:::.cal_spd(matrix(c(1, 2, 0, 1), 2), "Q"), "Q must be finite")
  expect_error(spDBL:::.cal_spd(matrix(1, 2, 3)), "square")
  expect_equal(crossprod(spDBL:::.cal_spd(diag(c(4, 9)))), diag(c(4, 9)))
  for (bad in list(-1, 1.5, c(1, 2), "1", NA_real_))
    expect_error(spDBL:::.cal_seed(bad), "nonnegative integer")
  # A local seed leaves no global RNG state behind when none existed before.
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  if (had) rm(".Random.seed", envir = globalenv())
  restore <- spDBL:::.cal_seed(3L)
  first <- stats::runif(1)
  restore()
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  restore <- spDBL:::.cal_seed(3L)
  expect_identical(stats::runif(1), first)
  restore()
})

test_that("the offline design GP checks bounds and handles fixed inputs", {
  input <- cbind(a = seq(0, 1, length.out = 6), b = 0.5)
  F_ls <- list(cbind(1, sin(3 * input[, 1])), cbind(1, cos(2 * input[, 1])))
  fit <- spDBL:::fit_offline_ar2(input, F_ls, lower = c(0, 0.5), upper = c(1, 0.5))
  expect_s3_class(fit, "offline_ar2_predictor")
  expect_identical(fit$active, 1L)
  # The constant intercept columns are predicted exactly.
  pred <- spDBL:::predict_offline_ar2(fit, c(0.37, 0.5))
  expect_equal(c(pred$F[[1]][1], pred$F[[2]][1]), c(1, 1))
  expect_length(pred$variance, 4L)
  expect_true(all(pred$variance >= 0))
  # Prediction at a training input recovers the training design (small nugget).
  at_train <- spDBL:::predict_offline_ar2(fit, input[3, ])
  expect_equal(at_train$F[[1]][2], F_ls[[1]][3, 2], tolerance = 0.05)
  fit_err <- function(...) spDBL:::fit_offline_ar2(...)
  expect_error(fit_err(input, F_ls, c(0, 1), c(1, 0)), "Invalid input bounds")
  expect_error(fit_err(input, F_ls, c(0.5, 0.5), c(0.5, 0.5)), "Invalid input bounds")
  expect_error(fit_err(input, F_ls, c(0, 0), c(1, 0)), "exceed bounds")
  expect_error(fit_err(input, F_ls, c(0.1, 0), c(1, 1)), "exceed bounds")
  expect_error(fit_err(input, list(F_ls[[1]][1:3, ]), c(0, 0), c(1, 1)), "finite training design")
  expect_error(fit_err(input, list(F_ls[[1]], cbind(F_ls[[1]], 1)), c(0, 0), c(1, 1)),
               "Design widths must agree")
  expect_error(spDBL:::predict_offline_ar2(fit, c(1.2, 0.5)), "eta exceeds bounds")
})

test_that("calibration inputs validate the episode mapping and draws", {
  f <- calibration_fixture()
  e <- f$emulator
  build <- function(...) {
    args <- list(input = e$input, Y = e$Y, F_ls = e$F_ls, posterior = e$posterior,
                 block_cols = e$block_cols, R_field = e$R_field, phi_para = 2,
                 lower = c(0, 0), upper = c(1, 1))
    new <- list(...)
    args[names(new)] <- new
    do.call(spDBL:::make_offline_calibration_input, args)
  }
  expect_error(build(block_cols = matrix(c(1, 1, 2, 3), 2)), "episode mapping")
  expect_error(build(phi_para = 0), "episode mapping")
  bad_post <- e$posterior
  bad_post$theta[[2]] <- bad_post$theta[[2]][, , 1, drop = FALSE]
  expect_error(build(posterior = bad_post), "Incompatible episode arrays")
  expect_error(spDBL:::predict_calibration(f$eta, 3L, e), "Invalid draw or locations")
  expect_error(spDBL:::predict_calibration(f$eta, 1L, e, c(1, 1)), "Invalid draw or locations")
  shrunk <- e
  shrunk$V_chol <- e$V_chol / 100
  expect_error(spDBL:::predict_calibration(f$eta, 1L, shrunk), "Nonpositive emulator")
})

test_that("calibration_prepare validates draws and names unnamed inputs", {
  em <- tiny_emulator()
  expect_error(calibration_prepare(em, c(0, 0), c(1, 1), n_draws = 0), "positive integer")
  unnamed <- em
  colnames(unnamed$setup$pde_para_train) <- NULL
  cal <- calibration_prepare(unnamed, c(0, 0), c(1, 1))
  expect_identical(cal$parameter_names, c("eta1", "eta2"))
  expect_identical(cal$time, 0:2)
  expect_identical(cal$AR_choice, 1L)
})

test_that("held-out validation reproduces training runs with known designs", {
  em <- tiny_emulator()
  cal <- calibration_prepare(em, c(0, 0), c(1, 1))
  input <- em$setup$pde_para_train
  responses <- lapply(1:3, function(t)
    outer(input[, 1] + t * input[, 2], 1:4, function(x, y) sin(x) + y / 10))
  rows <- c(2, 7)
  out <- calibration_validate(cal, input[rows, ], lapply(responses, function(y) y[rows, ]))
  expect_named(out, c("design", "field", "uncertainty"))
  expect_identical(out$uncertainty, cal$predictor$uncertainty)
  expect_identical(nrow(out$design), cal$E)
  expect_equal(out$design$time, rep(0:2, each = 2))
  # Design errors are summaries of the offline GP errors at the held-out runs.
  index <- data.frame(block_no = 1:2, U = 1, D = 2, L = c(1, 3), R = c(2, 4))
  F_true <- gen_F_ls_AR1_EP(lapply(responses, function(y) y[rows, ]), 3, 2, index)
  for (k in seq_len(cal$E)) {
    est <- rbind(spDBL:::predict_offline_ar2(cal$predictor, input[rows[1], ])$F[[k]],
                 spDBL:::predict_offline_ar2(cal$predictor, input[rows[2], ])$F[[k]])
    expect_equal(out$design$rmse[k], sqrt(mean((est - F_true[[k]])^2)))
    expect_equal(out$design$max_abs_error[k], max(abs(est - F_true[[k]])))
  }
  # With the true lagged design the kriging mean interpolates training runs.
  expect_identical(out$field$input, 1:2)
  expect_true(all(out$field$rmse_known_design < 1e-2))
  expect_true(all(out$field$rmse_offline_design >= 0))
  expect_error(calibration_validate(cal, input[rows, ], responses[1:2]), "do not match")
  expect_error(calibration_validate(em, input[rows, ], responses), "do not match")
  bad <- cal
  bad$block_cols <- matrix(c(1L, 3L, 2L, 4L), 2, byrow = TRUE)
  expect_error(calibration_validate(bad, input[rows, ], lapply(responses, function(y) y[rows, ])),
               "contiguous strip blocks")
  expect_error(spDBL:::validate_offline_ar2(cal, input[rows, ], F_true[1:2]), "incompatible")
})

test_that("AR2 calibration inputs validate with two-lag designs", {
  input <- cbind(a = seq(0.05, 0.95, length.out = 8),
                 b = c(0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4, 0.8))
  responses <- lapply(1:4, function(t)
    outer(input[, 1] + t * input[, 2], 1:4, function(x, y) sin(x) + y / 10))
  em <- emulator_learn(input, responses, Nx = 2L, Ny = 2L, AR_choice = 2L,
    covariance = "sigma2R", nsam = 2L, seed = 11L)
  cal <- calibration_prepare(em, c(0, 0), c(1, 1))
  expect_identical(cal$p, 4L)
  set.seed(12)
  out <- calibration_validate(cal, input[4, ], lapply(responses, function(y) y[4, , drop = FALSE]))
  expect_identical(nrow(out$field), 1L)
  expect_true(is.finite(out$field$rmse_offline_design))
})
