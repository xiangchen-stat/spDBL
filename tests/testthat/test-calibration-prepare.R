# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

test_that("the installed emulator API prepares PDE-free predictions", {
  local_mocked_bindings(gen_pde = function(...) stop("PDE solver must not run"),
                        .package = "spDBL")
  set.seed(10)
  before <- .Random.seed
  em <- tiny_emulator()
  expect_identical(before, .Random.seed)
  expect_s3_class(em, "spDBL_emulator")
  expect_identical(attr(em$res_ffbs, "sampling_method"), "marginal")
  expect_false(attr(em$res_ffbs, "joint_trajectories"))
  cal <- calibration_prepare(em, c(b = 0, a = 0), c(b = 1, a = 1))
  expect_s3_class(cal, "spDBL_calibration_input")
  expect_equal(cal$M, 3L)
  expect_identical(cal$sampling_method, "marginal")
  expect_false(cal$posterior$joint_trajectories)
  fresh <- calibration_prepare(em, c(0, 0), c(1, 1), n_draws = 4L)
  expect_equal(fresh$M, 4L)
  expect_identical(fresh$sampling_method, "marginal")
  expect_identical(fresh$F_ls, em$setup$F_ls)
  expect_identical(before, .Random.seed)
  pred <- emulator_predict(cal, c(0.4, 0.6))
  expect_length(pred, 3L)
  expect_identical(dim(pred[[1]]), c(1L, 4L))
  expect_true(all(is.finite(unlist(pred))))
  expect_error(calibration_prepare(em, c(a = 0, c = 0), c(a = 1, c = 1)), "names")
  expect_error(calibration_prepare(em, c(0.5, 0), c(1, 1)), "exceed")
  expect_error(calibration_prepare(tiny_emulator("iw"), c(0, 0), c(1, 1)), "sigma2R")
})

test_that("calibration preparation inherits or overrides coefficient sampling", {
  em <- tiny_emulator()
  before <- .Random.seed
  joint <- calibration_prepare(em, c(0, 0), c(1, 1), sampling = "joint")
  expect_true(joint$posterior$joint_trajectories)
  expect_identical(joint$sampling_method, "joint")
  expect_equal(joint$M, 3L)
  expect_identical(joint$F_ls, em$setup$F_ls)
  expect_identical(before, .Random.seed)
  em_joint <- tiny_emulator(sampling = "joint")
  inherited <- calibration_prepare(em_joint, c(0, 0), c(1, 1))
  expect_identical(inherited$sampling_method, "joint")
  expect_identical(inherited$posterior$theta, em_joint$res_ffbs[1:6])
})

test_that("custom designs retain wrapper dimensions", {
  em <- tiny_emulator()
  custom <- emulator_learn(em$setup$pde_para_train,
    lapply(1:3, function(t) do.call(cbind, em$setup$Y[((t - 1) * 2 + 1):(t * 2)])),
    2L, 2L, F_ls_train = em$setup$F_ls, covariance = "sigma2R", seed = 1L)
  new_F <- lapply(em$setup$F_ls, function(f) f[1, , drop = FALSE])
  predicted <- emulator_predict(custom, em$setup$pde_para_train[1, ], F_new_ls = new_F)
  expect_identical(dim(predicted[[1]]), c(1L, 4L))
})

test_that("a single episode works with scalar designs", {
  em <- emulator_learn(matrix(c(0, 1), 2), list(matrix(c(0.1, 0.3), 2)),
    1L, 1L, AR_choice = 1L, covariance = "sigma2R", nsam = 2L, seed = 2L)
  cal <- calibration_prepare(em, 0, 1)
  pred <- emulator_predict(cal, 0.4)
  expect_identical(dim(pred[[1]]), c(1L, 1L))
})

test_that("the bundled calibration design has independent provenance", {
  env <- new.env()
  data("dt_calibration", package = "spDBL", envir = env)
  d <- env$dt_calibration
  expect_identical(dim(d$training$input), c(40L, 5L))
  expect_identical(dim(d$field$z), c(26L, 36L))
  expect_equal(unname(diag(d$field$U)), rep(1, 36L))
  expect_false(d$truth$input_id %in% d$provenance$train_ids)
  expect_true(all(d$training$input[, 1] >= 2))
  expect_false(any(c("eta", "truth", "y") %in% names(d$field)))
})

test_that("calibration builds the field covariance only where it is needed", {
  input <- cbind(a = seq(0.05, 0.95, length.out = 8),
                 b = c(0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4, 0.8))
  responses <- lapply(1:3, function(t)
    outer(input[, 1] + t * input[, 2], 1:4, function(x, y) sin(x) + y / 10))
  learn <- function(...) emulator_learn(input, responses, Nx = 2L, Ny = 2L,
    AR_choice = 1L, covariance = "sigma2R", nsam = 3L, seed = 91L, ...)
  lazy <- learn()
  expect_null(lazy$setup$R_field)  # no S x S matrix in the fitted emulator
  s <- lazy$setup
  # The kernel emulator_learn() built and stored in spDBL 1.1.1.
  R_full <- s$gp_sigma2 * exp(-s$phi_sp * as.matrix(dist(s$ind_sp))) +
    diag(s$gp_tau2, s$N_sp)
  explicit <- learn(R_field = R_full)
  expect_identical(explicit$setup$R_field, R_full)
  cal_lazy <- calibration_prepare(lazy, c(0, 0), c(1, 1))
  cal_explicit <- calibration_prepare(explicit, c(0, 0), c(1, 1))
  expect_null(cal_lazy$R_field)
  expect_null(cal_explicit$field_kernel)
  for (loc in list(1:4, c(3L, 1L))) {
    expect_equal(spDBL:::.cal_R_field(cal_lazy, loc), R_full[loc, loc],
                 ignore_attr = TRUE)
    for (mode in c("full", "marginal"))
      expect_equal(spDBL:::predict_calibration(c(0.4, 0.6), 2L, cal_lazy, loc, mode),
                   spDBL:::predict_calibration(c(0.4, 0.6), 2L, cal_explicit, loc, mode),
                   ignore_attr = TRUE)
  }
  z <- matrix(0.5, 3, 2)
  expect_equal(spDBL:::make_stan_calibration_data(cal_lazy, 1L, z, diag(4), c(2L, 4L))$R_field,
               R_full[c(2, 4), c(2, 4)], ignore_attr = TRUE)
  expect_error(spDBL:::make_offline_calibration_input(cal_lazy$input, cal_lazy$Y,
    cal_lazy$F_ls, cal_lazy$posterior, cal_lazy$block_cols, NULL, 2, c(0, 0), c(1, 1)),
    "Supply R_field or a valid field_kernel")
})
