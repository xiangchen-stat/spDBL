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

test_that("the design GP defaults to a Matern 5/2 kernel per time step", {
  em <- tiny_emulator()
  cal <- calibration_prepare(em, c(0, 0), c(1, 1))
  expect_identical(cal$design_kernel, "matern52")
  expect_identical(cal$design_groups, "time")
  expect_identical(cal$predictor$kernel, "matern52")
  expect_identical(cal$predictor$groups, (seq_len(cal$E) - 1L) %/% cal$B + 1L)
  expect_identical(dim(cal$predictor$lengthscale), c(cal$T, 2L))
  expect_length(cal$predictor$K_chol, cal$T)
  # The earlier design GP: one squared-exponential kernel for all times.
  old <- calibration_prepare(em, c(0, 0), c(1, 1), design_kernel = "sqexp",
                             design_groups = "shared")
  expect_identical(old$predictor$groups, rep(1L, old$E))
  expect_identical(dim(old$predictor$lengthscale), c(1L, 2L))
  expect_identical(old$posterior, cal$posterior)
  shared <- calibration_prepare(em, c(0, 0), c(1, 1), design_groups = "shared")
  expect_identical(shared$predictor$kernel, "matern52")
  for (object in list(cal, old, shared)) {
    p <- spDBL:::predict_calibration(c(0.3, 0.6), 2L, object)
    expect_length(p$design_variance, object$T)
    expect_true(all(is.finite(p$mu)))
  }
  expect_error(calibration_prepare(em, c(0, 0), c(1, 1), design_kernel = "exp"),
               "should be one of")
  expect_error(calibration_prepare(em, c(0, 0), c(1, 1), design_groups = "site"),
               "should be one of")
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

test_that("fitted objects record their seed and spDBL version", {
  version <- as.character(utils::packageVersion("spDBL"))
  em <- tiny_emulator()
  expect_identical(em$setup$seed, 91L)
  expect_identical(em$setup$spDBL_version, version)
  reused <- calibration_prepare(em, c(0, 0), c(1, 1))
  expect_identical(reused$seed, 91L)
  expect_identical(reused$spDBL_version, version)
  fresh <- calibration_prepare(em, c(0, 0), c(1, 1), n_draws = 2L, seed = 7L)
  expect_identical(fresh$seed, 7L)
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
  # Cross blocks for site prediction match the square blocks entry by entry,
  # and both objects keep the grid coordinates for off-grid prediction.
  for (cal in list(cal_lazy, cal_explicit)) {
    expect_identical(spDBL:::.cal_grid_coords(cal), unname(as.matrix(s$ind_sp)))
    expect_identical(spDBL:::.cal_R_pairs(cal, c(1L, 2L, 4L), c(1L, 3L, 4L)),
      spDBL:::.cal_R_field(cal, 1:4)[cbind(c(1, 2, 4), c(1, 3, 4))])
    expect_identical(spDBL:::.cal_R_cross(cal, 1:4, 1:4), spDBL:::.cal_R_field(cal, 1:4))
    expect_identical(spDBL:::.cal_R_cross(cal, c(3L, 1L), c(4L, 2L)),
                     spDBL:::.cal_R_field(cal, 1:4)[c(3, 1), c(4, 2)])
  }
  z <- matrix(0.5, 3, 2)
  expect_equal(spDBL:::make_stan_calibration_data(cal_lazy, 1L, z, diag(4), c(2L, 4L))$R_field,
               R_full[c(2, 4), c(2, 4)], ignore_attr = TRUE)
  expect_error(spDBL:::make_offline_calibration_input(cal_lazy$input, cal_lazy$Y,
    cal_lazy$F_ls, cal_lazy$posterior, cal_lazy$block_cols, NULL, 2, c(0, 0), c(1, 1)),
    "Supply R_field or a valid field_kernel")
})

test_that("conditional covariances are factored even when singular", {
  x <- matrix(1, 2, 2)  # rank 1: chol() fails
  expect_error(chol(x))
  L <- spDBL:::.cal_psd_factor(x)
  expect_equal(L %*% t(L), x)
  expect_equal(spDBL:::.cal_psd_factor(diag(2)), diag(2))
  expect_error(spDBL:::.cal_psd_factor(diag(c(1, -1)), "Test covariance"),
               "Test covariance is not positive semidefinite")
})

test_that("points map to their containing cell or four surrounding cell centres", {
  grid <- cbind(rep(1:3, times = 2), rep(1:2, each = 3))  # 3 rows x 2 columns
  cell <- spDBL:::.cal_point_cells(grid, rbind(c(1.2, 1.4), c(2.5, 1.5), c(3.4, 1.49)), "cell")
  expect_identical(as.integer(cell$cells), c(1L, 6L, 3L))
  bil <- spDBL:::.cal_point_cells(grid, rbind(c(1.2, 1.4), c(3, 2), c(0.6, 1.5)), "bilinear")
  # (1.2, 1.4): corners (1,1), (2,1), (1,2), (2,2) = sites 1, 2, 4, 5.
  expect_equal(bil$weights[1, ], c(0.8 * 0.6, 0.2 * 0.6, 0.8 * 0.4, 0.2 * 0.4))
  expect_identical(bil$cells[1, ], c(1L, 2L, 4L, 5L))
  # A cell centre on the last row and column: all weight on site 6.
  expect_identical(bil$cells[2, !is.na(bil$cells[2, ])], 6L)
  expect_equal(sum(bil$weights[2, ]), 1)
  # Beyond the outer cell centres the point is clamped to them.
  expect_equal(bil$weights[3, ], c(0.5, 0, 0.5, 0))
  expect_identical(bil$cells[3, c(1, 3)], c(1L, 4L))
  expect_equal(rowSums(bil$weights), rep(1, 3))
  expect_error(spDBL:::.cal_point_cells(grid, rbind(c(3.6, 1)), "cell"), "within the grid's cells")
  expect_error(spDBL:::.cal_point_cells(grid[-2, ], rbind(c(1, 1)), "cell"), "full integer lattice")
  expect_error(spDBL:::.cal_point_cells(grid, rbind(c(NA, 1)), "cell"), "numeric matrix")
  one_row <- spDBL:::.cal_point_cells(cbind(1, 1:3), rbind(c(1.3, 2.25)), "bilinear")
  expect_equal(one_row$weights[1, ], c(0.75, 0, 0.25, 0))
})
