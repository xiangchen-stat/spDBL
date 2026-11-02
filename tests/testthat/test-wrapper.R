# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

wrapper_data <- function(nT = 3) {
  input <- cbind(a = seq(0.05, 0.95, length.out = 8),
                 b = c(0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4, 0.8))
  responses <- lapply(seq_len(nT), function(t)
    outer(input[, 1] + t * input[, 2], 1:4, function(x, y) sin(x) + y / 10))
  list(input = input, responses = responses)
}

# Episode-wise kriging written out with explicit matrices.
episode_kriging <- function(s, input_new, F_new, theta) {
  V <- exp_kernel_loop(s$pde_para_train, s$pde_para_train, s$phi_para, s$gp_sigma2) +
    diag(s$gp_tau2, s$n_train)
  J <- exp_kernel_loop(s$pde_para_train, input_new, s$phi_para, s$gp_sigma2)
  fields <- lapply(seq_len(s$nT_ori), function(t) matrix(NA_real_, nrow(input_new), s$N_sp))
  for (e in seq_len(s$nT)) {
    t <- (e - 1L) %/% s$n_b + 1L
    b <- (e - 1L) %% s$n_b + 1L
    pred <- F_new[[e]] %*% theta[[e]] +
      t(J) %*% solve(V, s$Y[[e]] - s$F_ls[[e]] %*% theta[[e]])
    fields[[t]][, s$block_cols[b, ]] <- pred
  }
  fields
}

test_that("prepare_data splits runs and reshapes outputs by time", {
  set.seed(71)
  dir <- tempfile("prep")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  para <- matrix(runif(10), 5, 2, dimnames = list(NULL, c("beta", "gamma")))
  utils::write.csv(para, file.path(dir, "para.csv"), row.names = FALSE)
  runs <- lapply(1:5, function(j) matrix(rnorm(12), 3, 4))
  files <- vapply(1:5, function(j) {
    path <- file.path(dir, paste0("run", j, ".csv"))
    utils::write.csv(runs[[j]], path, row.names = FALSE)
    path
  }, "")
  out <- prepare_data(file.path(dir, "para.csv"), files, prop_train = 0.6)
  expect_identical(c(out$n_train, out$n_test), c(3, 2))
  expect_equal(out$pde_para_train, para[1:3, ], tolerance = 1e-12)
  expect_equal(out$pde_para_test, para[4:5, ], tolerance = 1e-12)
  expect_named(out$dt_pde_train, c("T1", "T2", "T3"))
  for (t in 1:3) {
    expect_equal(unname(out$dt_pde_train[[t]]),
                 t(vapply(runs[1:3], function(r) r[t, ], numeric(4))), tolerance = 1e-12)
    expect_equal(unname(out$dt_pde_test[[t]]),
                 t(vapply(runs[4:5], function(r) r[t, ], numeric(4))), tolerance = 1e-12)
  }
  # A single test run stays a one-row matrix.
  one <- prepare_data(file.path(dir, "para.csv"), files, prop_train = 0.8)
  expect_identical(dim(one$pde_para_test), c(1L, 2L))
  expect_identical(colnames(one$pde_para_test), c("beta", "gamma"))
  expect_equal(one$pde_para_test, para[5, , drop = FALSE], tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_error(prepare_data(file.path(dir, "para.csv"), files, prop_train = 1),
               "at least one training run and one test run")
  expect_error(prepare_data(file.path(dir, "para.csv"), files, prop_train = 0.05),
               "at least one training run and one test run")
})

test_that("emulator_learn rejects invalid settings and data", {
  w <- wrapper_data()
  learn <- function(...) emulator_learn(w$input, w$responses, 2L, 2L, AR_choice = 1L,
                                        nsam = 2L, seed = 1L, ...)
  expect_error(learn(episode_window = 3), "Invalid grid")
  expect_error(emulator_learn(w$input, w$responses, 2L, 2L, AR_choice = 3), "Invalid grid")
  expect_error(learn(gp_tune = 0), "Invalid grid")
  expect_error(emulator_learn(w$input, lapply(w$responses, function(y) y[, 1:3]),
                              2L, 2L, AR_choice = 1L), "incompatible dimensions")
  expect_error(emulator_learn(w$input[1, , drop = FALSE],
                              lapply(w$responses, function(y) y[1, , drop = FALSE]),
                              2L, 2L, AR_choice = 1L), "incompatible dimensions")
  expect_error(emulator_learn(w$input[c(1, 1), ],
                              lapply(w$responses, function(y) y[1:2, ]),
                              2L, 2L, AR_choice = 1L), "two distinct training inputs")
  expect_error(learn(R_episode = diag(3)), "Spatial covariance dimensions")
  # The causal AR(2) design needs no third observation.
  expect_s3_class(emulator_learn(w$input, w$responses[1:2], 2L, 2L, AR_choice = 2L,
                                 nsam = 2L, seed = 1L), "spDBL_emulator")
  expect_error(learn(F_ls_train = list(diag(2))), "one compatible design matrix")
})

test_that("default AR designs predict training runs up to the nugget correction", {
  w <- wrapper_data()
  for (model in c("iw", "sigma2R")) {
    em <- emulator_learn(w$input, w$responses, 2L, 2L, AR_choice = 1L,
                         covariance = model, nsam = 2L, seed = 3L)
    s <- em$setup
    rows <- c(2, 5)
    pred <- emulator_predict(em, w$input[rows, ],
                             dt_pde_test = lapply(w$responses, function(y) y[rows, ]))
    expect_named(pred, c("T1", "T2", "T3"))
    # With J = V - tau2 I on training inputs, the kriging mean at a training
    # run equals Y_i - tau2 (V^{-1} r)_i, with r the smoothed residuals.
    V <- exp_kernel_loop(w$input, w$input, s$phi_para, s$gp_sigma2) +
      diag(s$gp_tau2, s$n_train)
    for (t in 1:3) {
      expected <- matrix(NA_real_, 2, 4)
      for (b in 1:2) {
        e <- (t - 1) * 2 + b
        r <- s$Y[[e]] - s$F_ls[[e]] %*% em$para_ffbs$bs$st[, , e]
        expected[, s$block_cols[b, ]] <- s$Y[[e]][rows, ] - s$gp_tau2 * solve(V, r)[rows, ]
      }
      expect_equal(pred[[t]], expected, tolerance = 1e-8, ignore_attr = TRUE)
      expect_equal(pred[[t]], w$responses[[t]][rows, ], tolerance = 1e-2, ignore_attr = TRUE)
    }
  }
  expect_error(emulator_predict(em, w$input[1, ]), "Supply F_new_ls or dt_pde_test")
})

test_that("AR2 emulators build default designs for new runs", {
  w <- wrapper_data(nT = 4)
  em <- emulator_learn(w$input, w$responses, 2L, 2L, AR_choice = 2L,
                       covariance = "sigma2R", nsam = 2L, seed = 4L)
  expect_identical(ncol(em$setup$F_ls[[1]]), 4L)
  pred <- emulator_predict(em, w$input[3, ],
                           dt_pde_test = lapply(w$responses, function(y) y[3, , drop = FALSE]))
  expect_length(pred, 4L)
  # Lags at t >= 3 are actual responses, so the training run is reproduced.
  for (t in 3:4) expect_equal(pred[[t]], w$responses[[t]][3, , drop = FALSE],
                              tolerance = 1e-2, ignore_attr = TRUE)
})

test_that("MC prediction applies the kriging mean to every coefficient draw", {
  w <- wrapper_data()
  em <- emulator_learn(w$input, w$responses, 2L, 2L, AR_choice = 1L,
                       covariance = "sigma2R", nsam = 3L, seed = 5L)
  s <- em$setup
  input_new <- rbind(c(0.33, 0.45), c(0.71, 0.25))
  F_new <- lapply(s$F_ls, function(f) f[1:2, , drop = FALSE] + 0.05)
  mc <- emulator_predict(em, input_new, F_new_ls = F_new, MC = TRUE)
  expect_named(mc, c("T1", "T2", "T3"))
  expect_identical(dim(mc$T1), c(2L, 4L, 3L))
  for (m in 1:3) {
    theta <- lapply(seq_len(s$nT), function(e) em$res_ffbs[[e]][, , m])
    expected <- episode_kriging(s, input_new, F_new, theta)
    for (t in 1:3) expect_equal(mc[[t]][, , m], expected[[t]], tolerance = 1e-10,
                                ignore_attr = TRUE)
  }
  smoothed <- lapply(seq_len(s$nT), function(e) em$para_ffbs$bs$st[, , e])
  exact <- emulator_predict(em, input_new, F_new_ls = F_new)
  expect_equal(exact, setNames(episode_kriging(s, input_new, F_new, smoothed),
                               c("T1", "T2", "T3")), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("prepared emulators return draw-wise or draw-averaged predictions", {
  w <- wrapper_data()
  em <- emulator_learn(w$input, w$responses, 2L, 2L, AR_choice = 1L,
                       covariance = "sigma2R", nsam = 3L, seed = 6L)
  cal <- calibration_prepare(em, c(0, 0), c(1, 1))
  new <- rbind(c(0.2, 0.4), c(0.6, 0.3))
  mc <- emulator_predict(cal, new, MC = TRUE)
  avg <- emulator_predict(cal, new)
  expect_identical(dim(mc[[1]]), c(2L, 4L, 3L))
  for (t in 1:3) {
    expect_equal(avg[[t]], apply(mc[[t]], 1:2, mean))
    for (m in 1:3) expect_equal(mc[[t]][2, , m],
      spDBL:::predict_calibration(new[2, ], m, cal)$mu[t, ])
  }
})

test_that("verbose emulator fitting reports filtering and sampling", {
  w <- wrapper_data()
  for (model in c("iw", "sigma2R")) {
    msgs <- capture_messages(emulator_learn(w$input, w$responses, 2L, 2L,
      AR_choice = 1L, covariance = model, nsam = 2L, seed = 7L, verbose = TRUE))
    expect_true(any(grepl("Finish FFBS", msgs)))
    expect_match(msgs[length(msgs)], "Sampled 2 marginal coefficient sets")
  }
})
