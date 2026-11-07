# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

test_that("backward smoothing matches exact conditioning on all observations", {
  d <- dlm_problem(21)
  ref <- batch_dlm(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0)
  iw <- FFBS(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$n0, d$D0, d$nT)
  ig <- FFBS_sigma2R(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0,
                     n0 = 3, D0 = 1, nT = d$nT, R = diag(d$q))
  id <- FFBS_I(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$nT)
  expect_identical(iw$ff, FF(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0,
                             d$n0, d$D0, d$nT))
  expect_identical(id$ff, FF_I(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$nT))
  for (fit in list(iw, ig, id)) {
    expect_identical(dim(fit$bs$st), c(d$p, d$q, d$nT))
    expect_identical(dim(fit$bs$St), c(d$p, d$p, d$nT))
    for (t in seq_len(d$nT)) {
      expect_equal(fit$bs$st[, , t], ref[[t]]$st, tolerance = 1e-9, ignore_attr = TRUE)
      expect_equal(fit$bs$St[, , t], ref[[t]]$St, tolerance = 1e-9, ignore_attr = TRUE)
    }
  }
})

test_that("the C++ backward step is the Rauch-Tung-Striebel recursion", {
  set.seed(22)
  p <- 3; q <- 2
  mt <- matrix(rnorm(p * q), p, q); Mt <- rand_spd(p)
  G <- diag(p) + 0.1 * matrix(rnorm(p * p), p, p)
  At1 <- G %*% Mt %*% t(G) + rand_spd(p, 0.3)
  at1 <- G %*% mt
  st1 <- matrix(rnorm(p * q), p, q); St1 <- rand_spd(p, 0.5)
  B <- Mt %*% t(G) %*% solve(At1)
  out <- spDBL:::BS_1step_cpp(mt, Mt, st1, St1, at1, At1, G)
  expect_equal(out$st, mt + B %*% (st1 - at1), tolerance = 1e-10)
  expect_equal(out$St, Mt - B %*% (At1 - St1) %*% t(B), tolerance = 1e-10)
})

test_that("FFBS wrappers report progress when verbose", {
  d <- dlm_problem(23, N = rep(2L, 6), time_varying = FALSE)
  msgs <- capture_messages(FFBS(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0,
    d$M0, d$n0, d$D0, d$nT, verbose = TRUE))
  expect_true(any(grepl("^FF: 3/6", msgs)))
  expect_true(any(grepl("^BS: 3/6", msgs)))
  expect_match(msgs[length(msgs)], "Finish FFBS")
  for (fit in list(
    function() FFBS_sigma2R(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0,
                            3, 1, d$nT, diag(d$q), verbose = TRUE),
    function() FFBS_I(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$nT,
                      verbose = TRUE))) {
    msgs <- capture_messages(fit())
    expect_match(msgs[length(msgs)], "Finish FFBS")
  }
})

test_that("identity-covariance sampler draws from the smoothed matrix normal", {
  d <- dlm_problem(24, N = c(3L, 3L), time_varying = FALSE)
  fit <- FFBS_I(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$nT)
  set.seed(25)
  nsam <- 4000L
  draws <- FFBS_sampling_I(nsam, fit, d$F_ls, d$G_ls, d$nT)
  expect_named(draws, c("T1", "T2"))
  for (t in seq_len(d$nT)) {
    expect_identical(dim(draws[[t]]), c(d$p, d$q, nsam))
    v <- t(matrix(draws[[t]], d$p * d$q))
    target <- kronecker(diag(d$q), fit$bs$St[, , t])
    se <- sqrt(diag(target) / nsam)
    expect_lt(max(abs(colMeans(v) - as.numeric(fit$bs$st[, , t])) / se), 4.5)
    expect_equal(stats::cov(v), target, tolerance = 0.1)
  }
  d3 <- dlm_problem(31, N = c(3L, 3L, 3L), time_varying = FALSE)
  fit3 <- FFBS_I(d3$Y, d3$F_ls, d3$G_ls, d3$W_ls, d3$V_ls, d3$m0, d3$M0, d3$nT)
  msgs <- capture_messages(FFBS_sampling_I(5L, fit3, rep(list(d3$F_ls), 3),
                                           d3$G_ls, d3$nT, verbose = TRUE))
  expect_length(msgs, 3L)
})

test_that("identity-covariance sampler accepts a time-varying transition list", {
  d <- dlm_problem(32, N = c(3L, 3L, 3L), time_varying = FALSE)
  G_list <- rep(list(d$G_ls), d$nT)
  fit <- FFBS_I(d$Y, d$F_ls, G_list, d$W_ls, d$V_ls, d$m0, d$M0, d$nT)
  expect_length(FFBS_sampling_I(2L, fit, d$F_ls, G_list, d$nT), d$nT)
})

test_that("IW and sigma2R samplers validate their arguments", {
  d <- dlm_problem(26, N = c(3L, 3L), time_varying = FALSE)
  iw <- FFBS(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$n0, d$D0, d$nT)
  ig <- FFBS_sigma2R(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, 3, 1,
                     d$nT, diag(d$q))
  expect_error(FFBS_sampling(0, iw, d$F_ls, d$G_ls, d$nT), "Invalid sample count")
  expect_error(FFBS_sampling(2, iw, d$F_ls, d$G_ls, d$nT + 1), "episode count")
  expect_error(FFBS_sampling(2, iw, d$F_ls, d$G_ls, d$nT, delta = 0.9), "delta = 1")
  expect_error(FFBS_sampling_sigma2R(2, ig, d$F_ls, d$G_ls, d$nT, diag(d$q),
                                     delta = 0.9), "delta = 1")
  expect_error(FFBS_sampling_sigma2R(2, ig, d$F_ls, d$G_ls, d$nT + 1, diag(d$q)),
               "nT does not match")
  expect_error(spDBL:::sample_calibration_emulator(ig, diag(d$q), d$G_ls, 0),
               "Invalid filtering results or sample count")
  expect_error(spDBL:::sample_calibration_emulator(ig, diag(d$q + 1), d$G_ls, 2),
               "wrong size")
  expect_message(FFBS_sampling(2, iw, d$F_ls, d$G_ls, d$nT, verbose = TRUE, seed = 1L),
                 "Sampled 2 marginal coefficient sets")
  expect_message(FFBS_sampling_sigma2R(2, ig, d$F_ls, d$G_ls, d$nT, diag(d$q),
                                       verbose = TRUE, seed = 1L),
                 "Sampled 2 marginal coefficient sets")
})

test_that("joint samplers accept time-varying transition lists", {
  d <- dlm_problem(27, N = c(3L, 3L, 3L), time_varying = FALSE)
  G_list <- rep(list(d$G_ls), d$nT)
  iw <- FFBS(d$Y, d$F_ls, G_list, d$W_ls, d$V_ls, d$m0, d$M0, d$n0, d$D0, d$nT)
  expect_identical(
    FFBS_sampling(3, iw, d$F_ls, G_list, d$nT, sampling = "joint", seed = 4L),
    FFBS_sampling(3, iw, d$F_ls, d$G_ls, d$nT, sampling = "joint", seed = 4L))
  ig <- FFBS_sigma2R(d$Y, d$F_ls, G_list, d$W_ls, d$V_ls, d$m0, d$M0, 3, 1,
                     d$nT, diag(d$q))
  expect_identical(
    FFBS_sampling_sigma2R(3, ig, d$F_ls, G_list, d$nT, diag(d$q),
                          sampling = "joint", seed = 4L),
    FFBS_sampling_sigma2R(3, ig, d$F_ls, d$G_ls, d$nT, diag(d$q),
                          sampling = "joint", seed = 4L))
})

# A small prediction problem with explicit kernel matrices.
prediction_problem <- function(seed, nsam = 3L) {
  set.seed(seed)
  N <- 5; N_new <- 2; d <- 2; p <- 2; q <- 3; nT <- 2
  input <- matrix(runif(N * d), N, d)
  input_new <- matrix(runif(N_new * d), N_new, d)
  list(N = N, N_new = N_new, p = p, q = q, nT = nT, input = input,
       input_new = input_new, phi = 1.7, sigma2 = 1.3, tau2 = 0.05,
       Y = lapply(1:nT, function(t) matrix(rnorm(N * q), N, q)),
       F_ls = lapply(1:nT, function(t) matrix(rnorm(N * p), N, p)),
       F_new = lapply(1:nT, function(t) matrix(rnorm(N_new * p), N_new, p)),
       theta = lapply(1:nT, function(t) array(rnorm(p * q * nsam), c(p, q, nsam))))
}

kriging_mean <- function(pp, theta, t, input_new = pp$input_new, F_new = pp$F_new[[t]]) {
  V <- exp_kernel_loop(pp$input, pp$input, pp$phi, pp$sigma2) + diag(pp$tau2, pp$N)
  J <- exp_kernel_loop(pp$input, input_new, pp$phi, pp$sigma2)
  F_new %*% theta + t(J) %*% solve(V, pp$Y[[t]] - pp$F_ls[[t]] %*% theta)
}

test_that("Monte Carlo prediction applies the exact predictor to every draw", {
  pp <- prediction_problem(28)
  mc <- FFBS_predict_MC(3L, pp$Y, pp$theta, pp$F_ls, pp$F_new, pp$input,
                        pp$input_new, pp$nT, phi_para = pp$phi,
                        gp_sigma2 = pp$sigma2, gp_tau2 = pp$tau2)
  expect_named(mc, c("T1", "T2"))
  for (t in seq_len(pp$nT)) {
    expect_equal(dim(mc[[t]]), c(pp$N_new, pp$q, 3))
    for (m in 1:3)
      expect_equal(mc[[t]][, , m], kriging_mean(pp, pp$theta[[t]][, , m], t),
                   tolerance = 1e-10)
  }
  one <- FFBS_predict_MC(2L, pp$Y, pp$theta, pp$F_ls,
                         lapply(pp$F_new, function(f) f[1, , drop = FALSE]),
                         pp$input, pp$input_new[1, ], pp$nT, phi_para = pp$phi,
                         gp_sigma2 = pp$sigma2, gp_tau2 = pp$tau2)
  expect_equal(one$T2[1, , 2], as.numeric(mc$T2[1, , 2]), tolerance = 1e-10)
})

test_that("exact prediction uses smoothed means and returns the kriging covariance", {
  pp <- prediction_problem(29)
  st <- array(c(pp$theta[[1]][, , 1], pp$theta[[2]][, , 1]), c(pp$p, pp$q, pp$nT))
  out <- FFBS_predict_exact(pp$Y, list(bs = list(st = st)), pp$F_ls, pp$F_new,
                            pp$input, pp$input_new, pp$nT, pp$phi, pp$sigma2,
                            pp$tau2, return_covariance = TRUE)
  for (t in seq_len(pp$nT))
    expect_equal(out$mean[[t]], kriging_mean(pp, st[, , t], t), tolerance = 1e-10,
                 ignore_attr = TRUE)
  V <- exp_kernel_loop(pp$input, pp$input, pp$phi, pp$sigma2) + diag(pp$tau2, pp$N)
  J <- exp_kernel_loop(pp$input, pp$input_new, pp$phi, pp$sigma2)
  V_new <- exp_kernel_loop(pp$input_new, pp$input_new, pp$phi, pp$sigma2) +
    diag(pp$tau2, pp$N_new)
  expect_equal(out$input_covariance, V_new - t(J) %*% solve(V, J), tolerance = 1e-10,
               ignore_attr = TRUE)
  # Designs may be a single shared matrix, a length-one list, or a vector.
  th <- list(st[, , 1], st[, , 1])
  shared <- FFBS_predict_exact(pp$Y, F_ls = pp$F_ls[[1]], F_new_ls = list(pp$F_new[[1]]),
    input = pp$input, input_new = pp$input_new, nT = pp$nT, phi_para = pp$phi,
    gp_sigma2 = pp$sigma2, gp_tau2 = pp$tau2, theta = th)
  pp_shared <- pp
  pp_shared$F_ls <- rep(pp$F_ls[1], 2)
  expect_equal(shared$T2, kriging_mean(pp_shared, st[, , 1], 2, F_new = pp$F_new[[1]]),
               tolerance = 1e-10, ignore_attr = TRUE)
  vec <- FFBS_predict_exact(pp$Y, F_ls = pp$F_ls, F_new_ls = as.numeric(pp$F_new[[1]][1, ]),
    input = pp$input, input_new = pp$input_new[1, ], nT = pp$nT, phi_para = pp$phi,
    gp_sigma2 = pp$sigma2, gp_tau2 = pp$tau2, theta = th)
  expect_equal(vec$T1, shared$T1[1, , drop = FALSE], tolerance = 1e-10)
})

test_that("exact prediction rejects inconsistent inputs", {
  pp <- prediction_problem(30)
  th <- lapply(pp$theta, function(x) x[, , 1])
  call <- function(...) {
    args <- list(Y = pp$Y, F_ls = pp$F_ls, F_new_ls = pp$F_new,
      input = pp$input, input_new = pp$input_new, nT = pp$nT, phi_para = pp$phi,
      theta = th)
    new <- list(...)
    args[names(new)] <- new
    do.call(FFBS_predict_exact, args)
  }
  expect_error(call(phi_para = -1), "Invalid prediction inputs")
  expect_error(call(nT = 3), "Invalid prediction inputs")
  expect_error(call(input_new = matrix(0.5, 1, 3)), "Invalid prediction inputs")
  expect_error(call(theta = NULL), "Supply smoothed coefficients")
  expect_error(call(theta = th[1]), "nT entries")
  expect_error(call(residual_weights = list(1)), "nT entries")
  expect_error(call(V_chol = diag(2)), "V_chol has incompatible dimensions")
  expect_error(call(F_ls = pp$F_ls[c(1, 1, 1)]), "one or nT entries")
  expect_error(call(F_new_ls = lapply(pp$F_new, cbind, 1)), "Incompatible design")
  expect_error(call(residual_weights = list(diag(2), diag(2))),
               "Residual weights have incompatible dimensions")
})
