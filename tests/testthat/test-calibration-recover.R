# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Joint Gaussian moments of (u_0, ..., u_T, z_1, ..., z_T) for the random-walk
# discrepancy u_t = u_{t-1} + N(0, tau2_t U), u_0 ~ N(m0, tau2_0 M0), and
# z_t = mu_t + y_t + u_t + N(0, tau2_t I) with y_t ~ N(0, Sigma) independent.
discrepancy_joint <- function(pred, tau2_all, U, m0, M0) {
  T <- nrow(pred$mu)
  S <- ncol(pred$mu)
  cov_u <- function(i, j) tau2_all[1] * M0 +
    if (min(i, j) > 0) sum(tau2_all[seq_len(min(i, j)) + 1]) * U else 0
  n <- (2 * T + 1) * S
  K <- matrix(0, n, n)
  blk <- function(k) (k - 1) * S + seq_len(S)
  # Blocks 1..T+1 are u_0..u_T; blocks T+2..2T+1 are z_1..z_T.
  time_of <- function(k) if (k <= T + 1) k - 1 else k - T - 1
  for (a in seq_len(2 * T + 1)) for (b in seq_len(2 * T + 1)) {
    K[blk(a), blk(b)] <- cov_u(time_of(a), time_of(b))
    if (a == b && a > T + 1)
      K[blk(a), blk(b)] <- K[blk(a), blk(b)] + pred$Sigma +
        design_cov(pred, time_of(a)) + diag(tau2_all[time_of(a) + 1], S)
  }
  mean <- c(rep(m0, T + 1), as.numeric(t(pred$mu)) + rep(m0, T))
  list(K = K, mean = mean, u = seq_len((T + 1) * S), z = (T + 1) * S + seq_len(T * S))
}

test_that("Stan data carry the emulator draw, residual weights and factors", {
  f <- calibration_fixture()
  e <- f$emulator
  data <- spDBL:::make_stan_calibration_data(e, 2L, f$z, f$U)
  V <- exp_kernel_loop(e$input, e$input, e$phi_para, e$gp_sigma2) + diag(e$gp_tau2, 8)
  for (k in seq_len(e$E)) {
    th <- e$posterior$theta[[k]][, , 2]
    expect_equal(data$theta[k, , ], th)
    expect_equal(data$residual_weights[k, , ], solve(V, e$Y[[k]] - e$F_ls[[k]] %*% th),
                 tolerance = 1e-8)
  }
  expect_equal(data$L_V %*% t(data$L_V), V, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(data$L_U %*% t(data$L_U), f$U, tolerance = 1e-10)
  expect_equal(data$sigma2, e$posterior$sigma2[2])
  expect_equal(c(data$N, data$T, data$S, data$E, data$P), c(8, 2, 4, 4, 4))
  sub <- spDBL:::make_stan_calibration_data(e, 1L, f$z[, c(1, 3)], f$U,
    locations = c(4, 2), m0 = c(1, 2), M0 = diag(2) * 2, discount = 0.9)
  expect_equal(sub$R_field, e$R_field[c(4, 2), c(4, 2)])
  expect_equal(sub$L_U %*% t(sub$L_U), f$U[c(4, 2), c(4, 2)], tolerance = 1e-10)
  expect_equal(sub$L_M0 %*% t(sub$L_M0), diag(2) * 2)
  expect_identical(as.integer(sub$locations), c(4L, 2L))
  bad <- function(...) spDBL:::make_stan_calibration_data(e, ...)
  expect_error(bad(1L, f$z[, 1:3], f$U), "Invalid calibration observations")
  expect_error(bad(3L, f$z, f$U), "Invalid calibration observations")
  expect_error(bad(1L, f$z, f$U, discount = 1.5), "Invalid calibration observations")
  expect_error(bad(1L, f$z, f$U, M0 = diag(3)), "Invalid calibration observations")
  expect_error(bad(1L, f$z, f$U[1:3, 1:3]), "Invalid calibration observations")
})

test_that("residual weights honour precomputed values and the draw cache", {
  f <- calibration_fixture()
  e <- f$emulator
  w2 <- spDBL:::.cal_residual_weights(e, 2L)
  expect_identical(e$residual_cache$draw, 2L)
  expect_identical(spDBL:::.cal_residual_weights(e, 2L), w2)
  uncached <- e
  uncached$residual_cache <- NULL
  expect_equal(spDBL:::.cal_residual_weights(uncached, 2L), w2)
  supplied <- e
  supplied$residual_weights <- list("first", "second")
  expect_identical(spDBL:::.cal_residual_weights(supplied, 2L), "second")
})

test_that("the noncentred log density is the prior times the Gaussian likelihood", {
  f <- calibration_fixture()
  e <- f$emulator
  loc <- c(3, 1)
  m0 <- c(0.1, -0.2)
  M0 <- matrix(c(1, 0.3, 0.3, 2), 2)
  precision <- c(3, 5, 4)
  u_std <- matrix(c(0.3, -0.5, 1.1, 0.2, -0.7, 0.4), 3)
  eta_unit <- c(0.25, 0.7)
  got <- spDBL:::calibration_log_density(eta_unit, precision, u_std, e, 2L,
    f$z[, loc], f$U, locations = loc, n0 = 10, d0 = 3, discount = 0.9, m0 = m0, M0 = M0)
  pred <- spDBL:::predict_calibration(eta_unit, 2L, e, loc)
  tau2 <- 1 / precision
  u0 <- m0 + sqrt(tau2[1]) * t(chol(M0)) %*% u_std[1, ]
  u1 <- u0 + sqrt(tau2[2]) * t(chol(f$U[loc, loc])) %*% u_std[2, ]
  u2 <- u1 + sqrt(tau2[3]) * t(chol(f$U[loc, loc])) %*% u_std[3, ]
  lik <- function(t, u) log_dmatnorm(matrix(f$z[t, loc]), matrix(pred$mu[t, ] + u),
                                     pred$Sigma + design_cov(pred, t) +
                                       diag(tau2[t + 1], 2), matrix(1))
  expected <- sum(dgamma(precision, 10 * 0.9^(0:2), 3 * 0.9^(0:2), log = TRUE)) +
    sum(dnorm(u_std, log = TRUE)) + lik(1, u1) + lik(2, u2)
  expect_equal(got, expected, tolerance = 1e-10)
})

test_that("the marginal density integrates the discrepancy out exactly", {
  f <- calibration_fixture()
  e <- f$emulator
  loc <- c(2, 4, 1)
  m0 <- c(0.2, 0, -0.1)
  M0 <- diag(c(1, 2, 0.5))
  precision <- c(2, 6, 3)
  got <- spDBL:::calibration_marginal_log_density(c(0.6, 0.2), precision, e, 1L,
    f$z[, loc], f$U, locations = loc, n0 = 12, d0 = 2, discount = 0.95, m0 = m0, M0 = M0)
  pred <- spDBL:::predict_calibration(c(0.6, 0.2), 1L, e, loc)
  joint <- discrepancy_joint(pred, 1 / precision, f$U[loc, loc], m0, M0)
  zvec <- as.numeric(t(f$z[, loc]))
  expected <- sum(dgamma(precision, 12 * 0.95^(0:2), 2 * 0.95^(0:2), log = TRUE)) +
    log_dmatnorm(matrix(zvec), matrix(joint$mean[joint$z]), joint$K[joint$z, joint$z], matrix(1))
  expect_equal(got, expected, tolerance = 1e-10)
  # Averaging the conditional likelihood over prior discrepancy draws agrees.
  set.seed(101)
  n <- 2000
  prior <- sum(dgamma(precision, 12 * 0.95^(0:2), 2 * 0.95^(0:2), log = TRUE))
  log_lik <- vapply(seq_len(n), function(i) {
    u_std <- matrix(rnorm(9), 3)
    spDBL:::calibration_log_density(c(0.6, 0.2), precision, u_std, e, 1L, f$z[, loc],
      f$U, locations = loc, n0 = 12, d0 = 2, discount = 0.95, m0 = m0, M0 = M0) -
      prior - sum(dnorm(u_std, log = TRUE))
  }, numeric(1))
  w <- exp(log_lik - max(log_lik))
  mc <- max(log_lik) + log(mean(w))
  se <- stats::sd(w) / (sqrt(n) * mean(w))  # delta-method s.e. of log(mean(w))
  expect_lt(abs(mc - (got - prior)), 4 * se)
})

test_that("backward discrepancy draws follow the exact Gaussian posterior", {
  f <- calibration_fixture()
  e <- f$emulator
  n <- 2000
  tau2_0 <- 0.5
  tau2 <- c(0.3, 0.2)
  m0 <- c(0.1, 0, -0.1, 0.2)
  set.seed(102)
  before <- .Random.seed
  out <- spDBL:::recover_calibration_discrepancy(matrix(f$eta, n, 2, byrow = TRUE),
    matrix(tau2, n, 2, byrow = TRUE), rep(tau2_0, n), rep(2L, n), e, f$z, f$U,
    m0 = m0, seed = 7L)
  expect_identical(.Random.seed, before)
  expect_identical(dim(out$u), c(as.integer(n), 2L, 4L))
  expect_identical(dim(out$u0), c(as.integer(n), 4L))
  pred <- spDBL:::predict_calibration(f$eta, 2L, e)
  j <- discrepancy_joint(pred, c(tau2_0, tau2), f$U, m0, diag(4))
  Kzz <- j$K[j$z, j$z]
  Kuz <- j$K[j$u, j$z]
  post_mean <- j$mean[j$u] + Kuz %*% solve(Kzz, as.numeric(t(f$z)) - j$mean[j$z])
  post_cov <- j$K[j$u, j$u] - Kuz %*% solve(Kzz, t(Kuz))
  draws <- cbind(out$u0, out$u[, 1, ], out$u[, 2, ])
  se <- sqrt(diag(post_cov) / n)
  expect_lt(max(abs(colMeans(draws) - post_mean) / se), 4.5)
  expect_equal(stats::cov(draws), post_cov, tolerance = 0.1)
})

test_that("predictive field draws condition on the discrepancy and noise", {
  f <- calibration_fixture()
  e <- f$emulator
  n <- 2000
  loc <- c(1, 3, 4)
  tau2 <- c(0.15, 0.25)
  u <- array(rep(c(0.1, -0.2, 0.05), each = n * 2), c(n, 2, 3))
  out <- spDBL:::recover_calibration_predictions(matrix(f$eta, n, 2, byrow = TRUE),
    matrix(tau2, n, 2, byrow = TRUE), u, rep(1L, n), e, f$z[, loc], loc, seed = 8L)
  pred <- spDBL:::predict_calibration(f$eta, 1L, e, loc)
  for (t in 1:2) {
    prior <- pred$Sigma + design_cov(pred, t)
    K <- prior %*% solve(prior + diag(tau2[t], 3))
    mean <- pred$mu[t, ] + K %*% (f$z[t, loc] - u[1, t, ] - pred$mu[t, ])
    cov <- prior - K %*% prior
    y <- out$y[, t, ]
    expect_lt(max(abs(colMeans(y) - mean) / sqrt(diag(cov) / n)), 4.5)
    expect_equal(stats::cov(y), cov, tolerance = 0.1, ignore_attr = TRUE)
    noise <- out$z_rep[, t, ] - y - u[, t, ]
    expect_equal(apply(noise, 2, stats::var), rep(tau2[t], 3), tolerance = 0.1)
  }
})

test_that("Stan output columns are reshaped by their indices", {
  cols <- c("u[2,1]", "lp__", "u[1,1]", "u[1,2]", "u[2,2]", "u[1,3]", "u[2,3]")
  samples <- matrix(0, 4, length(cols), dimnames = list(NULL, cols))
  for (k in seq_along(cols)) if (cols[k] != "lp__") {
    idx <- as.integer(strsplit(gsub("u\\[|\\]", "", cols[k]), ",")[[1]])
    samples[, k] <- 10 * idx[1] + idx[2] + (1:4) / 10
  }
  arr <- spDBL:::.cal_extract_array(samples, "u", c(2, 3))
  expect_identical(dim(arr), c(4L, 2L, 3L))
  for (i in 1:2) for (j in 1:3) expect_equal(arr[, i, j], 10 * i + j + (1:4) / 10)
})
