# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

test_that("FFBS predictions match direct Gaussian conditioning", {
  f <- calibration_fixture()
  e <- f$emulator
  V <- 1.1 * exp(-2 * as.matrix(dist(e$input))) + diag(1e-4, nrow(e$input))
  for (eta in list(f$eta, e$input[1, ])) for (m in 1:2) {
    p <- spDBL:::predict_calibration(eta, m, e)
    F <- spDBL:::predict_offline_ar2(e$predictor, eta)$F
    J <- 1.1 * exp(-2 * sqrt(rowSums(sweep(e$input, 2, eta, "-")^2)))
    expected <- matrix(0, e$T, e$S)
    for (episode in seq_len(e$E)) {
      t <- (episode - 1L) %/% e$B + 1L
      b <- (episode - 1L) %% e$B + 1L
      theta <- e$posterior$theta[[episode]][, , m]
      expected[t, e$block_cols[b, ]] <- F[[episode]] %*% theta +
        crossprod(J, solve(V, e$Y[[episode]] - e$F_ls[[episode]] %*% theta))
    }
    expect_equal(p$mu, expected, tolerance = 1e-10)
    expect_equal(p$Sigma,
      e$posterior$sigma2[m] * (1.1001 - drop(crossprod(J, solve(V, J)))) * e$R_field,
      tolerance = 1e-10)
    expect_true(all(eigen(p$Sigma, symmetric = TRUE)$values > 0))
  }
  p <- spDBL:::predict_calibration(f$eta, 1L, e)
  sites <- c(4L, 1L, 3L)
  subset <- spDBL:::predict_calibration(f$eta, 1L, e, sites)
  expect_equal(subset$mu, p$mu[, sites, drop = FALSE])
  expect_equal(subset$Sigma, p$Sigma[sites, sites])
  for (t in seq_len(e$T))
    expect_equal(design_cov(subset, t), design_cov(p, t)[sites, sites])
  means_only <- spDBL:::predict_calibration(f$eta, 1L, e, covariance = "none")
  expect_identical(means_only, list(mu = p$mu))
  pointwise <- spDBL:::predict_calibration(f$eta, 1L, e, sites, covariance = "marginal")
  expect_identical(pointwise$mu, subset$mu)
  expect_equal(pointwise$variance, t(vapply(seq_len(e$T), function(t)
    diag(spDBL:::.cal_time_covariance(subset, t)), numeric(3))), ignore_attr = TRUE)
  mapped <- e
  mapped$block_cols <- matrix(c(4L, 2L, 1L, 3L), 2, byrow = TRUE)
  remapped <- spDBL:::predict_calibration(f$eta, 1L, mapped)
  expect_equal(remapped$mu[, as.vector(t(mapped$block_cols))], p$mu)
  before <- .Random.seed
  expect_identical(spDBL:::predict_calibration(f$eta, 1L, e), p)
  expect_identical(before, .Random.seed)
})

test_that("design-GP covariance enters through each coefficient draw", {
  f <- calibration_fixture()
  e <- f$emulator
  # Independent calculation of the shared-kernel GP for the design columns:
  # unit-amplitude SE kernel on unit inputs, maximum-likelihood column
  # covariance Fc' K^{-1} Fc / N.
  x <- e$input  # the fixture's bounds are [0, 1], so inputs are already unit
  N <- nrow(x)
  ell <- e$predictor$lengthscale
  kern <- function(a, b) exp(-0.5 * outer(seq_len(nrow(a)), seq_len(nrow(b)),
    Vectorize(function(i, j) sum(((a[i, ] - b[j, ]) / ell)^2))))
  K <- kern(x, x) + diag(e$predictor$nugget, N)
  Fc <- lapply(e$F_ls, function(F) sweep(F, 2, colMeans(F)))
  column_cov <- function(i, j) crossprod(Fc[[i]], solve(K, Fc[[j]])) / N
  for (eta in list(f$eta, c(0.02, 0.97), e$input[5, ])) {
    k <- kern(x, matrix(eta, 1))
    v <- max(0, 1 - drop(crossprod(k, solve(K, k))))
    design <- spDBL:::predict_offline_ar2(e$predictor, eta)
    expect_equal(design$variance_factor, v, tolerance = 1e-8)
    expect_equal(design$variance, v * unlist(lapply(seq_len(e$E), function(i)
      diag(column_cov(i, i)))), tolerance = 1e-8)
    for (m in 1:2) {
      p <- spDBL:::predict_calibration(eta, m, e)
      expect_length(p$design_factors, e$T)
      for (t in seq_len(e$T)) {
        expected <- matrix(0, e$S, e$S)
        episodes <- (t - 1L) * e$B + seq_len(e$B)
        for (i in episodes) for (j in episodes) {
          bi <- e$block_cols[(i - 1L) %% e$B + 1L, ]
          bj <- e$block_cols[(j - 1L) %% e$B + 1L, ]
          expected[bi, bj] <- v * t(e$posterior$theta[[i]][, , m]) %*%
            column_cov(i, j) %*% e$posterior$theta[[j]][, , m]
        }
        expect_equal(design_cov(p, t), expected, tolerance = 1e-8, ignore_attr = TRUE)
      }
    }
  }
  # Column 1 of every design is the constant 1, which has no GP error.
  design <- spDBL:::predict_offline_ar2(e$predictor, f$eta)
  expect_equal(design$variance[seq(1, 4 * e$E, by = 4)], rep(0, e$E))
  # The variance is smaller at a training input than far from the design.
  near <- spDBL:::predict_offline_ar2(e$predictor, e$input[5, ])$variance_factor
  far <- spDBL:::predict_offline_ar2(e$predictor, c(0.02, 0.97))$variance_factor
  expect_lt(near, far)
  # Monte Carlo: the field mean is linear in the design rows, so drawing all
  # design rows jointly from the GP predictive gives the stated covariance.
  set.seed(83)
  L <- t(chol(solve(K)))
  Fc_all <- do.call(cbind, Fc)
  draws <- t(replicate(4000, {
    err <- sqrt(design$variance_factor / N) * drop(crossprod(Fc_all, L %*% rnorm(N)))
    F_new <- lapply(seq_len(e$E), function(i)
      matrix(design$F[[i]] + err[(i - 1) * 4 + 1:4], 1))
    th <- lapply(e$posterior$theta, function(a) a[, , 1])
    mean <- FFBS_predict_exact(Y = e$Y, F_ls = e$F_ls, F_new_ls = F_new,
      input = e$input, input_new = matrix(f$eta, 1), nT = e$E, phi_para = 2,
      theta = th)
    as.numeric(recover_from_EP_exact(mean, e$T, e$E, block_cols = e$block_cols)[[2]])
  }))
  p <- spDBL:::predict_calibration(f$eta, 1L, e)
  expect_equal(stats::cov(draws), design_cov(p, 2L), tolerance = 0.1,
               ignore_attr = TRUE)
})

test_that("the Matern 5/2 kernel matches its closed form", {
  k <- function(r) spDBL:::.cal_kernel(matrix(0, 1, 2), matrix(c(r, 0), 1),
                                       c(1, 1), "matern52")
  expect_equal(drop(k(0)), 1)
  for (r in c(0.1, 1 / sqrt(5), 1, 2.5)) {
    a <- sqrt(5) * r
    expect_equal(drop(k(r)), (1 + a + a^2 / 3) * exp(-a))
  }
  # Anisotropic: each coordinate is scaled by its own lengthscale.
  x <- matrix(c(0.2, 0.9, 0.5, 0.1), 2)
  y <- matrix(c(0.6, 0.3), 1)
  ell <- c(0.5, 2)
  r <- sqrt(rowSums(sweep(sweep(x, 2, y), 2, ell, "/")^2))
  expect_equal(drop(spDBL:::.cal_kernel(x, y, ell, "matern52")),
               (1 + sqrt(5) * r + 5 * r^2 / 3) * exp(-sqrt(5) * r))
  expect_equal(spDBL:::.cal_kernel(x, y, ell, "sqexp"), spDBL:::.cal_sqexp(x, y, ell))
})

test_that("a design GP per time step propagates each time's own variance", {
  f <- calibration_fixture(design_kernel = "matern52", design_groups = "time")
  e <- f$emulator
  expect_identical(e$design_kernel, "matern52")
  expect_identical(e$design_groups, "time")
  expect_identical(e$predictor$groups, c(1L, 1L, 2L, 2L))
  expect_identical(dim(e$predictor$lengthscale), c(2L, 2L))
  expect_false(isTRUE(all.equal(e$predictor$lengthscale[1, ],
                                e$predictor$lengthscale[2, ])))
  # Independent calculation, time by time: Matern 5/2 kernel on unit inputs
  # with that time's hyperparameters, maximum-likelihood column covariance
  # Fc' K_t^{-1} Fc / N over the time's episodes.
  x <- e$input
  N <- nrow(x)
  matern <- function(a, b, ell) {
    r <- sqrt(outer(seq_len(nrow(a)), seq_len(nrow(b)), Vectorize(function(i, j)
      sum(((a[i, ] - b[j, ]) / ell)^2))))
    (1 + sqrt(5) * r + 5 * r^2 / 3) * exp(-sqrt(5) * r)
  }
  Fc <- lapply(e$F_ls, function(F) sweep(F, 2, colMeans(F)))
  for (eta in list(f$eta, c(0.02, 0.97), e$input[5, ])) {
    design <- spDBL:::predict_offline_ar2(e$predictor, eta)
    v <- numeric(e$T)
    for (t in seq_len(e$T)) {
      ell <- e$predictor$lengthscale[t, ]
      K <- matern(x, x, ell) + diag(e$predictor$nugget[t], N)
      k <- matern(x, matrix(eta, 1), ell)
      v[t] <- max(0, 1 - drop(crossprod(k, solve(K, k))))
      episodes <- (t - 1L) * e$B + seq_len(e$B)
      for (i in episodes) {
        expect_equal(design$F[[i]], colMeans(e$F_ls[[i]]) +
          drop(crossprod(k, solve(K, Fc[[i]]))), tolerance = 1e-8)
        expect_equal(design$variance[(i - 1L) * 4L + 1:4],
          v[t] * diag(crossprod(Fc[[i]], solve(K, Fc[[i]])) / N), tolerance = 1e-8)
      }
      for (m in 1:2) {
        p <- spDBL:::predict_calibration(eta, m, e)
        expected <- matrix(0, e$S, e$S)
        for (i in episodes) for (j in episodes) {
          bi <- e$block_cols[(i - 1L) %% e$B + 1L, ]
          bj <- e$block_cols[(j - 1L) %% e$B + 1L, ]
          expected[bi, bj] <- v[t] * t(e$posterior$theta[[i]][, , m]) %*%
            (crossprod(Fc[[i]], solve(K, Fc[[j]])) / N) %*% e$posterior$theta[[j]][, , m]
        }
        expect_equal(design_cov(p, t), expected, tolerance = 1e-8, ignore_attr = TRUE)
        pointwise <- spDBL:::predict_calibration(eta, m, e, covariance = "marginal")
        expect_equal(pointwise$variance[t, ], diag(spDBL:::.cal_time_covariance(p, t)),
                     tolerance = 1e-10, ignore_attr = TRUE)
      }
    }
    expect_equal(design$variance_factor, v, tolerance = 1e-8)
    expect_equal(design$episode_factor, rep(v, each = e$B), tolerance = 1e-8)
    expect_equal(spDBL:::predict_calibration(eta, 1L, e)$design_variance, v,
                 tolerance = 1e-8)
  }
})

test_that("design GP groups are fitted as separate GPs", {
  f <- calibration_fixture()
  e <- f$emulator
  grouped <- spDBL:::fit_offline_ar2(e$input, e$F_ls, c(0, 0), c(1, 1),
                                     kernel = "matern52", groups = c(1, 1, 2, 2))
  for (g in 1:2) {
    alone <- spDBL:::fit_offline_ar2(e$input, e$F_ls[2 * g - 1:0], c(0, 0), c(1, 1),
                                     kernel = "matern52")
    cols <- (g - 1L) * 8L + 1:8
    expect_equal(grouped$lengthscale[g, ], drop(alone$lengthscale))
    expect_equal(grouped$nugget[g], alone$nugget)
    expect_equal(grouped$alpha[, cols], alone$alpha)
    expect_equal(grouped$column_variance[cols], alone$column_variance)
    expect_equal(spDBL:::predict_offline_ar2(grouped, f$eta)$F[2 * g - 1:0],
                 spDBL:::predict_offline_ar2(alone, f$eta)$F)
  }
  # A group of constant designs has no GP error and keeps fixed settings.
  constant <- rep(list(matrix(1, nrow(e$input), 4)), 2)
  mixed <- spDBL:::fit_offline_ar2(e$input, c(constant, e$F_ls[3:4]), c(0, 0), c(1, 1),
                                   kernel = "matern52", groups = c(1, 1, 2, 2))
  expect_true(is.na(mixed$optimizer_convergence[1]))
  expect_equal(mixed$lengthscale[1, ], c(0.5, 0.5))
  design <- spDBL:::predict_offline_ar2(mixed, f$eta)
  expect_equal(unlist(design$F[1:2]), rep(1, 8))
  expect_equal(design$variance[1:8], rep(0, 8))
  for (bad in list(c(1, 1, 3, 3), c(1, 2), c(0, 1, 1, 1), c(1, 1.5, 2, 2)))
    expect_error(spDBL:::fit_offline_ar2(e$input, e$F_ls, c(0, 0), c(1, 1),
                                         groups = bad), "groups must number")
  expect_error(spDBL:::fit_offline_ar2(e$input, e$F_ls, c(0, 0), c(1, 1),
                                       kernel = "linear"), "should be one of")
})

test_that("design predictors saved before kernel groups still predict", {
  f <- calibration_fixture()
  e <- f$emulator
  old <- e$predictor
  old$kernel <- NULL
  old$groups <- NULL
  old$lengthscale <- drop(old$lengthscale)
  old$K_chol <- old$K_chol[[1]]
  old$episode_factor <- NULL
  expect_equal(spDBL:::predict_offline_ar2(old, f$eta),
               spDBL:::predict_offline_ar2(e$predictor, f$eta))
  expect_equal(spDBL:::.cal_design_whitened(old), spDBL:::.cal_design_whitened(e$predictor))
})

test_that("design uncertainty can be switched off", {
  f <- calibration_fixture()
  e <- f$emulator
  off <- e
  off$design_uncertainty <- FALSE
  on <- spDBL:::predict_calibration(f$eta, 2L, e)
  plain <- spDBL:::predict_calibration(f$eta, 2L, off)
  expect_null(plain$design_factors)
  expect_identical(plain[c("mu", "Sigma")], on[c("mu", "Sigma")])
  expect_identical(spDBL:::.cal_time_covariance(plain, 2L), plain$Sigma)
  expect_equal(spDBL:::.cal_time_covariance(on, 2L),
               on$Sigma + on$design_variance[2] * crossprod(on$design_factors[[2]]))
  expect_equal(spDBL:::predict_calibration(f$eta, 2L, off, covariance = "marginal")$variance,
               matrix(diag(plain$Sigma), e$T, e$S, byrow = TRUE), ignore_attr = TRUE)
  expect_identical(spDBL:::.cal_design_factors(off, 1L, c(2L, 3L)), list())
  expect_identical(dim(spDBL:::.cal_design_factors(e, 1L, c(4L, 1L))[[2]]), c(8L, 2L))
  # Sites 1 and 2 form block 1, so block 2 contributes nothing.
  full <- spDBL:::.cal_design_factors(e, 1L)
  only_first <- spDBL:::.cal_design_factors(e, 1L, c(1L, 2L))
  for (t in seq_len(e$T)) expect_equal(only_first[[t]], full[[t]][, 1:2])
  # Objects saved before the option existed behave as if it were off.
  legacy <- e
  legacy$design_uncertainty <- NULL
  expect_null(spDBL:::predict_calibration(f$eta, 2L, legacy)$design_factors)
  expect_error(spDBL:::make_offline_calibration_input(e$input, e$Y, e$F_ls,
    e$posterior, e$block_cols, e$R_field, 2, c(0, 0), c(1, 1),
    design_uncertainty = NA), "TRUE or FALSE")
  expect_null(spDBL:::make_offline_calibration_input(e$input, e$Y, e$F_ls,
    e$posterior, e$block_cols, e$R_field, 2, c(0, 0), c(1, 1),
    design_uncertainty = FALSE)$design_whitened)
  # Cached factors are recomputed when missing.
  uncached <- e
  uncached$design_whitened <- NULL
  expect_equal(spDBL:::predict_calibration(f$eta, 2L, uncached), on)
  old_predictor <- e$predictor
  old_predictor$K_chol <- NULL
  old_predictor$column_variance <- NULL
  expect_equal(spDBL:::predict_offline_ar2(old_predictor, f$eta),
               spDBL:::predict_offline_ar2(e$predictor, f$eta))
})

test_that("backward coefficient draws retain cross-time covariance", {
  ff <- list(
    list(mt = matrix(0), Mt = matrix(2/3), at = matrix(0), At = matrix(2), nt = 1000, Dt = 1000),
    list(mt = matrix(0), Mt = matrix(5/8), at = matrix(0), At = matrix(5/3), nt = 1000, Dt = 1000))
  post <- spDBL:::sample_calibration_emulator(list(ff = ff), matrix(1),
    matrix(1), 6000L, seed = 1234L, sampling = "joint")
  expect_equal(cov(as.numeric(post$theta[[1]]), as.numeric(post$theta[[2]])),
    0.25, tolerance = 0.03)
})

test_that("the Kalman likelihood equals full Gaussian conditioning", {
  f <- calibration_fixture()
  # Preserve the non-unit discrepancy scale used in the development study.
  f$U <- 1.1001 * f$U
  e <- f$emulator
  p <- spDBL:::predict_calibration(f$eta, 1L, e)
  tau <- c(0.2, 0.3, 0.4)
  n <- e$S
  block <- function(t) ((t - 1L) * n + 1L):(t * n)
  covariance <- matrix(0, n * e$T, n * e$T)
  for (i in seq_len(e$T)) for (j in seq_len(e$T)) {
    covariance[block(i), block(j)] <- tau[1] * diag(n) +
      sum(tau[seq_len(min(i, j)) + 1L]) * f$U
    if (i == j) covariance[block(i), block(j)] <-
      covariance[block(i), block(j)] + p$Sigma + design_cov(p, i) +
      diag(tau[i + 1L], n)
  }
  L <- chol(covariance)
  residual <- as.numeric(t(f$z - p$mu))
  expected <- -n * e$T / 2 * log(2 * pi) - sum(log(diag(L))) -
    sum(forwardsolve(t(L), residual)^2) / 2
  filtered <- spDBL:::.cal_filter_discrepancy(p, tau, f$z, f$U, rep(0, n), diag(n))
  expect_equal(filtered$log_likelihood, expected, tolerance = 1e-10)
})
