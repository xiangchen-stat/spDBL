# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

waic_fixture <- function() {
  list(
    X = cbind(1, c(-0.5, 0.2, 0.7, 1.2)),
    Y = matrix(c(0.2, 0.7, -0.1, 1.1, -0.3, 0.4, 0.8, 0.6), 4, 2),
    Mn = matrix(c(0.18, 0.04, 0.04, 0.12), 2),
    mn = matrix(c(0.3, 0.1, -0.1, 0.4), 2),
    Psi_n = matrix(c(4, 1.2, 1.2, 3), 2),
    nu_n = 16)
}

test_that("analytical row densities agree with independent matrix-t densities", {
  f <- waic_fixture()
  rownames(f$Y) <- paste0("row", seq_len(nrow(f$Y)))
  score <- do.call(compute_analytical_waic, f)
  h <- diag(f$X %*% f$Mn %*% t(f$X))
  reference <- vapply(seq_len(nrow(f$Y)), function(i)
    mniw::dMT(f$Y[i, , drop = FALSE],
      Lambda = f$X[i, , drop = FALSE] %*% f$mn,
      SigmaR = matrix(1 + h[i]), SigmaC = f$Psi_n,
      nu = f$nu_n - ncol(f$Y) + 1, log = TRUE), numeric(1))
  expect_equal(unname(score$pointwise[, "lppd"]), reference, tolerance = 1e-12)
  expect_identical(rownames(score$pointwise), rownames(f$Y))
  expect_equal(score$lppd, sum(reference))
  expect_equal(score$WAIC, -2 * (score$lppd - score$p_waic))
  expect_equal(unname(colSums(score$pointwise)),
    unname(unlist(score[c("lppd", "p_waic", "elpd_waic", "WAIC")])))
})

test_that("analytical WAIC supports a univariate inverse-gamma posterior", {
  X <- matrix(c(0, 0.7, 1.2), 3, 1)
  Y <- matrix(c(0.2, -0.1, 0.9), 3, 1)
  score <- compute_analytical_waic(Y, X, matrix(0.2), matrix(0.4),
    matrix(3), 10)
  h <- as.numeric(X)^2 * 0.2
  d <- (as.numeric(Y) - as.numeric(X) * 0.4)^2 / 3
  scale <- sqrt((1 + h) * 3 / 10)
  reference <- stats::dt((as.numeric(Y) - as.numeric(X) * 0.4) / scale,
    df = 10, log = TRUE) - log(scale)
  expect_equal(score$pointwise[, "lppd"], reference, tolerance = 1e-12)

  # For q=1 and a zero design row, integrate over precision directly.
  # tau = 1 / sigma^2 follows Gamma(nu/2, rate=Psi/2).
  ll <- function(tau) 0.5 * (log(tau) - log(2 * pi) - 0.2^2 * tau)
  integrand <- function(tau, power)
    ll(tau)^power * stats::dgamma(tau, shape = 5, rate = 1.5)
  first <- stats::integrate(integrand, 0, Inf, power = 1)$value
  second <- stats::integrate(integrand, 0, Inf, power = 2)$value
  expect_equal(unname(score$pointwise[1, "p_waic"]), second - first^2,
    tolerance = 1e-8)
})

test_that("analytical scores agree with independent MNIW posterior draws", {
  f <- waic_fixture()
  set.seed(163)
  draws <- mniw::rMNIW(30000L, Lambda = f$mn, Sigma = f$Mn,
    Psi = f$Psi_n, nu = f$nu_n)
  before <- .Random.seed
  analytical <- do.call(compute_analytical_waic, f)
  sampled <- compute_sample_waic(f$Y, f$X, draws$X, draws$V)
  expect_identical(before, .Random.seed)
  expect_lt(max(abs(analytical$pointwise[, "lppd"] -
    sampled$pointwise[, "lppd"])), 0.025)
  expect_lt(max(abs(analytical$pointwise[, "p_waic"] -
    sampled$pointwise[, "p_waic"])), 0.05)
})

test_that("sample log likelihoods agree with direct normal densities", {
  f <- waic_fixture()
  set.seed(72)
  draws <- mniw::rMNIW(12L, Lambda = f$mn, Sigma = f$Mn,
    Psi = f$Psi_n, nu = f$nu_n)
  score <- compute_sample_waic(f$Y, f$X, draws$X, draws$V,
    return_log_lik = TRUE)
  reference <- vapply(seq_len(nrow(f$Y)), function(i)
    vapply(seq_len(12L), function(s)
      mniw::dMNorm(f$Y[i, , drop = FALSE],
        Lambda = f$X[i, , drop = FALSE] %*% draws$X[, , s],
        SigmaR = matrix(1), SigmaC = draws$V[, , s], log = TRUE),
      numeric(1)), numeric(12))
  expect_equal(score$log_lik, reference, tolerance = 1e-12)
  expect_equal(score$pointwise[, "p_waic"], apply(reference, 2, stats::var))
  expect_null(compute_sample_waic(f$Y, f$X, draws$X, draws$V)$log_lik)
})

test_that("sample scores match loo aggregation on identical likelihood draws", {
  skip_if_not_installed("loo")
  f <- waic_fixture()
  set.seed(73)
  draws <- mniw::rMNIW(100L, Lambda = f$mn, Sigma = f$Mn,
    Psi = f$Psi_n, nu = f$nu_n)
  score <- compute_sample_waic(f$Y, f$X, draws$X, draws$V,
    return_log_lik = TRUE)
  # loo diagnoses large row penalties; this test checks the aggregation only.
  reference <- suppressWarnings(loo::waic(score$log_lik))
  expect_equal(score$WAIC, unname(reference$estimates["waic", "Estimate"]),
    tolerance = 1e-12)
  expect_equal(score$p_waic, unname(reference$estimates["p_waic", "Estimate"]),
    tolerance = 1e-12)
  expect_equal(score$pointwise[, "elpd_waic"],
    unname(reference$pointwise[, "elpd_waic"]), tolerance = 1e-12)
})

test_that("singleton axes and very small likelihoods retain stable scores", {
  X <- Y <- matrix(1, 1, 1)
  B <- array(c(0, 0.2, 0.4), c(1, 1, 3))
  Sigma <- array(c(0.5, 1, 1.5), c(1, 1, 3))
  score <- compute_sample_waic(Y, X, B, Sigma, return_log_lik = TRUE)
  expect_identical(dim(score$pointwise), c(1L, 4L))
  expect_identical(dim(score$log_lik), c(3L, 1L))
  expect_equal(as.numeric(score$log_lik),
    stats::dnorm(1, c(0, 0.2, 0.4), sqrt(c(0.5, 1, 1.5)), log = TRUE))
  analytical <- compute_analytical_waic(Y, X, matrix(0.2), matrix(0.4),
    matrix(3), 10)
  expect_identical(dim(analytical$pointwise), c(1L, 4L))
  extreme <- compute_sample_waic(matrix(1000), X, B, Sigma)
  expect_true(all(is.finite(unlist(extreme))))

  # One predictor and three responses still require matrix-valued slices.
  multivariate <- compute_sample_waic(matrix(c(0, 0.1, 0.2), 1), X,
    array(0, c(1, 3, 3)), array(rep(diag(3), 3), c(3, 3, 3)),
    return_log_lik = TRUE)
  expect_identical(dim(multivariate$log_lik), c(3L, 1L))
  expect_equal(multivariate$p_waic, 0)
})

test_that("invalid matrices, posterior parameters and draws fail clearly", {
  f <- waic_fixture()
  check <- function(name, value) {
    args <- f
    args[name] <- list(value)
    do.call(compute_analytical_waic, args)
  }
  expect_error(check("Y", matrix(NA_real_, 4, 2)), "Y must")
  expect_error(check("Y", matrix(character(), 0, 2)), "Y must")
  expect_error(check("X", f$X[-1, , drop = FALSE]), "same number of rows")
  expect_error(check("Mn", -diag(2)), "Mn must be symmetric positive definite")
  expect_error(check("Mn", matrix(c(1, 0.1, 0.3, 1), 2)), "symmetric")
  expect_error(check("mn", matrix(0, 1, 2)), "mn must have dimensions")
  expect_error(check("Psi_n", matrix(1, 2, 2)), "Psi_n must be symmetric")
  for (nu in list(1, NA_real_, Inf, c(2, 3)))
    expect_error(check("nu_n", nu), "nu_n must")
  expect_true(is.finite(check("nu_n", 1.5)$WAIC))

  B <- array(0, c(2, 2, 3))
  Sigma <- array(rep(diag(2), 3), c(2, 2, 3))
  expect_error(compute_sample_waic(f$Y, f$X, B[, , 1], Sigma), "B_samples")
  expect_error(compute_sample_waic(f$Y, f$X, B[, , 1, drop = FALSE],
    Sigma[, , 1, drop = FALSE]), "S >= 2")
  expect_error(compute_sample_waic(f$Y, f$X, B, Sigma[, , 1]), "Sigma_samples")
  expect_error(compute_sample_waic(f$Y, f$X, B, Sigma, NA), "return_log_lik")
  Sigma[, , 2] <- -diag(2)
  expect_error(compute_sample_waic(f$Y, f$X, B, Sigma), "positive definite")
})
