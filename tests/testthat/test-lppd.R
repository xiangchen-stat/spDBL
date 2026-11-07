# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

lmvgamma <- function(a, q) q * (q - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(q)) / 2))

# Matrix-t: X | Sigma ~ MN(M, A, Sigma), Sigma ~ IW(nu, B) (Gupta and Nagar).
log_dmatt_iw <- function(X, M, A, nu, B) {
  n <- nrow(X)
  q <- ncol(X)
  Q <- B + t(X - M) %*% solve(A, X - M)
  logdet <- function(x) as.numeric(determinant(x)$modulus)
  -n * q / 2 * log(pi) + lmvgamma((nu + n) / 2, q) - lmvgamma(nu / 2, q) -
    q / 2 * logdet(A) + nu / 2 * logdet(B) - (nu + n) / 2 * logdet(Q)
}

# Marginal of X | s2 ~ MN(M, A, s2 R), s2 ~ IG(shape, rate) by quadrature.
log_dmatt_ig_numeric <- function(X, M, A, R, shape, rate) {
  integrand <- function(s2) vapply(s2, function(s) exp(
    log_dmatnorm(X, M, A, s * R) + shape * log(rate) - lgamma(shape) -
      (shape + 1) * log(s) - rate / s), numeric(1))
  log(stats::integrate(integrand, 0, Inf, rel.tol = 1e-10)$value)
}

test_that("dMTig is the inverse-gamma mixture of matrix normals", {
  set.seed(51)
  p <- 3; S <- 2
  Y <- matrix(rnorm(p * S), p, S)
  m <- matrix(rnorm(p * S, sd = 0.3), p, S)
  M <- rand_spd(p)
  R <- rand_spd(S)
  expected <- log_dmatt_ig_numeric(Y, m, M, R, shape = 4, rate = 2.5)
  expect_equal(dMTig(Y, m, M, nu = 4, d = 2.5, R = R), expected, tolerance = 1e-7)
  expect_equal(dMTig(Y, m, M, nu = 4, d = 2.5, R = R, log = FALSE), exp(expected),
               tolerance = 1e-7)
})

test_that("lppd_IG_1t integrates the scalar variance analytically", {
  set.seed(52)
  N <- 3; p <- 2; S <- 2
  Yt <- matrix(rnorm(N * S), N, S)
  Ft <- matrix(rnorm(N * p), N, p)
  Vt <- rand_spd(N)
  st <- matrix(rnorm(p * S), p, S)
  St <- rand_spd(p)
  R <- rand_spd(S)
  expected <- log_dmatt_ig_numeric(Yt, Ft %*% st, Ft %*% St %*% t(Ft) + Vt, R,
                                   shape = 5, rate = 1.5)
  expect_equal(lppd_IG_1t(Yt, Ft, Vt, st, St, nt = 5, Dt = 1.5, R = R), expected,
               tolerance = 1e-7)
})

test_that("lppd_IW_1t is the matrix-t predictive of the MNIW model", {
  set.seed(53)
  N <- 4; p <- 2; q <- 3; nt <- 8
  Yt <- matrix(rnorm(N * q), N, q)
  Ft <- matrix(rnorm(N * p), N, p)
  Vt <- rand_spd(N)
  st <- matrix(rnorm(p * q), p, q)
  St <- rand_spd(p)
  Dt <- rand_spd(q)
  expected <- log_dmatt_iw(Yt, Ft %*% st, Ft %*% St %*% t(Ft) + Vt, nt, Dt)
  expect_equal(lppd_IW_1t(Yt, Ft, Vt, st, St, nt, Dt), expected, tolerance = 1e-8)
  # Monte Carlo over Sigma ~ IW(nt, Dt) agrees with the closed form.
  set.seed(54)
  Sig <- matrixsampling::rinvwishart(4000, nt, Dt)
  dens <- vapply(seq_len(4000), function(i) exp(log_dmatnorm(Yt, Ft %*% st,
    Ft %*% St %*% t(Ft) + Vt, Sig[, , i])), numeric(1))
  expect_equal(log(mean(dens)), expected, tolerance = 0.05)
})

test_that("lppd_id_1t is a matrix normal with identity column covariance", {
  set.seed(55)
  N <- 3; p <- 2; S <- 2
  Yt <- matrix(rnorm(N * S), N, S)
  Ft <- matrix(rnorm(N * p), N, p)
  Vt <- rand_spd(N)
  st <- matrix(rnorm(p * S), p, S)
  St <- rand_spd(p)
  A <- Ft %*% St %*% t(Ft) + Vt
  expect_equal(lppd_id_1t(Yt, Ft, Vt, st, St),
               log_dmatnorm(Yt, Ft %*% st, A, diag(S)), tolerance = 1e-10)
})
