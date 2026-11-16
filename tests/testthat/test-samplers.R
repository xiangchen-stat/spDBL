# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Draws are checked against their exact first two moments with fixed seeds.
# Mean tolerances are ~4.5 Monte Carlo standard errors; covariance tolerances
# are relative (all.equal-style) and well above the sampling error at n = 4000.

vec_draws <- function(arr) t(matrix(arr, prod(dim(arr)[1:2])))

expect_mean_close <- function(draws, mean, cov, z = 4.5) {
  se <- sqrt(diag(cov) / nrow(draws))
  expect_lt(max(abs(colMeans(draws) - as.numeric(mean)) / se), z)
}

test_that("C++ and R matrix-normal generators have the requested moments", {
  set.seed(41)
  p <- 2; q <- 3; n <- 4000
  m <- matrix(rnorm(p * q), p, q)
  U <- rand_spd(p)
  V <- rand_spd(q)
  target <- kronecker(V, U)
  RU <- chol(U)
  RV <- chol(V)
  cpp <- vec_draws(simplify2array(lapply(seq_len(n), function(i) spDBL:::rmn_cpp(m, U, V))))
  chol_cpp <- vec_draws(simplify2array(lapply(seq_len(n), function(i)
    spDBL:::rmn_chol_cpp(m, RU, RV))))
  chol_r <- vec_draws(simplify2array(lapply(seq_len(n), function(i) rmn_chol(m, RU, RV))))
  more <- rmn_chol_more(n, m, RU, RV)
  expect_identical(dim(more), c(as.integer(p), as.integer(q), as.integer(n)))
  for (x in list(cpp, chol_cpp, chol_r, vec_draws(more))) {
    expect_mean_close(x, m, target)
    expect_equal(stats::cov(x), target, tolerance = 0.1)
  }
})

test_that("rmn_chol and rmn_chol_more use the same normal stream", {
  m <- matrix(0, 2, 2)
  RM <- chol(matrix(c(2, 0.5, 0.5, 1), 2))
  RS <- chol(matrix(c(1, 0.2, 0.2, 3), 2))
  set.seed(42)
  a <- rmn_chol(m, RM, RS)
  b <- rmn_chol(m, RM, RS)
  set.seed(42)
  both <- rmn_chol_more(2, m, RM, RS)
  expect_equal(both[, , 1], a)
  expect_equal(both[, , 2], b)
})

test_that("MNIW_sampler draws from the matrix-normal inverse-Wishart law", {
  set.seed(43)
  p <- 2; q <- 2; n <- 4000; v <- 9
  S <- rand_spd(q)
  C <- matrix(rnorm(p * q), p, q)
  Vb <- rand_spd(p)
  X <- matrix(rnorm(5 * p), 5, p)
  out <- MNIW_sampler(n, X, v, S, C, Vb)
  expect_named(out, c("sigma", "B"))
  expect_identical(dim(out$B), c(as.integer(p), as.integer(q), as.integer(n)))
  expect_identical(dim(out$sigma), c(as.integer(q), as.integer(q), as.integer(n)))
  E_sigma <- S / (v - q - 1)
  expect_equal(apply(out$sigma, 1:2, mean), E_sigma, tolerance = 0.06)
  B <- vec_draws(out$B)
  expect_mean_close(B, C, kronecker(E_sigma, Vb))
  expect_equal(stats::cov(B), kronecker(E_sigma, Vb), tolerance = 0.1)
})

test_that("MNIG_sampler draws from the matrix-normal inverse-gamma law", {
  set.seed(44)
  p <- 2; q <- 2; n <- 4000; v <- 7; rate <- 3
  C <- matrix(rnorm(p * q), p, q)
  Vb <- rand_spd(p)
  R <- rand_spd(q)
  out <- MNIG_sampler(n, matrix(1, 4, p), v, rate, C, Vb, R)
  expect_named(out, c("sigma2", "B"))
  E_sigma2 <- rate / (v - 1)
  sd_sigma2 <- E_sigma2 / sqrt(v - 2)
  expect_lt(abs(mean(out$sigma2) - E_sigma2) / (sd_sigma2 / sqrt(n)), 4.5)
  B <- vec_draws(out$B)
  expect_mean_close(B, C, E_sigma2 * kronecker(R, Vb))
  expect_equal(stats::cov(B), E_sigma2 * kronecker(R, Vb), tolerance = 0.1)
})

test_that("check_pds and make_pds symmetrize, repair, or reject matrices", {
  A <- matrix(c(2, 0.5, 0.5, 1), 2)
  expect_identical(check_pds(A), A)
  nudged <- A + matrix(c(0, 1e-3, 0, 0), 2)
  expect_equal(check_pds(nudged), (nudged + t(nudged)) / 2)
  singular <- matrix(1, 2, 2)
  expect_warning(fixed <- check_pds(singular), "not positive definite")
  expect_equal(fixed, singular + 3e-4 * diag(2))
  expect_error(check_pds(matrix(c(1, 5, 0, 1), 2)), "symmetry is outside of tolerance")
  expect_identical(make_pds(A), A)
  expect_warning(fixed <- make_pds(singular, eps = 0.01), "not positive definite")
  expect_equal(fixed, singular + 0.01 * diag(2))
  expect_equal(suppressWarnings(make_pds(matrix(c(1, 2, 0, 1), 2))),
               matrix(1, 2, 2) + 1e-4 * diag(2))
})
