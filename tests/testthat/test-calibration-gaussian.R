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
  mapped <- e
  mapped$block_cols <- matrix(c(4L, 2L, 1L, 3L), 2, byrow = TRUE)
  remapped <- spDBL:::predict_calibration(f$eta, 1L, mapped)
  expect_equal(remapped$mu[, as.vector(t(mapped$block_cols))], p$mu)
  before <- .Random.seed
  expect_identical(spDBL:::predict_calibration(f$eta, 1L, e), p)
  expect_identical(before, .Random.seed)
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
      covariance[block(i), block(j)] + p$Sigma + diag(tau[i + 1L], n)
  }
  L <- chol(covariance)
  residual <- as.numeric(t(f$z - p$mu))
  expected <- -n * e$T / 2 * log(2 * pi) - sum(log(diag(L))) -
    sum(forwardsolve(t(L), residual)^2) / 2
  filtered <- spDBL:::.cal_filter_discrepancy(p, tau, f$z, f$U, rep(0, n), diag(n))
  expect_equal(filtered$log_likelihood, expected, tolerance = 1e-10)
})
