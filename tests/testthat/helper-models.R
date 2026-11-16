# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Independent reference computations shared by the unit tests. None of these
# helpers call spDBL, so they can be used as ground truth.

# A well-conditioned random symmetric positive-definite matrix.
rand_spd <- function(n, scale = 1) {
  A <- matrix(stats::rnorm(n * n), n, n)
  scale * (crossprod(A) / n + diag(n))
}

rep_list <- function(x, n) if (is.list(x)) x else rep(list(x), n)

# Exact moments of the matrix-normal dynamic linear model
#   Theta_0 ~ MN(m0, M0, Sigma),
#   Theta_t = G_t Theta_{t-1} + MN(0, W_t, Sigma),
#   Y_t     = F_t Theta_t     + MN(0, V_t, Sigma),
# obtained by writing every state and observation as a linear map of the
# independent innovations and conditioning one stacked Gaussian (no recursion).
# For each t the result holds the one-step prediction (at, At), the filtered
# moments (mt, Mt), the smoothed moments (st, St) and the cumulative quadratic
# form quad = (Y_1:t - E Y_1:t)' Cov(Y_1:t)^{-1} (Y_1:t - E Y_1:t).
batch_dlm <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, nT = length(Y)) {
  F_ls <- rep_list(F_ls, nT)
  G_ls <- rep_list(G_ls, nT)
  W_ls <- rep_list(W_ls, nT)
  V_ls <- rep_list(V_ls, nT)
  p <- nrow(m0)
  N <- vapply(Y[seq_len(nT)], nrow, integer(1))
  nx <- p * (nT + 1L) + sum(N)
  gamma_rows <- function(t) p * t + seq_len(p)
  eps_rows <- function(t) p * (nT + 1L) + sum(N[seq_len(t - 1L)]) + seq_len(N[t])
  Omega <- matrix(0, nx, nx)
  Omega[seq_len(p), seq_len(p)] <- M0
  for (t in seq_len(nT)) {
    Omega[gamma_rows(t), gamma_rows(t)] <- W_ls[[t]]
    Omega[eps_rows(t), eps_rows(t)] <- V_ls[[t]]
  }
  A <- B <- mA <- mB <- vector("list", nT)
  prev <- cbind(diag(p), matrix(0, p, nx - p))
  prev_mean <- m0
  for (t in seq_len(nT)) {
    innov <- matrix(0, p, nx)
    innov[, gamma_rows(t)] <- diag(p)
    A[[t]] <- G_ls[[t]] %*% prev + innov
    mA[[t]] <- G_ls[[t]] %*% prev_mean
    noise <- matrix(0, N[t], nx)
    noise[, eps_rows(t)] <- diag(N[t])
    B[[t]] <- F_ls[[t]] %*% A[[t]] + noise
    mB[[t]] <- F_ls[[t]] %*% mA[[t]]
    prev <- A[[t]]
    prev_mean <- mA[[t]]
  }
  condition <- function(t, s) {
    prior_cov <- A[[t]] %*% Omega %*% t(A[[t]])
    if (s == 0) return(list(mean = mA[[t]], cov = prior_cov))
    Bs <- do.call(rbind, B[seq_len(s)])
    resid <- do.call(rbind, Y[seq_len(s)]) - do.call(rbind, mB[seq_len(s)])
    Cyy <- Bs %*% Omega %*% t(Bs)
    Cty <- A[[t]] %*% Omega %*% t(Bs)
    list(mean = mA[[t]] + Cty %*% solve(Cyy, resid),
         cov = prior_cov - Cty %*% solve(Cyy, t(Cty)),
         quad = t(resid) %*% solve(Cyy, resid))
  }
  lapply(seq_len(nT), function(t) {
    pred <- condition(t, t - 1L)
    filt <- condition(t, t)
    smooth <- condition(t, nT)
    list(at = pred$mean, At = pred$cov, mt = filt$mean, Mt = filt$cov,
         st = smooth$mean, St = smooth$cov, quad = filt$quad)
  })
}

# A small random DLM problem used by the filtering and smoothing tests.
dlm_problem <- function(seed, p = 2L, q = 2L, N = c(3L, 4L, 3L),
                        time_varying = TRUE) {
  set.seed(seed)
  nT <- length(N)
  F_ls <- lapply(N, function(n) matrix(stats::rnorm(n * p), n, p))
  G_ls <- lapply(seq_len(nT), function(t)
    diag(p) + 0.2 * matrix(stats::rnorm(p * p), p, p))
  W_ls <- lapply(seq_len(nT), function(t) rand_spd(p, 0.3))
  V_ls <- lapply(N, function(n) rand_spd(n, 0.5))
  if (!time_varying) {
    F_ls <- F_ls[[1]]
    G_ls <- G_ls[[1]]
    W_ls <- W_ls[[1]]
    V_ls <- V_ls[[1]]
  }
  list(Y = lapply(N, function(n) matrix(stats::rnorm(n * q), n, q)),
       F_ls = F_ls, G_ls = G_ls, W_ls = W_ls, V_ls = V_ls,
       m0 = matrix(stats::rnorm(p * q), p, q), M0 = rand_spd(p),
       n0 = q + 3, D0 = rand_spd(q), nT = nT, N = N, p = p, q = q)
}

# Exponential input kernel written out element by element.
exp_kernel_loop <- function(x, y, phi, sigma2) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  out <- matrix(0, nrow(x), nrow(y))
  for (i in seq_len(nrow(x))) for (j in seq_len(nrow(y)))
    out[i, j] <- sigma2 * exp(-phi * sqrt(sum((x[i, ] - y[j, ])^2)))
  out
}

# Log density of vec(X) ~ N(vec(M), SigmaC %x% SigmaR), computed directly.
log_dmatnorm <- function(X, M, SigmaR, SigmaC) {
  K <- kronecker(SigmaC, SigmaR)
  r <- as.numeric(X - M)
  L <- chol(K)
  z <- forwardsolve(t(L), r)
  -length(r) / 2 * log(2 * pi) - sum(log(diag(L))) - sum(z^2) / 2
}
