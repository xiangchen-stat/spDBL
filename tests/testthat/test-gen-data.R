# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

test_that("random matrix generators are scaled and positive definite", {
  set.seed(61)
  A <- gen_ran_matrix(3, 4)
  expect_identical(dim(A), c(3L, 4L))
  expect_equal(max(abs(A)), 1)
  # Entries stay in [-1, 1] even when the most negative entry dominates.
  for (s in 1:30) {
    set.seed(s)
    expect_equal(max(abs(gen_ran_matrix(3, 4))), 1)
  }
  P <- gen_pd_matrix(4)
  expect_equal(P, t(P))
  expect_equal(max(P), 1)
  expect_true(all(eigen(P, symmetric = TRUE)$values > 0))
})

test_that("kernel generators match elementwise definitions", {
  set.seed(62)
  loc <- matrix(runif(10), 5, 2)
  D <- as.matrix(dist(loc))
  expect_equal(gen_gp_kernel(loc, phi = 2, sigma2 = 1.5, tau2 = 0.1),
               exp_kernel_loop(loc, loc, 2, 1.5) + diag(0.1, 5), ignore_attr = TRUE)
  expect_equal(gen_exp_kernel(loc, phi = 2, sigma2 = 1.5, tau2 = 0.1),
               exp_kernel_loop(loc, loc, 2, 1.5) + diag(0.1, 5), ignore_attr = TRUE)
  expect_equal(gen_exp_kernel(loc, phi = 3), exp(-3 * D), ignore_attr = TRUE)
  expect_equal(gen_expsq_kernel(loc, phi = 2, sigma2 = 1.5, tau2 = 0.1),
               1.5 * exp(-2 * D^2) + diag(0.1, 5), ignore_attr = TRUE)
  expect_equal(gen_expsq_kernel(loc, phi = 2), exp(-2 * D^2), ignore_attr = TRUE)
})

test_that("gen_ffbs_data simulates the documented MNIW DLM structure", {
  set.seed(63)
  N <- 6; S <- 2; p <- 3; nT <- 4
  out <- gen_ffbs_data(N, S, p, nT)
  expect_named(out, c("Y", "para"))
  expect_named(out$Y, paste0("Y", seq_len(nT)))
  expect_true(all(vapply(out$Y, function(y) identical(dim(y), c(6L, 2L)), logical(1))))
  pa <- out$para
  expect_equal(pa$F0[, 1], rep(1, N))
  expect_equal(pa$G0, diag(p))
  expect_equal(pa$n0, S + 2)
  expect_equal(pa$dist, as.matrix(dist(pa$F0)), ignore_attr = TRUE)
  expect_equal(diag(pa$V0), rep(1.01, N), ignore_attr = TRUE)
  expect_true(all(eigen(pa$V0, symmetric = TRUE)$values > 0))
  expect_equal(crossprod(pa$RV0), pa$V0)
  expect_equal(crossprod(pa$RSigma), pa$Sigma)
  expect_true(all(eigen(pa$Sigma, symmetric = TRUE)$values > 0))
  expect_error(gen_ffbs_data(N, S, 1, nT), "p is less than 2")
})

test_that("gen_ffbs_csv writes every simulated matrix", {
  set.seed(64)
  dir <- tempfile("ffbscsv")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  prefix <- paste0(dir, "/")
  gen_ffbs_csv(N = 5, S = 2, p = 2, nT = 3, path = prefix)
  expected <- c("F0.csv", "G0.csv", "m0.csv", "MM0.csv", "RSigma.csv", "RV0.csv",
                "Sigma.csv", "V0.csv", "W0.csv", "Y1.csv", "Y2.csv", "Y3.csv")
  expect_setequal(list.files(dir), expected)
  read <- function(name) unname(as.matrix(utils::read.csv(file.path(dir, name))))
  expect_identical(dim(read("Y3.csv")), c(5L, 2L))
  expect_equal(crossprod(read("RSigma.csv")), read("Sigma.csv"), tolerance = 1e-10)
  expect_equal(read("G0.csv"), diag(2))
  expect_error(gen_ffbs_csv(5, 2, 1, 3, prefix), "p is less than 2")
})

test_that("AR design builders lag the responses", {
  Y <- lapply(1:4, function(t) matrix(t * 10 + 1:6, 2, 3))
  ones <- matrix(1, 2, 3)
  ar1 <- gen_F_ls_AR1(Y, 4)
  expect_identical(ar1[[1]], ones)
  for (t in 2:4) expect_identical(ar1[[t]], Y[[t - 1]])
  ar2 <- gen_F_ls_AR2(Y, 4)
  expect_length(ar2, 4L)
  expect_identical(ar2[[1]], cbind(ones, ones))
  expect_identical(ar2[[2]], cbind(ones, Y[[1]]))
  for (t in 3:4) expect_identical(ar2[[t]], cbind(Y[[t - 2]], Y[[t - 1]]))
  # Short series work too.
  expect_identical(gen_F_ls_AR2(Y[1], 1), list(cbind(ones, ones)))
  expect_length(gen_F_ls_AR2(Y[1:2], 2), 2L)
})

test_that("AR designs are causal and deterministic", {
  set.seed(67)
  Y <- lapply(1:5, function(t) matrix(rnorm(6), 2, 3))
  before <- .Random.seed
  designs <- list(gen_F_ls_AR1(Y, 5), gen_F_ls_AR2(Y, 5))
  expect_identical(.Random.seed, before)  # no random numbers are drawn
  for (t in 1:5) {
    # Changing the outputs from time t on leaves the designs up to t unchanged.
    changed <- Y
    for (k in t:5) changed[[k]] <- changed[[k]] + 100
    for (d in 1:2) {
      redone <- list(gen_F_ls_AR1, gen_F_ls_AR2)[[d]](changed, 5)
      expect_identical(redone[seq_len(t)], designs[[d]][seq_len(t)])
    }
  }
})

test_that("episode designs slice the full designs block by block", {
  set.seed(66)
  Y <- lapply(1:3, function(t) matrix(rnorm(2 * 6), 2, 6))
  ind <- generate.grid.rowsnake(2, 6, 2, 2)
  ar1 <- gen_F_ls_AR1_EP(Y, 3, nrow(ind), ind)
  ar2 <- gen_F_ls_AR2_EP(Y, 3, nrow(ind), ind)
  expect_length(ar1, 9L)
  expect_length(ar2, 9L)
  full1 <- gen_F_ls_AR1(Y, 3)
  for (t in 1:3) for (j in 1:3) {
    cols <- ind$L[j]:ind$R[j]
    expect_identical(ar1[[(t - 1) * 3 + j]], full1[[t]][, cols, drop = FALSE])
    if (t == 3) expect_identical(ar2[[(t - 1) * 3 + j]],
                                 cbind(Y[[1]][, cols], Y[[2]][, cols]))
  }
})

test_that("recover_from_EP_exact reassembles episode blocks", {
  fields <- lapply(1:2, function(t) matrix(t * 100 + 1:12, 2, 6))
  ind <- generate.grid.rowsnake(2, 6, 2, 2)
  blocks <- unlist(lapply(fields, function(f) lapply(1:3, function(j)
    f[, ind$L[j]:ind$R[j], drop = FALSE])), recursive = FALSE)
  rec <- recover_from_EP_exact(blocks, nT_ori = 2, nT = 6)
  expect_named(rec, c("T1", "T2"))
  expect_equal(rec$T1, fields[[1]], ignore_attr = TRUE)
  expect_equal(rec$T2, fields[[2]], ignore_attr = TRUE)
  # Blocks stored out of grid order are placed by block_cols.
  perm <- matrix(c(5, 6, 1, 2, 3, 4), 3, byrow = TRUE)
  shuffled <- unlist(lapply(fields, function(f) lapply(1:3, function(j)
    f[, perm[j, ], drop = FALSE])), recursive = FALSE)
  rec <- recover_from_EP_exact(shuffled, 2, 6, block_cols = perm)
  expect_equal(rec$T2, fields[[2]])
  expect_error(recover_from_EP_exact(blocks, 4, 6), "divisible")
  expect_error(recover_from_EP_exact(blocks, 2, 5), "divisible")
  expect_error(recover_from_EP_exact(blocks, 2, 6, block_cols = perm[1:2, ]),
               "permutation")
  expect_error(recover_from_EP_exact(blocks, 2, 6, block_cols = matrix(c(1, 1:5), 3)),
               "permutation")
})

test_that("recover_from_EP_MC reassembles every Monte Carlo draw", {
  set.seed(67)
  nsam <- 3
  fields <- lapply(1:2, function(t) array(rnorm(2 * 6 * nsam), c(2, 6, nsam)))
  blocks <- unlist(lapply(fields, function(f) lapply(1:3, function(j)
    f[, (2 * j - 1):(2 * j), , drop = FALSE])), recursive = FALSE)
  rec <- recover_from_EP_MC(blocks, nT_ori = 2, nT = 6, nsam = nsam)
  expect_named(rec, c("T1", "T2"))
  expect_equal(rec$T1, fields[[1]], ignore_attr = TRUE)
  expect_equal(rec$T2, fields[[2]], ignore_attr = TRUE)
  expect_error(recover_from_EP_MC(blocks, 2, 5, nsam), "doesn't match nT")
})
