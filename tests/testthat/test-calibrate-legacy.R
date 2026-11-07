# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

test_that("legacy helpers are deprecated and warn once per session", {
  legacy <- c("update_y_eta", "update_y_eta_one", "update_muSigma_eta_one",
              "sample_y_eta_one", "cal_Bt_bt", "scale_uniform",
              "scale_back_uniform", "expit", "cal_jacobian_logit_uniform")
  ns <- asNamespace("spDBL")
  for (f in legacy) expect_true(any(grepl(sprintf('.legacy_deprecated("%s")', f),
    deparse(body(get(f, ns))), fixed = TRUE)), info = f)
  seen <- get(".spdbl_deprecated", ns)
  rm(list = ls(seen), envir = seen)
  old <- options(spDBL.warn_deprecated = NULL)
  on.exit({ options(old); rm(list = ls(seen), envir = seen) }, add = TRUE)
  expect_warning(expit(0), "expit\\(\\) is deprecated since spDBL 1.2.0",
                 class = "deprecatedWarning")
  expect_no_warning(expit(1))                      # once per session
  expect_warning(scale_uniform(0.5, 0, 1), "scale_uniform")
  rm(list = ls(seen), envir = seen)
  options(spDBL.warn_deprecated = FALSE)
  expect_no_warning(expit(0))                      # silenced by the option
})

# The remaining tests exercise the deprecated helpers themselves.
legacy_options <- options(spDBL.warn_deprecated = FALSE)

test_that("scaling, logistic and Cholesky helpers are exact transforms", {
  x <- c(2, 3.5, 5)
  u <- scale_uniform(x, 2, 5)
  expect_equal(u, c(0, 0.5, 1))
  expect_equal(scale_back_uniform(u, 2, 5), x)
  expect_equal(scale_uniform(c(1, 1), c(1, 0), c(1, 2)), c(0, 0.5))
  expect_error(scale_uniform(6, 2, 5), "out of range")
  expect_equal(expit(c(-3, 0, 2)), stats::plogis(c(-3, 0, 2)))
  A <- matrix(c(4, 1, 1, 3), 2)
  expect_equal(inv_chol(A), solve(A))
})

test_that("the logit-uniform Jacobian is the derivative of the logit map", {
  eta <- c(0.3, 2.5)
  low <- c(0, 1)
  high <- c(1, 4)
  h <- 1e-6
  logit_map <- function(e) stats::qlogis((e - low) / (high - low))
  derivs <- (logit_map(eta + h) - logit_map(eta - h)) / (2 * h)
  expect_equal(cal_jacobian_logit_uniform(eta, low, high), prod(derivs), tolerance = 1e-6)
  expect_equal(cal_jacobian_logit_uniform(eta, low, high, log = TRUE),
               log(prod(derivs)), tolerance = 1e-6)
})

test_that("the Jacobian skips fixed coordinates, not coordinates equal to zero", {
  low <- c(-1, 0, 0.5)
  high <- c(1, 0, 0.5)
  # Coordinate 1 is free and equals zero; coordinates 2 and 3 are fixed.
  expect_equal(cal_jacobian_logit_uniform(c(0, 0, 0.5), low, high), 2)
  expect_equal(cal_jacobian_logit_uniform(c(0, 0, 0.5), low, high, log = TRUE), log(2))
})

test_that("gen_Jt is the exponential cross-covariance to new inputs", {
  set.seed(91)
  input <- matrix(runif(10), 5, 2)
  new <- matrix(runif(4), 2, 2)
  expect_equal(gen_Jt(input, new, phi_para = 2, gp_sigma2 = 1.3),
               exp_kernel_loop(input, new, 2, 1.3), ignore_attr = TRUE)
  expect_equal(gen_Jt(input, new[1, ], phi_para = 2),
               as.numeric(exp_kernel_loop(input, new[1, , drop = FALSE], 2, 1.1)),
               ignore_attr = TRUE)
})

test_that("gen_Jt leaves the nugget out when a new input equals a training input", {
  input <- matrix(c(0, 1, 0.5, 0.2), 2)
  J <- gen_Jt(input, input[1, ], phi_para = 2, gp_sigma2 = 1.3, gp_tau2 = 0.5)
  expect_equal(J, c(1.3, 1.3 * exp(-2 * sqrt(sum((input[1, ] - input[2, ])^2)))),
               ignore_attr = TRUE)
})

test_that("cal_Bt_bt is the Gaussian posterior of y given z - u", {
  set.seed(92)
  S <- 3
  Sigma <- rand_spd(S)
  mu <- rnorm(S)
  zu <- rnorm(S)
  tau2 <- 0.4
  out <- cal_Bt_bt(mu, solve(Sigma), tau2, zu, diag(S))
  K <- Sigma %*% solve(Sigma + tau2 * diag(S))
  expect_equal(out$Bt, Sigma - K %*% Sigma, tolerance = 1e-10)
  expect_equal(as.numeric(out$Btbt), as.numeric(mu + K %*% (zu - mu)), tolerance = 1e-10)
})

test_that("simulated discrepancies have the stated increments and noise", {
  set.seed(93)
  s <- 1000
  ycal <- matrix(rnorm(3 * s), 3, s)
  para <- list(b = 0.9, n0 = 20, d0 = 2, u0 = rep(0.5, s))
  U <- 0.5 * diag(s)
  out <- gen_calibrate_data(ycal, para, U)
  expect_named(out, c("u", "z", "tau2"))
  expect_identical(dim(out$z), dim(ycal))
  expect_true(all(out$tau2 > 0))
  ratio <- function(x, v) stats::var(x) / v
  incr <- rbind(out$u[1, ] - para$u0, diff(out$u))
  for (t in 1:3) {
    expect_equal(ratio(out$z[t, ] - ycal[t, ] - out$u[t, ], out$tau2[t]), 1, tolerance = 0.2)
    expect_equal(ratio(incr[t, ], 0.5 * out$tau2[t]), 1, tolerance = 0.2)
  }
  unc <- gen_calibrate_data_uncorr(ycal, para, U)
  expect_named(unc, c("z", "tau2"))
  # z_t - y_t - u0 accumulates every increment plus the current noise.
  for (t in 1:3) {
    v <- 0.5 * sum(unc$tau2[1:t]) + unc$tau2[t]
    expect_equal(ratio(unc$z[t, ] - ycal[t, ] - para$u0, v), 1, tolerance = 0.2)
  }
})

test_that("gen_prior_u_tau2 draws the discrepancy prior", {
  set.seed(94)
  s <- 1000
  out <- gen_prior_u_tau2(n0 = 20, d0 = 2, b = 0.95, m0 = rep(1, s), M0 = diag(s),
                          U = 0.5 * diag(s), nT_cal = 3)
  expect_named(out, c("u0", "u", "tau2_0", "tau2"))
  expect_identical(dim(out$u), c(3L, as.integer(s)))
  expect_equal(stats::var(as.numeric(out$u0) - 1) / out$tau2_0, 1, tolerance = 0.2)
  incr <- rbind(out$u[1, ] - out$u0, diff(out$u))
  for (t in 1:3)
    expect_equal(stats::var(incr[t, ]) / (0.5 * out$tau2[t]), 1, tolerance = 0.2)
})

sir_parms <- function(eta, Nx = 4, Ny = 4, N = 100) {
  list(eta = eta, Nx = Nx, Ny = Ny, N = N,
       Gridx = ReacTran::setup.grid.1D(x.up = 0, x.down = Nx, N = Nx),
       Gridy = ReacTran::setup.grid.1D(x.up = 0, x.down = Ny, N = Ny))
}

test_that("SIR reduces to the reaction terms without diffusion", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("ReacTran")
  set.seed(95)
  S <- runif(16, 50, 100); I <- runif(16, 0, 10); R <- runif(16, 0, 5)
  d <- SIR(0, c(S, I, R), sir_parms(c(0.3, 0.1, 0, 0, 0)))[[1]]
  expect_length(d, 48L)
  expect_equal(d[1:16], -0.3 * S * I / 100)
  expect_equal(d[17:32], 0.3 * S * I / 100 - 0.1 * I)
  expect_equal(d[33:48], 0.1 * I)
  expect_equal(sum(d), 0)
})

test_that("SIR diffusion leaves uniform interiors unchanged and drains boundaries", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("ReacTran")
  y <- c(rep(5, 16), rep(0, 32))
  d <- SIR(0, y, sir_parms(c(0, 0, 0.3, 0.2, 0.1)))[[1]]
  dS <- matrix(d[1:16], 4, 4)
  expect_equal(dS[2:3, 2:3], matrix(0, 2, 2))
  expect_true(all(dS[c(1, 4), ] < 0) && all(dS[, c(1, 4)] < 0))
  expect_equal(d[17:48], rep(0, 32))
})

test_that("gen_pde returns log infections on the grid", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("ReacTran")
  frozen <- gen_pde(c(0, 0, 0, 0, 0), Nx = 4, Ny = 4, N = 100, nT_ori = 3)
  expect_identical(dim(frozen), c(3L, 16L))
  expect_true(all(is.finite(frozen)) && all(frozen >= 0))
  # Zero rates give a constant solution.
  expect_equal(frozen[3, ], frozen[1, ], ignore_attr = TRUE)
  expect_equal(max(frozen), log(2))
  moving <- gen_pde(c(0.5, 0.1, 0.05, 0.05, 0.05), Nx = 4, Ny = 4, N = 100, nT_ori = 3)
  expect_equal(moving[1, ], frozen[1, ])
  expect_gt(sum(moving[3, ]), sum(moving[1, ]))
})

test_that("gen_pde sizes the solver work space by the grid", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("ReacTran")
  requested <- numeric()
  local_mocked_bindings(ode.2D = function(y, ..., lrw) {
    requested <<- c(requested, lrw)
    cbind(0, matrix(0, 3, length(y)))
  }, .package = "deSolve")
  gen_pde(c(0.5, 0.1, 0.05, 0.05, 0.05), Nx = 4, Ny = 4, N = 100, nT_ori = 3)
  gen_pde(c(0.5, 0.1, 0.05, 0.05, 0.05), Nx = 20, Ny = 30, N = 100, nT_ori = 3)
  expect_equal(requested, c(1e5, 500 * 3 * 20 * 30))
})

# A tiny calibration problem for the legacy one-draw updates. gen_pde is mocked
# so the expected means can be computed from a known simulator output.
legacy_problem <- function(AR_choice = 1, nsam = 2) {
  set.seed(96)
  N <- 5; Nx <- 2; Ny <- 3; N_sp <- 6; bncol <- 2; n_b <- 3; nT_ori <- 4
  p <- AR_choice * bncol
  input <- matrix(runif(N * 2), N, 2)
  # gp_sigma2 + gp_tau2 = 1 and a common variance draw keep the expected values
  # free of the legacy conventions (prior variance fixed at 1, mean vs. selected
  # sigma^2 draw) that the newer calibration code handles differently.
  phi <- 2; gp_tau2 <- 1e-3; gp_sigma2 <- 1 - gp_tau2
  V <- exp_kernel_loop(input, input, phi, gp_sigma2) + diag(gp_tau2, N)
  res_ffbs <- lapply(seq_len(nT_ori * n_b), function(e)
    array(rnorm(p * bncol * nsam), c(p, bncol, nsam)))
  res_ffbs$Sigma <- rep(0.25, nsam)
  list(eta = c(0.3, 0.1, 0.01, 0.01, 0.01), Nx = Nx, Ny = Ny, N_people = 100,
       nT_ori = nT_ori, nT = nT_ori * n_b, nT_cal = 3, N_sp = N_sp, N_sp_train = 4,
       AR_choice = AR_choice, bncol = bncol, res_ffbs = res_ffbs, n_b_cal = n_b,
       ind_block_cal = generate.grid.rowsnake(1, N_sp, 1, bncol),
       t_VinvY_FTheta_ls = lapply(1:3, function(t) array(rnorm(N_sp * N * nsam),
                                                         c(N_sp, N, nsam))),
       Vinv = solve(V), R_cal_train = rand_spd(4, 0.5), input = input,
       input_cal = c(0.4, 0.6), phi_para = phi, gp_sigma2 = gp_sigma2,
       gp_tau2 = gp_tau2, ind_train_cal = c(1, 3, 4, 6),
       y_eta = matrix(runif(nT_ori * N_sp), nT_ori, N_sp))
}

legacy_expected <- function(lp, draw) {
  J <- exp_kernel_loop(lp$input, matrix(lp$input_cal, 1), lp$phi_para, lp$gp_sigma2)
  mu <- matrix(NA_real_, lp$nT_cal, lp$N_sp)
  for (t in seq_len(lp$nT_cal)) {
    # Causal designs: a lag before the first time is a row of ones.
    lag <- function(k) if (t - k >= 1) lp$y_eta[t - k, , drop = FALSE] else
      matrix(1, 1, lp$N_sp)
    lagged <- if (lp$AR_choice == 1) lag(1) else cbind(lag(2), lag(1))
    for (b in seq_len(lp$n_b_cal)) {
      cols <- (2 * b - 1):(2 * b)
      design <- if (lp$AR_choice == 1) lagged[, cols, drop = FALSE] else
        lagged[, c(cols, lp$N_sp + cols), drop = FALSE]
      mu[t, cols] <- design %*% lp$res_ffbs[[(t - 1) * lp$n_b_cal + b]][, , draw]
    }
    mu[t, ] <- mu[t, ] + lp$t_VinvY_FTheta_ls[[t]][, , draw] %*% J
  }
  # Kriging variance of the unit-variance input GP at input_cal.
  list(mu = mu[, lp$ind_train_cal, drop = FALSE],
       c = as.numeric(lp$gp_sigma2 + lp$gp_tau2 - t(J) %*% lp$Vinv %*% J))
}

call_legacy <- function(fun, lp, ...) {
  args <- lp[setdiff(names(lp), "y_eta")]
  do.call(fun, c(args, list(...)))
}

test_that("legacy y(eta) updates combine the PDE design with kriging residuals", {
  lp <- legacy_problem()
  local_mocked_bindings(gen_pde = function(...) lp$y_eta, .package = "spDBL")
  ex1 <- legacy_expected(lp, 1)
  ex2 <- legacy_expected(lp, 2)
  lp_tiny <- lp
  lp_tiny$res_ffbs$Sigma <- c(1e-14, 1e-14)
  set.seed(97)
  all <- call_legacy(update_y_eta, lp_tiny, nsam = 2)
  expect_identical(dim(all$y), c(3L, 4L, 2L))
  expect_equal(all$mu[, , 1], ex1$mu, tolerance = 1e-10)
  expect_equal(all$mu[, , 2], ex2$mu, tolerance = 1e-10)
  expect_equal(all$sigma2_cal, 1e-14 * ex1$c, tolerance = 1e-10)
  expect_equal(all$y, all$mu, tolerance = 1e-5)
  one <- call_legacy(update_y_eta_one, lp, ind_sam = 2)
  expect_equal(one$mu, ex2$mu, tolerance = 1e-10)
  expect_equal(one$Sigma, 0.25 * ex2$c * lp$R_cal_train, tolerance = 1e-10)
  expect_identical(dim(one$y), c(3L, 4L))
  ms <- call_legacy(update_muSigma_eta_one, lp, ind_sam = 1)
  expect_equal(ms$mu, ex1$mu, tolerance = 1e-10)
  expect_equal(ms$Sigma, 0.25 * ex1$c * lp$R_cal_train, tolerance = 1e-10)
  y <- sample_y_eta_one(1, 3, 4, ms$mu, 1e-14 * diag(4))
  expect_equal(y, ms$mu, tolerance = 1e-5)
})

test_that("legacy updates use the input-GP prior variance and the selected draw", {
  lp <- legacy_problem()
  lp$gp_sigma2 <- 1.1
  lp$gp_tau2 <- 1e-4
  lp$Vinv <- solve(exp_kernel_loop(lp$input, lp$input, lp$phi_para, 1.1) +
                     diag(1e-4, nrow(lp$input)))
  lp$res_ffbs$Sigma <- c(0.2, 0.3)
  local_mocked_bindings(gen_pde = function(...) lp$y_eta, .package = "spDBL")
  ex <- legacy_expected(lp, 2)
  expect_equal(call_legacy(update_y_eta_one, lp, ind_sam = 2)$Sigma,
               0.3 * ex$c * lp$R_cal_train, tolerance = 1e-10)
  expect_equal(call_legacy(update_muSigma_eta_one, lp, ind_sam = 2)$Sigma,
               0.3 * ex$c * lp$R_cal_train, tolerance = 1e-10)
  expect_equal(call_legacy(update_y_eta, lp, nsam = 2)$sigma2_cal, 0.25 * ex$c,
               tolerance = 1e-10)
  # At a training input the variance is about the nugget. With the prior
  # variance fixed at 1 it was 1 - 1.1 < 0 and sampling failed.
  lp$input_cal <- lp$input[2, ]
  near <- legacy_expected(lp, 1)
  expect_gt(near$c, 0)
  expect_lt(near$c, 1e-3)
  expect_equal(call_legacy(update_muSigma_eta_one, lp, ind_sam = 1)$Sigma,
               0.2 * near$c * lp$R_cal_train, tolerance = 1e-8)
  expect_identical(dim(call_legacy(update_y_eta_one, lp, ind_sam = 1)$y), c(3L, 4L))
})

test_that("legacy AR2 updates use two lags once both are observed", {
  lp <- legacy_problem(AR_choice = 2)
  local_mocked_bindings(gen_pde = function(...) lp$y_eta, .package = "spDBL")
  ex <- legacy_expected(lp, 1)
  for (out in list(call_legacy(update_y_eta, lp, nsam = 2)$mu[, , 1],
                   call_legacy(update_y_eta_one, lp, ind_sam = 1)$mu,
                   call_legacy(update_muSigma_eta_one, lp, ind_sam = 1)$mu))
    expect_equal(out, ex$mu, tolerance = 1e-10)
})

options(legacy_options)
