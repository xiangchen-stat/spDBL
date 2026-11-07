# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Two times x two blocks (E = 4 episodes, S = 4 sites). The design GP is the
# squared-exponential kernel shared by all episodes unless design_kernel and
# design_groups say otherwise ("matern52", "time": one kernel per time step).
# With design_groups = "time" the last design column oscillates faster at
# time 2, so the two times' kernels get different hyperparameters.
calibration_fixture <- function(n_draws = 2L, design_kernel = "sqexp",
                                design_groups = "shared") {
  set.seed(81)
  input <- cbind(a = seq(0.05, 0.95, length.out = 8),
                 b = c(0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4, 0.8))
  frequency <- function(e) if (design_groups == "time") 1 + 4 * (e > 2) else 1
  F <- lapply(1:4, function(e) cbind(1, input[, 1], input[, 2],
                                   sin(frequency(e) * input[, 1]) + e / 10))
  theta <- lapply(1:4, function(e)
    array(rnorm(4 * 2 * n_draws, sd = 0.15), c(4, 2, n_draws)))
  Y <- lapply(1:4, function(e) F[[e]] %*% theta[[e]][, , 1] +
                cbind(input[, 1]^2, cos(input[, 2])))
  posterior <- list(theta = theta, sigma2 = seq(0.08, 0.12, length.out = n_draws),
                    joint_trajectories = TRUE)
  R <- 0.2^as.matrix(dist(matrix(1:4)))
  emulator <- spDBL:::make_offline_calibration_input(input, Y, F, posterior,
    matrix(1:4, 2, byrow = TRUE), R, phi_para = 2,
    lower = c(0, 0), upper = c(1, 1), design_kernel = design_kernel,
    design_groups = design_groups)
  emulator$parameter_names <- c("a", "b")
  list(emulator = emulator,
    z = matrix(c(0.2, 0.3, 0.1, 0.4, 0.7, 0.5, 0.6, 0.2), 2),
    U = R, eta = c(0.4, 0.6))
}

# Design-row covariance of a calibration prediction at time t (dense).
design_cov <- function(pred, t) spDBL:::.cal_time_covariance(pred, t) - pred$Sigma

tiny_emulator <- function(covariance = "sigma2R", sampling = "marginal") {
  input <- cbind(a = seq(0.05, 0.95, length.out = 8),
                 b = c(0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4, 0.8))
  responses <- lapply(1:3, function(t)
    outer(input[, 1] + t * input[, 2], 1:4, function(x, y) sin(x) + y / 10))
  emulator_learn(input, responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
    covariance = covariance, nsam = 3L, seed = 91L, sampling = sampling)
}
