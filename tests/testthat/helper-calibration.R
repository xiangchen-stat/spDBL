calibration_fixture <- function(n_draws = 2L) {
  set.seed(81)
  input <- cbind(a = seq(0.05, 0.95, length.out = 8),
                 b = c(0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4, 0.8))
  F <- lapply(1:4, function(e) cbind(1, input[, 1], input[, 2],
                                   sin(input[, 1]) + e / 10))
  theta <- lapply(1:4, function(e)
    array(rnorm(4 * 2 * n_draws, sd = 0.15), c(4, 2, n_draws)))
  Y <- lapply(1:4, function(e) F[[e]] %*% theta[[e]][, , 1] +
                cbind(input[, 1]^2, cos(input[, 2])))
  posterior <- list(theta = theta, sigma2 = seq(0.08, 0.12, length.out = n_draws),
                    joint_trajectories = TRUE)
  R <- 0.2^as.matrix(dist(matrix(1:4)))
  emulator <- spDBL:::make_offline_calibration_input(input, Y, F, posterior,
    matrix(1:4, 2, byrow = TRUE), R, phi_para = 2,
    lower = c(0, 0), upper = c(1, 1))
  emulator$parameter_names <- c("a", "b")
  list(emulator = emulator,
    z = matrix(c(0.2, 0.3, 0.1, 0.4, 0.7, 0.5, 0.6, 0.2), 2),
    U = R, eta = c(0.4, 0.6))
}

tiny_emulator <- function(covariance = "sigma2R", sampling = "marginal") {
  input <- cbind(a = seq(0.05, 0.95, length.out = 8),
                 b = c(0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4, 0.8))
  responses <- lapply(1:3, function(t)
    outer(input[, 1] + t * input[, 2], 1:4, function(x, y) sin(x) + y / 10))
  emulator_learn(input, responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
    covariance = covariance, nsam = 3L, seed = 91L, sampling = sampling)
}
