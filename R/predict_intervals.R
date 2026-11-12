# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Pointwise predictive distributions behind predict(..., interval = "prediction").
#
# Known autoregressive rows (an ordinary emulator). For episode e, with the
# smoothed state Theta_e | Sigma ~ MN(s_e, S_e, Sigma), the input-GP cross
# covariance J and A = V^{-1} J, a new run's block is
#   y_new = (F_new - A' F_e) Theta_e + A' Y_e + eps,  eps ~ MN(0, q, Sigma),
#   q = gp_sigma2 + gp_tau2 - J' V^{-1} J,
# so y_new | Sigma ~ MN(mean, u, Sigma) with u = h S_e h' + q. Integrating Sigma:
#   inverse-Wishart, Sigma ~ IW(n_E, D_E): t_{n_E - C + 1}(mean, u D_cc / (n_E - C + 1));
#   sigma2R, sigma2 ~ IG(n_E, D_E):        t_{2 n_E}(mean, u R_cc D_E / n_E),
# with n_E, D_E the final forward-filter quantities, as in the samplers.
#
# Design-GP rows (a prepared calibration object). Given emulator draw m the
# field at a new input is Gaussian, with the mean and pointwise variance used
# by calibration_fit() (including the propagated design-GP variance when
# calibration_prepare(design_uncertainty = TRUE)). The predictive is the
# equal-weight mixture over the M draws; it is summarised by its mean and
# variance (law of total variance), with Gaussian intervals.

.check_level <- function(level) {
  if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must be a single number in (0, 1)")
  level
}

.time_list <- function(a) {
  out <- lapply(seq_len(dim(a)[3]), function(t) matrix(a[, , t], dim(a)[1], dim(a)[2]))
  names(out) <- paste0("T", seq_along(out))
  out
}

.list_array <- function(x) {
  x <- lapply(x, as.matrix)
  array(unlist(x), c(nrow(x[[1]]), ncol(x[[1]]), length(x)))
}

# Default autoregressive designs for new runs whose outputs are known, built
# as emulator_predict() builds them.
.emulator_default_designs <- function(emulator, input_new, dt_pde_test) {
  if (is.null(dt_pde_test))
    stop("Supply F_new_ls or dt_pde_test, or use a calibration_prepare() result")
  s <- emulator$setup
  n_test <- nrow(input_new)
  index <- generate.grid.rowsnake(n_test, s$fncol, n_test, s$bncol)
  if (s$AR_choice == 2) gen_F_ls_AR2_EP(dt_pde_test, s$nT_ori, s$n_b, index) else
    gen_F_ls_AR1_EP(dt_pde_test, s$nT_ori, s$n_b, index)
}

.predictive_known <- function(emulator, input_new, F_new_ls, level) {
  s <- emulator$setup
  ff <- emulator$para_ffbs$ff
  St <- emulator$para_ffbs$bs$St
  E <- s$nT
  C <- s$bncol
  input <- as.matrix(s$pde_para_train)
  N <- nrow(input)
  dd <- as.matrix(stats::dist(rbind(input, input_new)))
  J <- s$gp_sigma2 * exp(-s$phi_para * dd[seq_len(N), N + seq_len(nrow(input_new)), drop = FALSE])
  V_chol <- chol(s$V)
  A <- backsolve(V_chol, forwardsolve(t(V_chol), J))
  q <- s$gp_sigma2 + s$gp_tau2 - colSums(J * A)
  nt <- ff[[E]]$nt
  if (identical(s$covariance, "iw")) {
    df <- nt - C + 1
    column_scale <- diag(as.matrix(ff[[E]]$Dt)) / df
  } else {
    df <- 2 * nt
    column_scale <- diag(as.matrix(s$R_episode)) * as.numeric(ff[[E]]$Dt) / nt
  }
  p <- dim(St)[1]
  scale2 <- lapply(seq_len(E), function(e) {
    H <- as.matrix(F_new_ls[[e]]) - crossprod(A, as.matrix(s$F_ls[[e]]))
    u <- rowSums((H %*% matrix(St[, , e], p, p)) * H) + q
    outer(u, column_scale)
  })
  scale <- sqrt(.list_array(recover_from_EP_exact(scale2, s$nT_ori, E,
                                                  block_cols = s$block_cols)))
  mean <- .list_array(emulator_predict(emulator, input_new, F_new_ls = F_new_ls))
  half <- stats::qt((1 + level) / 2, df) * scale
  list(mean = .time_list(mean), sd = .time_list(scale * sqrt(df / (df - 2))),
       lower = .time_list(mean - half), upper = .time_list(mean + half),
       level = level, distribution = "student_t", df = df)
}

.predictive_design_gp <- function(cal, input_new, level) {
  n <- nrow(input_new)
  M <- cal$M
  S <- cal$S
  T <- cal$T
  # Work that depends only on the new inputs, done once: the design GP's rows
  # and variance factors, and the input GP's cross-covariance J and
  # conditional variance (as in FFBS_predict_exact() and predict_calibration()).
  designs <- lapply(seq_len(n), function(i) predict_offline_ar2(cal$predictor, input_new[i, ]))
  F_new <- lapply(seq_len(cal$E), function(e)
    do.call(rbind, lapply(designs, function(z) z$F[[e]])))
  time_factor <- matrix(vapply(designs, function(z) .cal_time_factor(z, cal), numeric(T)),
                        n, T, byrow = TRUE)
  input <- as.matrix(cal$input)
  N <- nrow(input)
  distances <- as.matrix(stats::dist(rbind(input, input_new)))
  J <- cal$gp_sigma2 * exp(-cal$phi_para * distances[seq_len(N), N + seq_len(n), drop = FALSE])
  conditional <- cal$gp_sigma2 + cal$gp_tau2 - colSums(forwardsolve(t(cal$V_chol), J)^2)
  if (any(conditional <= 0)) stop("Nonpositive emulator conditional variance")
  R_diag <- .cal_R_field_diag(cal, seq_len(S))
  first <- second <- array(0, c(n, S, T))
  design_uncertainty <- isTRUE(cal$design_uncertainty)
  if (design_uncertainty) {
    W <- cal$design_whitened
    if (is.null(W)) W <- .cal_design_whitened(cal$predictor)
  }
  # Time by time, episode by episode, so that each episode's draws are read
  # (or regenerated) once; per draw the arithmetic is that of
  # .cal_residual_weights(), .cal_design_factors() and predict_calibration().
  for (t in seq_len(T)) {
    mu <- array(0, c(n, S, M))
    design_sq <- matrix(0, M, S)
    for (b in seq_len(cal$B)) {
      e <- (t - 1L) * cal$B + b
      cols <- cal$block_cols[b, ]
      draws <- .cal_theta(cal, e)
      for (m in seq_len(M)) {
        theta <- matrix(draws[, , m], cal$p, cal$C)
        weights <- if (!is.null(cal$residual_weights)) cal$residual_weights[[m]][[e]] else
          backsolve(cal$V_chol, forwardsolve(t(cal$V_chol), cal$Y[[e]] - cal$F_ls[[e]] %*% theta))
        mu[, cols, m] <- F_new[[e]] %*% theta + crossprod(J, weights)
        if (design_uncertainty)
          design_sq[m, cols] <- colSums((W[, (e - 1L) * cal$p + seq_len(cal$p)] %*% theta /
                                         sqrt(N))^2)
      }
    }
    for (m in seq_len(M)) {
      variance <- outer(cal$posterior$sigma2[m] * conditional, R_diag)
      if (design_uncertainty) variance <- variance + outer(time_factor[, t], design_sq[m, ])
      first[, , t] <- first[, , t] + mu[, , m] / M
      second[, , t] <- second[, , t] + (variance + mu[, , m]^2) / M
    }
  }
  sd <- sqrt(pmax(second - first^2, 0))
  half <- stats::qnorm((1 + level) / 2) * sd
  list(mean = .time_list(first), sd = .time_list(sd), lower = .time_list(first - half),
       upper = .time_list(first + half), level = level,
       distribution = "gaussian_mixture_moments", n_draws = M,
       design_uncertainty = isTRUE(cal$design_uncertainty))
}
