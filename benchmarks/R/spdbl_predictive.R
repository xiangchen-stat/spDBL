# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Pointwise predictive distributions for spDBL emulators.
#
# Workaround: emulator_predict() returns predictive means (MC = FALSE) or
# coefficient-conditioned mean draws (MC = TRUE) only. It has no predictive
# variance or interval output, so interval coverage and CRPS cannot be
# computed from the public API alone. The functions below derive the
# marginal predictive distribution from objects that emulator_learn() and
# calibration_prepare() store, or from the package's internal predictor.
#
# Known AR rows (ordinary emulator). For episode e, with smoothed state
# Theta_e | Sigma ~ MN(s_e, S_e, Sigma), input-GP cross covariance J and
# A = V^{-1} J, a new run's block is
#   y_new = (F_new - A' F_e) Theta_e + A' Y_e + eps,  eps ~ MN(0, q, Sigma),
#   q = gp_sigma2 + gp_tau2 - J' V^{-1} J,
# so y_new | Sigma ~ MN(mean, u, Sigma) with u = h S_e h' + q. Integrating
#   IW model:      Sigma ~ IW(n_E, D_E)        -> t_{n_E - C + 1}(mean, u D_cc / (n_E - C + 1))
#   sigma2R model: sigma2 ~ IG(n_E, D_E), Sigma = sigma2 R -> t_{2 n_E}(mean, u R_cc D_E / n_E)
# where n_E, D_E are the final forward-filter quantities (the samplers in
# FFBS_sampling / sample_calibration_emulator use the same parameterisation).
#
# Design-GP rows (prepared calibration object). Conditional on draw m the
# field is Gaussian; its mean and pointwise variance come from the package's
# internal predict_calibration(covariance = "marginal"), the same quantities
# calibration_fit() uses. The variance is sigma2_m q(eta) R_field[s, s] plus,
# when calibration_prepare(design_uncertainty = TRUE) (the default since
# 1.2.0), the propagated design-GP variance.
# The predictive is the equal-weight mixture over draws; we summarise it by
# its mean and variance (law of total variance) and use a moment-matched
# Gaussian for intervals and CRPS.

.input_gp_terms <- function(input, input_new, phi, gp_sigma2, gp_tau2, V_chol) {
  input <- as.matrix(input)
  input_new <- as.matrix(input_new)
  N <- nrow(input)
  dd <- as.matrix(stats::dist(rbind(input, input_new)))
  J <- gp_sigma2 * exp(-phi * dd[seq_len(N), N + seq_len(nrow(input_new)), drop = FALSE])
  A <- backsolve(V_chol, forwardsolve(t(V_chol), J))
  list(A = A, q = gp_sigma2 + gp_tau2 - colSums(J * A))
}

# Default AR design rows for held-out runs, built exactly as emulator_predict()
# does internally. spDBL 1.1.1's AR(2) initialisation drew random numbers, so
# a seed is set; since 1.2.0 the designs are deterministic and it has no effect.
spdbl_known_designs <- function(emulator, dt_test, seed) {
  s <- emulator$setup
  n_test <- nrow(dt_test[[1]])
  index <- spDBL::generate.grid.rowsnake(n_test, s$fncol, n_test, s$bncol)
  set.seed(seed)
  if (s$AR_choice == 2)
    spDBL::gen_F_ls_AR2_EP(dt_test, s$nT_ori, s$n_b, index) else
    spDBL::gen_F_ls_AR1_EP(dt_test, s$nT_ori, s$n_b, index)
}

# Returns arrays [n_test, S, T] of mean, sd, 95% lower and upper bounds.
spdbl_known_predictive <- function(emulator, input_new, F_new_ls) {
  mean_list <- spDBL::emulator_predict(emulator, input_new, F_new_ls = F_new_ls)
  s <- emulator$setup
  ff <- emulator$para_ffbs$ff
  St <- emulator$para_ffbs$bs$St
  E <- s$nT
  C <- s$bncol
  g <- .input_gp_terms(s$pde_para_train, input_new, s$phi_para, s$gp_sigma2,
                       s$gp_tau2, chol(s$V))
  nt <- ff[[E]]$nt
  if (s$covariance == "iw") {
    df <- nt - C + 1
    column_scale <- diag(as.matrix(ff[[E]]$Dt)) / df
  } else {
    df <- 2 * nt
    column_scale <- diag(as.matrix(s$R_episode)) * as.numeric(ff[[E]]$Dt) / nt
  }
  p <- dim(St)[1]
  scale2 <- lapply(seq_len(E), function(e) {
    H <- as.matrix(F_new_ls[[e]]) - crossprod(g$A, as.matrix(s$F_ls[[e]]))
    u <- rowSums((H %*% matrix(St[, , e], p, p)) * H) + g$q
    outer(u, column_scale)
  })
  scale2 <- spDBL::recover_from_EP_exact(scale2, s$nT_ori, E, block_cols = s$block_cols)
  mean <- list_to_array(mean_list)
  scale <- sqrt(list_to_array(scale2))
  half <- stats::qt(0.975, df) * scale
  list(mean = mean, sd = scale * sqrt(df / (df - 2)), lower = mean - half,
       upper = mean + half, df = df)
}

spdbl_designgp_predictive <- function(cal, input_new) {
  input_new <- as.matrix(input_new)
  n <- nrow(input_new)
  M <- cal$M
  first <- second <- array(0, c(n, cal$S, cal$T))
  # Draw-outer order reuses the package's one-draw residual cache.
  for (m in seq_len(M)) for (i in seq_len(n)) {
    p <- spDBL:::predict_calibration(input_new[i, ], m, cal, covariance = "marginal")
    first[i, , ] <- first[i, , ] + t(p$mu) / M
    second[i, , ] <- second[i, , ] + t(p$variance + p$mu^2) / M
  }
  sd <- sqrt(pmax(second - first^2, 0))
  list(mean = first, sd = sd, lower = first - stats::qnorm(0.975) * sd,
       upper = first + stats::qnorm(0.975) * sd, n_draws = M,
       design_uncertainty = isTRUE(cal$design_uncertainty))
}

# Active coordinates rescaled to [0, 1] by the box; fixed coordinates -> 0.
unit_inputs <- function(X, lower, upper) {
  X <- as.matrix(X)
  active <- upper > lower
  out <- X
  out[, active] <- sweep(sweep(X[, active, drop = FALSE], 2, lower[active], "-"),
                         2, (upper - lower)[active], "/")
  out[, !active] <- 0
  out
}
