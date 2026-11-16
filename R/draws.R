# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Coefficient draws, stored or regenerated on demand.
#
# By default the samplers store every episode's p x C x nsam coefficient draws,
# which dominates the memory of a fitted emulator on large grids. With
# keep_draws = FALSE they keep instead the random-number state at the start of
# each episode's draws (2.5 KB per episode), and .regenerate_episode() replays
# that state with the samplers' own arithmetic, so regenerated draws are
# identical to stored ones. This needs marginal sampling: a jointly sampled
# episode depends on its successor's draws.

.lazy_draws <- function(type, Sigma, states, p, C, nsam, L_R = NULL)
  structure(list(Sigma = Sigma, draws = list(type = type, states = states, p = p,
    C = C, nsam = nsam, E = length(states), L_R = L_R)),
    class = "spDBL_lazy_draws", joint_trajectories = FALSE,
    sampling_method = "marginal")

#' @export
`[[.spDBL_lazy_draws` <- function(x, i, ...) {
  if (is.numeric(i))
    stop("Coefficient draws are not stored (keep_draws = FALSE); refit with ",
         "keep_draws = TRUE to index them", call. = FALSE)
  .subset2(x, i)
}

#' @export
length.spDBL_lazy_draws <- function(x) x$draws$E

#' @export
print.spDBL_lazy_draws <- function(x, ...) {
  cat(sprintf("spDBL coefficient draws, regenerated on demand: %d episodes x %d draws (%d x %d)\n",
    x$draws$E, x$draws$nsam, x$draws$p, x$draws$C))
  invisible(x)
}

# The smoother moments the draws are regenerated from.
.draw_moments <- function(para_ffbs) {
  ff <- para_ffbs$ff
  E <- length(ff) - as.integer(!is.null(ff$prior))
  list(st = para_ffbs$bs$st, St = para_ffbs$bs$St, mtE = ff[[E]]$mt, MtE = ff[[E]]$Mt)
}

.save_rng <- function() {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  function() {
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }
}

.rng_state <- function() get(".Random.seed", envir = globalenv())

# Episode e's p x C x nsam draws, exactly as FFBS_sampling(),
# FFBS_sampling_sigma2R() and sample_calibration_emulator() drew them.
.regenerate_episode <- function(lazy, moments, e) {
  d <- lazy$draws
  p <- d$p
  C <- d$C
  restore <- .save_rng()
  on.exit(restore(), add = TRUE)
  assign(".Random.seed", d$states[[e]], envir = globalenv())
  if (e == d$E) {
    mean <- moments$mtE
    L <- t(.cal_spd(moments$MtE))
  } else {
    mean <- matrix(moments$st[, , e], p, C)
    H <- matrix(moments$St[, , e], p, p)
    L <- t(.cal_spd((H + t(H)) / 2))
  }
  theta <- array(NA_real_, c(p, C, d$nsam))
  if (identical(d$type, "sigma2R")) {
    L_R <- d$L_R
    for (m in seq_len(d$nsam))
      theta[, , m] <- mean + sqrt(lazy$Sigma[m]) * L %*% matrix(stats::rnorm(p * C), p, C) %*% t(L_R)
  } else {
    for (m in seq_len(d$nsam)) {
      L_Sigma <- t(.cal_spd(matrix(lazy$Sigma[, , m], C, C)))
      theta[, , m] <- mean + L %*% matrix(stats::rnorm(p * C), p, C) %*% t(L_Sigma)
    }
  }
  theta
}

# Episode e's p x C x M draws of a prepared calibration object.
.cal_theta <- function(cal, e) {
  theta <- cal$posterior$theta
  if (inherits(theta, "spDBL_lazy_draws"))
    .regenerate_episode(theta, cal$posterior$moments, e) else theta[[e]]
}

# Draw m of every episode, as a list of p x C matrices. With regenerated
# draws this regenerates every episode, so loop over episodes instead where
# many draws are needed on a large grid.
.cal_theta_draw <- function(cal, m)
  lapply(seq_len(cal$E), function(e) matrix(.cal_theta(cal, e)[, , m], cal$p, cal$C))

# emulator_predict(MC = TRUE) for an emulator with regenerated draws: the same
# arithmetic as FFBS_predict_MC(), looping over episodes so that each
# episode's draws are regenerated once.
.predict_mc_regenerated <- function(emulator, Y, F_ls, F_new_ls, input, input_new,
                                    nT, phi_para, gp_sigma2, gp_tau2) {
  if (is.null(dim(input_new))) input_new <- matrix(input_new, nrow = 1)
  input <- as.matrix(input)
  input_new <- as.matrix(input_new)
  N <- nrow(input)
  n <- nrow(input_new)
  nsam <- emulator$res_ffbs$draws$nsam
  moments <- .draw_moments(emulator$para_ffbs)
  V <- gp_sigma2 * exp(-phi_para * as.matrix(stats::dist(input))) + diag(gp_tau2, N)
  R_V <- chol(V)
  distances <- as.matrix(stats::dist(rbind(input, input_new)))
  Jt <- gp_sigma2 * exp(-phi_para * distances[seq_len(N), N + seq_len(n), drop = FALSE])
  design_at <- function(design, i, rows) {
    if (is.list(design)) design <- design[[if (length(design) == 1L) 1L else i]]
    if (is.null(dim(design))) design <- matrix(design, nrow = rows)
    as.matrix(design)
  }
  out <- vector("list", nT)
  for (e in seq_len(nT)) {
    draws <- .regenerate_episode(emulator$res_ffbs, moments, e)
    Ft <- design_at(F_ls, e, N)
    Ft_new <- design_at(F_new_ls, e, n)
    Yt <- as.matrix(Y[[e]])
    out[[e]] <- array(NA_real_, c(n, ncol(Yt), nsam))
    for (m in seq_len(nsam)) {
      th <- matrix(draws[, , m], dim(draws)[1], dim(draws)[2])
      weights <- backsolve(R_V, forwardsolve(t(R_V), Yt - Ft %*% th))
      out[[e]][, , m] <- Ft_new %*% th + crossprod(Jt, weights)
    }
  }
  names(out) <- paste0("T", seq_len(nT))
  out
}
