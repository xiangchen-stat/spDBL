# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

make_stan_calibration_data <- function(emulator, draw, z, U,
    locations = seq_len(emulator$S), n0 = 30, d0 = 2, discount = 1,
    m0 = rep(0, length(locations)), M0 = diag(length(locations))) {
  z <- as.matrix(z)
  S <- length(locations)
  if (!identical(dim(z), c(emulator$T, S)) || any(!is.finite(z)) ||
      anyDuplicated(locations) || any(!locations %in% seq_len(emulator$S)) ||
      draw < 1 || draw > emulator$M || length(draw) != 1 ||
      length(m0) != S || n0 <= 0 || d0 <= 0 || discount <= 0 || discount > 1 ||
      !identical(dim(U), c(emulator$S, emulator$S)) ||
      !identical(dim(M0), c(S, S))) stop("Invalid calibration observations or priors")
  theta <- array(0, c(emulator$E, emulator$p, emulator$C))
  residual <- array(0, c(emulator$E, nrow(emulator$input), emulator$C))
  weights <- .cal_residual_weights(emulator, draw)
  for (e in seq_len(emulator$E)) {
    theta[e, , ] <- .cal_theta(emulator, e)[, , draw]
    residual[e, , ] <- weights[[e]]
  }
  f <- emulator$predictor
  spec <- .cal_predictor_spec(f)
  N <- nrow(emulator$input)
  # Design-row covariance at each time before scaling by the GP variance
  # factor, which Stan computes at eta. Zero when design uncertainty is off.
  factors <- .cal_design_factors(emulator, draw, locations)
  design_covariance <- array(0, c(emulator$T, S, S))
  for (t in seq_along(factors)) design_covariance[t, , ] <- crossprod(factors[[t]])
  # One kernel (lengthscales and lower Cholesky factor) per group; the
  # episodes of one time step share a group.
  predictor_L <- array(0, c(spec$G, N, N))
  for (g in seq_len(spec$G)) predictor_L[g, , ] <- t(.cal_predictor_chol(f, g))
  time_group <- spec$groups[(seq_len(emulator$T) - 1L) * emulator$B + 1L]
  # Explicit dimensions keep length-one vectors from becoming Stan scalars.
  vector_data <- function(x) array(as.numeric(x), length(x))
  index_data <- function(x) array(as.integer(x), length(x))
  list(N = N, D_full = ncol(emulator$input), D = length(f$active),
    T = emulator$T, S_full = emulator$S, S = S, B = emulator$B, C = emulator$C,
    P = emulator$p, E = emulator$E, active = index_data(f$active),
    eta_lower = vector_data(f$lower), eta_upper = vector_data(f$upper),
    input = emulator$input, predictor_x = f$x,
    G = spec$G, time_group = index_data(time_group),
    design_kernel = if (spec$kernel == "matern52") 2L else 1L,
    lengthscale = spec$lengthscale,
    predictor_center = vector_data(f$center), predictor_alpha = f$alpha,
    predictor_L = predictor_L, design_covariance = design_covariance,
    theta = theta, residual_weights = residual, block_cols = emulator$block_cols,
    locations = index_data(locations), phi_para = emulator$phi_para,
    gp_sigma2 = emulator$gp_sigma2, gp_tau2 = emulator$gp_tau2,
    sigma2 = emulator$posterior$sigma2[draw], L_V = t(emulator$V_chol),
    R_field = .cal_R_field(emulator, locations),
    L_U = t(.cal_spd(U[locations, locations, drop = FALSE], "U")),
    L_M0 = t(.cal_spd(M0, "M0")), m0 = vector_data(m0), z = z,
    n0 = n0, d0 = d0, discount = discount)
}

# Independent R density for checking Stan's target in its noncentered variables.
calibration_log_density <- function(eta_unit, precision, u_std, emulator,
    emulator_draw, z, U, locations = seq_len(emulator$S), n0 = 30, d0 = 2,
    discount = 1, m0 = rep(0, length(locations)), M0 = diag(length(locations))) {
  eta <- emulator$predictor$lower
  active <- emulator$predictor$active
  eta[active] <- eta[active] + eta_unit *
    (emulator$predictor$upper[active] - eta[active])
  pred <- predict_calibration(eta, emulator_draw, emulator, locations)
  tau2 <- 1 / precision
  u <- matrix(0, emulator$T + 1, length(locations))
  u[1, ] <- m0 + sqrt(tau2[1]) * t(chol(M0)) %*% u_std[1, ]
  L_U <- t(chol(U[locations, locations, drop = FALSE]))
  for (t in seq_len(emulator$T))
    u[t + 1, ] <- u[t, ] + sqrt(tau2[t + 1]) * L_U %*% u_std[t + 1, ]
  log_p <- sum(stats::dgamma(precision, shape = discount^(0:emulator$T) * n0,
                      rate = discount^(0:emulator$T) * d0, log = TRUE)) +
    sum(stats::dnorm(u_std, log = TRUE))
  for (t in seq_len(emulator$T)) {
    covariance <- .cal_time_covariance(pred, t) + diag(tau2[t + 1], length(locations))
    L <- chol(covariance)
    residual <- forwardsolve(t(L), z[t, ] - pred$mu[t, ] - u[t + 1, ])
    log_p <- log_p - length(locations) / 2 * log(2 * pi) -
      sum(log(diag(L))) - sum(residual^2) / 2
  }
  log_p
}

# Exact Gaussian filtering after integrating out the independent emulator field.
.cal_filter_discrepancy <- function(pred, tau2_all, z, U, m0, M0) {
  T <- nrow(z)
  S <- ncol(z)
  means <- matrix(0, T + 1L, S)
  means[1, ] <- m0
  covariances <- vector("list", T + 1L)
  priors <- vector("list", T)
  covariances[[1]] <- tau2_all[1] * M0
  log_likelihood <- 0
  for (t in seq_len(T)) {
    P <- covariances[[t]] + tau2_all[t + 1L] * U
    observation_covariance <- P + .cal_time_covariance(pred, t) + diag(tau2_all[t + 1L], S)
    L <- t(.cal_spd(observation_covariance))
    residual <- z[t, ] - pred$mu[t, ] - means[t, ]
    whitened <- forwardsolve(L, residual)
    whitened_prior <- forwardsolve(L, P)
    means[t + 1L, ] <- means[t, ] + crossprod(whitened_prior, whitened)
    C <- P - crossprod(whitened_prior)
    covariances[[t + 1L]] <- (C + t(C)) / 2
    priors[[t]] <- P
    log_likelihood <- log_likelihood - S / 2 * log(2 * pi) -
      sum(log(diag(L))) - sum(whitened^2) / 2
  }
  list(mean = means, covariance = covariances, prior_covariance = priors,
       log_likelihood = log_likelihood)
}

calibration_marginal_log_density <- function(eta_unit, precision, emulator,
    emulator_draw, z, U, locations = seq_len(emulator$S), n0 = 30, d0 = 2,
    discount = 1, m0 = rep(0, length(locations)), M0 = diag(length(locations))) {
  eta <- emulator$predictor$lower
  active <- emulator$predictor$active
  eta[active] <- eta[active] + eta_unit *
    (emulator$predictor$upper[active] - eta[active])
  pred <- predict_calibration(eta, emulator_draw, emulator, locations)
  filtered <- .cal_filter_discrepancy(pred, 1 / precision, as.matrix(z),
    U[locations, locations, drop = FALSE], m0, M0)
  sum(stats::dgamma(precision, shape = discount^(0:emulator$T) * n0,
    rate = discount^(0:emulator$T) * d0, log = TRUE)) + filtered$log_likelihood
}

# Joint backward draws of u0,...,uT, conditional on the marginal Stan samples.
# This is Gaussian conditioning, not a second MCMC approximation.
recover_calibration_discrepancy <- function(eta, tau2, tau2_0, emulator_draw,
    emulator, z, U, locations = seq_len(emulator$S),
    m0 = rep(0, length(locations)), M0 = diag(length(locations)), seed = 1234) {
  restore_seed <- .cal_seed(seed)
  on.exit(restore_seed(), add = TRUE)
  n <- nrow(eta)
  T <- emulator$T
  S <- length(locations)
  u <- array(NA_real_, c(n, T, S))
  u0 <- matrix(NA_real_, n, S)
  U_observed <- U[locations, locations, drop = FALSE]
  for (i in seq_len(n)) {
    pred <- predict_calibration(eta[i, ], emulator_draw[i], emulator, locations)
    filtered <- .cal_filter_discrepancy(pred, c(tau2_0[i], tau2[i, ]),
                                        as.matrix(z), U_observed, m0, M0)
    states <- matrix(NA_real_, T + 1L, S)
    states[T + 1L, ] <- filtered$mean[T + 1L, ] +
      t(.cal_spd(filtered$covariance[[T + 1L]])) %*% stats::rnorm(S)
    for (t in seq.int(T, 1L)) {
      P <- filtered$prior_covariance[[t]]
      C <- filtered$covariance[[t]]
      B <- C %*% chol2inv(.cal_spd(P))
      H <- C - B %*% P %*% t(B)
      mean <- filtered$mean[t, ] + B %*% (states[t + 1L, ] - filtered$mean[t, ])
      states[t, ] <- mean + t(.cal_spd((H + t(H)) / 2)) %*% stats::rnorm(S)
    }
    u[i, , ] <- states[-1, , drop = FALSE]
    u0[i, ] <- states[1, ]
  }
  list(u = u, u0 = u0)
}

.cal_extract_array <- function(samples, name, dims) {
  indices <- do.call(expand.grid, lapply(dims, seq_len))
  columns <- apply(indices, 1, function(i) paste0(name, "[", paste(i, collapse = ","), "]"))
  array(samples[, columns, drop = FALSE], c(nrow(samples), dims))
}

# Draw y conditional on eta/u/tau, then independent pointwise observation noise.
recover_calibration_predictions <- function(eta, tau2, u, emulator_draw,
    emulator, z, locations = seq_len(emulator$S), seed = 1234) {
  restore_seed <- .cal_seed(seed)
  on.exit(restore_seed(), add = TRUE)
  n <- nrow(eta)
  T <- emulator$T
  S <- length(locations)
  y <- z_rep <- array(NA_real_, c(n, T, S))
  for (i in seq_len(n)) {
    pred <- predict_calibration(eta[i, ], emulator_draw[i], emulator, locations)
    for (t in seq_len(T)) {
      prior <- .cal_time_covariance(pred, t)
      covariance <- prior + diag(tau2[i, t], S)
      K <- prior %*% chol2inv(chol(covariance))
      mean <- pred$mu[t, ] + K %*% (z[t, ] - u[i, t, ] - pred$mu[t, ])
      variance <- prior - K %*% prior
      y[i, t, ] <- as.numeric(mean + t(.cal_spd((variance + t(variance)) / 2)) %*% stats::rnorm(S))
      z_rep[i, t, ] <- y[i, t, ] + u[i, t, ] + stats::rnorm(S, sd = sqrt(tau2[i, t]))
    }
  }
  list(y = y, z_rep = z_rep)
}

# Draws of the emulator field y, the discrepancy u and new observations z at
# prediction targets, at the observation times. Each draw conditions on the
# same posterior draw of y and u at the observed sites (composition sampling:
# Appendix B of the JMLR MLOSS paper, with the calibration's own field
# covariance).
# - Targets: each is a weighted sum of up to `ncol(cells)` grid cells (cells,
#   weights: targets x slots; NA cells have weight 0). An on-grid site is one
#   cell with weight 1; an off-grid point is its containing cell or the
#   bilinear interpolation of the surrounding cell centres.
# - y: at each time, y at the observed sites and the cells is Gaussian with
#   the mean and covariance of predict_calibration(), design-row term
#   included. Observed cells contribute their drawn values; the rest of each
#   target is drawn given y at the observed sites.
# - u: u0 at the targets has prior N(m0, tau2_0 M0), independent of u0 at the
#   observed sites, and each increment at the targets is drawn given the
#   increment at the observed sites under the discrepancy covariance (U_oo at
#   the observed sites, U_oP to the targets, U_PP at the targets: a matrix
#   when joint, the diagonal otherwise). This extends the prior without
#   changing the calibration posterior.
# - z = y + u + N(0, tau2_t) at each target.
# joint = TRUE draws the targets jointly; FALSE draws each from its own
# conditional distribution (targets independent; each u path stays a path).
.cal_target_draws <- function(eta, tau2, tau2_0, u, u0, y, emulator_draw,
    emulator, locations, cells, weights, U_oo, U_oP, U_PP, m0, M0, joint, seed) {
  restore_seed <- .cal_seed(seed)
  on.exit(restore_seed(), add = TRUE)
  n <- nrow(eta)
  T <- emulator$T
  K <- nrow(cells)
  nslot <- ncol(cells)
  draw <- if (joint) function(L, x) as.numeric(L %*% x) else function(L, x) L * x
  # Discrepancy: increments at the targets given those at the observed sites.
  L_U <- .cal_spd(U_oo, "U")
  A_U <- forwardsolve(t(L_U), U_oP)
  gamma_u <- backsolve(L_U, A_U)
  L_u <- if (joint) .cal_psd_factor(U_PP - crossprod(A_U),
                                    "Conditional discrepancy covariance") else
    sqrt(pmax(U_PP - colSums(A_U^2), 0))
  L_0 <- if (joint) .cal_psd_factor(M0, "M0") else sqrt(diag(M0))
  # Emulator field: split each target into observed cells (known given the
  # draw) and unobserved cells (drawn), slot by slot.
  pos <- matrix(match(cells, locations), K)
  w_obs <- ifelse(!is.na(pos), weights, 0)
  w_new <- ifelse(is.na(pos) & !is.na(cells), weights, 0)
  new_cells <- unique(cells[w_new != 0])
  idx <- matrix(match(cells, new_cells), K)
  idx[w_new == 0] <- NA
  o <- seq_along(locations)
  s <- length(locations) + seq_along(new_cells)
  # Sum over slots of f(slot) * weight, with f giving one column per target.
  combine <- function(f, rows) {
    out <- matrix(0, rows, K)
    for (m in seq_len(nslot)) {
      use <- which(w_new[, m] != 0)
      if (length(use)) out[, use] <- out[, use] + sweep(f(idx[use, m]), 2, w_new[use, m], "*")
    }
    out
  }
  R_oo <- .cal_R_field(emulator, locations)
  R_oP <- if (length(new_cells))
    combine(function(i) .cal_R_cross(emulator, locations, new_cells[i]), length(o)) else
    matrix(0, length(o), K)
  R_PP <- if (joint) matrix(0, K, K) else numeric(K)
  for (m in seq_len(nslot)) for (k in seq_len(nslot)) {
    a <- which(w_new[, m] != 0)
    b <- which(w_new[, k] != 0)
    if (!length(a) || !length(b)) next
    if (joint) {
      R_PP[a, b] <- R_PP[a, b] + outer(w_new[a, m], w_new[b, k]) *
        .cal_R_cross(emulator, new_cells[idx[a, m]], new_cells[idx[b, k]])
    } else {
      both <- intersect(a, b)
      R_PP[both] <- R_PP[both] + w_new[both, m] * w_new[both, k] *
        .cal_R_pairs(emulator, new_cells[idx[both, m]], new_cells[idx[both, k]])
    }
  }
  y_s <- u_s <- z_s <- array(NA_real_, c(n, T, K))
  for (i in seq_len(n)) {
    current <- m0 + sqrt(tau2_0[i]) * draw(L_0, stats::rnorm(K))
    previous <- u0[i, ]
    for (t in seq_len(T)) {
      current <- current + as.numeric(crossprod(gamma_u, u[i, t, ] - previous)) +
        sqrt(tau2[i, t]) * draw(L_u, stats::rnorm(K))
      previous <- u[i, t, ]
      u_s[i, t, ] <- current
    }
    pred <- predict_calibration(eta[i, ], emulator_draw[i], emulator,
                                c(locations, new_cells), covariance = "factors")
    for (t in seq_len(T)) {
      known <- as.numeric(rowSums(w_obs * matrix(y[i, t, ][pos], K), na.rm = TRUE))
      mu_P <- as.numeric(combine(function(j) matrix(pred$mu[t, s[j]], 1), 1))
      S_oo <- pred$scale * R_oo
      S_oP <- pred$scale * R_oP
      S_PP <- pred$scale * R_PP
      if (length(pred$design_factors)) {
        Z <- pred$design_factors[[t]]
        v <- pred$design_variance[t]
        Z_P <- combine(function(j) Z[, s[j], drop = FALSE], nrow(Z))
        S_oo <- S_oo + v * crossprod(Z[, o, drop = FALSE])
        S_oP <- S_oP + v * crossprod(Z[, o, drop = FALSE], Z_P)
        S_PP <- S_PP + v * (if (joint) crossprod(Z_P) else colSums(Z_P^2))
      }
      L <- t(.cal_spd(S_oo, "Emulator field covariance"))
      A <- forwardsolve(L, S_oP)
      mean <- known + mu_P + as.numeric(crossprod(A, forwardsolve(L, y[i, t, ] - pred$mu[t, o])))
      L_y <- if (joint) .cal_psd_factor(S_PP - crossprod(A),
                                        "Conditional emulator field covariance") else
        sqrt(pmax(S_PP - colSums(A^2), 0))
      y_s[i, t, ] <- mean + draw(L_y, stats::rnorm(K))
      z_s[i, t, ] <- y_s[i, t, ] + u_s[i, t, ] + stats::rnorm(K, sd = sqrt(tau2[i, t]))
    }
  }
  list(y = y_s, u = u_s, z = z_s)
}

# Grid sites without observations: each target is one cell, and the
# discrepancy covariance comes from the full-grid matrix U.
recover_site_predictions <- function(eta, tau2, tau2_0, u, u0, y,
    emulator_draw, emulator, U, locations, sites,
    m0_sites = rep(0, length(sites)), M0_sites = diag(length(sites)),
    joint = FALSE, seed = 1234) {
  .cal_target_draws(eta, tau2, tau2_0, u, u0, y, emulator_draw, emulator,
    locations, cells = matrix(as.integer(sites)), weights = matrix(1, length(sites)),
    U_oo = U[locations, locations, drop = FALSE],
    U_oP = U[locations, sites, drop = FALSE],
    U_PP = if (joint) U[sites, sites, drop = FALSE] else diag(U)[sites],
    m0 = m0_sites, M0 = M0_sites, joint = joint, seed = seed)
}

# Points anywhere on the grid's domain, in grid coordinates (cell centres at
# the grid's integer coordinates). field = "cell" takes the emulated field of
# the cell containing each point; "bilinear" interpolates the four
# surrounding cell centres (at the edge, the nearest ones). The discrepancy
# covariance at the observed sites is the fitted U; kernel(a, b), C at the
# rows of two coordinate matrices, gives it at and to the points themselves.
recover_point_predictions <- function(eta, tau2, tau2_0, u, u0, y,
    emulator_draw, emulator, U, kernel, locations, coords, field = c("cell", "bilinear"),
    m0_points = rep(0, nrow(coords)), M0_points = diag(nrow(coords)),
    joint = FALSE, seed = 1234) {
  field <- match.arg(field)
  coords <- as.matrix(coords)
  grid <- .cal_grid_coords(emulator)
  map <- .cal_point_cells(grid, coords, field)
  K <- nrow(coords)
  grid_o <- grid[locations, , drop = FALSE]
  .cal_target_draws(eta, tau2, tau2_0, u, u0, y, emulator_draw, emulator,
    locations, cells = map$cells, weights = map$weights,
    U_oo = U[locations, locations, drop = FALSE],
    U_oP = .cal_kernel_matrix(kernel, grid_o, coords),
    U_PP = if (joint) .cal_kernel_matrix(kernel, coords, coords) else
      vapply(seq_len(K), function(k)
        .cal_kernel_matrix(kernel, coords[k, , drop = FALSE], coords[k, , drop = FALSE])[1, 1],
        numeric(1)),
    m0 = m0_points, M0 = M0_points, joint = joint, seed = seed)
}

# Each fit has one fixed emulator draw. Equal pooling approximates Eq. (29),
# rather than fitting a mixture likelihood that reweights the emulator draws.
fit_calibration_rstan <- function(emulator, z, U,
    model_file = NULL, draws = seq_len(emulator$M),
    locations = seq_len(emulator$S), chains = 4, iter = 2000, warmup = 1000,
    draws_per_emulator = 200, seed = 1234, cores = 1,
    keep_predictions = TRUE, n0 = 30, d0 = 2, discount = 1,
    m0 = rep(0, length(locations)), M0 = diag(length(locations)),
    control = list(adapt_delta = 0.95, max_treedepth = 12), compiled_model = NULL,
    marginal_discrepancy = TRUE, refresh = 100, checkpoint_dir = NULL,
    keep_discrepancy = TRUE, init = "midpoint", mode_starts = 20L,
    parallel_draws = 1L) {
  if (!requireNamespace("rstan", quietly = TRUE)) stop("Install rstan first")
  if (warmup >= iter || chains < 1 || draws_per_emulator < 1 ||
      draws_per_emulator > chains * (iter - warmup) || !length(draws) ||
      anyDuplicated(draws) || any(!draws %in% seq_len(emulator$M))) stop("Invalid sampling settings")
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  model <- if (is.null(compiled_model))
    .cal_stan_model(marginal_discrepancy, cache = FALSE) else compiled_model
  # Fits run at the same time (forked) run their chains one after another, so
  # that no process pool is started inside a forked process. rstan's draws do
  # not depend on the number of cores, so results are unchanged.
  forked <- parallel_draws > 1L && .Platform$OS.type != "windows" && length(draws) > 1L
  chain_cores <- if (forked) 1L else cores
  if (forked && cores > 1L)
    message("parallel_draws > 1: each conditional fit runs its chains in sequence ",
            "(cores = ", cores, " is not used); ", parallel_draws, " fits run at once")
  if (forked) .warn_threaded_blas()
  # One conditional fit per emulator draw. Every random number a fit uses comes
  # from seeds derived from seed and j, so fits are independent and can run
  # in any order or at the same time.
  fit_one <- function(j) {
    m <- draws[j]
    message("Calibration conditional on emulator draw ", m, " (", j, "/", length(draws), ")")
    data <- make_stan_calibration_data(emulator, m, z, U, locations, n0, d0,
                                       discount, m0, M0)
    # Starting eta on the unit scale, one row per chain. init is "mode" (set
    # after the checkpoint check below), "midpoint", "random" (uniform in the
    # box, seeded per conditional fit), or a unit vector from .cal_init_spec().
    mode_init <- identical(init, "mode")
    precision_start <- matrix(n0 / d0, chains, data$T + 1)
    eta_start <- if (mode_init) NULL else if (identical(init, "random")) {
      set.seed(seed + j)
      matrix(stats::runif(chains * data$D), chains, data$D)
    } else matrix(if (identical(init, "midpoint")) 0.5 else init,
                  chains, data$D, byrow = TRUE)
    # A fit saved by an earlier run with identical data and settings is
    # reused, so an interrupted calibration resumes where it stopped.
    checkpoint <- if (!is.null(checkpoint_dir))
      file.path(checkpoint_dir, sprintf("conditional_%03d.rds", m))
    fingerprint <- if (!is.null(checkpoint)) rlang::hash(list(data = data,
      eta_start = if (mode_init) list("mode", mode_starts) else eta_start,
      chains = chains, iter = iter, warmup = warmup,
      seed = seed + j, draws_per_emulator = draws_per_emulator, control = control,
      marginal_discrepancy = marginal_discrepancy,
      program = tryCatch(model@model_code, error = function(e) NULL)))
    if (!is.null(checkpoint) && file.exists(checkpoint)) {
      saved <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
      if (identical(saved$fingerprint, fingerprint)) {
        message("  reusing the saved fit ", basename(checkpoint))
        return(list(samples = saved$samples, diagnostics = saved$diagnostics))
      }
      message("  ", basename(checkpoint), " was saved with other data or settings; refitting")
    }
    # init = "mode": L-BFGS (rstan::optimizing) from mode_starts random starts;
    # the chains start, slightly jittered, at the optimum with the highest
    # log posterior density. Chains started far from it can otherwise settle
    # in local modes of negligible posterior mass.
    mode_search <- NULL
    if (mode_init) {
      set.seed(seed + j)
      starts <- matrix(stats::runif(mode_starts * data$D), mode_starts, data$D)
      optima <- lapply(seq_len(mode_starts), function(k) {
        values <- list(eta_unit = array(starts[k, ], data$D),
                       precision = rep(n0 / d0, data$T + 1))
        if (!marginal_discrepancy) values$u_std <- matrix(0, data$T + 1, data$S)
        tryCatch(rstan::optimizing(model, data = data, init = values,
            seed = seed + 1000L * j + k, algorithm = "LBFGS", iter = 2000L,
            as_vector = FALSE, refresh = 0L),
          error = function(e) NULL)
      })
      values <- vapply(optima, function(o)
        if (is.null(o) || !identical(as.integer(o$return_code), 0L)) -Inf else o$value, numeric(1))
      if (!any(is.finite(values)))
        stop("The mode search failed from every start for emulator draw ", m)
      best <- optima[[which.max(values)]]
      eta_mode <- pmin(pmax(as.numeric(best$par$eta_unit), 1e-6), 1 - 1e-6)
      precision_mode <- as.numeric(best$par$precision)
      eta_start <- t(vapply(seq_len(chains), function(ch)
        stats::plogis(stats::qlogis(eta_mode) + stats::rnorm(data$D, sd = 0.05)), numeric(data$D)))
      eta_start <- matrix(eta_start, chains, data$D)
      precision_start <- t(vapply(seq_len(chains), function(ch)
        precision_mode * exp(stats::rnorm(data$T + 1, sd = 0.05)), numeric(data$T + 1)))
      precision_start <- matrix(precision_start, chains, data$T + 1)
      mode_search <- list(log_density = values, eta_unit = eta_mode,
        near_best = sum(values > max(values) - 1))
    }
    # local() keeps this fit's values; the loop reassigns eta_start and data.
    start <- local({
      eta <- eta_start
      precision <- precision_start
      D <- data$D
      T <- data$T
      S <- data$S
      function(chain_id = 1) {
        values <- list(eta_unit = array(eta[chain_id, ], D),
                       precision = precision[chain_id, ])
        if (!marginal_discrepancy) values$u_std <- matrix(0, T + 1, S)
        values
      }
    })
    fit <- rstan::sampling(model, data = data, init = start, chains = chains,
      iter = iter, warmup = warmup, seed = seed + j, cores = chain_cores, refresh = refresh,
      pars = c("eta", "tau2", if (!marginal_discrepancy) "u"), control = control)
    raw <- rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE)
    if (length(dim(raw)) != 3L || !length(raw))
      stop("Stan returned no calibration samples for emulator draw ", m)
    samples <- matrix(raw, nrow = dim(raw)[1] * dim(raw)[2],
                      dimnames = list(NULL, dimnames(raw)[[3]]))
    set.seed(seed + j)
    selected <- sample(seq_len(nrow(samples)), draws_per_emulator)
    eta <- .cal_extract_array(samples, "eta", data$D_full)[selected, , drop = FALSE]
    tau <- .cal_extract_array(samples, "tau2", data$T + 1)[selected, , drop = FALSE]
    if (!marginal_discrepancy) {
      states <- .cal_extract_array(samples, "u", c(data$T + 1, data$S))
      states <- array(states[selected, , , drop = FALSE], c(length(selected), data$T + 1, data$S))
    }
    if (dim(raw)[1] < 2L) {
      # One transition per chain cannot support Rhat or effective sample size.
      parameter_names <- c(paste0("eta[", data$active, "]"),
        grep("^tau2\\[", colnames(samples), value = TRUE),
        if (!marginal_discrepancy) grep("^u\\[", colnames(samples), value = TRUE))
      summary <- cbind(mean = colMeans(samples[, parameter_names, drop = FALSE]),
        n_eff = NA_real_, Rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_)
    } else {
      summary <- rstan::summary(fit, pars = c(paste0("eta[", data$active, "]"),
        "tau2", if (!marginal_discrepancy) "u"))$summary
      # Rank-normalised split Rhat and bulk and tail effective sample sizes
      # (Vehtari et al., 2021) replace the classic diagnostics that
      # rstan::summary() reports.
      summary <- cbind(summary[, setdiff(colnames(summary), "Rhat"), drop = FALSE],
        .cal_convergence(raw[, , rownames(summary), drop = FALSE]))
    }
    sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    diagnostics_j <- list(emulator_draw = m, summary = summary, mode_search = mode_search,
      divergences = sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1))),
      max_treedepth_hits = sum(vapply(sampler, function(x)
        sum(x[, "treedepth__"] >= if (is.null(control$max_treedepth)) 10 else control$max_treedepth), numeric(1))),
      max_rhat = .cal_extreme(summary[, "Rhat"], max),
      min_n_eff = .cal_extreme(summary[, "ess_bulk"], min),
      min_ess_tail = .cal_extreme(summary[, "ess_tail"], min))
    samples_j <- list(eta = eta, tau2 = tau[, -1, drop = FALSE], tau2_0 = tau[, 1],
      chain = rep(seq_len(dim(raw)[2]), each = dim(raw)[1])[selected],
      iteration = rep(seq_len(dim(raw)[1]), dim(raw)[2])[selected], emulator_draw = rep(m, length(selected)))
    if (!marginal_discrepancy) {
      samples_j$u <- array(states[, -1, , drop = FALSE], c(length(selected), data$T, data$S))
      samples_j$u0 <- matrix(states[, 1, ], length(selected), data$S)
    }
    if (!is.null(checkpoint)) saveRDS(list(samples = samples_j,
      diagnostics = diagnostics_j, marginal_discrepancy = marginal_discrepancy,
      chains = chains, iter = iter, warmup = warmup, init = init,
      fingerprint = fingerprint), checkpoint)
    list(samples = samples_j, diagnostics = diagnostics_j)
  }
  # parallel_draws > 1 runs that many fits at once in forked processes (not
  # on Windows, where fits run one after another); each fit's chains use
  # `cores` processes of their own.
  fits <- if (forked) {
    parallel::mclapply(seq_along(draws), fit_one, mc.cores = parallel_draws,
                       mc.preschedule = FALSE)
  } else lapply(seq_along(draws), fit_one)
  failed <- vapply(fits, function(f) inherits(f, "try-error") || is.null(f), logical(1))
  if (any(failed)) {
    reason <- fits[[which(failed)[1]]]
    stop("Calibration failed for emulator draw ", draws[which(failed)[1]], ": ",
         if (inherits(reason, "try-error")) conditionMessage(attr(reason, "condition")) else
           "the worker process returned nothing (out of memory?)", call. = FALSE)
  }
  out <- lapply(fits, `[[`, "samples")
  diagnostics <- lapply(fits, `[[`, "diagnostics")
  # Warnings are raised here, so they also reach the caller from forked fits.
  for (d in diagnostics)
    if (isTRUE(d$divergences > 0) || (!is.na(d$max_rhat) && d$max_rhat > 1.01))
      warning("Emulator draw ", d$emulator_draw,
              " needs additional convergence checks; inspect diagnostics")
  bind_arrays <- function(name, dims) {
    # rbind flattened trajectories to keep the first dimension as posterior draws.
    flat <- do.call(rbind, lapply(out, function(o) matrix(o[[name]], nrow(o$eta))))
    array(flat, c(nrow(flat), dims))
  }
  result <- list(mcmc_eta = do.call(rbind, lapply(out, `[[`, "eta")),
    mcmc_tau2 = do.call(rbind, lapply(out, `[[`, "tau2")),
    tau2_0 = unlist(lapply(out, `[[`, "tau2_0")),
    emulator_draw = unlist(lapply(out, `[[`, "emulator_draw")),
    chain = unlist(lapply(out, `[[`, "chain")), iteration = unlist(lapply(out, `[[`, "iteration")),
    diagnostics = diagnostics, locations = locations,
    predictor_uncertainty = emulator$predictor$uncertainty,
    marginal_discrepancy = marginal_discrepancy)
  if (marginal_discrepancy && (keep_discrepancy || keep_predictions)) {
    message("Recover discrepancy trajectories by Gaussian backward sampling")
    states <- recover_calibration_discrepancy(result$mcmc_eta, result$mcmc_tau2,
      result$tau2_0, result$emulator_draw, emulator, z, U, locations, m0, M0, seed + 1000L)
    result$mcmc_u <- states$u
    result$u0 <- states$u0
  } else if (!marginal_discrepancy && (keep_discrepancy || keep_predictions)) {
    result$mcmc_u <- bind_arrays("u", c(emulator$T, length(locations)))
    result$u0 <- do.call(rbind, lapply(out, `[[`, "u0"))
  }
  if (keep_predictions) {
    predictions <- recover_calibration_predictions(result$mcmc_eta, result$mcmc_tau2,
      result$mcmc_u, result$emulator_draw, emulator, z, locations, seed)
    result$mcmc_y <- predictions$y
    result$mcmc_z <- predictions$z_rep
  }
  result
}

# Per-parameter convergence diagnostics from an iterations x chains x
# parameters array of post-warmup draws.
.cal_convergence <- function(draws) {
  one <- function(k) {
    x <- matrix(draws[, , k], dim(draws)[1], dim(draws)[2])
    c(Rhat = rstan::Rhat(x), ess_bulk = rstan::ess_bulk(x),
      ess_tail = rstan::ess_tail(x))
  }
  out <- t(vapply(seq_len(dim(draws)[3]), one, numeric(3)))
  rownames(out) <- dimnames(draws)[[3]]
  out
}

.cal_extreme <- function(x, f) if (any(is.finite(x))) f(x[is.finite(x)]) else NA_real_


# Forking after an OpenMP-threaded BLAS (oneMKL, OpenBLAS, BLIS) has run can
# leave the forked processes waiting forever on the thread pool's lock.
.warn_threaded_blas <- function(blas = paste(La_library(), extSoftVersion()[["BLAS"]]),
                                omp_threads = Sys.getenv("OMP_NUM_THREADS")) {
  if (grepl("mkl|openblas|blis", blas, ignore.case = TRUE) && !identical(omp_threads, "1"))
    warning("parallel_draws forks R, and this R uses a multithreaded BLAS (",
            basename(strsplit(blas, " ")[[1]][1]), "); forked fits can hang. Start R with ",
            "OMP_NUM_THREADS=1 (and MKL_NUM_THREADS=1 or OPENBLAS_NUM_THREADS=1).",
            call. = FALSE)
  invisible(NULL)
}
