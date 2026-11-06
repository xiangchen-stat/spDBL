# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

.cal_spd <- function(x, name = "covariance") {
  x <- as.matrix(x)
  if (nrow(x) != ncol(x) || any(!is.finite(x)) ||
      max(abs(x - t(x))) > 1e-7 * max(1, max(abs(x)))) {
    stop(name, " must be finite, square, and symmetric")
  }
  chol((x + t(x)) / 2)
}

# Unit-amplitude ARD kernels of the design GP: squared exponential, or
# Matern 5/2 k(r) = (1 + a + a^2 / 3) exp(-a) with a = sqrt(5) r, where r is
# the lengthscale-scaled distance.
.cal_kernel <- function(x, y, lengthscale, kernel = "sqexp") {
  x <- sweep(as.matrix(x), 2, lengthscale, "/")
  y <- sweep(as.matrix(y), 2, lengthscale, "/")
  d2 <- pmax(outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y), 0)
  if (identical(kernel, "matern52")) {
    a <- sqrt(5 * d2)
    return((1 + a + a^2 / 3) * exp(-a))
  }
  exp(-0.5 * d2)
}

.cal_sqexp <- function(x, y, lengthscale) .cal_kernel(x, y, lengthscale, "sqexp")

# Kernel, episode groups and per-group hyperparameters of a design
# predictor; predictors saved before 1.2.0 have one squared-exponential
# kernel shared by all episodes.
.cal_predictor_spec <- function(object) {
  groups <- if (is.null(object$groups)) rep(1L, object$n_episode) else object$groups
  G <- max(groups)
  list(kernel = if (is.null(object$kernel)) "sqexp" else object$kernel,
       groups = groups, G = G,
       lengthscale = matrix(object$lengthscale, G),
       nugget = rep_len(object$nugget, G))
}

# Design columns (across episodes) of hyperparameter group g.
.cal_group_columns <- function(groups, g, p)
  as.vector(outer(seq_len(p), (which(groups == g) - 1L) * p, "+"))

# Separate centered/scaled responses for every design column. Columns in the
# same group of episodes share one ARD kernel and its hyperparameters, fitted
# by maximum likelihood; groups are independent GPs. The default (one group)
# is the shared kernel of spDBL <= 1.2.0; calibration_prepare() uses one
# group per time step. Calibration uses the GP mean and, unless disabled, its
# predictive covariance (see .cal_design_factors()).
fit_offline_ar2 <- function(input, F_ls, lower, upper, kernel = c("sqexp", "matern52"),
                            groups = NULL) {
  kernel <- match.arg(kernel)
  input <- as.matrix(input)
  active <- which(upper > lower)
  if (!length(active) || length(lower) != ncol(input) ||
      length(upper) != ncol(input) || any(upper < lower) ||
      any(!is.finite(c(input, lower, upper)))) stop("Invalid input bounds")
  if (any(sweep(input, 2, lower, "-") < 0) ||
      any(sweep(input, 2, upper, "-") > 0)) stop("Training inputs exceed bounds")
  if (!length(F_ls) || !all(vapply(F_ls, function(f)
      is.matrix(f) && nrow(f) == nrow(input) && all(is.finite(f)), logical(1)))) {
    stop("F_ls must contain finite training design matrices")
  }
  widths <- vapply(F_ls, ncol, integer(1))
  if (length(unique(widths)) != 1L) stop("Design widths must agree")
  if (is.null(groups)) groups <- rep(1L, length(F_ls))
  if (length(groups) != length(F_ls) || any(!is.finite(groups)) ||
      any(groups != as.integer(groups)) ||
      !identical(sort(unique(as.integer(groups))), seq_len(max(groups))))
    stop("groups must number the episodes' kernel groups 1, ..., G")
  groups <- as.integer(groups)
  x <- sweep(sweep(input[, active, drop = FALSE], 2, lower[active], "-"),
             2, upper[active] - lower[active], "/")
  N <- nrow(x)
  p <- widths[1]
  targets <- do.call(cbind, F_ls)
  center <- colMeans(targets)
  scales <- apply(targets, 2, stats::sd)
  varying <- is.finite(scales) & scales > 1e-10
  scales[!varying] <- 1
  standardized <- sweep(sweep(targets, 2, center, "-"), 2, scales, "/")
  G <- max(groups)
  starts <- c(0.2, 0.5, 1)
  lengthscale <- matrix(NA_real_, G, length(active))
  nugget <- convergence <- numeric(G)
  K_chol <- vector("list", G)
  alpha <- matrix(0, N, ncol(targets))
  loo_error <- matrix(0, N, ncol(targets))
  for (g in seq_len(G)) {
    cols <- .cal_group_columns(groups, g, p)
    use <- cols[varying[cols]]
    if (length(use)) {
      # Constant responses add no information about kernel hyperparameters.
      response_gram <- tcrossprod(standardized[, use, drop = FALSE])
      n_response <- length(use)
      objective <- function(par) {
        K <- .cal_kernel(x, x, exp(par[seq_along(active)]), kernel)
        diag(K) <- diag(K) + exp(par[length(par)])
        L <- chol(K)
        (n_response * 2 * sum(log(diag(L))) + sum(chol2inv(L) * response_gram)) /
          n_response
      }
      fits <- lapply(starts, function(s) stats::optim(
        log(c(rep(s, length(active)), 1e-4)), objective, method = "L-BFGS-B",
        lower = log(c(rep(0.03, length(active)), 1e-8)),
        upper = log(c(rep(5, length(active)), 0.2))))
      best <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]
      lengthscale[g, ] <- exp(best$par[seq_along(active)])
      nugget[g] <- exp(utils::tail(best$par, 1))
      convergence[g] <- best$convergence
    } else {
      # A group of constant columns (e.g. the first AR designs) has no GP
      # error; fixed hyperparameters keep its (zero-variance) factor defined.
      lengthscale[g, ] <- starts[2]
      nugget[g] <- 1e-4
      convergence[g] <- NA_integer_
    }
    K <- .cal_kernel(x, x, lengthscale[g, ], kernel)
    diag(K) <- diag(K) + nugget[g]
    K_chol[[g]] <- chol(K)
    Kinv <- chol2inv(K_chol[[g]])
    alpha[, cols] <- Kinv %*% sweep(standardized[, cols, drop = FALSE], 2, scales[cols], "*")
    loo_error[, cols] <- sweep(alpha[, cols, drop = FALSE], 1, diag(Kinv), "/")
  }
  out <- structure(list(x = x, active = active, lower = lower, upper = upper,
    kernel = kernel, groups = groups, center = center, alpha = alpha,
    K_chol = K_chol, lengthscale = lengthscale, nugget = nugget, scales = scales,
    n_episode = length(F_ls), p = p,
    loo_rmse = sqrt(mean(loo_error^2)), optimizer_convergence = convergence,
    uncertainty = .cal_uncertainty_note(TRUE)),
    class = "offline_ar2_predictor")
  out$column_variance <- colSums(.cal_design_whitened(out)^2) / N
  out
}

.cal_uncertainty_note <- function(design_uncertainty) {
  if (design_uncertainty)
    "Design-row GP variance is propagated within each time, independently across times" else
    "Conditional on the offline GP mean; predictor uncertainty is not integrated"
}

# Upper Cholesky factor of the design-GP kernel matrix, recomputed for objects
# saved before it was stored.
# Upper Cholesky factor of the kernel matrix of group g (recomputed for
# predictors saved without it).
.cal_predictor_chol <- function(object, g = 1L) {
  if (!is.null(object$K_chol))
    return(if (is.list(object$K_chol)) object$K_chol[[g]] else object$K_chol)
  spec <- .cal_predictor_spec(object)
  K <- .cal_kernel(object$x, object$x, spec$lengthscale[g, ], spec$kernel)
  diag(K) <- diag(K) + spec$nugget[g]
  chol(K)
}

# Whitened centred training designs W = L_K^{-1} (F - center), N x (E * p),
# with each group's columns whitened by its own kernel matrix. Within a
# group the plug-in (maximum likelihood) covariance of the design columns is
# W'W / N, and the predictive covariance of the group's design rows at eta is
# v_g(eta) W'W / N, where v_g(eta) is the group's latent variance factor.
.cal_design_whitened <- function(object) {
  spec <- .cal_predictor_spec(object)
  W <- matrix(0, nrow(object$alpha), ncol(object$alpha))
  for (g in seq_len(spec$G)) {
    cols <- .cal_group_columns(spec$groups, g, object$p)
    W[, cols] <- .cal_predictor_chol(object, g) %*% object$alpha[, cols, drop = FALSE]
  }
  W
}

# Design rows at eta: $F (one row per episode), the pointwise variance of
# every design column, the latent variance factor of every group
# ($variance_factor, length G) and of every episode ($episode_factor).
predict_offline_ar2 <- function(object, eta) {
  if (length(eta) != length(object$lower) || any(!is.finite(eta)) ||
      any(eta < object$lower) || any(eta > object$upper)) stop("eta exceeds bounds")
  spec <- .cal_predictor_spec(object)
  a <- object$active
  x <- matrix((eta[a] - object$lower[a]) / (object$upper[a] - object$lower[a]), 1)
  mean <- object$center
  variance_factor <- numeric(spec$G)
  for (g in seq_len(spec$G)) {
    cols <- .cal_group_columns(spec$groups, g, object$p)
    k <- .cal_kernel(object$x, x, spec$lengthscale[g, ], spec$kernel)
    mean[cols] <- mean[cols] + as.numeric(crossprod(k, object$alpha[, cols, drop = FALSE]))
    # Latent GP variance on the standardized scale (unit amplitude).
    whitened <- backsolve(.cal_predictor_chol(object, g), k, transpose = TRUE)
    variance_factor[g] <- max(0, 1 - sum(whitened^2))
  }
  episode_factor <- variance_factor[spec$groups]
  column_variance <- if (is.null(object$column_variance))
    colSums(.cal_design_whitened(object)^2) / nrow(object$x) else object$column_variance
  list(F = lapply(seq_len(object$n_episode), function(e)
    mean[((e - 1) * object$p + 1):(e * object$p)]),
    variance = rep(episode_factor, each = object$p) * column_variance,
    variance_factor = variance_factor, episode_factor = episode_factor)
}

# Read explicit sampling metadata, including objects saved before this option.
.cal_sampling_method <- function(sampling = NULL, joint = NULL) {
  if (identical(sampling, "independent_smoothed_marginals")) sampling <- "marginal"
  if (is.null(sampling)) {
    if (isTRUE(joint)) return("joint")
    if (identical(joint, FALSE)) return("marginal")
    stop("Emulator samples need sampling-method metadata")
  }
  sampling <- match.arg(sampling, c("marginal", "joint"))
  if (!is.null(joint) && !identical(joint, sampling == "joint"))
    stop("Emulator sampling metadata disagree")
  sampling
}

# Both methods share the variance draw and the terminal coefficient draw.
# Marginal sampling uses Eq. (9); joint sampling conditions on each successor.
sample_calibration_emulator <- function(para_ffbs, R_episode, G_ls, nsam,
    seed = NULL, sampling = c("marginal", "joint")) {
  sampling <- match.arg(sampling)
  ff <- para_ffbs$ff
  E <- length(ff) - as.integer(!is.null(ff$prior))
  if (E < 1 || length(nsam) != 1L || !is.finite(nsam) || nsam < 1 ||
      nsam != as.integer(nsam)) stop("Invalid filtering results or sample count")
  restore_seed <- .cal_seed(seed)
  on.exit(restore_seed(), add = TRUE)
  L_R <- t(.cal_spd(R_episode, "R_episode"))
  p <- nrow(ff[[E]]$mt)
  C <- ncol(ff[[E]]$mt)
  if (nrow(L_R) != C) stop("R_episode has the wrong size")
  sigma2 <- 1 / stats::rgamma(nsam, shape = ff[[E]]$nt, rate = ff[[E]]$Dt)
  theta <- lapply(seq_len(E), function(e) array(NA_real_, c(p, C, nsam)))
  matrix_normal <- function(mean, L, variance)
    mean + sqrt(variance) * L %*% matrix(stats::rnorm(p * C), p, C) %*% t(L_R)
  L <- t(.cal_spd(ff[[E]]$Mt))
  for (m in seq_len(nsam)) theta[[E]][, , m] <- matrix_normal(ff[[E]]$mt, L, sigma2[m])
  if (E > 1) for (e in seq.int(E - 1, 1)) {
    if (sampling == "marginal") {
      mean <- matrix(para_ffbs$bs$st[, , e], p, C)
      H <- matrix(para_ffbs$bs$St[, , e], p, p)
    } else {
      G <- if (is.list(G_ls)) G_ls[[e + 1]] else G_ls
      Ainv <- chol2inv(.cal_spd(ff[[e + 1]]$At))
      B <- ff[[e]]$Mt %*% t(G) %*% Ainv
      H <- ff[[e]]$Mt - B %*% ff[[e + 1]]$At %*% t(B)
    }
    L <- t(.cal_spd((H + t(H)) / 2))
    for (m in seq_len(nsam)) {
      if (sampling == "joint")
        mean <- ff[[e]]$mt + B %*% (theta[[e + 1]][, , m] - ff[[e + 1]]$at)
      theta[[e]][, , m] <- matrix_normal(mean, L, sigma2[m])
    }
  }
  list(theta = theta, sigma2 = sigma2,
       joint_trajectories = sampling == "joint", sampling_method = sampling)
}

# Y and F_ls have one entry per episode, theta[[e]] is p x C x M.
# block_cols[b, ] gives the original spatial indices of episode block b.
# R_field is the full-grid covariance, or NULL with field_kernel =
# list(coords, phi, sigma2, tau2) describing the exponential kernel that
# .cal_R_field() builds at the observed sites when it is needed.
# design_kernel and design_groups choose the design GP: "sqexp"/"shared" is
# one squared-exponential kernel for all episodes (spDBL <= 1.2.0 before the
# change of default); calibration_prepare() passes "matern52"/"time", one
# Matern 5/2 kernel per time step.
make_offline_calibration_input <- function(input, Y, F_ls, posterior, block_cols,
    R_field, phi_para, lower, upper, gp_sigma2 = 1.1, gp_tau2 = 1e-4,
    design_uncertainty = TRUE, field_kernel = NULL,
    design_kernel = c("sqexp", "matern52"), design_groups = c("shared", "time")) {
  design_kernel <- match.arg(design_kernel)
  design_groups <- match.arg(design_groups)
  input <- as.matrix(input)
  block_cols <- as.matrix(block_cols)
  E <- length(Y)
  B <- nrow(block_cols)
  C <- ncol(block_cols)
  if (is.null(R_field)) {
    k <- field_kernel
    if (!is.list(k) || !is.matrix(k$coords) || any(!is.finite(k$coords)) ||
        length(c(k$phi, k$sigma2, k$tau2)) != 3L ||
        any(!is.finite(c(k$phi, k$sigma2, k$tau2))) ||
        k$phi <= 0 || k$sigma2 <= 0 || k$tau2 < 0)
      stop("Supply R_field or a valid field_kernel")
    S <- nrow(k$coords)
  } else S <- nrow(R_field)
  p <- ncol(F_ls[[1]])
  M <- length(posterior$sigma2)
  sampling <- .cal_sampling_method(posterior$sampling_method,
                                   posterior$joint_trajectories)
  posterior$sampling_method <- sampling
  posterior$joint_trajectories <- sampling == "joint"
  if (E < 1 || E %% B != 0 || length(F_ls) != E || length(posterior$theta) != E ||
      !identical(sort(as.integer(block_cols)), seq_len(S)) ||
      M < 1 || any(!is.finite(posterior$sigma2)) || any(posterior$sigma2 <= 0) ||
      !is.finite(phi_para) || phi_para <= 0 || gp_sigma2 <= 0 || gp_tau2 <= 0) {
    stop("Invalid emulator input or episode mapping")
  }
  if (!isTRUE(design_uncertainty) && !isFALSE(design_uncertainty))
    stop("design_uncertainty must be TRUE or FALSE")
  if (!is.null(R_field)) .cal_spd(R_field, "R_field")
  for (e in seq_len(E)) {
    if (!identical(dim(posterior$theta[[e]]), c(p, C, M)) ||
        any(!is.finite(posterior$theta[[e]])) ||
        !identical(dim(Y[[e]]), c(nrow(input), C)) ||
        !identical(dim(F_ls[[e]]), c(nrow(input), p)))
      stop("Incompatible episode arrays")
  }
  dist <- as.matrix(stats::dist(input))
  V <- gp_sigma2 * exp(-phi_para * dist) + diag(gp_tau2, nrow(input))
  V_chol <- .cal_spd(V, "input kernel")
  # Episodes are time-major: episode e belongs to time (e - 1) %/% B + 1.
  groups <- if (design_groups == "time") (seq_len(E) - 1L) %/% B + 1L else rep(1L, E)
  predictor <- fit_offline_ar2(input, F_ls, lower, upper, kernel = design_kernel,
                               groups = groups)
  predictor$uncertainty <- .cal_uncertainty_note(design_uncertainty)
  structure(list(version = 1L, input = input, Y = Y, F_ls = F_ls,
    posterior = posterior, block_cols = block_cols, R_field = R_field,
    field_kernel = if (is.null(R_field)) field_kernel,
    phi_para = phi_para, gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2,
    V_chol = V_chol, residual_cache = new.env(parent = emptyenv()),
    predictor = predictor, design_uncertainty = design_uncertainty,
    design_kernel = design_kernel, design_groups = design_groups,
    design_whitened = if (design_uncertainty) .cal_design_whitened(predictor),
    T = E %/% B, S = S, E = E, B = B, C = C, p = p, M = M,
    sampling_method = sampling),
    class = c("spDBL_calibration_input", "offline_calibration_emulator"))
}

# Design-row covariance factors for one coefficient draw. The design rows
# of all episodes at one time share one GP (one kernel group per time, or one
# group for all times), so their errors at eta are jointly Gaussian: the row
# of episode e has error W_e' xi * sqrt(v_t(eta) / N), xi ~ N(0, I_N), with
# W_e the episode's columns of .cal_design_whitened() and v_t(eta) the
# variance factor of time t's group. The outputs F_e Theta_e of block b at
# time t therefore have error sqrt(v_t(eta)) * Z_t' xi, where Z_t (N x S)
# holds W_e Theta_e / sqrt(N) in the block's columns. The design-row
# covariance at time t is v_t(eta) * Z_t' Z_t; it couples sites in different
# blocks. Times are treated as independent.
# Returns T matrices over `locations`, or an empty list when the option is off.
.cal_design_factors <- function(emulator, draw, locations = seq_len(emulator$S)) {
  if (!isTRUE(emulator$design_uncertainty)) return(list())
  W <- emulator$design_whitened
  if (is.null(W)) W <- .cal_design_whitened(emulator$predictor)
  N <- nrow(W)
  factors <- lapply(seq_len(emulator$T), function(t) matrix(0, N, length(locations)))
  for (e in seq_len(emulator$E)) {
    t <- (e - 1L) %/% emulator$B + 1L
    b <- (e - 1L) %% emulator$B + 1L
    position <- match(emulator$block_cols[b, ], locations)
    observed <- !is.na(position)
    if (!any(observed)) next
    theta <- matrix(emulator$posterior$theta[[e]][, , draw], emulator$p, emulator$C)
    factors[[t]][, position[observed]] <- W[, (e - 1L) * emulator$p + seq_len(emulator$p)] %*%
      theta[, observed, drop = FALSE] / sqrt(N)
  }
  factors
}

# Full-grid emulator covariance at `locations`: the stored R_field, or the
# exponential kernel that emulator_learn() used to store, built only there.
.cal_R_field <- function(emulator, locations = seq_len(emulator$S)) {
  if (!is.null(emulator$R_field))
    return(emulator$R_field[locations, locations, drop = FALSE])
  k <- emulator$field_kernel
  x <- k$coords[locations, , drop = FALSE]
  unname(k$sigma2 * exp(-k$phi * as.matrix(stats::dist(x))) +
           diag(k$tau2, length(locations)))
}

.cal_R_field_diag <- function(emulator, locations = seq_len(emulator$S)) {
  if (!is.null(emulator$R_field)) return(diag(emulator$R_field)[locations])
  rep(emulator$field_kernel$sigma2 + emulator$field_kernel$tau2, length(locations))
}

# Field covariance at time t: the emulator term shared by all times plus the
# design-row term, if any.
.cal_time_covariance <- function(pred, t) {
  if (!length(pred$design_factors)) return(pred$Sigma)
  pred$Sigma + pred$design_variance[t] * crossprod(pred$design_factors[[t]])
}

# Latent variance factor of the design GP at each time, from
# predict_offline_ar2(); the episodes of one time share a kernel group.
.cal_time_factor <- function(design, emulator)
  design$episode_factor[(seq_len(emulator$T) - 1L) * emulator$B + 1L]

# Field mean and covariance at eta for one coefficient draw.
# covariance = "full": Sigma is the emulator covariance shared by all times;
#   the design-row term at time t is design_variance[t] * Z_t' Z_t with Z_t
#   from design_factors (absent when design uncertainty is off). Use
#   .cal_time_covariance() for the full matrix at time t.
# covariance = "marginal": variance is the T x S matrix of pointwise
#   variances, computed without forming S x S matrices.
# covariance = "none": the mean only.
predict_calibration <- function(eta, emulator_draw, emulator,
    locations = seq_len(emulator$S), covariance = c("full", "marginal", "none")) {
  covariance <- match.arg(covariance)
  m <- as.integer(emulator_draw)
  if (length(m) != 1 || is.na(m) || m < 1 || m > emulator$M ||
      !length(locations) || anyDuplicated(locations) ||
      any(!locations %in% seq_len(emulator$S))) stop("Invalid draw or locations")
  design <- predict_offline_ar2(emulator$predictor, eta)
  F_new_ls <- lapply(design$F, function(f) matrix(f, nrow = 1))
  theta <- lapply(emulator$posterior$theta, function(th)
    matrix(th[, , m], emulator$p, emulator$C))
  prediction <- FFBS_predict_exact(Y = emulator$Y, F_ls = emulator$F_ls,
    F_new_ls = F_new_ls, input = emulator$input, input_new = eta,
    nT = emulator$E, phi_para = emulator$phi_para,
    gp_sigma2 = emulator$gp_sigma2, gp_tau2 = emulator$gp_tau2,
    theta = theta, V_chol = emulator$V_chol,
    residual_weights = .cal_residual_weights(emulator, m), return_covariance = TRUE)
  factor <- prediction$input_covariance[1, 1]
  if (factor <= 0) stop("Nonpositive emulator conditional variance")
  fields <- recover_from_EP_exact(prediction$mean, nT_ori = emulator$T,
    nT = emulator$E, block_cols = emulator$block_cols)
  mu <- unname(do.call(rbind, fields))
  out <- list(mu = mu[, locations, drop = FALSE])
  if (covariance == "none") return(out)
  scale <- emulator$posterior$sigma2[m] * factor
  factors <- .cal_design_factors(emulator, m, locations)
  time_factor <- .cal_time_factor(design, emulator)
  if (covariance == "marginal") {
    variance <- matrix(scale * .cal_R_field_diag(emulator, locations),
                       emulator$T, length(locations), byrow = TRUE)
    for (t in seq_along(factors)) variance[t, ] <- variance[t, ] +
      time_factor[t] * colSums(factors[[t]]^2)
    out$variance <- variance
    return(out)
  }
  out$Sigma <- scale * .cal_R_field(emulator, locations)
  if (length(factors)) {
    out$design_variance <- time_factor
    out$design_factors <- factors
  }
  out
}

validate_offline_ar2 <- function(emulator, input_new, F_new_ls) {
  input_new <- as.matrix(input_new)
  if (length(F_new_ls) != emulator$E || ncol(input_new) != ncol(emulator$input) ||
      !all(vapply(F_new_ls, function(f) identical(dim(f), c(nrow(input_new), emulator$p)), logical(1)))) {
    stop("Held-out inputs/designs have incompatible dimensions")
  }
  pred <- lapply(seq_len(nrow(input_new)), function(i)
    predict_offline_ar2(emulator$predictor, input_new[i, ])$F)
  do.call(rbind, lapply(seq_len(emulator$E), function(e) {
    truth <- F_new_ls[[e]]
    estimate <- do.call(rbind, lapply(pred, function(f) f[[e]]))
    error <- estimate - truth
    data.frame(episode = e, time = (e - 1) %/% emulator$B,
      rmse = sqrt(mean(error^2)), max_abs_error = max(abs(error)),
      relative_rmse = sqrt(mean(error^2)) / max(1e-10, sqrt(mean(truth^2))))
  }))
}


# Cache only shared factors. Draw-specific residual coefficients are temporary.
.cal_residual_weights <- function(emulator, draw) {
  if (!is.null(emulator$residual_weights))
    return(emulator$residual_weights[[draw]])
  cache <- emulator$residual_cache
  if (!is.null(cache) && identical(cache$draw, draw)) return(cache$weights)
  weights <- lapply(seq_len(emulator$E), function(e) {
    theta <- matrix(emulator$posterior$theta[[e]][, , draw], emulator$p, emulator$C)
    residual <- emulator$Y[[e]] - emulator$F_ls[[e]] %*% theta
    backsolve(emulator$V_chol, forwardsolve(t(emulator$V_chol), residual))
  })
  if (!is.null(cache)) {
    cache$draw <- draw
    cache$weights <- weights
  }
  weights
}

# A supplied seed is local to the call; NULL follows R's normal RNG behavior.
.cal_seed <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed))
    stop("seed must be NULL or a nonnegative integer")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
    invisible(NULL)
  }
}

#' Prepare a fitted emulator for PDE-free calibration
#'
#' Fits a Gaussian process (the design GP) to the recorded autoregressive
#' design matrices, so that design rows can be predicted at untried inputs.
#' Reuses the fitted emulator's coefficient samples or draws additional samples
#' from saved filtering results, without refitting the dynamic emulator.
#'
#' @param emulator A sigma2R object from [emulator_learn()].
#' @param lower,upper Finite parameter bounds in training-column order.
#'   Equal bounds specify fixed coordinates. Named bounds are matched to
#'   the training parameter names.
#' @param n_draws Number of fresh coefficient sample sets. NULL reuses existing
#'   samples when the selected sampling method agrees with the fitted emulator.
#' @param seed Optional local seed used when drawing fresh samples.
#' @param sampling NULL inherits the fitted emulator's method (normally
#'   "marginal"). "marginal" draws independent smoothed episode coefficients
#'   conditional on a shared variance; "joint" draws coherent backward paths.
#'   Changing the method resamples the existing number of draws unless n_draws
#'   is supplied.
#' @param design_uncertainty If TRUE (the default), the predictive covariance
#'   of the AR design GP is added to the field covariance used by
#'   [calibration_fit()] and [calibration_predict()]. FALSE uses the GP mean
#'   only, as spDBL 1.1.1 did; its intervals are too narrow.
#' @param design_kernel Kernel of the design GP: `"matern52"` (Matern 5/2,
#'   the default) or `"sqexp"` (squared exponential). Both are anisotropic,
#'   with one lengthscale per free parameter on the unit scale.
#' @param design_groups `"time"` (the default) fits separate kernel
#'   hyperparameters, by maximum likelihood, for the design columns of each
#'   time step; `"shared"` fits one set for all times.
#'   `design_kernel = "sqexp", design_groups = "shared"` reproduces the design
#'   GP of earlier spDBL versions.
#' @return A spDBL_calibration_input object containing the training data,
#'   coefficient samples, offline design predictor, and shared GP factors.
#'   At most one draw's residual weights are cached at a time.
#' @details The AR design GP supplies the design rows at a new input. At each
#'   time step, all design columns share one kernel (with `design_groups =
#'   "time"`; one kernel for all times with `"shared"`), so their prediction
#'   errors at an input \eqn{\eta} are jointly Gaussian with covariance
#'   \eqn{C_{F,t}(\eta) = v_t(\eta)\hat\Sigma_{F,t}}: \eqn{v_t(\eta)} is the
#'   GP's latent variance factor and \eqn{\hat\Sigma_{F,t}} is the
#'   maximum-likelihood covariance between the design columns of time
#'   \eqn{t}. With `design_uncertainty = TRUE` this is propagated through each
#'   coefficient draw: at each time, the field covariance gains
#'   \eqn{\Theta^\top C_{F,t}(\eta) \Theta} over the episodes of that time,
#'   including covariance between spatial blocks. Design errors are treated as
#'   independent across times and independent of the emulator's own
#'   prediction error. The emulator's conditional covariance and coefficient
#'   uncertainty are always retained. Calibration support is limited to the
#'   scalar sigma2R model on the original time/grid.
#'   Objects prepared by spDBL 1.1.1 or earlier have no design variance.
#' @seealso [calibration_validate()], [calibration_fit()], [emulator_predict()]
#' @examples
#' data("dt_calibration")
#' d <- dt_calibration
#' em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
#'                      covariance = "sigma2R", episode_window = 2, nsam = 3,
#'                      seed = 1)
#' cal <- calibration_prepare(em, d$lower, d$upper)
#' class(cal)
#' c(times = cal$T, locations = cal$S, draws = cal$M)
#' cal$parameter_names
#' @export
calibration_prepare <- function(emulator, lower, upper, n_draws = NULL,
                                seed = 1234, sampling = NULL,
                                design_uncertainty = TRUE,
                                design_kernel = c("matern52", "sqexp"),
                                design_groups = c("time", "shared")) {
  design_kernel <- match.arg(design_kernel)
  design_groups <- match.arg(design_groups)
  s <- emulator$setup
  if (is.null(s$covariance) || s$covariance != "sigma2R")
    stop("Calibration requires a sigma2R emulator; use emulator_learn(covariance = 'sigma2R')")
  parameter_names <- colnames(s$pde_para_train)
  align_bounds <- function(x) {
    if (!is.null(names(x))) {
      if (is.null(parameter_names) || anyDuplicated(names(x)) ||
          !setequal(names(x), parameter_names))
        stop("Named bounds must match the training parameter names")
      x <- x[parameter_names]
    }
    as.numeric(x)
  }
  lower <- align_bounds(lower)
  upper <- align_bounds(upper)
  existing_method <- .cal_sampling_method(
    attr(emulator$res_ffbs, "sampling_method"),
    attr(emulator$res_ffbs, "joint_trajectories"))
  if (is.null(sampling)) sampling <- existing_method
  sampling <- match.arg(sampling, c("marginal", "joint"))
  reuse_draws <- is.null(n_draws) && sampling == existing_method
  if (reuse_draws) {
    posterior <- list(theta = emulator$res_ffbs[seq_len(s$nT)],
      sigma2 = emulator$res_ffbs$Sigma, joint_trajectories = sampling == "joint",
      sampling_method = sampling)
  } else {
    if (is.null(n_draws)) n_draws <- length(emulator$res_ffbs$Sigma)
    if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1 ||
        n_draws != as.integer(n_draws)) stop("n_draws must be a positive integer")
    posterior <- sample_calibration_emulator(emulator$para_ffbs,
      s$R_episode, s$G_ls, n_draws, seed, sampling)
  }
  out <- make_offline_calibration_input(s$pde_para_train, s$Y, s$F_ls,
    posterior, s$block_cols, s$R_field, s$phi_para, lower, upper,
    s$gp_sigma2, s$gp_tau2, design_uncertainty,
    field_kernel = list(coords = as.matrix(s$ind_sp), phi = s$phi_sp,
                        sigma2 = s$gp_sigma2, tau2 = s$gp_tau2),
    design_kernel = design_kernel, design_groups = design_groups)
  out$parameter_names <- if (is.null(parameter_names))
    paste0("eta", seq_len(ncol(out$input))) else parameter_names
  out$AR_choice <- s$AR_choice
  # The seed behind the emulator draws: the emulator's own when they are
  # reused, this call's when they are drawn afresh.
  out$seed <- if (reuse_draws) s$seed else seed
  out$spDBL_version <- .spdbl_version()
  out$time <- seq_len(out$T) - 1L
  out$grid <- s$ind_sp
  out
}

#' Validate offline calibration predictions against held-out simulations
#'
#' Held-out trajectories are used only for this accuracy assessment. They are
#' not added to either GP's training data.
#' @param object A prepared object from [calibration_prepare()].
#' @param input_new Matrix of held-out mechanistic inputs.
#' @param responses Time-indexed list of held-out response matrices.
#' @param seed Not used since spDBL 1.2.0, whose AR designs are deterministic;
#'   kept for compatibility.
#' @return A list with episode design errors and per-input field RMSE using
#'   known versus GP-predicted autoregressive rows.
#' @examples
#' data("dt_calibration")
#' d <- dt_calibration
#' em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
#'                      covariance = "sigma2R", episode_window = 2, nsam = 3,
#'                      seed = 1)
#' cal <- calibration_prepare(em, d$lower, d$upper)
#' held_out <- lapply(d$validation$responses, function(y) y[1:3, , drop = FALSE])
#' val <- calibration_validate(cal, d$validation$input[1:3, ], held_out)
#' val$field
#' @export
calibration_validate <- function(object, input_new, responses, seed = 1234) {
  restore_seed <- .cal_seed(seed)
  on.exit(restore_seed(), add = TRUE)
  if (is.null(dim(input_new))) input_new <- matrix(input_new, nrow = 1)
  input_new <- as.matrix(input_new)
  n <- nrow(input_new)
  if (!inherits(object, "spDBL_calibration_input") ||
      length(responses) != object$T || !n ||
      !all(vapply(responses, function(y)
        identical(dim(y), c(n, object$S)) && all(is.finite(y)), logical(1))))
    stop("Held-out responses do not match the prepared emulator")
  index <- data.frame(block_no = seq_len(object$B), U = 1L, D = n,
    L = object$block_cols[, 1], R = object$block_cols[, object$C])
  # The current strip blocks are contiguous in original-grid order.
  if (any(apply(object$block_cols, 1, function(x) any(diff(x) != 1))))
    stop("Validation of AR designs requires contiguous strip blocks")
  F_new <- if (object$AR_choice == 2)
    gen_F_ls_AR2_EP(responses, object$T, object$B, index) else
    gen_F_ls_AR1_EP(responses, object$T, object$B, index)
  design <- validate_offline_ar2(object, input_new, F_new)
  F_hat <- lapply(seq_len(n), function(i)
    predict_offline_ar2(object$predictor, input_new[i, ])$F)
  offline <- known <- rep(list(matrix(0, object$T, object$S)), n)
  # Draw-outer order reuses the cached residual weights of each draw.
  for (m in seq_len(object$M)) for (i in seq_len(n)) {
    pred <- predict_calibration(input_new[i, ], m, object, covariance = "none")$mu
    correction <- matrix(0, object$T, object$S)
    for (e in seq_len(object$E)) {
      t <- (e - 1L) %/% object$B + 1L
      b <- (e - 1L) %% object$B + 1L
      th <- matrix(object$posterior$theta[[e]][, , m], object$p, object$C)
      correction[t, object$block_cols[b, ]] <- (F_new[[e]][i, ] - F_hat[[i]][[e]]) %*% th
    }
    offline[[i]] <- offline[[i]] + pred / object$M
    known[[i]] <- known[[i]] + (pred + correction) / object$M
  }
  field <- do.call(rbind, lapply(seq_len(n), function(i) {
    truth <- do.call(rbind, lapply(responses, function(y) y[i, , drop = FALSE]))
    data.frame(input = i, rmse_known_design = sqrt(mean((known[[i]] - truth)^2)),
      rmse_offline_design = sqrt(mean((offline[[i]] - truth)^2)))
  }))
  list(design = design, field = field,
       uncertainty = object$predictor$uncertainty)
}
