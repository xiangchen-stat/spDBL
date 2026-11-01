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
    theta[e, , ] <- emulator$posterior$theta[[e]][, , draw]
    residual[e, , ] <- weights[[e]]
  }
  f <- emulator$predictor
  # Explicit dimensions keep length-one vectors from becoming Stan scalars.
  vector_data <- function(x) array(as.numeric(x), length(x))
  index_data <- function(x) array(as.integer(x), length(x))
  list(N = nrow(emulator$input), D_full = ncol(emulator$input), D = length(f$active),
    T = emulator$T, S_full = emulator$S, S = S, B = emulator$B, C = emulator$C,
    P = emulator$p, E = emulator$E, active = index_data(f$active),
    eta_lower = vector_data(f$lower), eta_upper = vector_data(f$upper),
    input = emulator$input, predictor_x = f$x, lengthscale = vector_data(f$lengthscale),
    predictor_center = vector_data(f$center), predictor_alpha = f$alpha,
    theta = theta, residual_weights = residual, block_cols = emulator$block_cols,
    locations = index_data(locations), phi_para = emulator$phi_para,
    gp_sigma2 = emulator$gp_sigma2, gp_tau2 = emulator$gp_tau2,
    sigma2 = emulator$posterior$sigma2[draw], L_V = t(emulator$V_chol),
    R_field = emulator$R_field[locations, locations, drop = FALSE],
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
    covariance <- pred$Sigma + diag(tau2[t + 1], length(locations))
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
    observation_covariance <- P + pred$Sigma + diag(tau2_all[t + 1L], S)
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
      covariance <- pred$Sigma + diag(tau2[i, t], S)
      K <- pred$Sigma %*% chol2inv(chol(covariance))
      mean <- pred$mu[t, ] + K %*% (z[t, ] - u[i, t, ] - pred$mu[t, ])
      variance <- pred$Sigma - K %*% pred$Sigma
      y[i, t, ] <- as.numeric(mean + t(.cal_spd((variance + t(variance)) / 2)) %*% stats::rnorm(S))
      z_rep[i, t, ] <- y[i, t, ] + u[i, t, ] + stats::rnorm(S, sd = sqrt(tau2[i, t]))
    }
  }
  list(y = y, z_rep = z_rep)
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
    keep_discrepancy = TRUE) {
  if (!requireNamespace("rstan", quietly = TRUE)) stop("Install rstan first")
  if (warmup >= iter || chains < 1 || draws_per_emulator < 1 ||
      draws_per_emulator > chains * (iter - warmup) || !length(draws) ||
      anyDuplicated(draws) || any(!draws %in% seq_len(emulator$M))) stop("Invalid sampling settings")
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  model <- if (is.null(compiled_model))
    .cal_stan_model(marginal_discrepancy, cache = FALSE) else compiled_model
  out <- vector("list", length(draws))
  diagnostics <- vector("list", length(draws))
  for (j in seq_along(draws)) {
    m <- draws[j]
    message("Calibration conditional on emulator draw ", m, " (", j, "/", length(draws), ")")
    data <- make_stan_calibration_data(emulator, m, z, U, locations, n0, d0,
                                       discount, m0, M0)
    init <- function() {
      values <- list(eta_unit = array(rep(0.5, data$D), data$D),
                     precision = rep(n0 / d0, data$T + 1))
      if (!marginal_discrepancy) values$u_std <- matrix(0, data$T + 1, data$S)
      values
    }
    fit <- rstan::sampling(model, data = data, init = init, chains = chains,
      iter = iter, warmup = warmup, seed = seed + j, cores = cores, refresh = refresh,
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
        n_eff = NA_real_, Rhat = NA_real_)
    } else {
      summary <- rstan::summary(fit, pars = c(paste0("eta[", data$active, "]"),
        "tau2", if (!marginal_discrepancy) "u"))$summary
    }
    sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    diagnostics[[j]] <- list(emulator_draw = m, summary = summary,
      divergences = sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1))),
      max_treedepth_hits = sum(vapply(sampler, function(x)
        sum(x[, "treedepth__"] >= if (is.null(control$max_treedepth)) 10 else control$max_treedepth), numeric(1))),
      max_rhat = if (any(is.finite(summary[, "Rhat"]))) max(summary[, "Rhat"], na.rm = TRUE) else NA_real_,
      min_n_eff = if (any(is.finite(summary[, "n_eff"])))
        min(summary[is.finite(summary[, "n_eff"]), "n_eff"]) else NA_real_)
    if (diagnostics[[j]]$divergences > 0 ||
        (!is.na(diagnostics[[j]]$max_rhat) && diagnostics[[j]]$max_rhat > 1.01))
      warning("Emulator draw ", m, " needs additional convergence checks; inspect diagnostics")
    out[[j]] <- list(eta = eta, tau2 = tau[, -1, drop = FALSE], tau2_0 = tau[, 1],
      chain = rep(seq_len(dim(raw)[2]), each = dim(raw)[1])[selected],
      iteration = rep(seq_len(dim(raw)[1]), dim(raw)[2])[selected], emulator_draw = rep(m, length(selected)))
    if (!marginal_discrepancy) {
      out[[j]]$u <- array(states[, -1, , drop = FALSE], c(length(selected), data$T, data$S))
      out[[j]]$u0 <- matrix(states[, 1, ], length(selected), data$S)
    }
    if (!is.null(checkpoint_dir)) saveRDS(list(samples = out[[j]],
      diagnostics = diagnostics[[j]], marginal_discrepancy = marginal_discrepancy,
      chains = chains, iter = iter, warmup = warmup),
      file.path(checkpoint_dir, sprintf("conditional_%03d.rds", m)))
  }
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
