.cal_spd <- function(x, name = "covariance") {
  x <- as.matrix(x)
  if (nrow(x) != ncol(x) || any(!is.finite(x)) ||
      max(abs(x - t(x))) > 1e-7 * max(1, max(abs(x)))) {
    stop(name, " must be finite, square, and symmetric")
  }
  chol((x + t(x)) / 2)
}

.cal_sqexp <- function(x, y, lengthscale) {
  x <- sweep(as.matrix(x), 2, lengthscale, "/")
  y <- sweep(as.matrix(y), 2, lengthscale, "/")
  d2 <- outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)
  exp(-0.5 * pmax(d2, 0))
}

# One shared ARD kernel, separate centered/scaled responses for every design
# column. All hyperparameters are fitted offline; calibration uses its mean.
fit_offline_ar2 <- function(input, F_ls, lower, upper) {
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
  x <- sweep(sweep(input[, active, drop = FALSE], 2, lower[active], "-"),
             2, upper[active] - lower[active], "/")
  targets <- do.call(cbind, F_ls)
  center <- colMeans(targets)
  scales <- apply(targets, 2, stats::sd)
  varying <- is.finite(scales) & scales > 1e-10
  scales[!varying] <- 1
  standardized <- sweep(sweep(targets, 2, center, "-"), 2, scales, "/")
  # Constant responses add no information about kernel hyperparameters.
  response_gram <- tcrossprod(standardized[, varying, drop = FALSE])
  n_response <- sum(varying)
  objective <- function(par) {
    K <- .cal_sqexp(x, x, exp(par[seq_along(active)]))
    diag(K) <- diag(K) + exp(par[length(par)])
    L <- chol(K)
    (n_response * 2 * sum(log(diag(L))) + sum(chol2inv(L) * response_gram)) /
      max(1, n_response)
  }
  starts <- c(0.2, 0.5, 1)
  fits <- lapply(starts, function(s) stats::optim(
    log(c(rep(s, length(active)), 1e-4)), objective, method = "L-BFGS-B",
    lower = log(c(rep(0.03, length(active)), 1e-8)),
    upper = log(c(rep(5, length(active)), 0.2))))
  best <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]
  lengthscale <- exp(best$par[seq_along(active)])
  nugget <- exp(utils::tail(best$par, 1))
  K <- .cal_sqexp(x, x, lengthscale)
  diag(K) <- diag(K) + nugget
  Kinv <- chol2inv(chol(K))
  alpha <- Kinv %*% sweep(standardized, 2, scales, "*")
  loo_error <- sweep(alpha, 1, diag(Kinv), "/")
  structure(list(x = x, active = active, lower = lower, upper = upper,
    center = center, alpha = alpha, Kinv = Kinv, lengthscale = lengthscale,
    nugget = nugget, scales = scales, n_episode = length(F_ls), p = widths[1],
    loo_rmse = sqrt(mean(loo_error^2)), optimizer_convergence = best$convergence,
    uncertainty = "Conditional on the offline GP mean; predictor uncertainty is not integrated"),
    class = "offline_ar2_predictor")
}

predict_offline_ar2 <- function(object, eta) {
  if (length(eta) != length(object$lower) || any(!is.finite(eta)) ||
      any(eta < object$lower) || any(eta > object$upper)) stop("eta exceeds bounds")
  a <- object$active
  x <- matrix((eta[a] - object$lower[a]) / (object$upper[a] - object$lower[a]), 1)
  k <- .cal_sqexp(object$x, x, object$lengthscale)
  mean <- as.numeric(object$center + crossprod(k, object$alpha))
  variance_factor <- max(0, as.numeric(1 - crossprod(k, object$Kinv %*% k)))
  list(F = lapply(seq_len(object$n_episode), function(e)
    mean[((e - 1) * object$p + 1):(e * object$p)]),
    variance = variance_factor * object$scales^2)
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
make_offline_calibration_input <- function(input, Y, F_ls, posterior, block_cols,
    R_field, phi_para, lower, upper, gp_sigma2 = 1.1, gp_tau2 = 1e-4) {
  input <- as.matrix(input)
  block_cols <- as.matrix(block_cols)
  E <- length(Y)
  B <- nrow(block_cols)
  C <- ncol(block_cols)
  S <- nrow(R_field)
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
  .cal_spd(R_field, "R_field")
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
  structure(list(version = 1L, input = input, Y = Y, F_ls = F_ls,
    posterior = posterior, block_cols = block_cols, R_field = R_field,
    phi_para = phi_para, gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2,
    V_chol = V_chol, residual_cache = new.env(parent = emptyenv()),
    predictor = fit_offline_ar2(input, F_ls, lower, upper),
    T = E %/% B, S = S, E = E, B = B, C = C, p = p, M = M,
    sampling_method = sampling),
    class = c("spDBL_calibration_input", "offline_calibration_emulator"))
}

predict_calibration <- function(eta, emulator_draw, emulator, locations = seq_len(emulator$S)) {
  m <- as.integer(emulator_draw)
  if (length(m) != 1 || is.na(m) || m < 1 || m > emulator$M ||
      !length(locations) || anyDuplicated(locations) ||
      any(!locations %in% seq_len(emulator$S))) stop("Invalid draw or locations")
  F_new_ls <- lapply(predict_offline_ar2(emulator$predictor, eta)$F,
                     function(f) matrix(f, nrow = 1))
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
  list(mu = mu[, locations, drop = FALSE], Sigma =
    emulator$posterior$sigma2[m] * factor * emulator$R_field[locations, locations, drop = FALSE])
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
#' Fits a shared GP mean to the recorded autoregressive design matrices.
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
#' @return A spDBL_calibration_input object containing the training data,
#'   coefficient samples, offline design predictor, and shared GP factors.
#'   At most one draw's residual weights are cached at a time.
#' @details The AR design GP is used through its fitted mean; interpolation
#'   uncertainty is not integrated into calibration. The emulator's conditional
#'   covariance and coefficient uncertainty are retained. Initial calibration
#'   support is limited to the scalar sigma2R model on the original time/grid.
#' @seealso [calibration_validate()], [calibration_fit()], [emulator_predict()]
#' @export
calibration_prepare <- function(emulator, lower, upper, n_draws = NULL,
                                seed = 1234, sampling = NULL) {
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
  if (is.null(n_draws) && sampling == existing_method) {
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
    s$gp_sigma2, s$gp_tau2)
  out$parameter_names <- if (is.null(parameter_names))
    paste0("eta", seq_len(ncol(out$input))) else parameter_names
  out$AR_choice <- s$AR_choice
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
#' @param seed Local seed for the special AR2 validation initialization.
#' @return A list with episode design errors and per-input field RMSE using
#'   known versus GP-predicted autoregressive rows.
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
  field <- do.call(rbind, lapply(seq_len(n), function(i) {
    F_hat <- predict_offline_ar2(object$predictor, input_new[i, ])$F
    offline <- known <- matrix(0, object$T, object$S)
    for (m in seq_len(object$M)) {
      pred <- predict_calibration(input_new[i, ], m, object)$mu
      correction <- matrix(0, object$T, object$S)
      for (e in seq_len(object$E)) {
        t <- (e - 1L) %/% object$B + 1L
        b <- (e - 1L) %% object$B + 1L
        th <- matrix(object$posterior$theta[[e]][, , m], object$p, object$C)
        correction[t, object$block_cols[b, ]] <- (F_new[[e]][i, ] - F_hat[[e]]) %*% th
      }
      offline <- offline + pred / object$M
      known <- known + (pred + correction) / object$M
    }
    truth <- do.call(rbind, lapply(responses, function(y) y[i, , drop = FALSE]))
    data.frame(input = i, rmse_known_design = sqrt(mean((known - truth)^2)),
      rmse_offline_design = sqrt(mean((offline - truth)^2)))
  }))
  list(design = design, field = field,
       uncertainty = object$predictor$uncertainty)
}
