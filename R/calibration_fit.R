# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

.spdbl_stan_cache <- new.env(parent = emptyenv())

.cal_stan_model <- function(marginal = TRUE, cache = TRUE, cache_dir = NULL) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Install the suggested package 'rstan' to fit calibration models")
  filename <- if (marginal) "calibration_marginal.stan" else "calibration_offline.stan"
  path <- system.file("stan", filename, package = "spDBL", mustWork = TRUE)
  key <- paste(unname(tools::md5sum(path)), getRversion(),
    utils::packageVersion("rstan"), utils::packageVersion("StanHeaders"),
    R.version$platform, "O1", sep = "-")
  if (cache && exists(key, .spdbl_stan_cache, inherits = FALSE))
    return(get(key, .spdbl_stan_cache, inherits = FALSE))
  if (is.null(cache_dir)) cache_dir <- tools::R_user_dir("spDBL", "cache")
  cache_file <- file.path(cache_dir, paste0(key, ".rds"))
  if (cache && file.exists(cache_file)) {
    model <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(model) && methods::is(model, "stanmodel")) {
      assign(key, model, .spdbl_stan_cache)
      return(model)
    }
  }
  # stanc's O1 optimisations: about 10% faster gradients, same results.
  model <- rstan::stan_model(file = path, auto_write = FALSE, save_dso = TRUE,
                             allow_optimizations = TRUE)
  if (cache) {
    assign(key, model, .spdbl_stan_cache)
    if (dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE) || dir.exists(cache_dir))
      tryCatch(saveRDS(model, cache_file), error = function(e)
        warning("Could not persist the Stan cache: ", conditionMessage(e)))
  }
  model
}

# Validate init and convert parameter values to the unit scale of the free
# coordinates. Returns "mode", "midpoint", "random", or a numeric vector.
.cal_init_spec <- function(init, object) {
  f <- object$predictor
  if (is.character(init) && length(init) == 1L && init %in% c("mode", "midpoint", "random"))
    return(init)
  if (!is.numeric(init) || length(init) != length(f$lower) || any(!is.finite(init)))
    stop("init must be \"mode\", \"midpoint\", \"random\", or one value per parameter")
  if (!is.null(names(init))) {
    if (anyDuplicated(names(init)) || !setequal(names(init), object$parameter_names))
      stop("Named init values must match the parameter names")
    init <- init[object$parameter_names]
  }
  a <- f$active
  unit <- (as.numeric(init[a]) - f$lower[a]) / (f$upper[a] - f$lower[a])
  if (any(unit <= 0 | unit >= 1))
    stop("init must lie strictly inside the bounds of the free parameters")
  unit
}

#' Fit modular calibration conditional on emulator posterior draws
#'
#' Targets one conditional posterior per fixed emulator coefficient set and pools equal
#' numbers of samples. The emulator is not updated by field observations.
#' Gaussian emulator fields are integrated out. By default discrepancy states
#' are also integrated out using a Kalman filter.
#'
#' @param object Prepared emulator from [calibration_prepare()].
#' @param z Finite time-by-observed-location matrix of field observations.
#' @param U Full-grid symmetric positive-definite discrepancy covariance scale.
#'   Correlation matrices are supported, but the diagonal need not equal one.
#'   Discrepancy increments have covariance tau2[t] * U; observation errors
#'   retain covariance tau2[t] * I.
#' @param locations Ordered original-grid indices corresponding to columns of z.
#' @param backend Currently "rstan".
#' @param draws Emulator draw indices to use, without duplicates.
#' @param chains,iter,warmup Number of chains, total iterations per chain, and
#'   warmup iterations for each emulator-conditional fit. `iter` has no
#'   default and must be given. For inference use several chains and enough
#'   iterations, for example `chains = 4, iter = 2000`; `iter = 1` (with no
#'   warmup) is a quick execution check, not inference. `chains` defaults to
#'   one and `warmup` to `floor(iter / 2)`.
#' @param draws_per_emulator Equal number of retained samples per conditional.
#'   Defaults to the smaller of 200 and the available post-warmup samples.
#' @param n0,d0 Positive inverse-gamma shape and scale hyperparameters.
#' @param discount Calibration variance discount factor in (0, 1].
#' @param m0,M0 Initial discrepancy mean and covariance scale at observed sites.
#' @param seed Nonnegative integer RNG seed, local to the function call.
#' @param cores Positive number of parallel chains. Defaults to one.
#' @param parallel_draws Number of conditional fits (one per emulator draw)
#'   run at the same time, in forked processes. Each forked fit runs its
#'   chains one after another (rstan's parallel chains can deadlock inside a
#'   forked process), so `cores` is then not used and `parallel_draws`
#'   processes run at once. Results do not depend on it: every fit draws its
#'   random numbers from its own seeds, and rstan's draws do not depend on the
#'   number of cores. Forking is unavailable on Windows, where fits run one
#'   after another. Defaults to one.
#' @param marginal_discrepancy Whether to integrate out discrepancy states.
#' @param keep_predictions,keep_discrepancy Whether to recover and retain latent
#'   fields, replicated observations, and discrepancy states.
#' @param control Named list passed to rstan sampling. With zero warmup the
#'   defaults disable adaptation and use stepsize = 0.001. With warmup they
#'   use adapt_delta = 0.95. Both use max_treedepth = 12.
#' @param refresh Stan progress reporting interval.
#' @param checkpoint_dir Optional directory for completed conditional fits.
#'   Each fit is saved as it completes; a fit already saved there with the same
#'   data, settings, seed and Stan program is reused instead of refitted, so an
#'   interrupted calibration resumes where it stopped.
#' @param cache,cache_dir Whether to reuse compiled models and optional user
#'   cache directory. Installed package files are never modified.
#' @param init Starting values of the free parameters for every chain.
#'   `"mode"` (the default) runs L-BFGS (`rstan::optimizing()`) from
#'   `mode_starts` random starts in each conditional fit and starts every chain,
#'   slightly jittered, at the optimum with the highest log posterior density,
#'   with the precisions of that optimum. Chains started far from the main mode
#'   can settle in local modes of negligible posterior mass, which this avoids;
#'   the optima found are recorded in `diagnostics[[j]]$mode_search`.
#'   `"random"` draws each chain's start uniformly within the bounds,
#'   reproducibly from `seed`, so that \eqn{\hat R} can reveal chains that
#'   settle in different modes; `"midpoint"` starts every chain at the centre of
#'   the bounds; or give a numeric vector of parameter values in training-column
#'   order (named vectors are matched by name), strictly inside the bounds of
#'   the free parameters. Values given for fixed parameters are ignored. Each
#'   conditional fit starts afresh; except for `"mode"`, precisions start at
#'   n0 / d0.
#' @param mode_starts Number of random starts of the mode search for
#'   `init = "mode"`.
#' @param compiled_model Optional compatible compiled Stan model for advanced use.
#' @return An ordinary serializable list of class spDBL_calibration containing
#'   samples, per-conditional diagnostics, settings, and prediction context.
#'   No compiled Stan model is stored in the result.
#' @details Requires the suggested package rstan and a C++ toolchain on first
#'   compilation. One iteration per emulator draw is an execution experiment,
#'   not converged posterior inference; its samples stay near `init`. When the
#'   midpoint of the bounds is close to the true value, as in `dt_calibration`,
#'   use `init = "random"` to check that results do not depend on the start. Missing observations and new time/grid
#'   locations are not supported in this initial interface.
#' @seealso [calibration_predict()], [calibration_validate()]
#' @examplesIf interactive() && requireNamespace("rstan", quietly = TRUE)
#' data("dt_calibration")
#' d <- dt_calibration
#' em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
#'                      covariance = "sigma2R", episode_window = 2, nsam = 3,
#'                      seed = 1)
#' cal <- calibration_prepare(em, d$lower, d$upper)
#' fit <- calibration_fit(cal, d$field$z, d$field$U,
#'                        locations = d$field$locations, draws = 1:2,
#'                        chains = 1, iter = 200, seed = 1)
#' summary(fit)
#' @export
calibration_fit <- function(object, z, U, locations = seq_len(object$S),
    backend = "rstan", draws = seq_len(object$M), chains = 1L,
    iter, warmup = iter %/% 2L,
    draws_per_emulator = min(200L, chains * (iter - warmup)),
    n0 = 30, d0 = 2, discount = 1,
    m0 = rep(0, length(locations)), M0 = diag(length(locations)),
    seed = 1234L, cores = 1L, marginal_discrepancy = TRUE,
    keep_predictions = FALSE, keep_discrepancy = TRUE,
    control = if (warmup == 0L)
      list(adapt_engaged = FALSE, stepsize = 0.001, max_treedepth = 12L) else
      list(adapt_delta = 0.95, max_treedepth = 12L), refresh = 100L,
    checkpoint_dir = NULL, cache = TRUE, cache_dir = NULL, compiled_model = NULL,
    init = "mode", mode_starts = 20L, parallel_draws = 1L) {
  if (missing(iter))
    stop("Set iter, the number of iterations per chain: for example ",
         "chains = 4, iter = 2000 for inference, or iter = 1 for a quick execution check")
  if (is.null(seed)) stop("Calibration fitting requires an integer seed")
  restore_seed <- .cal_seed(seed)
  on.exit(restore_seed(), add = TRUE)
  if (!inherits(object, "spDBL_calibration_input")) stop("Use calibration_prepare() first")
  if (!identical(backend, "rstan")) stop("Only the rstan backend is currently supported")
  counts <- c(chains, iter, warmup, draws_per_emulator, cores, parallel_draws)
  if (length(parallel_draws) != 1L || any(!is.finite(counts)) || any(counts != as.integer(counts)) ||
      min(chains, iter, draws_per_emulator, cores, parallel_draws) < 1 || warmup < 0 || warmup >= iter ||
      draws_per_emulator > chains * (iter - warmup) || !length(draws) ||
      anyDuplicated(draws) || any(!draws %in% seq_len(object$M)))
    stop("Invalid calibration sampling settings")
  z <- as.matrix(z)
  if (any(!is.finite(c(n0, d0, discount, m0))) ||
      !identical(dim(U), c(object$S, object$S)) ||
      any(!is.finite(U)))
    stop("U must be a finite full-grid covariance scale matrix; check priors")
  .cal_spd(U, "U")
  init_unit <- .cal_init_spec(init, object)
  if (length(mode_starts) != 1L || !is.finite(mode_starts) || mode_starts < 1 ||
      mode_starts != as.integer(mode_starts)) stop("mode_starts must be a positive integer")
  invisible(make_stan_calibration_data(object, draws[1], z, U, locations,
    n0, d0, discount, m0, M0))
  if (is.null(compiled_model))
    compiled_model <- .cal_stan_model(marginal_discrepancy, cache, cache_dir)
  result <- fit_calibration_rstan(object, z, U, draws = draws,
    locations = locations, chains = chains, iter = iter, warmup = warmup,
    draws_per_emulator = draws_per_emulator, seed = seed, cores = cores,
    keep_predictions = keep_predictions, n0 = n0, d0 = d0, discount = discount,
    m0 = m0, M0 = M0, control = control, compiled_model = compiled_model,
    marginal_discrepancy = marginal_discrepancy, refresh = refresh,
    checkpoint_dir = checkpoint_dir, keep_discrepancy = keep_discrepancy,
    init = init_unit, mode_starts = as.integer(mode_starts),
    parallel_draws = as.integer(parallel_draws))
  colnames(result$mcmc_eta) <- object$parameter_names
  result$execution_only <- warmup == 0L || iter - warmup < 2L
  result$settings <- list(backend = backend, draws = draws, chains = chains,
    iter = iter, warmup = warmup, draws_per_emulator = draws_per_emulator,
    seed = seed, control = control, init = init, mode_starts = mode_starts,
    parallel_draws = parallel_draws,
    chain_cores = if (parallel_draws > 1L && .Platform$OS.type != "windows" &&
                      length(draws) > 1L) 1L else as.integer(cores),
    emulator_sampling = .cal_sampling_method(object$posterior$sampling_method,
                                             object$posterior$joint_trajectories),
    spDBL_version = .spdbl_version())
  result$emulator <- object
  result$field <- list(z = z, U = U, locations = locations, n0 = n0, d0 = d0,
    discount = discount, m0 = m0, M0 = M0)
  result$call <- match.call()
  class(result) <- "spDBL_calibration"
  result
}

#' Recover fields at the calibrated observation times and locations
#'
#' Recovers discrepancy trajectories when needed, then draws latent emulator
#' fields conditional on the observations and fresh replicated observations.
#' This is not forecasting at new times or new spatial sites.
#' @param object Result from [calibration_fit()].
#' @param seed Local RNG seed for Gaussian recovery.
#' @return A list containing y, u, and z_rep arrays, each with dimensions
#'   posterior draws by observation times by observed locations.
#' @examplesIf interactive() && requireNamespace("rstan", quietly = TRUE)
#' data("dt_calibration")
#' d <- dt_calibration
#' em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
#'                      covariance = "sigma2R", episode_window = 2, nsam = 3,
#'                      seed = 1)
#' cal <- calibration_prepare(em, d$lower, d$upper)
#' fit <- calibration_fit(cal, d$field$z, d$field$U,
#'                        locations = d$field$locations, draws = 1:2,
#'                        chains = 1, iter = 200, seed = 1)
#' pred <- calibration_predict(fit)
#' dim(pred$y)
#' @export
calibration_predict <- function(object, seed = 1234L) {
  if (!inherits(object, "spDBL_calibration")) stop("Supply a calibration_fit() result")
  restore_seed <- .cal_seed(seed)
  on.exit(restore_seed(), add = TRUE)
  f <- object$field
  u <- object$mcmc_u
  if (is.null(u)) u <- recover_calibration_discrepancy(object$mcmc_eta,
    object$mcmc_tau2, object$tau2_0, object$emulator_draw, object$emulator,
    f$z, f$U, f$locations, f$m0, f$M0, seed)$u
  pred <- recover_calibration_predictions(object$mcmc_eta, object$mcmc_tau2,
    u, object$emulator_draw, object$emulator, f$z, f$locations, seed + 1L)
  list(y = pred$y, u = u, z_rep = pred$z_rep)
}

#' Summarize modular calibration samples
#' @param object A calibration result.
#' @param ... Unused.
#' @return A list of parameter summaries, conditional diagnostics, and settings.
#' @examplesIf interactive() && requireNamespace("rstan", quietly = TRUE)
#' data("dt_calibration")
#' d <- dt_calibration
#' em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
#'                      covariance = "sigma2R", episode_window = 2, nsam = 3,
#'                      seed = 1)
#' cal <- calibration_prepare(em, d$lower, d$upper)
#' fit <- calibration_fit(cal, d$field$z, d$field$U,
#'                        locations = d$field$locations, draws = 1:2,
#'                        chains = 1, iter = 200, seed = 1)
#' s <- summary(fit)
#' s$parameters
#' @export
summary.spDBL_calibration <- function(object, ...) {
  q <- t(apply(object$mcmc_eta, 2, stats::quantile, c(0.025, 0.5, 0.975)))
  parameters <- data.frame(parameter = colnames(object$mcmc_eta),
    mean = colMeans(object$mcmc_eta), lower = q[, 1], median = q[, 2], upper = q[, 3],
    row.names = NULL)
  diagnostics <- do.call(rbind, lapply(object$diagnostics, function(d)
    data.frame(emulator_draw = d$emulator_draw, divergences = d$divergences,
      max_treedepth_hits = d$max_treedepth_hits, max_rhat = d$max_rhat,
      min_n_eff = d$min_n_eff,
      min_ess_tail = if (is.null(d$min_ess_tail)) NA_real_ else d$min_ess_tail)))
  structure(list(parameters = parameters, diagnostics = diagnostics,
    n_samples = nrow(object$mcmc_eta), settings = object$settings,
    execution_only = object$execution_only,
    emulator_sampling = .cal_sampling_method(object$emulator$posterior$sampling_method,
                                             object$emulator$posterior$joint_trajectories),
    predictor_uncertainty = object$predictor_uncertainty),
    class = "summary.spDBL_calibration")
}

#' @export
print.summary.spDBL_calibration <- function(x, ...) {
  cat(if (isTRUE(x$execution_only)) "Calibration execution experiment\n" else "Modular calibration\n")
  cat(x$n_samples, "retained samples from", nrow(x$diagnostics), "emulator conditionals\n")
  cat("Emulator coefficient sampling:", x$emulator_sampling, "\n")
  print(x$parameters, row.names = FALSE)
  cat("Divergences:", sum(x$diagnostics$divergences),
      "; tree-depth hits:", sum(x$diagnostics$max_treedepth_hits), "\n")
  cat("Rhat/ESS apply to each conditional fit, not to the pooled mixture.\n")
  if (any(is.finite(x$diagnostics$max_rhat))) {
    cat(sprintf("Largest conditional Rhat: %.3f; smallest conditional bulk ESS: %.0f",
        max(x$diagnostics$max_rhat, na.rm = TRUE), min(x$diagnostics$min_n_eff, na.rm = TRUE)))
    tail_ess <- x$diagnostics$min_ess_tail
    if (any(is.finite(tail_ess))) cat(sprintf(", tail ESS: %.0f", min(tail_ess, na.rm = TRUE)))
    cat("\n")
  }
  if (isTRUE(x$execution_only)) cat("One-transition/no-warmup samples do not establish convergence.\n")
  cat(x$predictor_uncertainty, "\n")
  invisible(x)
}

#' @export
print.spDBL_calibration <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Plot calibrated mechanistic parameter samples
#' @param x Calibration result.
#' @param ... Graphical arguments passed to hist.
#' @return Invisibly returns x.
#' @examplesIf interactive() && requireNamespace("rstan", quietly = TRUE)
#' data("dt_calibration")
#' d <- dt_calibration
#' em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
#'                      covariance = "sigma2R", episode_window = 2, nsam = 3,
#'                      seed = 1)
#' cal <- calibration_prepare(em, d$lower, d$upper)
#' fit <- calibration_fit(cal, d$field$z, d$field$U,
#'                        locations = d$field$locations, draws = 1:2,
#'                        chains = 1, iter = 200, seed = 1)
#' plot(fit)
#' @export
plot.spDBL_calibration <- function(x, ...) {
  active <- x$emulator$predictor$active
  old <- graphics::par(mfrow = c(1, length(active)))
  on.exit(graphics::par(old), add = TRUE)
  for (d in active)
    graphics::hist(x$mcmc_eta[, d], xlab = colnames(x$mcmc_eta)[d],
      main = if (isTRUE(x$execution_only)) "Execution samples" else "Calibration posterior",
      col = "lightblue", border = "white", ...)
  invisible(x)
}
