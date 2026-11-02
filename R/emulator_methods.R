# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Settings of a fitted emulator, tolerating objects saved by older versions.
.emulator_settings <- function(object) {
  s <- object$setup
  sampling <- s$sampling_method
  if (is.null(sampling)) sampling <- attr(object$res_ffbs, "sampling_method")
  if (is.null(sampling)) sampling <- "unknown"
  data.frame(covariance = if (is.null(s$covariance)) "iw" else s$covariance,
    sampling = sampling, n_train = s$n_train, n_inputs = ncol(s$pde_para_train),
    Nx = s$Nx, Ny = s$Ny, n_sites = s$N_sp, n_times = s$nT_ori,
    n_blocks = s$n_b, block_sites = s$bncol, n_episodes = s$nT,
    AR_order = s$AR_choice, design_columns = ncol(s$F_ls[[1]]),
    n_draws = s$nsam)
}

#' Print, summarize and predict from a fitted emulator
#'
#' Methods for the result of [emulator_learn()]. `print()` shows the model
#' settings. `summary()` adds the in-sample fit of the smoothed coefficients
#' and the posterior mean of the residual variance. `predict()` calls
#' [emulator_predict()].
#'
#' @param x,object An `spDBL_emulator` from [emulator_learn()].
#' @param newdata Matrix (or vector, for one run) of new mechanistic inputs.
#' @param dt_pde_test,F_new_ls,MC Passed to [emulator_predict()].
#' @param ... Unused.
#' @return `print()` returns `x` invisibly. `summary()` returns an object of
#'   class `summary.spDBL_emulator`, a list with
#'   \describe{
#'     \item{settings}{One-row data frame of model settings.}
#'     \item{rmse_by_time}{Root mean squared error of the smoothed fit
#'       \eqn{F_e \hat\Theta_e} to the training responses, at each time.}
#'     \item{rmse}{The same over all times.}
#'     \item{response_sd}{Standard deviation of the training responses.}
#'     \item{residual_variance}{Posterior mean of the residual variance,
#'       averaged over the sites of a block: \eqn{D/(n-C-1)} for the
#'       inverse-Wishart model and \eqn{D/(n-1)} times the mean diagonal of
#'       `R_episode` for the \eqn{\sigma^2 R} model, from the final filter.}
#'   }
#'   `predict()` returns what [emulator_predict()] returns.
#' @seealso [emulator_learn()], [emulator_predict()]
#' @examples
#' data("dt_emulation")
#' d <- dt_emulation
#' em <- emulator_learn(d$pde_para_train, d$dt_pde_train, Nx = 10, Ny = 10,
#'                      nsam = 5, seed = 1)
#' em
#' summary(em)
#' pred <- predict(em, d$pde_para_test, dt_pde_test = d$dt_pde_test)
#' dim(pred[[3]])
#' @name spDBL_emulator-methods
NULL

#' @rdname spDBL_emulator-methods
#' @export
print.spDBL_emulator <- function(x, ...) {
  s <- .emulator_settings(x)
  cat("spDBL emulator (", if (s$covariance == "iw") "inverse-Wishart" else
        "sigma2 R", " covariance, ", s$sampling, " coefficient sampling)\n", sep = "")
  cat(sprintf("  Training runs: %d with %d inputs\n", s$n_train, s$n_inputs))
  cat(sprintf("  Grid: %d x %d (%d sites) over %d times\n", s$Nx, s$Ny,
              s$n_sites, s$n_times))
  cat(sprintf("  Episodes: %d = %d times x %d blocks of %d sites\n",
              s$n_episodes, s$n_times, s$n_blocks, s$block_sites))
  cat(sprintf("  Design: AR(%d), %d columns\n", s$AR_order, s$design_columns))
  cat(sprintf("  Stored posterior draws: %d\n", s$n_draws))
  invisible(x)
}

#' @rdname spDBL_emulator-methods
#' @export
summary.spDBL_emulator <- function(object, ...) {
  s <- object$setup
  settings <- .emulator_settings(object)
  st <- object$para_ffbs$bs$st
  p <- dim(st)[1]
  C <- dim(st)[2]
  sq_error <- vapply(seq_len(s$nT), function(e) {
    fitted <- s$F_ls[[e]] %*% matrix(st[, , e], p, C)
    sum((s$Y[[e]] - fitted)^2)
  }, numeric(1))
  time <- (seq_len(s$nT) - 1L) %/% s$n_b + 1L
  cells <- s$n_train * C * s$n_b
  last <- object$para_ffbs$ff[[s$nT]]
  residual_variance <- if (settings$covariance == "iw") {
    if (last$nt > C + 1) mean(diag(as.matrix(last$Dt))) / (last$nt - C - 1) else NA_real_
  } else {
    if (last$nt > 1) as.numeric(last$Dt) / (last$nt - 1) *
      mean(diag(as.matrix(s$R_episode))) else NA_real_
  }
  structure(list(settings = settings,
    rmse_by_time = sqrt(as.numeric(tapply(sq_error, time, sum)) / cells),
    rmse = sqrt(sum(sq_error) / (cells * s$nT_ori)),
    response_sd = stats::sd(unlist(s$Y)),
    residual_variance = residual_variance),
    class = "summary.spDBL_emulator")
}

#' @export
print.summary.spDBL_emulator <- function(x, ...) {
  s <- x$settings
  cat("spDBL emulator summary\n")
  cat(sprintf("  Covariance: %s; sampling: %s; AR order: %d\n", s$covariance,
              s$sampling, s$AR_order))
  cat(sprintf("  %d runs x %d sites x %d times; %d episodes; %d stored draws\n",
              s$n_train, s$n_sites, s$n_times, s$n_episodes, s$n_draws))
  cat(sprintf("  In-sample RMSE of the smoothed fit: %.4g (response SD %.4g)\n",
              x$rmse, x$response_sd))
  cat(sprintf("  Posterior mean residual variance per site: %.4g\n",
              x$residual_variance))
  invisible(x)
}

#' @rdname spDBL_emulator-methods
#' @export
predict.spDBL_emulator <- function(object, newdata, dt_pde_test = NULL,
                                   F_new_ls = "default", MC = FALSE, ...) {
  emulator_predict(object, newdata, dt_pde_test = dt_pde_test,
                   F_new_ls = F_new_ls, MC = MC)
}

#' Print and predict from a prepared calibration object
#'
#' Methods for the result of [calibration_prepare()]. `print()` shows the
#' parameter bounds, the design GP and the emulator draws. `predict()` calls
#' [emulator_predict()], which needs no new simulator runs.
#'
#' @param x,object An `spDBL_calibration_input` from [calibration_prepare()].
#' @param newdata Matrix (or vector, for one input) of mechanistic inputs.
#' @param MC Passed to [emulator_predict()].
#' @param ... Unused.
#' @return `print()` returns `x` invisibly; `predict()` returns what
#'   [emulator_predict()] returns.
#' @seealso [calibration_prepare()], [emulator_predict()]
#' @examples
#' data("dt_calibration")
#' d <- dt_calibration
#' em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
#'   covariance = "sigma2R", episode_window = 2, nsam = 3, seed = 1)
#' cal <- calibration_prepare(em, d$lower, d$upper)
#' cal
#' pred <- predict(cal, d$validation$input[1:2, ])
#' dim(pred[[5]])
#' @name spDBL_calibration_input-methods
NULL

#' @rdname spDBL_calibration_input-methods
#' @export
print.spDBL_calibration_input <- function(x, ...) {
  f <- x$predictor
  names <- x$parameter_names
  if (is.null(names)) names <- paste0("eta", seq_along(f$lower))
  free <- seq_along(f$lower) %in% f$active
  sampling <- if (is.null(x$sampling_method)) "unknown" else x$sampling_method
  cat(sprintf("spDBL calibration input: %d emulator draws (%s sampling), %d sites, %d times\n",
              x$M, sampling, x$S, x$T))
  bounds <- data.frame(parameter = names, lower = f$lower, upper = f$upper,
    status = ifelse(free, "free", "fixed"), lengthscale = NA_real_)
  bounds$lengthscale[free] <- signif(f$lengthscale, 3)
  print(bounds, row.names = FALSE)
  cat(sprintf("Design GP: leave-one-out RMSE %.3g, nugget %.3g (lengthscales on the unit scale)\n",
              f$loo_rmse, f$nugget))
  cat(f$uncertainty, "\n")
  invisible(x)
}

#' @rdname spDBL_calibration_input-methods
#' @export
predict.spDBL_calibration_input <- function(object, newdata, MC = FALSE, ...) {
  emulator_predict(object, newdata, MC = MC)
}
