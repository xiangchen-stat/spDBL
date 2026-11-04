# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

#' Analytical WAIC for an MNIW multivariate regression posterior
#'
#' Computes the log pointwise predictive density and the posterior-variance
#' form of WAIC without posterior simulation. Each row of the response matrix
#' is one multivariate observation.
#'
#' @param Y Finite numeric response matrix with dimensions n by q.
#' @param X Finite numeric design matrix with dimensions n by p. Include an
#'   intercept column explicitly when required.
#' @param Mn Symmetric positive-definite p by p posterior row covariance of
#'   the regression coefficients, conditional on the response covariance.
#' @param mn Finite numeric p by q posterior coefficient mean. Supply the
#'   actual mean, not its precision-weighted parameterization.
#' @param Psi_n Symmetric positive-definite q by q inverse-Wishart posterior
#'   scale matrix.
#' @param nu_n Finite scalar inverse-Wishart posterior degrees of freedom,
#'   strictly greater than q - 1.
#'
#' @details
#' The likelihood and posterior are
#' \deqn{Y_i \mid B, \Sigma \sim N_q(x_i B, \Sigma),}
#' independently across rows, and
#' \deqn{B \mid \Sigma, D \sim MN(m_n, M_n, \Sigma),
#'       \quad \Sigma \mid D \sim IW(\nu_n, \Psi_n).}
#' The inverse-Wishart density is proportional to
#' \eqn{|\Sigma|^{-(\nu_n+q+1)/2}
#' \exp\{-\mathrm{tr}(\Psi_n\Sigma^{-1})/2\}}.
#' In particular, \code{mn} is the coefficient mean itself. When using
#' [MNIW_R()] with \code{H = diag(n)}, pass \code{Vbnew}, \code{Cnew},
#' \code{Snew}, and \code{vnew} as \code{Mn}, \code{mn}, \code{Psi_n}, and
#' \code{nu_n}, respectively.
#'
#' Let \eqn{h_i=x_i M_n x_i^T} and
#' \eqn{d_i=(y_i-x_i m_n)\Psi_n^{-1}(y_i-x_i m_n)^T}.
#' The row predictive density is multivariate Student t with
#' \eqn{\nu_n-q+1} degrees of freedom, mean \eqn{x_i m_n}, and scale
#' \eqn{(1+h_i)\Psi_n/(\nu_n-q+1)}. The WAIC penalty for row i is
#' \deqn{p_i = \tfrac14\left\{
#' \sum_{j=1}^q\psi_1((\nu_n-j+1)/2) + 2q h_i^2 +
#' 4h_i\nu_n d_i + 2\nu_n d_i^2 - 4d_i\right\}.}
#' Here \eqn{\psi_1} is the trigamma function. WAIC is
#' \eqn{-2\sum_i(\mathrm{lppd}_i-p_i)}.
#'
#' Use the observations that fitted the posterior to compute training-data
#' WAIC. Scoring other observations evaluates the supplied posterior at those
#' observations; it does not refit the model. These functions assume independent
#' likelihood rows and do not directly compute WAIC for correlated spatial GP
#' rows or a full state-space emulator. Compare models using the same response
#' data and the same definition of an observation. Lower WAIC indicates higher
#' estimated predictive accuracy among the compared models.
#'
#' @return A list containing the original scalar summaries \code{WAIC},
#'   \code{lppd}, and \code{p_waic}; their difference \code{elpd_waic}; and an
#'   n by 4 \code{pointwise} matrix with columns \code{lppd}, \code{p_waic},
#'   \code{elpd_waic}, and \code{WAIC}.
#' @author Sudipto Banerjee (original derivation and development implementation).
#' @references
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' Journal of Machine Learning Research, 11, 3571--3594.
#' @seealso [compute_sample_waic()], [MNIW_R()], [lppd_IW_1t()]
#' @examples
#' X <- cbind(1, c(-1, 0, 1, 2))
#' Y <- cbind(c(0.1, 0.4, 1.2, 1.7), c(1.2, 0.8, 0.5, 0.1))
#' post <- MNIW_R(X, Y, H = diag(nrow(Y)), v = 4, S = diag(2),
#'                C = matrix(0, 2, 2), Vb = diag(2))
#' score <- compute_analytical_waic(Y, X, post$Vbnew, post$Cnew,
#'                                  post$Snew, post$vnew)
#' score$WAIC
#' score$pointwise
#' @export
compute_analytical_waic <- function(Y, X, Mn, mn, Psi_n, nu_n) {
  dims <- .mniw_waic_data(Y, X)
  U_Mn <- .mniw_waic_chol(Mn, dims$p, "Mn")
  .mniw_waic_matrix(mn, "mn", c(dims$p, dims$q))
  U_Psi <- .mniw_waic_chol(Psi_n, dims$q, "Psi_n")
  if (!is.numeric(nu_n) || length(nu_n) != 1L || !is.finite(nu_n) ||
      nu_n <= dims$q - 1)
    stop("nu_n must be a finite scalar greater than q - 1", call. = FALSE)

  # Norms of X times a covariance factor keep h nonnegative numerically.
  h <- rowSums(tcrossprod(X, U_Mn)^2)
  residual <- Y - X %*% mn
  whitened <- residual %*% backsolve(U_Psi, diag(dims$q))
  d <- rowSums(whitened^2)
  log_det_Psi <- 2 * sum(log(diag(U_Psi)))
  lppd <- lgamma((nu_n + 1) / 2) -
    lgamma((nu_n - dims$q + 1) / 2) - dims$q / 2 * log(pi) -
    (dims$q * log1p(h) + log_det_Psi) / 2 -
    (nu_n + 1) / 2 * log1p(d / (1 + h))
  var_log_det <- sum(trigamma((nu_n - seq_len(dims$q) + 1) / 2))
  penalty <- (var_log_det + 2 * dims$q * h^2 +
    4 * h * nu_n * d + 2 * nu_n * d^2 - 4 * d) / 4
  if (any(!is.finite(lppd)) || any(!is.finite(penalty)))
    stop("WAIC calculation overflowed; rescale the data or posterior parameters", call. = FALSE)
  .mniw_waic_result(lppd, penalty, rownames(Y))
}

#' WAIC from paired coefficient and covariance draws
#'
#' Evaluates row-wise multivariate normal log likelihoods for paired posterior
#' draws, using log-sum-exp aggregation for the predictive density and sample
#' variances for the WAIC penalty.
#'
#' @inheritParams compute_analytical_waic
#' @param B_samples Finite numeric array with dimensions p by q by S,
#'   containing posterior coefficient draws. At least two draws are required.
#' @param Sigma_samples Finite numeric array with dimensions q by q by S,
#'   containing symmetric positive-definite posterior response covariance
#'   draws. Slice s must be paired with slice s of \code{B_samples}.
#' @param return_log_lik Logical scalar. If \code{TRUE}, also return the S by n
#'   \code{log_lik} matrix, suitable for \code{loo::waic()} or \code{loo::loo()}.
#'   The loo package is optional and is not called by this function.
#'
#' @details Each response row is one observation, with likelihood
#'   \eqn{Y_i\mid B,\Sigma\sim N_q(x_i B,\Sigma)} independently across rows.
#'   The draws can come from any posterior with this likelihood; they need not
#'   follow an MNIW law. The analytical comparison requires an MNIW posterior.
#'   Arrays retain three dimensions even when p, q, or n equals one.
#'   The penalty uses sample variance with denominator S - 1, matching
#'   \code{loo::waic()} for a matrix of log likelihoods. Correlated draws can
#'   increase Monte Carlo error; this function does not diagnose convergence.
#' @return The same summaries and pointwise matrix as
#'   [compute_analytical_waic()], with an additional \code{log_lik} matrix when
#'   \code{return_log_lik = TRUE}.
#' @author Sudipto Banerjee (original development implementation).
#' @seealso [compute_analytical_waic()]
#' @examples
#' X <- matrix(1, 2, 1)
#' Y <- matrix(c(0.1, 0.3), 2, 1)
#' B <- array(c(0, 0.2, 0.4), c(1, 1, 3))
#' Sigma <- array(c(0.5, 1, 1.5), c(1, 1, 3))
#' score <- compute_sample_waic(Y, X, B, Sigma, return_log_lik = TRUE)
#' score$WAIC
#' dim(score$log_lik)
#' @export
compute_sample_waic <- function(Y, X, B_samples, Sigma_samples,
                                return_log_lik = FALSE) {
  dims <- .mniw_waic_data(Y, X)
  if (!is.numeric(B_samples) || length(dim(B_samples)) != 3L ||
      any(!is.finite(B_samples)) ||
      !identical(dim(B_samples)[1:2], c(dims$p, dims$q)) ||
      dim(B_samples)[3L] < 2L)
    stop("B_samples must be a finite p by q by S array with S >= 2", call. = FALSE)
  S <- dim(B_samples)[3L]
  if (!is.numeric(Sigma_samples) ||
      !identical(dim(Sigma_samples), c(dims$q, dims$q, S)) ||
      any(!is.finite(Sigma_samples)))
    stop("Sigma_samples must be a finite q by q by S array paired with B_samples", call. = FALSE)
  if (!is.logical(return_log_lik) || length(return_log_lik) != 1L ||
      is.na(return_log_lik))
    stop("return_log_lik must be TRUE or FALSE", call. = FALSE)

  I_q <- diag(dims$q)
  log_lik <- do.call(rbind, lapply(seq_len(S), function(s) {
    beta <- matrix(B_samples[, , s], dims$p, dims$q)
    sigma <- matrix(Sigma_samples[, , s], dims$q, dims$q)
    U <- .mniw_waic_chol(sigma, dims$q,
      paste0("Sigma_samples[, , ", s, "]"))
    whitened <- (Y - X %*% beta) %*% backsolve(U, I_q)
    -dims$q / 2 * log(2 * pi) - sum(log(diag(U))) -
      rowSums(whitened^2) / 2
  }))
  if (any(!is.finite(log_lik)))
    stop("Log-likelihood calculation overflowed; rescale the data or draws", call. = FALSE)
  colnames(log_lik) <- rownames(Y)
  max_ll <- apply(log_lik, 2L, max)
  lppd <- max_ll + log(colMeans(exp(sweep(log_lik, 2L, max_ll, "-"))))
  penalty <- apply(log_lik, 2L, stats::var)
  result <- .mniw_waic_result(lppd, penalty, rownames(Y))
  if (return_log_lik) result$log_lik <- log_lik
  result
}

.mniw_waic_matrix <- function(x, name, expected_dim = NULL) {
  if (!is.matrix(x) || !is.numeric(x) || any(dim(x) < 1L) ||
      any(!is.finite(x)))
    stop(name, " must be a nonempty finite numeric matrix", call. = FALSE)
  if (!is.null(expected_dim) && !identical(dim(x), expected_dim))
    stop(name, " must have dimensions ", paste(expected_dim, collapse = " by "), call. = FALSE)
  invisible(NULL)
}

.mniw_waic_data <- function(Y, X) {
  .mniw_waic_matrix(Y, "Y")
  .mniw_waic_matrix(X, "X")
  if (nrow(Y) != nrow(X))
    stop("Y and X must have the same number of rows", call. = FALSE)
  list(n = nrow(Y), p = ncol(X), q = ncol(Y))
}

.mniw_waic_chol <- function(x, n, name) {
  .mniw_waic_matrix(x, name, c(n, n))
  if (!isSymmetric(x, tol = 1e-10, check.attributes = FALSE))
    stop(name, " must be symmetric positive definite", call. = FALSE)
  tryCatch(chol(x), error = function(e)
    stop(name, " must be symmetric positive definite", call. = FALSE))
}

.mniw_waic_result <- function(lppd, penalty, labels) {
  elpd <- lppd - penalty
  pointwise <- cbind(lppd = lppd, p_waic = penalty,
    elpd_waic = elpd, WAIC = -2 * elpd)
  rownames(pointwise) <- labels
  list(WAIC = sum(pointwise[, "WAIC"]), lppd = sum(lppd),
    p_waic = sum(penalty), elpd_waic = sum(elpd), pointwise = pointwise)
}
