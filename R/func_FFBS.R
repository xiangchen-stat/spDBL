# MNIW----

#' Forward Filter Backward Sampler (MNIW model)
#'
#' Runs the complete FFBS algorithm under the Matrix Normal Inverse Wishart
#' (MNIW) model: first applies the forward filter (\code{\link{FF}}) and then
#' the backward sampler (\code{\link{BS}}).
#'
#' @param Y List of length \code{nT}. Each element is the \eqn{N \times q}
#'   data matrix at the corresponding time step.
#' @param F_ls Covariate matrix or list of matrices (see \code{\link{FF}}).
#' @param G_ls State transition matrix or list (see \code{\link{FF}}).
#' @param W_ls State noise left-covariance matrix or list (see \code{\link{FF}}).
#' @param V_ls Observation noise left-covariance matrix or list
#'   (see \code{\link{FF}}).
#' @param m0 Numeric matrix. Prior state mean (\eqn{p \times q}).
#' @param M0 Numeric matrix. Prior state left-covariance (\eqn{p \times p}).
#' @param n0 Numeric scalar. Prior degrees of freedom of the inverse-Wishart.
#' @param D0 Numeric matrix. Prior scale matrix of the inverse-Wishart.
#' @param nT Integer. Number of time steps.
#' @param delta Numeric scalar. Right-variance discount factor. Defaults to
#'   \code{1.0}.
#' @param verbose Logical. If \code{TRUE}, emits progress messages. Defaults
#'   to \code{FALSE}.
#'
#' @return A named list with:
#'   \describe{
#'     \item{ff}{Output of \code{\link{FF}}: filtered distributions for each
#'       time step.}
#'     \item{bs}{Output of \code{\link{BS}}: smoothed state means (\code{st})
#'       and left-covariances (\code{St}).}
#'   }
#'
#' @seealso \code{\link{FF}}, \code{\link{BS}}, \code{\link{FFBS_sampling}}
#' @export
FFBS <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, delta = 1.0,
                 verbose = FALSE){
  # TODO FF_cpp, BS_cpp
  out <- list()
  res_ff <- FF(Y = Y, F_ls = F_ls, G_ls = G_ls,
              W_ls = W_ls, V_ls = V_ls,
              m0 = m0, M0 = M0,
              n0 = n0, D0 = D0,
              nT = nT, delta = delta, verbose = verbose)
  res_bs <- BS(res_ff, G_ls, nT = nT, delta = delta, verbose = verbose)
  if(verbose){
    message("Finish FFBS")
  }
  out$ff <- res_ff
  out$bs <- res_bs
  return(out)
}


#' Draw posterior samples from FFBS output (MNIW model)
#'
#' Given the filtered and smoothed distributions from \code{\link{FFBS}},
#' draws \code{nsam} posterior sample sets of the state matrices
#' \eqn{\Theta_1, \ldots, \Theta_T} and the covariance \eqn{\Sigma} using the
#' MNIW sampler.
#'
#' @param nsam Integer. Number of posterior samples to draw.
#' @param para_ffbs List. Output of \code{\link{FFBS}}, containing elements
#'   \code{ff} and \code{bs}.
#' @param F_ls Covariate matrix or list of matrices (see \code{\link{FF}}).
#' @param G_ls State transition matrix or list (see \code{\link{FF}}).
#' @param nT Integer. Number of time steps.
#' @param delta Numeric scalar. Discount factor. Defaults to \code{1}.
#' @param verbose Logical. If \code{TRUE}, emits progress messages. Defaults
#'   to \code{FALSE}.
#' @param sampling "marginal" (default) samples the smoothed marginal at each
#'   episode conditional on the same covariance draw. "joint" samples backward
#'   conditionals to obtain coherent trajectories.
#' @param seed Optional local RNG seed; NULL follows the session RNG.
#'
#' @return A named list of length \code{nT + 1}. Elements \code{"T1"} through
#'   \code{"T<nT>"} are arrays of dimension \code{c(p, q, nsam)} containing
#'   posterior samples of \eqn{\Theta_t}. The element \code{Sigma} is an array
#'   of dimension \code{c(q, q, nsam)} containing posterior samples of
#'   \eqn{\Sigma}.
#'   Attributes sampling_method and joint_trajectories identify the chosen law.
#'
#' @seealso \code{\link{FFBS}}, \code{\link{MNIW_sampler}}
#' @export
FFBS_sampling <- function(nsam, para_ffbs, F_ls, G_ls, nT, delta = 1,
    verbose = FALSE, sampling = c("marginal", "joint"), seed = NULL) {
  sampling <- match.arg(sampling)
  restore_seed <- .cal_seed(seed)
  on.exit(restore_seed(), add = TRUE)
  ff <- para_ffbs$ff
  E <- length(ff) - as.integer(!is.null(ff$prior))
  if (length(nsam) != 1L || !is.finite(nsam) || nsam < 1 ||
      nsam != as.integer(nsam) || length(nT) != 1L || nT != E || E < 1L)
    stop("Invalid sample count or episode count")
  if (delta != 1) stop("The IW samplers currently require delta = 1")
  p <- nrow(ff[[E]]$mt)
  C <- ncol(ff[[E]]$mt)
  Sigma <- matrixsampling::rinvwishart(nsam, nu = ff[[E]]$nt,
    Omega = check_pds(ff[[E]]$Dt), epsilon = 0, checkSymmetry = FALSE)
  Sigma <- array(Sigma, c(C, C, nsam))
  L_Sigma <- lapply(seq_len(nsam), function(m)
    t(.cal_spd(matrix(Sigma[, , m], C, C))))
  theta <- lapply(seq_len(E), function(e) array(NA_real_, c(p, C, nsam)))
  matrix_normal <- function(mean, L, m)
    mean + L %*% matrix(stats::rnorm(p * C), p, C) %*% t(L_Sigma[[m]])
  L <- t(.cal_spd(ff[[E]]$Mt))
  for (m in seq_len(nsam))
    theta[[E]][, , m] <- matrix_normal(ff[[E]]$mt, L, m)
  if (E > 1L) for (e in seq.int(E - 1L, 1L)) {
    if (sampling == "marginal") {
      mean <- matrix(para_ffbs$bs$st[, , e], p, C)
      H <- matrix(para_ffbs$bs$St[, , e], p, p)
    } else {
      G <- if (is.list(G_ls)) G_ls[[e + 1L]] else G_ls
      B <- ff[[e]]$Mt %*% t(G) %*% chol2inv(.cal_spd(ff[[e + 1L]]$At))
      H <- ff[[e]]$Mt - B %*% ff[[e + 1L]]$At %*% t(B)
    }
    L <- t(.cal_spd((H + t(H)) / 2))
    for (m in seq_len(nsam)) {
      if (sampling == "joint")
        mean <- ff[[e]]$mt + B %*% (theta[[e + 1L]][, , m] - ff[[e + 1L]]$at)
      theta[[e]][, , m] <- matrix_normal(mean, L, m)
    }
  }
  names(theta) <- paste0("T", seq_len(E))
  theta$Sigma <- Sigma
  attr(theta, "joint_trajectories") <- sampling == "joint"
  attr(theta, "sampling_method") <- sampling
  if (verbose) message("Sampled ", nsam, " ", sampling, " coefficient sets")
  theta
}



#' Monte Carlo prediction using FFBS output (MNIW model)
#'
#' Estimates posterior predictive means at new spatial locations using Monte
#' Carlo integration over the posterior samples of the state \eqn{\Theta_t}.
#' Uses an exponential GP kernel to compute the cross-covariance between
#' observed and new locations.
#'
#' @param nsam Integer. Number of posterior samples to average over.
#' @param Y List of length \code{nT}. Observed data matrices.
#' @param res_ffbs List. Posterior samples of \eqn{\Theta_t}, as returned by
#'   \code{\link{FFBS_sampling}}.
#' @param F_ls Covariate matrix or list for observed locations.
#' @param F_new_ls Covariate matrix or list for new (prediction) locations.
#' @param input Numeric matrix. Coordinates of observed locations
#'   (\eqn{N \times d}).
#' @param input_new Numeric matrix or vector. Coordinates of new locations.
#' @param nT Integer. Number of time steps.
#' @param phi_para Numeric scalar. Range parameter of the exponential kernel.
#' @param gp_sigma2 Numeric scalar. Variance parameter of the GP kernel.
#'   Defaults to \code{1.1}.
#' @param gp_tau2 Numeric scalar. Nugget variance of the GP kernel. Defaults to
#'   \code{1e-4}.
#' @param delta Numeric scalar. Discount factor. Defaults to \code{1.0}.
#'
#' @return A named list of length \code{nT}. Each element \code{"T<t>"} is an
#'   array of dimension \code{c(N_new, q, nsam)} with posterior predictive mean
#'   samples at the new locations.
#'
#' @seealso \code{\link{FFBS_predict_exact}}, \code{\link{gen_exp_kernel}}
#' @export
FFBS_predict_MC <- function(nsam, Y, res_ffbs, F_ls, F_new_ls, input, input_new,
    nT, phi_para, gp_sigma2 = 1.1, gp_tau2 = 1e-4, delta = 1.0) {
  if (is.null(dim(input_new))) input_new <- matrix(input_new, nrow = 1)
  input <- as.matrix(input)
  V <- gp_sigma2 * exp(-phi_para * as.matrix(stats::dist(input))) +
    diag(gp_tau2, nrow(input))
  R_V <- chol(V)
  out <- lapply(seq_len(nT), function(e)
    array(NA_real_, c(nrow(input_new), ncol(Y[[e]]), nsam)))
  for (m in seq_len(nsam)) {
    theta <- lapply(seq_len(nT), function(e) {
      th <- res_ffbs[[e]]
      matrix(th[, , m], dim(th)[1], dim(th)[2])
    })
    pred <- FFBS_predict_exact(Y, F_ls = F_ls, F_new_ls = F_new_ls,
      input = input, input_new = input_new, nT = nT, phi_para = phi_para,
      gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2, theta = theta, V_chol = R_V)
    for (e in seq_len(nT)) out[[e]][, , m] <- pred[[e]]
  }
  names(out) <- paste0("T", seq_len(nT))
  out
}

#' Exact posterior predictive mean using FFBS smoothed states (MNIW model)
#'
#' Computes the analytical posterior predictive mean at new spatial locations
#' using the smoothed state means \eqn{s_t} from \code{\link{FFBS}} and an
#' exponential GP kernel for the cross-covariance.
#'
#' @param Y List of length \code{nT}. Observed data matrices.
#' @param para_ffbs List. Output of \code{\link{FFBS}}, containing \code{ff}
#'   and \code{bs}.
#' @param F_ls Covariate matrix or list for observed locations.
#' @param F_new_ls Covariate matrix or list for new locations.
#' @param input Numeric matrix. Coordinates of observed locations.
#' @param input_new Numeric matrix or vector. Coordinates of new locations.
#' @param nT Integer. Number of time steps.
#' @param phi_para Numeric scalar. Range parameter of the exponential kernel.
#' @param gp_sigma2 Numeric scalar. Variance parameter of the GP kernel.
#'   Defaults to \code{1.1}.
#' @param gp_tau2 Numeric scalar. Nugget variance of the GP kernel. Defaults to
#'   \code{1e-4}.
#' @param delta Numeric scalar. Discount factor. Defaults to \code{1.0}.
#'
#' @return A named list of length \code{nT}. Each element \code{"T<t>"} is a
#'   matrix of dimension \code{c(N_new, q)} with the posterior predictive mean
#'   at the new locations.
#'
#' @seealso \code{\link{FFBS_predict_MC}}, \code{\link{gen_exp_kernel}}
#' @param theta Optional list of fixed coefficient matrices, one per episode.
#' @param V_chol Optional cached upper Cholesky factor of the input covariance.
#' @param residual_weights Optional cached inverse-covariance-weighted residuals
#'   for the supplied coefficient trajectory.
#' @param return_covariance If TRUE, return a list with mean and
#'   input_covariance before scaling by the emulator variance and spatial covariance.
#' @export
FFBS_predict_exact <- function(Y, para_ffbs = NULL, F_ls, F_new_ls, input, input_new,
                               nT, phi_para, gp_sigma2 = 1.1, gp_tau2 = 10^(-4),
                               delta = 1.0, theta = NULL, V_chol = NULL,
                               residual_weights = NULL, return_covariance = FALSE){
  input <- as.matrix(input)
  if (is.null(dim(input_new))) input_new <- matrix(input_new, nrow = 1)
  input_new <- as.matrix(input_new)
  N <- nrow(input)
  N_new <- nrow(input_new)
  if (length(nT) != 1L || !is.finite(nT) || nT < 1 || nT != as.integer(nT) ||
      length(Y) != nT || N < 1 || N_new < 1 || ncol(input) != ncol(input_new) ||
      any(!is.finite(c(input, input_new))) ||
      length(phi_para) != 1L || !is.finite(phi_para) || phi_para <= 0 ||
      length(gp_sigma2) != 1L || !is.finite(gp_sigma2) || gp_sigma2 <= 0 ||
      length(gp_tau2) != 1L || !is.finite(gp_tau2) || gp_tau2 < 0) {
    stop("Invalid prediction inputs or kernel parameters")
  }
  if (is.null(theta)) {
    st <- para_ffbs$bs$st
    if (length(dim(st)) != 3L || dim(st)[3] != nT)
      stop("Supply smoothed coefficients in para_ffbs or a theta trajectory")
    theta <- lapply(seq_len(nT), function(i) matrix(st[, , i], dim(st)[1], dim(st)[2]))
  }
  if (!is.list(theta) || length(theta) != nT ||
      (!is.null(residual_weights) && length(residual_weights) != nT))
    stop("Coefficient and residual trajectories must have nT entries")

  # The new-run nugget is independent of training-run noise, including when
  # input_new coincides with a training input. J contains no nugget.
  distances <- as.matrix(stats::dist(rbind(input, input_new)))
  Jt <- gp_sigma2 * exp(-phi_para * distances[seq_len(N), N + seq_len(N_new), drop = FALSE])
  if (is.null(V_chol)) {
    Vt <- gp_sigma2 * exp(-phi_para * distances[seq_len(N), seq_len(N), drop = FALSE]) +
      diag(gp_tau2, N)
    V_chol <- chol(Vt)
  }
  if (!identical(dim(V_chol), c(N, N))) stop("V_chol has incompatible dimensions")
  design_at <- function(design, i, n) {
    if (is.list(design)) {
      if (!length(design) %in% c(1L, nT)) stop("Design lists must have one or nT entries")
      design <- design[[if (length(design) == 1L) 1L else i]]
    }
    if (is.null(dim(design))) design <- matrix(design, nrow = n)
    as.matrix(design)
  }
  post_mean_ls <- vector("list", nT)
  for (i in seq_len(nT)) {
    Ft <- design_at(F_ls, i, N)
    Ft_new <- design_at(F_new_ls, i, N_new)
    Yt <- as.matrix(Y[[i]])
    th <- as.matrix(theta[[i]])
    if (nrow(Ft) != N || nrow(Ft_new) != N_new || nrow(Yt) != N ||
        ncol(Ft) != nrow(th) || ncol(Ft_new) != nrow(th) || ncol(Yt) != ncol(th))
      stop("Incompatible design, response, or coefficient dimensions at episode ", i)
    weights <- if (is.null(residual_weights))
      backsolve(V_chol, forwardsolve(t(V_chol), Yt - Ft %*% th)) else residual_weights[[i]]
    if (!identical(dim(weights), c(N, ncol(Yt))))
      stop("Residual weights have incompatible dimensions at episode ", i)
    post_mean_ls[[i]] <- Ft_new %*% th + crossprod(Jt, weights)
  }
  names(post_mean_ls) <- paste0("T", seq_len(nT))
  if (!return_covariance) return(post_mean_ls)
  V_new <- gp_sigma2 * exp(-phi_para *
    distances[N + seq_len(N_new), N + seq_len(N_new), drop = FALSE]) + diag(gp_tau2, N_new)
  covariance <- V_new - crossprod(forwardsolve(t(V_chol), Jt))
  list(mean = post_mean_ls, input_covariance = (covariance + t(covariance)) / 2)
}

# sigma2R----

#' Forward Filter Backward Sampler (scalar sigma-squared-times-R model)
#'
#' Runs the complete FFBS algorithm under the model where the right covariance
#' is \eqn{\sigma^2 R} with \eqn{\sigma^2 \sim IG(n_0, d_0)}.
#'
#' @param Y List of length \code{nT}. Data matrices.
#' @param F_ls Covariate matrix or list (see \code{\link{FF}}).
#' @param G_ls State transition matrix or list (see \code{\link{FF}}).
#' @param W_ls State noise left-covariance matrix or list (see \code{\link{FF}}).
#' @param V_ls Observation noise left-covariance matrix or list
#'   (see \code{\link{FF}}).
#' @param m0 Numeric matrix. Prior state mean.
#' @param M0 Numeric matrix. Prior state left-covariance.
#' @param n0 Numeric scalar. Prior shape of \eqn{\sigma^2}.
#' @param D0 Numeric scalar. Prior rate of \eqn{\sigma^2}.
#' @param nT Integer. Number of time steps.
#' @param R Numeric matrix. Fixed spatial correlation matrix.
#' @param delta Numeric scalar. Discount factor. Defaults to \code{1.0}.
#' @param verbose Logical. If \code{TRUE}, emits progress messages. Defaults
#'   to \code{FALSE}.
#'
#' @return A named list with \code{ff} (output of \code{\link{FF_sigma2R}}) and
#'   \code{bs} (output of \code{\link{BS}}).
#'
#' @seealso \code{\link{FF_sigma2R}}, \code{\link{BS}},
#'   \code{\link{FFBS_sampling_sigma2R}}
#' @export
FFBS_sigma2R <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, R,
                         delta = 1.0, verbose = FALSE){
  # TODO FF_cpp, BS_cpp
  out <- list()
  res_ff <- FF_sigma2R(Y = Y, F_ls = F_ls, G_ls = G_ls,
                        W_ls = W_ls, V_ls = V_ls,
                        m0 = m0, M0 = M0,
                        n0 = n0, D0 = D0,
                        nT = nT, R = R, delta = delta, verbose = verbose)
  res_bs <- BS(res_ff, G_ls, nT = nT, delta = delta, verbose = verbose)
  if(verbose){
    message("Finish FFBS")
  }
  out$ff <- res_ff
  out$bs <- res_bs
  return(out)
}


#' Draw posterior samples from FFBS output (scalar sigma-squared-times-R model)
#'
#' Draws \code{nsam} posterior samples of the state \eqn{\Theta_t} and the
#' scalar variance \eqn{\sigma^2} using the MNIG sampler, given the FFBS output
#' under the scalar-sigma model.
#'
#' @param nsam Integer. Number of posterior samples.
#' @param para_ffbs List. Output of \code{\link{FFBS_sigma2R}}.
#' @param F_ls Covariate matrix or list (see \code{\link{FF}}).
#' @param G_ls State transition matrix or list (see \code{\link{FF}}).
#' @param nT Integer. Number of time steps.
#' @param R Numeric matrix. Fixed spatial correlation matrix.
#' @param delta Numeric scalar. Discount factor. Defaults to \code{1}.
#' @param verbose Logical. If \code{TRUE}, emits progress messages. Defaults
#'   to \code{FALSE}.
#'
#' @return A named list of length \code{nT + 1}. Elements \code{"T1"} through
#'   \code{"T<nT>"} are arrays of dimension \code{c(p, q, nsam)} with state
#'   samples. The element \code{Sigma} is a vector of length \code{nsam} with
#'   samples of \eqn{\sigma^2}.
#'   Attributes sampling_method and joint_trajectories identify the chosen law.
#'
#' @seealso \code{\link{FFBS_sigma2R}}, \code{\link{MNIG_sampler}}
#' @param seed Optional local seed. NULL uses the session RNG.
#' @param sampling "marginal" (default) draws smoothed episode marginals
#'   conditional on a shared variance. "joint" preserves cross-episode
#'   dependence by sampling backward conditionals.
#' @export
FFBS_sampling_sigma2R <- function(nsam, para_ffbs, F_ls, G_ls, nT, R,
    delta = 1, verbose = FALSE, seed = NULL, sampling = c("marginal", "joint")) {
  sampling <- match.arg(sampling)
  if (delta != 1) stop("The sigma2R samplers currently require delta = 1")
  if (nT != length(para_ffbs$ff) - as.integer(!is.null(para_ffbs$ff$prior)))
    stop("nT does not match the filtering results")
  sampled <- sample_calibration_emulator(para_ffbs, R, G_ls, nsam, seed, sampling)
  out <- sampled$theta
  names(out) <- paste0("T", seq_len(nT))
  out$Sigma <- sampled$sigma2
  attr(out, "joint_trajectories") <- sampled$joint_trajectories
  attr(out, "sampling_method") <- sampling
  if (verbose) message("Sampled ", nsam, " ", sampling, " coefficient sets")
  out
}


# I-----

#' Forward Filter Backward Sampler (identity right-covariance)
#'
#' Runs the complete FFBS algorithm with the right covariance fixed at the
#' identity matrix (no inverse-Wishart update for \eqn{\Sigma}).
#'
#' @param Y List of length \code{nT}. Data matrices.
#' @param F_ls Covariate matrix or list (see \code{\link{FF}}).
#' @param G_ls State transition matrix or list (see \code{\link{FF}}).
#' @param W_ls State noise left-covariance matrix or list (see \code{\link{FF}}).
#' @param V_ls Observation noise left-covariance matrix or list
#'   (see \code{\link{FF}}).
#' @param m0 Numeric matrix. Prior state mean.
#' @param M0 Numeric matrix. Prior state left-covariance.
#' @param nT Integer. Number of time steps.
#' @param delta Numeric scalar. Discount factor. Defaults to \code{1.0}.
#' @param verbose Logical. If \code{TRUE}, emits progress messages. Defaults
#'   to \code{FALSE}.
#'
#' @return A named list with \code{ff} (output of \code{\link{FF_I}}) and
#'   \code{bs} (output of \code{\link{BS}}).
#'
#' @seealso \code{\link{FF_I}}, \code{\link{BS}}, \code{\link{FFBS_sampling_I}}
#' @export
FFBS_I <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, nT, delta = 1.0,
                   verbose = FALSE){
  # TODO FF_cpp, BS_cpp
  out <- list()
  res_ff <- FF_I(Y = Y, F_ls = F_ls, G_ls = G_ls,
                                W_ls = W_ls, V_ls = V_ls,
                                m0 = m0, M0 = M0,
                                nT = nT, delta = delta, verbose = verbose)
  res_bs <- BS(res_ff, G_ls, nT = nT, delta = delta, verbose = verbose)
  if(verbose){
    message("Finish FFBS")
  }
  out$ff <- res_ff
  out$bs <- res_bs
  return(out)
}

#' Draw posterior samples from FFBS output (identity right-covariance)
#'
#' Draws \code{nsam} posterior samples of the state \eqn{\Theta_t} from the
#' smoothed distributions produced by \code{\link{FFBS_I}}, assuming an
#' identity right-covariance matrix.
#'
#' @param nsam Integer. Number of posterior samples.
#' @param para_ffbs List. Output of \code{\link{FFBS_I}}.
#' @param F_ls Covariate matrix or list (see \code{\link{FF}}).
#' @param G_ls State transition matrix or list (see \code{\link{FF}}).
#' @param nT Integer. Number of time steps.
#' @param delta Numeric scalar. Discount factor. Defaults to \code{1}.
#' @param verbose Logical. If \code{TRUE}, emits progress messages. Defaults
#'   to \code{FALSE}.
#'
#' @return A named list of length \code{nT}. Each element \code{"T<t>"} is an
#'   array of dimension \code{c(p, q, nsam)} with posterior samples of the
#'   state at time \eqn{t}.
#'
#' @seealso \code{\link{FFBS_I}}
#' @export
FFBS_sampling_I <- function(nsam, para_ffbs, F_ls, G_ls, nT, delta = 1,
                            verbose = FALSE){
  out <- list()
  res_ff <- para_ffbs$ff
  res_bs <- para_ffbs$bs

  # Firstly, draw sigma and B at t = nT
  # initialize matirces at t = nT
  Gt1 <- G_ls
  if(is.list(F_ls)){
    Ft <- as.matrix(F_ls[[nT]])
  } else{
    Ft <- F_ls
  }

  dim_array <- c(dim(res_bs$st[,,nT]), nsam)
  IS <- diag(dim_array[2])

  # backward nT
  for (i in nT:1) {
    if(verbose && i %% round(0.5*(nT-1)) == 0){
      message("Sampling: ", i, "/", nT)
    }
    # Check if G_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(G_ls)){
      Gt1 <- as.matrix(G_ls[[i+1]])
    }

    out[[i]] <- array(dim = dim_array)
    para_nt <- res_ff[[i]]
    st <- res_bs$st[,,i]
    St <- check_pds(res_bs$St[,,i])

    # get nsam samples for different Sigma
    for (j in 1:nsam) {
      out[[i]][,,j] <- rmn_cpp(m = st, U = St, V = IS)
    }
    if(i == 1){
      names(out) <- paste0("T", 1:nT)
    }
  }

  return(out)
}
