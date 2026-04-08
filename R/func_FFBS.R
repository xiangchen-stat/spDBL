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
FFBS <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, delta = 1.0){
  # TODO FF_cpp, BS_cpp
  out <- list()
  res_ff <- FF(Y = Y, F_ls = F_ls, G_ls = G_ls,
              W_ls = W_ls, V_ls = V_ls,
              m0 = m0, M0 = M0,
              n0 = n0, D0 = D0,
              nT = nT, delta = delta)
  res_bs <- BS(res_ff, G_ls, nT = nT, delta = delta)
  print("Finish FFBS")
  out$ff <- res_ff
  out$bs <- res_bs
  return(out)
}


#' Draw posterior samples from FFBS output (MNIW model)
#'
#' Given the filtered and smoothed distributions from \code{\link{FFBS}},
#' draws \code{nsam} joint posterior samples of the state matrices
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
#'
#' @return A named list of length \code{nT + 1}. Elements \code{"T1"} through
#'   \code{"T<nT>"} are arrays of dimension \code{c(p, q, nsam)} containing
#'   posterior samples of \eqn{\Theta_t}. The element \code{Sigma} is an array
#'   of dimension \code{c(q, q, nsam)} containing posterior samples of
#'   \eqn{\Sigma}.
#'
#' @seealso \code{\link{FFBS}}, \code{\link{MNIW_sampler}}
#' @export
FFBS_sampling <- function(nsam, para_ffbs, F_ls, G_ls, nT, delta = 1){
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

  para_nt <- res_ff[[nT]]
  para_nt$Dt <- check_pds(para_nt$Dt)
  res_bs$St[,,nT] <- check_pds(res_bs$St[,,nT])
  post_nT <- MNIW_sampler(nsam = nsam, X = Ft, v = para_nt$nt,
                          S = para_nt$Dt, C = res_bs$st[,,nT], Vb = res_bs$St[,,nT])
  out[[nT]] <- post_nT$B
  dim_array <- dim(post_nT$B)

  # chol decomp Sigma
  RSigma <- array(dim = dim(post_nT$sigma))
  for (j in 1:nsam) {
    RSigma[,,j] <- chol(post_nT$sigma[,,j])
  }

  # then backward t
  for (i in (nT-1):1) {
    if(i %% round(0.1*(nT-1)) == 0){
      print(paste("Sampling:", i, "/", nT))
      print(Sys.time())
    }
    # Check if G_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(G_ls)){
      Gt1 <- as.matrix(G_ls[[i+1]])
    }

    out[[i]] <- array(dim = dim_array)
    para_nt <- res_ff[[i]]
    para_nt$Dt <- check_pds(para_nt$Dt)
    st <- res_bs$st[,,i]
    St <- check_pds(res_bs$St[,,i])
    RSt <- chol(St)

    # get nsam samples for different Sigma
    for (j in 1:nsam) {
      out[[i]][,,j] <- rmn_chol_cpp(m = st, RM = RSt, RSigma = RSigma[,,j])
    }
    if(i == 1){
      names(out) <- paste0("T", 1:nT)
    }
  }

  out$Sigma <- post_nT$sigma
  return(out)
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
                                     nT, phi_para, gp_sigma2 = 1.1, gp_tau2 = 10^(-4),
                            delta = 1.0){
  out <- list()
  post_mean_ls <- list()
  N <- dim(input)[1]

  # deal with one input prediction
  if(is.vector(input_new)){
    input_new <- t(input_new)
  }

  # get spatial kernel
  input_full <- rbind(input, input_new)
  dist_input_full <- as.matrix(stats::dist(input_full, method = "euclidean", diag = T, upper = T))
  V_full <- gen_exp_kernel(loc = input_full, phi = phi_para, sigma2 = gp_sigma2, tau2 = gp_tau2) # exponential kernel
  Vt <- V_full[1:N, 1:N]
  Jt <- V_full[1:N, (N+1):dim(V_full)[2]]
  # prepare matrix
  Vt_chol <- chol(Vt)
  Vtinv <- chol2inv(Vt_chol)
  tJtVtinv <- crossprod(Jt, Vtinv)

  # initialize Ft, Ft_new
  Ft <- F_ls
  Ft_new <- F_new_ls

  for (i in 1:nT) {
    # get Ft's
    if(is.list(F_ls) && length(F_ls) != 1){
      Ft <- F_ls[[i]]
    }
    if(is.list(F_new_ls) && length(F_new_ls) != 1){
      Ft_new <- F_new_ls[[i]]
    }
    Yt <- Y[[i]]
    post_mean <- array(dim = c(dim(input_new)[1], dim(Yt)[2], nsam))

    for(j in 1:nsam){
      # Get Thetat
      Thetat <- res_ffbs[[i]][,,j]
      post_mean[,,j] <- Ft_new %*% Thetat + tJtVtinv %*% (Yt - Ft %*% Thetat)

    }
    post_mean_ls[[i]] <- post_mean
  }
  names(post_mean_ls) <- paste("T", seq(1:nT), sep = "")
  out <- post_mean_ls
  return(out)
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
#' @export
FFBS_predict_exact <- function(Y, para_ffbs, F_ls, F_new_ls, input, input_new,
                               nT, phi_para, gp_sigma2 = 1.1, gp_tau2 = 10^(-4),
                               delta = 1.0){
  out <- list()
  post_mean_ls <- list()
  N <- dim(input)[1]

  # deal with one input prediction
  if(is.vector(input_new)){
    input_new <- t(input_new)
  }

  # get spatial kernel
  input_full <- rbind(input, input_new)
  dist_input_full <- as.matrix(stats::dist(input_full, method = "euclidean", diag = T, upper = T))
  V_full <- gen_exp_kernel(loc = input_full, phi = phi_para, sigma2 = gp_sigma2, tau2 = gp_tau2) # exponential kernel
  Vt <- V_full[1:N, 1:N]
  Jt <- V_full[1:N, (N+1):dim(V_full)[2]]
  # prepare matrix
  Vt_chol <- chol(Vt)
  Vtinv <- chol2inv(Vt_chol)
  tJtVtinv <- crossprod(Jt, Vtinv)

  # initialize Ft, Ft_new
  Ft <- F_ls
  Ft_new <- F_new_ls

  for (i in 1:nT) {
    # get Ft's
    if(is.list(F_ls) && length(F_ls) != 1){
      Ft <- F_ls[[i]]
    }
    if(is.list(F_new_ls) && length(F_new_ls) != 1){
      Ft_new <- F_new_ls[[i]]
    }
    Yt <- Y[[i]]

    # Get st
    st <- para_ffbs$bs$st[,,i]
    post_mean <- Ft_new %*% st + tJtVtinv %*% (Yt - Ft %*% st)

    post_mean_ls[[i]] <- post_mean
  }
  names(post_mean_ls) <- paste("T", seq(1:nT), sep = "")
  out <- post_mean_ls
  return(out)
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
#'
#' @return A named list with \code{ff} (output of \code{\link{FF_sigma2R}}) and
#'   \code{bs} (output of \code{\link{BS}}).
#'
#' @seealso \code{\link{FF_sigma2R}}, \code{\link{BS}},
#'   \code{\link{FFBS_sampling_sigma2R}}
#' @export
FFBS_sigma2R <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, R, delta = 1.0){
  # TODO FF_cpp, BS_cpp
  out <- list()
  res_ff <- FF_sigma2R(Y = Y, F_ls = F_ls, G_ls = G_ls,
                        W_ls = W_ls, V_ls = V_ls,
                        m0 = m0, M0 = M0,
                        n0 = n0, D0 = D0,
                        nT = nT, R = R, delta = delta)
  res_bs <- BS(res_ff, G_ls, nT = nT, delta = delta)
  print("Finish FFBS")
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
#'
#' @return A named list of length \code{nT + 1}. Elements \code{"T1"} through
#'   \code{"T<nT>"} are arrays of dimension \code{c(p, q, nsam)} with state
#'   samples. The element \code{Sigma} is a vector of length \code{nsam} with
#'   samples of \eqn{\sigma^2}.
#'
#' @seealso \code{\link{FFBS_sigma2R}}, \code{\link{MNIG_sampler}}
#' @export
FFBS_sampling_sigma2R <- function(nsam, para_ffbs, F_ls, G_ls, nT, R, delta = 1){
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

  para_nt <- res_ff[[nT]]
  res_bs$St[,,nT] <- check_pds(res_bs$St[,,nT])
  post_nT <- MNIG_sampler(nsam = nsam, X = Ft, v = para_nt$nt, S = para_nt$Dt,
                          C = res_bs$st[,,nT], Vb = res_bs$St[,,nT], R = R)
  out[[nT]] <- post_nT$B
  dim_array <- dim(post_nT$B)

  # then backward t
  for (i in (nT-1):1) {
    if(i %% round(0.1*(nT-1)) == 0){
      print(paste("Sampling:", i, "/", nT))
      print(Sys.time())
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
      Sigma_j <- post_nT$sigma[j]
      out[[i]][,,j] <- rmn_cpp(m = st, U = St, V = Sigma_j * R)
    }
    if(i == 1){
      names(out) <- paste0("T", 1:nT)
    }
  }

  out$Sigma <- post_nT$sigma
  return(out)
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
#'
#' @return A named list with \code{ff} (output of \code{\link{FF_I}}) and
#'   \code{bs} (output of \code{\link{BS}}).
#'
#' @seealso \code{\link{FF_I}}, \code{\link{BS}}, \code{\link{FFBS_sampling_I}}
#' @export
FFBS_I <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, nT, delta = 1.0){
  # TODO FF_cpp, BS_cpp
  out <- list()
  res_ff <- FF_I(Y = Y, F_ls = F_ls, G_ls = G_ls,
                                W_ls = W_ls, V_ls = V_ls,
                                m0 = m0, M0 = M0,
                                nT = nT, delta = delta)
  res_bs <- BS(res_ff, G_ls, nT = nT, delta = delta)
  print("Finish FFBS")
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
#'
#' @return A named list of length \code{nT}. Each element \code{"T<t>"} is an
#'   array of dimension \code{c(p, q, nsam)} with posterior samples of the
#'   state at time \eqn{t}.
#'
#' @seealso \code{\link{FFBS_I}}
#' @export
FFBS_sampling_I <- function(nsam, para_ffbs, F_ls, G_ls, nT, delta = 1){
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
    if(i %% round(0.1*(nT-1)) == 0){
      print(paste("Sampling:", i, "/", nT))
      print(Sys.time())
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
