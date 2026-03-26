# MNIW----
#' Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#'
#' @param Y The data matrix.
#' @param F_ls The matrix of covariates F_t.
#' @param G_ls G_t, the beta transition matrix.
#' @param m0 mean of beta_{t-1} | D_{t-1}
#' @param M0 leF_ls-covariance matrix of beta_{t-1} | D_{t-1}
#' @param W_ls leF_ls-covariance matrix of the noise parameter of beta_{t-1} | D_{t-1}
#' @param V_ls leF_ls-covariance matrix of the noise parameter of Y
#' @param n0 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | D_{t-1}
#' @param D0 the scale matrix of the right-covariance matrix Sigma | D_{t-1}
#' @param delta The right-variance matrix discount factor
#' @returns The mean and covariance matrices m_t and C_t of one filtering step, the updated inverse-Wishart parameters a_t and B_t for the right-covariance matrix, plus other relevant parameters.
#' @export
FFBS <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, delta = 1.0){
  # TODO FF_cpp, BS_cpp
  out <- list()
  res_ff <- FF(Y = Y, F_ls = F_ls, G_ls = G_ls,
              W_ls = W_ls, V_ls = V_ls,
              m0 = m0, M0 = M0,
              n0 = n0, D0 = D0,
              nT = nT, delta = delta)
  # print("Finish FF")
  res_bs <- BS(res_ff, G_ls, nT = nT, delta = delta)
  print("Finish FFBS")
  out$ff <- res_ff
  out$bs <- res_bs
  return(out)
}


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



#' prediction using monte carlo method to estimate post mean:st, to get mean of Y_new
#'
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
      # post_mean <- Ft_new %*% Thetat + t(Jt) %*% solve(Vt) %*% (Yt - Ft %*% Thetat)
      post_mean[,,j] <- Ft_new %*% Thetat + tJtVtinv %*% (Yt - Ft %*% Thetat)

    }
    post_mean_ls[[i]] <- post_mean
  }
  names(post_mean_ls) <- paste("T", seq(1:nT), sep = "")
  out <- post_mean_ls
  return(out)
}

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
#' Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#'
#' @param Y The data matrix.
#' @param F_ls The matrix of covariates F_t.
#' @param G_ls G_t, the beta transition matrix.
#' @param m0 mean of beta_{t-1} | D_{t-1}
#' @param M0 leF_ls-covariance matrix of beta_{t-1} | D_{t-1}
#' @param W_ls leF_ls-covariance matrix of the noise parameter of beta_{t-1} | D_{t-1}
#' @param V_ls leF_ls-covariance matrix of the noise parameter of Y
#' @param n0 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | D_{t-1}
#' @param D0 the scale matrix of the right-covariance matrix Sigma | D_{t-1}
#' @param delta The right-variance matrix discount factor
#' @returns The mean and covariance matrices m_t and C_t of one filtering step, the updated inverse-Wishart parameters a_t and B_t for the right-covariance matrix, plus other relevant parameters.
#' @export
FFBS_sigma2R <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, R, delta = 1.0){
  # TODO FF_cpp, BS_cpp
  out <- list()
  res_ff <- FF_sigma2R(Y = Y, F_ls = F_ls, G_ls = G_ls,
                        W_ls = W_ls, V_ls = V_ls,
                        m0 = m0, M0 = M0,
                        n0 = n0, D0 = D0,
                        nT = nT, R = R, delta = delta)
  # print("Finish FF")
  res_bs <- BS(res_ff, G_ls, nT = nT, delta = delta)
  print("Finish FFBS")
  out$ff <- res_ff
  out$bs <- res_bs
  return(out)
}


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
#' @export
FFBS_I <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, nT, delta = 1.0){
  # TODO FF_cpp, BS_cpp
  out <- list()
  res_ff <- FF_I(Y = Y, F_ls = F_ls, G_ls = G_ls,
                                W_ls = W_ls, V_ls = V_ls,
                                m0 = m0, M0 = M0,
                                nT = nT, delta = delta)
  # print("Finish FF")
  res_bs <- BS(res_ff, G_ls, nT = nT, delta = delta)
  print("Finish FFBS")
  out$ff <- res_ff
  out$bs <- res_bs
  return(out)
}

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

