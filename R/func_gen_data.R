#' Generate a random matrix with entries scaled to \eqn{[-1, 1]}
#'
#' Draws a matrix of independent standard normal entries and divides by the
#' absolute maximum, so all entries lie in \eqn{[-1, 1]}.
#'
#' @param nrow Integer. Number of rows.
#' @param ncol Integer. Number of columns.
#'
#' @return A numeric matrix of dimension \code{c(nrow, ncol)}.
#'
#' @export
gen_ran_matrix <- function(nrow, ncol){
  S <- matrix(rnorm(nrow * ncol), nrow  = nrow, ncol  = ncol)
  S <- S / abs(max(S))
  return(S)
}

#' Generate a random positive definite matrix
#'
#' Constructs a positive definite matrix of dimension \code{dim} by forming
#' \eqn{A^\top A} from a random normal matrix \eqn{A} and then rescaling so the
#' maximum entry is 1.
#'
#' @param dim Integer. Number of rows (and columns) of the output matrix.
#'
#' @return A symmetric positive definite numeric matrix of dimension
#'   \code{c(dim, dim)}.
#'
#' @export
gen_pd_matrix <- function(dim){
  tempS <- matrix(rnorm(dim * dim), nrow  = dim, ncol  = dim)
  S <- crossprod(tempS)
  S <- S / abs(max(S))
  return(S)
}

#' Draw one sample from a matrix-normal distribution (Cholesky parameterisation)
#'
#' Samples \eqn{\Theta \sim MN(m, U, \Sigma)} using the upper-triangular
#' Cholesky factors \code{RM} and \code{RSigma} of \eqn{U} and \eqn{\Sigma}
#' respectively: \eqn{\Theta = m + \text{RM}^\top Z \, \text{RSigma}}, where
#' \eqn{Z} has i.i.d. standard normal entries.
#'
#' @param m Numeric matrix. Mean matrix (\eqn{p \times S}).
#' @param RM Numeric matrix. Upper-triangular Cholesky factor of the row
#'   covariance \eqn{U} (\eqn{p \times p}).
#' @param RSigma Numeric matrix. Upper-triangular Cholesky factor of the column
#'   covariance \eqn{\Sigma} (\eqn{S \times S}).
#'
#' @return A numeric matrix of dimension \code{c(p, S)}.
#'
#' @seealso \code{\link{rmn_chol_more}}
#' @export
rmn_chol <- function(m, RM, RSigma){
  p <- nrow(m)
  S <- ncol(m)
  vec <- rnorm(p * S)
  Z <- matrix(vec, nrow = p, ncol = S)
  theta <- m + crossprod(RM, Z) %*% RSigma
  return(theta)
}


#' Draw multiple samples from a matrix-normal distribution (Cholesky parameterisation)
#'
#' Draws \code{nsam} independent samples from \eqn{MN(m, U, \Sigma)} using
#' the Cholesky factors of \eqn{U} and \eqn{\Sigma}.
#'
#' @param nsam Integer. Number of samples to draw.
#' @param m Numeric matrix. Mean matrix (\eqn{p \times S}).
#' @param RM Numeric matrix. Upper-triangular Cholesky factor of the row
#'   covariance \eqn{U} (\eqn{p \times p}).
#' @param RSigma Numeric matrix. Upper-triangular Cholesky factor of the column
#'   covariance \eqn{\Sigma} (\eqn{S \times S}).
#'
#' @return An array of dimension \code{c(p, S, nsam)}.
#'
#' @seealso \code{\link{rmn_chol}}
#' @export
rmn_chol_more <- function(nsam, m, RM, RSigma){
  p <- nrow(m)
  S <- ncol(m)
  out <- array(dim = c(p, S, nsam))
  for (i in 1:nsam) {
    vec <- rnorm(p * S)
    Z <- matrix(vec, nrow = p, ncol = S)
    theta <- m + crossprod(RM, Z) %*% RSigma
    out[,,i] <- theta
  }

  return(out)
}


#' Generate synthetic FFBS data in memory
#'
#' Simulates data from the MNIW dynamic linear model for \code{nT} time steps
#' and returns the generated observations and model parameters as in-memory
#' objects.
#'
#' @param N Integer. Number of spatial locations (rows of \eqn{Y_t}).
#' @param S Integer. Number of response variables (columns of \eqn{Y_t}).
#' @param p Integer. Dimension of the state vector (must be \eqn{\geq 2}).
#' @param nT Integer. Number of time steps.
#'
#' @return A named list with:
#'   \describe{
#'     \item{Y}{Named list of length \code{nT} (\code{"Y1"}, \ldots,
#'       \code{"Y<nT>"}), each element an \eqn{N \times S} data matrix.}
#'     \item{para}{Named list of model parameters: \code{loc}, \code{dist},
#'       \code{m0}, \code{M0}, \code{G0}, \code{F0}, \code{Sigma},
#'       \code{RSigma}, \code{V0}, \code{RV0}, \code{W0}, \code{n0},
#'       \code{D0}.}
#'   }
#'
#' @seealso \code{\link{gen_ffbs_csv}}
#' @export
gen_ffbs_data <- function(N, S, p, nT){
  out <- list()
  para <- list()
  Y_out <- list()
  m0 <- gen_ran_matrix(p, S)
  M0 <- gen_pd_matrix(p)
  n0 <- S + 2
  D0 <- gen_pd_matrix(S)
  Sigmapre <- rinvwishart(1, nu = n0, Omega = D0, epsilon = 0, checkSymmetry = TRUE)
  Sigma <- matrix(Sigmapre, ncol = S, nrow = S)
  G0 <- diag(1, nrow = p)
  # generate F0
  F0 <- matrix(nrow = N, ncol = p)
  if(p > 1){
    F0[,1] <- 1
    for (i in 2:p) {
      F0[,i] <- rnorm(n = N)
    }
  } else{
    stop("p is less than 2")
  }
  # generate V0 as Gaussian kernel
  loc <- F0 # location
  dist <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  phi <- 0.4 * 3 / max(dist)
  sigma2 <- 1
  tau2 <- 0.01
  V0 <- gen_gp_kernel(loc = F0, phi = phi, sigma2 = sigma2, tau2 = tau2)
  W0 <- gen_pd_matrix(p)
  RM0 <- chol(M0)
  RSigma <- chol(Sigma)
  RV0 <- chol(V0)
  zero_NS <- matrix(0, nrow = N, ncol = S)
  zero_pS <- matrix(0, nrow = p, ncol = S)

  # generate and save Y's
  for(i in 1:nT){
    if(i == 1){
      Theta = rmn_chol(m0, RM0, RSigma)
      E <- rmn_chol(zero_NS, RV0, RSigma)
      Y <- F0 %*% Theta + E
    }
    else{
      Gamma <- rmn_chol(zero_pS, W0, RSigma)
      Theta <- G0 %*% Theta + Gamma
      E <- rmn_chol(zero_NS, RV0, RSigma)
      Y <- F0 %*% Theta + E
    }
    Yl <- list(Y)
    names(Yl) <- paste("Y", i, sep = "")
    Y_out <- append(Y_out, Yl)
  }


  # save other matrices
  para$loc <- loc
  para$dist <- dist
  para$m0 <- m0
  para$M0 <- M0
  para$G0 <- G0
  para$F0 <- F0
  para$Sigma <- Sigma
  para$RSigma <- RSigma
  para$V0 <- V0
  para$RV0 <- RV0
  para$W0 <- W0
  para$n0 <- n0
  para$D0 <- D0

  out <- list(Y = Y_out, para = para)
  return(out)
}


#' Generate synthetic FFBS data and write to CSV files
#'
#' Simulates data from the MNIW dynamic linear model for \code{nT} time steps
#' and writes the observations \eqn{Y_1, \ldots, Y_T} and all model parameters
#' to CSV files in \code{path}.
#'
#' @param N Integer. Number of spatial locations (rows of \eqn{Y_t}).
#' @param S Integer. Number of response variables (columns of \eqn{Y_t}).
#' @param p Integer. Dimension of the state vector (must be \eqn{\geq 2}).
#' @param nT Integer. Number of time steps.
#' @param path Character. Directory path (with trailing \code{/}) where CSV
#'   files are written.
#'
#' @return Invisibly returns \code{NULL}. Files written: \code{m0.csv},
#'   \code{MM0.csv}, \code{G0.csv}, \code{F0.csv}, \code{Sigma.csv},
#'   \code{RSigma.csv}, \code{V0.csv}, \code{RV0.csv}, \code{W0.csv}, and
#'   \code{Y1.csv} through \code{Y<nT>.csv}.
#'
#' @seealso \code{\link{gen_ffbs_data}}
#' @export
gen_ffbs_csv <- function(N, S, p, nT, path){
  m0 <- gen_ran_matrix(p, S)
  M0 <- gen_pd_matrix(p)
  Sigma <- gen_pd_matrix(S)
  G0 <- diag(1, nrow = p)
  # generate F0
  F0 <- matrix(nrow = N, ncol = p)
  if(p > 1){
    F0[,1] <- 1
    for (i in 2:p) {
      F0[,i] <- rnorm(n = N)
    }
  } else{
    stop("p is less than 2")
  }
  # generate V0 as Gaussian kernel
  loc <- F0 # location
  dist <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  phi <- 0.4 * 3 / max(dist)
  sigma2 <- 1
  tau2 <- 0.01
  V0 <- gen_gp_kernel(loc = F0, phi = phi, sigma2 = sigma2, tau2 = tau2)
  W0 <- gen_pd_matrix(p)
  RM0 <- chol(M0)
  RSigma <- chol(Sigma)
  RV0 <- chol(V0)
  zero_NS <- matrix(0, nrow = N, ncol = S)
  zero_pS <- matrix(0, nrow = p, ncol = S)

  # save big matrices
  write.csv(x = m0, file = paste0(path, "m0.csv"), row.names=FALSE)
  write.csv(x = M0, file = paste0(path, "MM0.csv"), row.names=FALSE)
  write.csv(x = G0, file = paste0(path, "G0.csv"), row.names=FALSE)
  write.csv(x = F0, file = paste0(path, "F0.csv"), row.names=FALSE)
  write.csv(x = Sigma, file = paste0(path, "Sigma.csv"), row.names=FALSE)
  write.csv(x = RSigma, file = paste0(path, "RSigma.csv"), row.names=FALSE)
  write.csv(x = V0, file = paste0(path, "V0.csv"), row.names=FALSE)
  write.csv(x = RV0, file = paste0(path, "RV0.csv"), row.names=FALSE)
  write.csv(x = W0, file = paste0(path, "W0.csv"), row.names=FALSE)

  # generate and save Y's
  for(i in 1:nT){
    if(i == 1){
      Theta = rmn_chol(m0, RM0, RSigma)
      E <- rmn_chol(zero_NS, RV0, RSigma)
      Y <- F0 %*% Theta + E
    }
    else{
      Gamma <- rmn_chol(zero_pS, W0, RSigma)
      Theta <- G0 %*% Theta + Gamma
      E <- rmn_chol(zero_NS, RV0, RSigma)
      Y <- F0 %*% Theta + E
    }
    write.csv(x = Y, file = paste0(path, "Y", i, ".csv"), row.names=FALSE)
  }

}


#' Compute a Gaussian Process covariance kernel matrix
#'
#' Evaluates the isotropic GP exponential kernel
#' \eqn{k(s, s') = \sigma^2 \exp(-\phi \|s - s'\|) + \tau^2 \mathbf{1}_{s=s'}}
#' for all pairs of rows in \code{loc}.
#'
#' @param loc Numeric matrix. Spatial locations (\eqn{n \times d}).
#' @param phi Numeric scalar. Range (decay) parameter of the exponential kernel.
#' @param sigma2 Numeric scalar. Marginal variance parameter.
#' @param tau2 Numeric scalar. Nugget (white-noise) variance.
#'
#' @return A symmetric positive definite matrix of dimension \code{c(n, n)}.
#'
#' @seealso \code{\link{gen_exp_kernel}}, \code{\link{gen_expsq_kernel}}
#' @export
gen_gp_kernel <- function(loc, phi, sigma2, tau2){
  dist_k <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  n <- dim(dist_k)[1]
  var_gp <- sigma2 * exp(-phi * dist_k) + tau2 * (dist_k == 0)
  return(var_gp)
}

#' Compute a squared-exponential (Gaussian) GP kernel matrix
#'
#' Evaluates the squared-exponential kernel
#' \eqn{k(s, s') = \sigma^2 \exp(-\phi \|s - s'\|^2) + \tau^2 \mathbf{1}_{s=s'}}
#' for all pairs of rows in \code{loc}.
#'
#' @param loc Numeric matrix. Spatial locations (\eqn{n \times d}).
#' @param phi Numeric scalar. Inverse length-scale parameter.
#' @param sigma2 Numeric scalar. Marginal variance. Defaults to \code{1}.
#' @param tau2 Numeric scalar. Nugget variance. Defaults to \code{0}.
#'
#' @return A symmetric positive definite matrix of dimension \code{c(n, n)}.
#'
#' @seealso \code{\link{gen_gp_kernel}}, \code{\link{gen_exp_kernel}}
#' @export
gen_expsq_kernel <- function(loc, phi, sigma2 = 1, tau2 = 0){
  dist_k <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  n <- dim(dist_k)[1]
  var_gp <- sigma2 * exp(-phi * dist_k^2) + tau2 * (dist_k == 0)
  return(var_gp)
}

#' Compute an exponential GP kernel matrix
#'
#' Evaluates the exponential kernel
#' \eqn{k(s, s') = \sigma^2 \exp(-\phi \|s - s'\|) + \tau^2 \mathbf{1}_{s=s'}}
#' for all pairs of rows in \code{loc}.
#'
#' @param loc Numeric matrix. Spatial locations (\eqn{n \times d}).
#' @param phi Numeric scalar. Range (decay) parameter.
#' @param sigma2 Numeric scalar. Marginal variance. Defaults to \code{1}.
#' @param tau2 Numeric scalar. Nugget variance. Defaults to \code{0}.
#'
#' @return A symmetric positive definite matrix of dimension \code{c(n, n)}.
#'
#' @seealso \code{\link{gen_gp_kernel}}, \code{\link{gen_expsq_kernel}}
#' @export
gen_exp_kernel <- function(loc, phi, sigma2 = 1, tau2 = 0){
  dist_k <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  n <- dim(dist_k)[1]
  var_gp <- sigma2 * exp(-phi * dist_k) + tau2 * (dist_k == 0)
  return(var_gp)
}


#' Build AR(1) covariate list from a list of response matrices
#'
#' Constructs the time-varying covariate list \eqn{F_t} for a first-order
#' autoregressive model: \eqn{F_t = Y_{t-1}} (with \eqn{F_1 = Y_1}).
#'
#' @param Y List of length \code{nT}. Each element is the \eqn{N \times S}
#'   data matrix at time \eqn{t}.
#' @param nT Integer. Number of time steps.
#'
#' @return A list of length \code{nT} where element \eqn{t} is the covariate
#'   matrix \eqn{F_t}.
#'
#' @seealso \code{\link{gen_F_ls_AR2}}, \code{\link{gen_F_ls_AR1_EP}}
#' @export
gen_F_ls_AR1 <- function(Y, nT){
  F_ls <- list()
  for (i in 1:nT) {
    if(i == 1){
      Ft <- Y[[1]]
      F_ls[[i]] <- Ft
    } else{
      Ft <- Y[[i - 1]]
      F_ls[[i]] <- Ft
    }
  }
  return(F_ls)
}

#' Build AR(2) covariate list from a list of response matrices
#'
#' Constructs the time-varying covariate list \eqn{F_t} for a second-order
#' autoregressive model by column-binding the two most recent lags:
#' \eqn{F_t = [Y_{t-2}, Y_{t-1}]} (with special handling for \eqn{t = 1, 2}).
#'
#' @param Y List of length \code{nT}. Each element is the \eqn{N \times S}
#'   data matrix at time \eqn{t}.
#' @param nT Integer. Number of time steps.
#'
#' @return A list of length \code{nT} where element \eqn{t} is the
#'   \eqn{N \times 2S} covariate matrix \eqn{F_t}.
#'
#' @seealso \code{\link{gen_F_ls_AR1}}, \code{\link{gen_F_ls_AR2_EP}}
#' @export
gen_F_ls_AR2 <- function(Y, nT){
  F_ls <- list()
  for (i in 1:nT) {
    if(i == 1){
      Ft <- cbind(Y[[1]],
                  0.01 * gen_ran_matrix(nrow = nrow(Y[[1]]), ncol = ncol(Y[[1]])))
      F_ls[[i]] <- Ft
    } else if(i == 2){
      Ft <- cbind(Y[[1]],
                  1/2 * (Y[[1]] + Y[[3]]))
      F_ls[[i]] <- Ft
    } else{
      Ft <- cbind(Y[[i - 2]], Y[[i - 1]])
      F_ls[[i]] <- Ft
    }
  }
  return(F_ls)
}

#' Build AR(1) covariate list for the episode-block model
#'
#' Generates a flattened covariate list for the episode-partition (EP) model
#' with AR(1) lags. The output list is indexed by \eqn{(t-1) \times n_b + j},
#' where \eqn{j} is the block index and \eqn{t} is the time index.
#'
#' @param Y List of length \code{nT}. Each element is the full spatial data
#'   matrix at time \eqn{t}.
#' @param nT Integer. Number of time steps.
#' @param n_b Integer. Number of spatial blocks.
#' @param ind Data frame of block indices as returned by
#'   \code{\link{generate_grid}}.
#'
#' @return A list of length \code{nT * n_b} with block-wise AR(1) covariate
#'   matrices.
#'
#' @seealso \code{\link{gen_F_ls_AR1}}, \code{\link{gen_F_ls_AR2_EP}}
#' @export
gen_F_ls_AR1_EP <- function(Y, nT, n_b, ind){
  out <- list()
  F_ls_full <- gen_F_ls_AR1(Y = Y, nT = nT)
  for (i in 1:nT) {
    temp_F <- F_ls_full[[i]]
    for (j in 1:n_b) {
      out[[(i-1)*n_b + j]] <- temp_F[ind[j,2]:ind[j,3], ind[j,4]:ind[j,5]]
    }
  }
  return(out)
}

#' Build AR(2) covariate list for the episode-block model
#'
#' Generates a flattened covariate list for the episode-partition (EP) model
#' with AR(2) lags. Both lag-1 and lag-2 column blocks are extracted and
#' column-bound for each spatial block.
#'
#' @param Y List of length \code{nT}. Each element is the full spatial data
#'   matrix at time \eqn{t}.
#' @param nT Integer. Number of time steps.
#' @param n_b Integer. Number of spatial blocks.
#' @param ind Data frame of block indices as returned by
#'   \code{\link{generate_grid}}.
#'
#' @return A list of length \code{nT * n_b} with block-wise AR(2) covariate
#'   matrices (each of width \eqn{2 \times \text{bncol}}).
#'
#' @seealso \code{\link{gen_F_ls_AR2}}, \code{\link{gen_F_ls_AR1_EP}}
#' @export
gen_F_ls_AR2_EP <- function(Y, nT, n_b, ind){
  out <- list()
  S <- ncol(Y[[1]])
  F_ls_full <- gen_F_ls_AR2(Y = Y, nT = nT)
  for (i in 1:nT) {
    temp_F <- F_ls_full[[i]]
    for (j in 1:n_b) {
      out[[(i-1)*n_b + j]] <- temp_F[ind[j,2]:ind[j,3],
                                     c(ind[j,4]:ind[j,5], (S + (ind[j,4]:ind[j,5])))]
    }
  }
  return(out)
}


#' Recover episode-partitioned data to original time dimension (exact)
#'
#' Reassembles a flattened episode-partitioned list (as produced by the EP
#' model) back into a list indexed by the original \code{nT_ori} time steps.
#' Assumes the number of blocks per season \code{n_block = nT / nT_ori} is an
#' integer.
#'
#' @param dat_EP List of length \code{nT}. Flattened episode-partitioned data.
#' @param nT_ori Integer. Number of original (non-partitioned) time steps.
#' @param nT Integer. Total number of episode time steps
#'   (\code{nT = nT_ori * n_block}).
#'
#' @return A named list of length \code{nT_ori} (\code{"T1"}, \ldots,
#'   \code{"T<nT_ori>"}), each element a matrix covering one original time step.
#'
#' @seealso \code{\link{recover_from_EP_MC}}
#' @export
recover_from_EP_exact <- function(dat_EP, nT_ori, nT){
  if(length(dat_EP) != nT){
    stop("dimension of dat_EP doesn't match nT")
  }

  n_block <- nT / nT_ori
  out <- list()

  # cbind n_block of EP blocks to recover season
  for (i in 1:nT_ori) {
    out[[i]] <- as.matrix(as.data.frame(dat_EP[(1+(i-1)*n_block):(i*n_block)]))
  }

  names(out) <- paste0("T", seq(1:nT_ori))
  return(out)
}

#' Recover episode-partitioned posterior samples to original time dimension
#'
#' Reassembles a flattened episode-partitioned list of posterior sample arrays
#' back into a list indexed by the original \code{nT_ori} time steps.
#'
#' @param dat_EP List of length \code{nT}. Each element is an array of
#'   dimension \code{c(p, bncol, nsam)}.
#' @param nT_ori Integer. Number of original time steps.
#' @param nT Integer. Total number of episode time steps.
#' @param nsam Integer. Number of posterior samples.
#'
#' @return A named list of length \code{nT_ori}. Each element is an array of
#'   dimension \code{c(p, bncol * n_block, nsam)}.
#'
#' @seealso \code{\link{recover_from_EP_exact}}
#' @export
recover_from_EP_MC <- function(dat_EP, nT_ori, nT, nsam){
  if(length(dat_EP) != nT){
    stop("dimension of dat_EP doesn't match nT")
  }

  n_block <- nT / nT_ori
  dim_block <- dim(dat_EP[[1]])
  dim_no_block <- dim_block
  dim_no_block[2] <- dim_block[2] * n_block
  n_input <- dim(dat_EP)
  out <- list()

  # transform list to array to simplify modification
  dat_ary <- array(dim = c(dim_block, nT))
  for (i in 1:nT) {
    dat_ary[,,,i] <- dat_EP[[i]]
  }

  # cbind n_block of EP blocks to recover season
  for (i in 1:nT_ori) {
    out[[i]] <- array(dim = dim_no_block)
    for (j in 1:nsam) {
      out[[i]][,,j] <- as.matrix(as.data.frame(dat_ary[,,j,(1+(i-1)*n_block):(i*n_block)]))
    }
  }

  names(out) <- paste0("T", seq(1:nT_ori))
  return(out)
}
