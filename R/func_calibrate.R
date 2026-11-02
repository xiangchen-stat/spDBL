# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Legacy research helpers are deprecated in 1.2.0. Each warns once per session;
# options(spDBL.warn_deprecated = FALSE) silences the warnings.
.spdbl_deprecated <- new.env(parent = emptyenv())
.legacy_deprecated <- function(name) {
  if (isFALSE(getOption("spDBL.warn_deprecated", TRUE)) ||
      exists(name, envir = .spdbl_deprecated, inherits = FALSE)) return(invisible())
  assign(name, TRUE, envir = .spdbl_deprecated)
  warning(structure(class = c("deprecatedWarning", "warning", "condition"),
    list(message = paste0(name, "() is deprecated since spDBL 1.2.0 and will ",
      "be removed in a future release. It comes from the original research ",
      "code; new work should use calibration_prepare() and calibration_fit()."),
      call = NULL)))
}

#' Scale a vector to the unit interval via a uniform transformation
#'
#' Maps \code{x} to \eqn{(x - \text{low}) / (\text{high} - \text{low})}.
#' Raises an error if \code{x} is outside \code{[low, high]}. Sets
#' \code{NaN} results to 0.
#'
#' @param x Numeric vector or matrix. Values to scale.
#' @param low Numeric. Lower bound of the original range.
#' @param high Numeric. Upper bound of the original range.
#'
#' @return Numeric vector or matrix with values in \eqn{[0, 1]}.
#'
#' @seealso \code{\link{scale_back_uniform}}
#' @examples
#' x <- c(2, 3.5, 5)
#' u <- scale_uniform(x, low = 2, high = 5)
#' u
#' all.equal(scale_back_uniform(u, low = 2, high = 5), x)
#'
#' @section Lifecycle:
#' Deprecated since spDBL 1.2.0 and will be removed in a future release.
#' The first call in a session warns; set
#' \code{options(spDBL.warn_deprecated = FALSE)} to silence it.
#' @export
scale_uniform <- function(x, low, high){
  .legacy_deprecated("scale_uniform")
  if (all(x >= low) && all(x <= high)) {
    res <- (x - low) / (high - low)
  } else{
    stop("x is out of range.")
  }
  res[is.nan(res)] <- 0
  return(res)
}

#' Invert a uniform scaling transformation
#'
#' Maps \code{x} from the unit interval back to the original range
#' \eqn{[\text{low}, \text{high}]}: \eqn{x \cdot (\text{high} - \text{low}) + \text{low}}.
#'
#' @param x Numeric vector or matrix. Values in \eqn{[0, 1]}.
#' @param low Numeric. Lower bound of the target range.
#' @param high Numeric. Upper bound of the target range.
#'
#' @return Numeric vector or matrix rescaled to \eqn{[\text{low}, \text{high}]}.
#'
#' @seealso \code{\link{scale_uniform}}
#' @examples
#' u <- c(0, 0.25, 1)
#' x <- scale_back_uniform(u, low = -1, high = 3)
#' x
#' all.equal(scale_uniform(x, low = -1, high = 3), u)
#'
#' @section Lifecycle:
#' Deprecated since spDBL 1.2.0 and will be removed in a future release.
#' The first call in a session warns; set
#' \code{options(spDBL.warn_deprecated = FALSE)} to silence it.
#' @export
scale_back_uniform <- function(x, low, high){
  .legacy_deprecated("scale_back_uniform")
  res <- x * (high - low) + low
  return(res)
}

#' Logistic (expit) function
#'
#' Computes the logistic sigmoid \eqn{1 / (1 + e^{-x})}.
#'
#' @param x Numeric vector or matrix.
#'
#' @return Numeric vector or matrix with values in \eqn{(0, 1)}.
#'
#' @examples
#' expit(c(-2, 0, 2))
#' all.equal(expit(1.5), stats::plogis(1.5))
#'
#' @section Lifecycle:
#' Deprecated since spDBL 1.2.0 and will be removed in a future release.
#' The first call in a session warns; set
#' \code{options(spDBL.warn_deprecated = FALSE)} to silence it.
#' @export
expit <- function(x){
  .legacy_deprecated("expit")
  return((1 + exp(1)^(-x))^(-1))
}

#' Invert a matrix via its Cholesky factorisation
#'
#' Computes \eqn{X^{-1}} using \code{chol} followed by \code{chol2inv} for
#' numerical stability.
#'
#' @param X Numeric symmetric positive definite matrix.
#'
#' @return The inverse of \code{X}.
#'
#' @examples
#' A <- matrix(c(4, 1, 1, 3), 2)
#' inv_chol(A)
#' all.equal(inv_chol(A), solve(A))
#'
#' @export
inv_chol <- function(X){
  RX <- chol(X)
  invX <- chol2inv(RX)
  return(invX)
}

#' Log absolute Jacobian of the logit-uniform transformation
#'
#' Computes the log (or raw) absolute Jacobian of the transformation from the
#' constrained parameter \code{eta} (bounded in
#' \code{[eta_limit_low, eta_limit_high]}) to the unconstrained logit scale.
#' Fixed coordinates (equal lower and upper bounds) are excluded from the
#' product. Versions up to 1.1.1 excluded zero entries of \code{eta}
#' instead, which wrongly dropped a free coordinate equal to zero.
#'
#' @param eta Numeric vector. Constrained parameter values.
#' @param eta_limit_low Numeric vector. Lower bounds for each element of
#'   \code{eta}.
#' @param eta_limit_high Numeric vector. Upper bounds for each element of
#'   \code{eta}.
#' @param log Logical. If \code{TRUE}, returns the log Jacobian. Defaults to
#'   \code{FALSE}.
#'
#' @return Numeric scalar (log absolute Jacobian or absolute Jacobian).
#'
#' @examples
#' eta <- c(0.3, 2.5)
#' low <- c(0, 1)
#' high <- c(1, 4)
#' cal_jacobian_logit_uniform(eta, low, high)
#' cal_jacobian_logit_uniform(eta, low, high, log = TRUE)
#' # A fixed coordinate (low == high) does not contribute
#' all.equal(cal_jacobian_logit_uniform(c(0.3, 2), c(0, 2), c(1, 2)), 1 / (0.3 * 0.7))
#'
#' @section Lifecycle:
#' Deprecated since spDBL 1.2.0 and will be removed in a future release.
#' The first call in a session warns; set
#' \code{options(spDBL.warn_deprecated = FALSE)} to silence it.
#' @export
cal_jacobian_logit_uniform <- function(eta, eta_limit_low, eta_limit_high, log = FALSE){
  .legacy_deprecated("cal_jacobian_logit_uniform")
  free <- eta_limit_high > eta_limit_low
  res <- abs(prod(1/(eta[free] - eta_limit_low[free]) + 1/(eta_limit_high[free] - eta[free])))
  if (log == TRUE) {
    res <- log(res)
  }
  return(res)
}

#' Generate synthetic calibration data with correlated discrepancy
#'
#' Simulates observed calibration data \eqn{z_t} by adding a spatially
#' correlated random discrepancy \eqn{u_t} and an independent noise term to
#' the computer model output \eqn{y_t}. The variance \eqn{\tau^2_t} decays
#' over time via the discount factor \eqn{b}.
#'
#' @param ycal_mat Numeric matrix of dimension \code{c(t_cal, s_cal)}. Computer
#'   model output at calibration locations and time steps.
#' @param para_gen_cal Named list with elements:
#'   \describe{
#'     \item{b}{Numeric scalar. Discount factor in \eqn{(0, 1)}.}
#'     \item{n0}{Numeric scalar. Prior shape of \eqn{\tau^2_0}.}
#'     \item{d0}{Numeric scalar. Prior scale of \eqn{\tau^2_0}.}
#'     \item{u0}{Numeric vector of length \code{s_cal}. Initial discrepancy.}
#'   }
#' @param U_gen Numeric matrix. Spatial covariance matrix for the discrepancy
#'   (\code{s_cal x s_cal}).
#'
#' @return A named list with:
#'   \describe{
#'     \item{u}{Matrix of dimension \code{c(t_cal, s_cal)} of discrepancy
#'       realisations.}
#'     \item{z}{Matrix of dimension \code{c(t_cal, s_cal)} of simulated
#'       observations.}
#'     \item{tau2}{Numeric vector of length \code{t_cal} with sampled
#'       variances.}
#'   }
#'
#' @seealso \code{\link{gen_prior_u_tau2}}
#' @examples
#' set.seed(1)
#' ycal <- matrix(rnorm(3 * 4), nrow = 3, ncol = 4)
#' para <- list(b = 0.9, n0 = 20, d0 = 2, u0 = rep(0, 4))
#' sim <- gen_calibrate_data(ycal, para, U_gen = 0.5 * diag(4))
#' names(sim)
#' dim(sim$z)
#' sim$tau2
#'
#' @export
gen_calibrate_data <- function(ycal_mat, para_gen_cal, U_gen){
  res <- list()
  tau2 <- c()
  t_cal <- dim(ycal_mat)[1]
  s_cal <- dim(ycal_mat)[2]
  u <- matrix(nrow = dim(ycal_mat)[1], ncol = dim(ycal_mat)[2])
  z <- matrix(nrow = dim(ycal_mat)[1], ncol = dim(ycal_mat)[2])
  for (t in 1:t_cal) {
    y_t <- ycal_mat[t,]
    alpha <- para_gen_cal$b^t * para_gen_cal$n0
    beta <- para_gen_cal$b^t * para_gen_cal$d0
    tau2_t <- LaplacesDemon::rinvgamma(n = 1, shape = alpha, scale = beta)
    tau2 <- c(tau2, tau2_t)
    e_u <- mniw::rmNorm(n = 1, mu = rep(0, s_cal), Sigma = tau2_t * U_gen)
    if (t == 1) {
      u_t <- para_gen_cal$u0 + e_u
      u[t,] <- u_t
    } else{
      u_t <- u_t + e_u
      u[t,] <- u_t
    }
    e_z <- rnorm(n = s_cal, mean = 0, sd = sqrt(tau2_t))
    z[t,] <- y_t + u_t + e_z
  }
  res$u <- u
  res$z <- z
  res$tau2 <- tau2
  return(res)
}

#' Generate synthetic calibration data with uncorrelated discrepancy
#'
#' Like \code{\link{gen_calibrate_data}} but uses independent (uncorrelated)
#' discrepancy increments drawn from a univariate normal with variance
#' \eqn{\tau^2_t \cdot U_{11}}.
#'
#' @param ycal_mat Numeric matrix. Computer model output (\code{t_cal x s_cal}).
#' @param para_gen_cal Named list (same structure as in
#'   \code{\link{gen_calibrate_data}}).
#' @param U_gen Numeric matrix. Only the \code{[1,1]} entry is used as the
#'   marginal variance scale.
#'
#' @return A named list with:
#'   \describe{
#'     \item{z}{Matrix of dimension \code{c(t_cal, s_cal)} of simulated
#'       observations.}
#'     \item{tau2}{Numeric vector of length \code{t_cal}.}
#'   }
#'
#' @seealso \code{\link{gen_calibrate_data}}
#' @examples
#' set.seed(1)
#' ycal <- matrix(rnorm(3 * 4), nrow = 3, ncol = 4)
#' para <- list(b = 0.9, n0 = 20, d0 = 2, u0 = rep(0, 4))
#' sim <- gen_calibrate_data_uncorr(ycal, para, U_gen = 0.5 * diag(4))
#' names(sim)
#' dim(sim$z)
#' sim$tau2
#'
#' @export
gen_calibrate_data_uncorr <- function(ycal_mat, para_gen_cal, U_gen){
  res <- list()
  tau2 <- c()
  t_cal <- dim(ycal_mat)[1]
  s_cal <- dim(ycal_mat)[2]
  z <- matrix(nrow = dim(ycal_mat)[1], ncol = dim(ycal_mat)[2])
  for (t in 1:t_cal) {
    y_t <- ycal_mat[t,]
    alpha <- para_gen_cal$b^t * para_gen_cal$n0
    beta <- para_gen_cal$b^t * para_gen_cal$d0
    tau2_t <- LaplacesDemon::rinvgamma(n = 1, shape = alpha, scale = beta)
    tau2 <- c(tau2, tau2_t)
    e_u <- rnorm(n = s_cal, mean = 0, sd = sqrt(tau2_t * U_gen[1,1]))
    if (t == 1) {
      u_t <- para_gen_cal$u0 + e_u
    } else{
      u_t <- u_t + e_u
    }
    e_z <- rnorm(n = s_cal, mean = 0, sd = sqrt(tau2_t))
    z[t,] <- y_t + u_t + e_z
  }
  res$z <- z
  res$tau2 <- tau2
  return(res)
}

#' Sample prior discrepancy trajectory and variance sequence
#'
#' Simulates a prior trajectory \eqn{u_1, \ldots, u_T} of the spatial
#' discrepancy field and the associated inverse-gamma variance sequence
#' \eqn{\tau^2_1, \ldots, \tau^2_T} from the dynamic calibration model.
#'
#' @param n0 Numeric scalar. Prior shape of \eqn{\tau^2_0}.
#' @param d0 Numeric scalar. Prior scale of \eqn{\tau^2_0}.
#' @param b Numeric scalar. Discount factor in \eqn{(0, 1)}.
#' @param m0 Numeric vector. Prior mean of the initial discrepancy \eqn{u_0}.
#' @param M0 Numeric matrix. Prior covariance of \eqn{u_0}.
#' @param U Numeric matrix. Spatial covariance matrix for the discrepancy
#'   increments.
#' @param nT_cal Integer. Number of calibration time steps.
#'
#' @return A named list with:
#'   \describe{
#'     \item{u0}{Numeric vector. Sampled initial discrepancy.}
#'     \item{u}{Matrix of dimension \code{c(nT_cal, s_cal)} of discrepancy
#'       realisations.}
#'     \item{tau2_0}{Numeric scalar. Sampled \eqn{\tau^2_0}.}
#'     \item{tau2}{Numeric vector of length \code{nT_cal}.}
#'   }
#'
#' @seealso \code{\link{gen_calibrate_data}}
#' @examples
#' set.seed(1)
#' prior <- gen_prior_u_tau2(n0 = 20, d0 = 2, b = 0.95, m0 = rep(0, 4),
#'                           M0 = diag(4), U = 0.5 * diag(4), nT_cal = 3)
#' names(prior)
#' dim(prior$u)
#' prior$tau2
#'
#' @export
gen_prior_u_tau2 <- function(n0, d0, b, m0, M0, U, nT_cal){
  res <- list()
  s_cal <- length(m0)
  tau2 <- c()
  u <- array(dim = c(nT_cal, s_cal))
  zero_s <- rep(0, s_cal)
  # t = 0
  tau2_t <- LaplacesDemon::rinvgamma(n = 1, shape = n0, scale = d0)
  tau2_0 <- tau2_t
  u_t <- mniw::rmNorm(n = 1, mu = m0, Sigma = tau2_t * M0)
  u0 <- u_t
  # t = 1:T
  for (t in 1:nT_cal) {
    alpha <- b^t * n0
    beta <- b^t * d0
    tau2_t <- LaplacesDemon::rinvgamma(n = 1, shape = alpha, scale = beta)
    tau2 <- c(tau2, tau2_t)
    e_u <- mniw::rmNorm(n = 1, mu = zero_s, Sigma = tau2_t * U)
    u_t <- u_t + e_u
    u[t,] <- u_t
  }

  res$u0 <- u0
  res$u <- u
  res$tau2_0 <- tau2_0
  res$tau2 <- tau2
  return(res)
}

# gen y | Y_1:T

#' Right-hand side of the spatially extended SIR ODE
#'
#' Defines the right-hand side of a spatial SIR (Susceptible-Infected-Recovered)
#' partial differential equation suitable for use with \code{deSolve::ode.2D}.
#' Diffusion is implemented via \code{ReacTran::tran.2D} with zero-flux
#' boundary conditions.
#'
#' @param t Numeric. Current time (passed by the ODE solver).
#' @param y Numeric vector. Flattened state vector of length \eqn{3 N_x N_y}:
#'   \eqn{S} values, then \eqn{I} values, then \eqn{R} values.
#' @param parms Named list of parameters:
#'   \describe{
#'     \item{Nx}{Integer. Number of grid cells in the x-direction.}
#'     \item{Ny}{Integer. Number of grid cells in the y-direction.}
#'     \item{N}{Numeric. Total population size.}
#'     \item{eta}{Numeric vector of length 5:
#'       \code{eta[1]} transmission rate,
#'       \code{eta[2]} recovery rate,
#'       \code{eta[3]} diffusion coefficient for S,
#'       \code{eta[4]} diffusion coefficient for I,
#'       \code{eta[5]} diffusion coefficient for R.}
#'     \item{Gridx}{1D grid object for x (from \code{ReacTran::setup.grid.1D}).}
#'     \item{Gridy}{1D grid object for y.}
#'   }
#'
#' @return A list with one element: a numeric vector of derivatives of the same
#'   length as \code{y}, as required by \code{deSolve}.
#'
#' @seealso \code{\link{gen_pde}}
#' @examples
#' set.seed(1)
#' parms <- list(eta = c(0.3, 0.1, 0, 0, 0), Nx = 2, Ny = 2, N = 100,
#'               Gridx = ReacTran::setup.grid.1D(x.up = 0, x.down = 2, N = 2),
#'               Gridy = ReacTran::setup.grid.1D(x.up = 0, x.down = 2, N = 2))
#' y <- c(runif(4, 50, 100), runif(4, 0, 10), rep(0, 4))
#' d <- SIR(t = 0, y = y, parms = parms)[[1]]
#' length(d)
#' all.equal(sum(d), 0)  # without diffusion, S + I + R is conserved
#'
#' @export
SIR <- function(t, y, parms) {
  Nx <- parms$Nx
  Ny <- parms$Ny
  N <- parms$N
  eta1 <- parms$eta[1]
  eta2 <- parms$eta[2]
  alpha1 <- parms$eta[3]
  alpha2 <- parms$eta[4]
  alpha3 <- parms$eta[5]
  Gridx <- parms$Gridx
  Gridy <- parms$Gridy

  S <- matrix(nrow = Nx, ncol = Ny, data = y[1:(Nx*Ny)])
  I <- matrix(nrow = Nx, ncol = Ny, data = y[(Nx*Ny+1) : (2*Nx*Ny)])
  R <- matrix(nrow = Nx, ncol = Ny, data = y[(2*Nx*Ny+1) : (3*Nx*Ny)])
  #NOTE: tran.2D as called here is equivalent to computing the Laplacian of S, times alpha1;
  #	with dC being the result.
  dS <- -eta1*S*I/N +
    ReacTran::tran.2D(C = S, D.x = alpha1, D.y = alpha1,
                      dx = Gridx, dy= Gridy,
                      C.x.up = 0,
                      C.y.up = 0,
                      C.x.down = 0,
                      C.y.down = 0)$dC
  dI <- eta1*S*I/N - eta2*I +
    ReacTran::tran.2D(C = I, D.x = alpha2, D.y = alpha2,
                      dx = Gridx, dy = Gridy,
                      C.x.up = 0,
                      C.y.up = 0,
                      C.x.down = 0,
                      C.y.down = 0)$dC
  dR <- eta2*I +
    ReacTran::tran.2D(C = R, D.x = alpha3, D.y = alpha3,
                      dx = Gridx, dy = Gridy,
                      C.x.up = 0,
                      C.y.up = 0,
                      C.x.down = 0,
                      C.y.down = 0)$dC
  list(c(dS, dI, dR))
}

#' Simulate a spatially extended SIR PDE model
#'
#' Numerically integrates the spatial SIR ODE (\code{\link{SIR}}) on a
#' two-dimensional grid using \code{deSolve::ode.2D} and returns the log(1 + I)
#' values of the infected compartment.
#'
#' @param eta Numeric vector of length 5. PDE parameters:
#'   transmission rate, recovery rate, and diffusion coefficients for S, I, R.
#' @param Nx Integer. Number of grid cells in the x-direction (rows).
#' @param Ny Integer. Number of grid cells in the y-direction (columns).
#' @param N Numeric. Total population size.
#' @param nT_ori Integer. Number of time steps to simulate (times
#'   \eqn{0, 1, \ldots, \text{nT\_ori} - 1}).
#'
#' @return Numeric matrix of dimension \code{c(nT_ori, Nx * Ny)} containing
#'   \eqn{\log(1 + I_t)} for each time step and grid cell.
#'
#' @seealso \code{\link{SIR}}, \code{\link{update_y_eta}}
#' @examples
#' y <- gen_pde(eta = c(0.5, 0.1, 0.05, 0.05, 0.05), Nx = 4, Ny = 4, N = 100,
#'              nT_ori = 3)
#' dim(y)  # nT_ori x (Nx * Ny)
#' round(rowSums(y), 3)  # total log(1 + I) grows over time
#'
#' @export
gen_pde <- function(eta, Nx, Ny, N, nT_ori){
  # NOTE: Nx is length of x-dimension, or rows, Ny is the length of the y-dimension, or columns.
  # set up parameters
  Nx <- Nx
  Ny <- Ny
  N = N
  times <- seq(0, nT_ori - 1)

  # generate grid
  locs <- expand.grid(y=1:Ny, x=1:Nx)
  locs <- locs[,c("x", "y")]
  locs = as.numeric(rownames(locs))
  S = length(locs)
  Gridx <- ReacTran::setup.grid.1D(x.up = 0, x.down = Nx, N = Nx)
  Gridy <- ReacTran::setup.grid.1D(x.up = 0, x.down = Ny, N = Ny)

  # set up initial values for SIR
  S0 <- matrix(nrow = Nx, ncol = Ny, data = N)
  I0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
  I0[Nx/2,Ny/2] <- 1 # one infected person at the midpoint
  I0[round(Nx/10), round(Ny/10)] <- N * 0.05
  R0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
  yini <- c(S0, I0, R0)

  # pass parms to SIR within ode.2D
  parms <- list(eta = eta, Gridx = Gridx, Gridy = Gridy, Nx = Nx, Ny = Ny, N = N,
                y = yini, t = times)

  # Generate SIR data from the SIR ODE. lsodes needs about 50-110 doubles of
  # work space per state on 4 x 4 to 40 x 40 grids (slowly rising with size);
  # 500 per state leaves a wide margin without a fixed 1.6 GB allocation.
  n_state <- 3 * Nx * Ny
  out <- deSolve::ode.2D(y = yini, times = times, func = SIR, parms = parms,
                         nspec = 3, dimens = c(Nx, Ny),
                         lrw = max(1e5, 500 * n_state), names=c("X1", "X2", "X3"))[,-1]
  lgout <- log(out + 1)
  res <- lgout[,(1 + Nx * Ny):(2 * Nx * Ny)]
  return(res)
}

#' Compute the cross-covariance matrix between observed and new locations
#'
#' Evaluates the exponential GP kernel between observed and new locations,
#' giving the \eqn{N \times N_{\text{new}}} cross-covariance block \eqn{J_t}.
#' The nugget does not enter a cross-covariance, even when a new location
#' coincides with an observed one; this matches \code{\link{FFBS_predict_exact}}.
#'
#' @param input Numeric matrix. Coordinates of observed locations
#'   (\eqn{N \times d}).
#' @param input_new Numeric matrix or vector. Coordinates of new locations
#'   (\eqn{N_{\text{new}} \times d}).
#' @param phi_para Numeric scalar. Range parameter of the exponential kernel.
#' @param gp_sigma2 Numeric scalar. Marginal variance of the kernel. Defaults
#'   to \code{1.1}.
#' @param gp_tau2 Numeric scalar. Nugget variance. Not used; kept for
#'   compatibility. Defaults to \code{1e-4}.
#'
#' @return Numeric matrix of dimension \code{c(N, N_new)}.
#'
#' @seealso \code{\link{gen_exp_kernel}}
#' @examples
#' set.seed(1)
#' input <- matrix(runif(10), nrow = 5, ncol = 2)
#' J <- gen_Jt(input, input_new = c(0.5, 0.5), phi_para = 2)
#' length(J)
#' # At an observed input the cross-covariance equals gp_sigma2
#' gen_Jt(input, input[1, ], phi_para = 2)[1]
#'
#' @export
gen_Jt <- function(input, input_new, phi_para, gp_sigma2 = 1.1, gp_tau2 = 10^(-4)){
  # deal with one input prediction
  if(is.vector(input_new)){
    input_new <- t(input_new)
  }
  N <- dim(input)[1]
  distances <- as.matrix(stats::dist(rbind(input, input_new)))
  Jt <- gp_sigma2 * exp(-phi_para * distances[1:N, (N+1):ncol(distances)])
  return(Jt)
}

#' Update the likelihood of observations given PDE parameters (Monte Carlo)
#'
#' For a given candidate PDE parameter vector \code{eta}, simulates the
#' computer model output, computes the posterior predictive mean and covariance
#' at calibration locations using all \code{nsam} FFBS samples, and draws new
#' observations from the resulting predictive distribution.
#'
#' @param eta Numeric vector. PDE parameters (length 5).
#' @param nsam Integer. Number of posterior samples.
#' @param Nx Integer. Grid cells in x-direction.
#' @param Ny Integer. Grid cells in y-direction.
#' @param N_people Numeric. Total population size.
#' @param nT_ori Integer. Number of original time steps.
#' @param nT Integer. Number of episode time steps.
#' @param nT_cal Integer. Number of calibration time steps.
#' @param N_sp Integer. Total number of spatial locations.
#' @param N_sp_train Integer. Number of training spatial locations.
#' @param AR_choice Integer. AR order (1 or 2).
#' @param bncol Integer. Number of columns per block.
#' @param res_ffbs List. FFBS posterior samples.
#' @param n_b_cal Integer. Number of calibration blocks.
#' @param ind_block_cal Data frame. Block indices for calibration.
#' @param t_VinvY_FTheta_ls List. Precomputed \eqn{V^{-1}(Y - F\Theta)} terms.
#' @param Vinv Numeric matrix. Inverse of the observation covariance.
#' @param R_cal_train Numeric matrix. Correlation matrix for training locations.
#' @param input Numeric matrix. Observed location coordinates.
#' @param input_cal Numeric matrix. Calibration location coordinates.
#' @param phi_para Numeric scalar. GP range parameter.
#' @param gp_sigma2 Numeric scalar. GP marginal variance.
#' @param gp_tau2 Numeric scalar. GP nugget variance.
#' @param ind_train_cal Integer vector. Indices of training locations within
#'   calibration locations.
#'
#' @return A named list with:
#'   \describe{
#'     \item{y}{Array \code{c(nT_cal, N_sp_train, nsam)} of predictive samples.}
#'     \item{mu}{Array \code{c(nT_cal, N_sp_train, nsam)} of predictive means.}
#'     \item{sigma2_cal}{Numeric scalar. Predictive variance constant.}
#'     \item{Sigma}{Numeric matrix. Predictive covariance matrix.}
#'   }
#'
#' @seealso \code{\link{update_y_eta_one}}, \code{\link{gen_pde}}
#' @examples
#' set.seed(1)
#' input <- matrix(runif(10), nrow = 5, ncol = 2)
#' res_ffbs <- lapply(1:12, function(e) array(rnorm(8), c(2, 2, 2)))
#' res_ffbs$Sigma <- c(0.2, 0.3)
#' args <- list(eta = c(0.3, 0.1, 0.01, 0.01, 0.01), Nx = 2, Ny = 3, N_people = 100,
#'   nT_ori = 4, nT = 12, nT_cal = 3, N_sp = 6, N_sp_train = 4, AR_choice = 1,
#'   bncol = 2, res_ffbs = res_ffbs, n_b_cal = 3,
#'   ind_block_cal = generate_grid(1, 6, 1, 2),
#'   t_VinvY_FTheta_ls = lapply(1:3, function(t) array(rnorm(60), c(6, 5, 2))),
#'   Vinv = solve(gen_exp_kernel(input, phi = 2, sigma2 = 1.1, tau2 = 1e-4)),
#'   R_cal_train = diag(4), input = input, input_cal = c(0.4, 0.6), phi_para = 2,
#'   gp_sigma2 = 1.1, gp_tau2 = 1e-4, ind_train_cal = c(1, 3, 4, 6))
#' out <- do.call(update_y_eta, c(args, nsam = 2))
#' dim(out$y)  # nT_cal x N_sp_train x nsam
#' out$sigma2_cal
#'
#' @section Lifecycle:
#' Deprecated since spDBL 1.2.0 and will be removed in a future release.
#' The first call in a session warns; set
#' \code{options(spDBL.warn_deprecated = FALSE)} to silence it.
#' @export
update_y_eta <- function(eta, nsam, Nx, Ny, N_people, nT_ori, nT, nT_cal, N_sp, N_sp_train, AR_choice, bncol,
                     res_ffbs, n_b_cal, ind_block_cal, t_VinvY_FTheta_ls, Vinv, R_cal_train, input, input_cal, phi_para,
                     gp_sigma2, gp_tau2,
                     ind_train_cal){
  .legacy_deprecated("update_y_eta")
  res <- list()
  # 1. predict y_eta using PDE
  yt_eta_tt <- gen_pde(eta = eta, Nx = Nx, Ny = Ny, N = N_people,
                       nT_ori = nT_ori)
  yt_eta_tt_ls <- list()
  for (i in 1:dim(yt_eta_tt)[1]) {
    yt_eta_tt_ls[[i]] <- t(yt_eta_tt[i,])
  }

  # 2. calculate post mean and variance
  # 2.1 generate F_eta
  if(AR_choice == 1){
    F_eta_ls <- gen_F_ls_AR1_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  } else if(AR_choice == 2){
    F_eta_ls <- gen_F_ls_AR2_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  }

  F_eta_Theta_EP <- list()
  for (t in 1:nT) {
    F_eta_Theta_EP[[t]] <- array(dim = c(1, bncol, nsam))
    F_eta_ls_t <- F_eta_ls[[t]]
    for (i in 1:nsam) {
      F_eta_Theta_EP[[t]][,,i] <- F_eta_ls_t %*% res_ffbs[[t]][,,i]
    }
  }
  F_eta_Theta <- recover_from_EP_MC(dat_EP = F_eta_Theta_EP, nT_ori = nT_ori, nT = nT, nsam = nsam)

  # 2.2 generate J_eta
  J_eta <- gen_Jt(input = input, input_new = input_cal, phi_para = phi_para,
                  gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2)
  # 2.3 calculate mu
  mu_arry <- array(dim = c(nT_cal, N_sp, nsam))
  for (t in 1:nT_cal) {
    t_VinvY_FTheta_t <- t_VinvY_FTheta_ls[[t]]
    for (i in 1:nsam) {
      mu_arry[t,,i] <- F_eta_Theta[[t]][,,i] + as.vector(t_VinvY_FTheta_t[,,i] %*% J_eta)
    }
  }

  mu_arry_train <- array(dim = c(nT_cal, N_sp_train, nsam))
  for (t in 1:nT_cal) {
    mu_arry_train[t,,] <- mu_arry[t,ind_train_cal,]
  }
  # 2.4 calculate Sigma constant
  # Conditional variance of the new run, as in FFBS_predict_exact().
  Sigma_eta_c <- gp_sigma2 + gp_tau2 - crossprod(J_eta, Vinv) %*% J_eta
  sigma2_cal <- as.numeric(mean(res_ffbs$Sigma) * Sigma_eta_c)
  Sigma_eta_train <- sigma2_cal * R_cal_train
  # 3. sampling
  y_arry <- array(dim = c(nT_cal, N_sp_train, nsam))
  for (t in 1:nT_cal) {
    for (i in 1:nsam) {
      y_arry[t,,i] <- mniw::rmNorm(n = 1, mu = mu_arry_train[t,,i],
                                    Sigma = Sigma_eta_train)
    }
  }

  res$y <- y_arry
  res$mu <- mu_arry_train
  res$sigma2_cal <- sigma2_cal
  res$Sigma <- Sigma_eta_train
  return(res)
}

#' Update the likelihood of observations given PDE parameters (single sample)
#'
#' Like \code{\link{update_y_eta}} but uses only one FFBS sample (indexed by
#' \code{ind_sam}) instead of averaging over all samples, making it faster for
#' MCMC proposals.
#'
#' @param eta Numeric vector. PDE parameters (length 5).
#' @param ind_sam Integer. Index of the FFBS sample to use.
#' @param Nx Integer. Grid cells in x-direction.
#' @param Ny Integer. Grid cells in y-direction.
#' @param N_people Numeric. Total population size.
#' @param nT_ori Integer. Number of original time steps.
#' @param nT Integer. Number of episode time steps.
#' @param nT_cal Integer. Number of calibration time steps.
#' @param N_sp Integer. Total number of spatial locations.
#' @param N_sp_train Integer. Number of training spatial locations.
#' @param AR_choice Integer. AR order (1 or 2).
#' @param bncol Integer. Number of columns per block.
#' @param res_ffbs List. FFBS posterior samples.
#' @param n_b_cal Integer. Number of calibration blocks.
#' @param ind_block_cal Data frame. Block indices for calibration.
#' @param t_VinvY_FTheta_ls List. Precomputed \eqn{V^{-1}(Y - F\Theta)} terms.
#' @param Vinv Numeric matrix. Inverse of the observation covariance.
#' @param R_cal_train Numeric matrix. Correlation matrix for training locations.
#' @param input Numeric matrix. Observed location coordinates.
#' @param input_cal Numeric matrix. Calibration location coordinates.
#' @param phi_para Numeric scalar. GP range parameter.
#' @param gp_sigma2 Numeric scalar. GP marginal variance.
#' @param gp_tau2 Numeric scalar. GP nugget variance.
#' @param ind_train_cal Integer vector. Indices of training locations.
#'
#' @return A named list with \code{y}, \code{mu}, and \code{Sigma}.
#'
#' @seealso \code{\link{update_y_eta}}
#' @examples
#' set.seed(1)
#' input <- matrix(runif(10), nrow = 5, ncol = 2)
#' res_ffbs <- lapply(1:12, function(e) array(rnorm(8), c(2, 2, 2)))
#' res_ffbs$Sigma <- c(0.2, 0.3)
#' args <- list(eta = c(0.3, 0.1, 0.01, 0.01, 0.01), Nx = 2, Ny = 3, N_people = 100,
#'   nT_ori = 4, nT = 12, nT_cal = 3, N_sp = 6, N_sp_train = 4, AR_choice = 1,
#'   bncol = 2, res_ffbs = res_ffbs, n_b_cal = 3,
#'   ind_block_cal = generate_grid(1, 6, 1, 2),
#'   t_VinvY_FTheta_ls = lapply(1:3, function(t) array(rnorm(60), c(6, 5, 2))),
#'   Vinv = solve(gen_exp_kernel(input, phi = 2, sigma2 = 1.1, tau2 = 1e-4)),
#'   R_cal_train = diag(4), input = input, input_cal = c(0.4, 0.6), phi_para = 2,
#'   gp_sigma2 = 1.1, gp_tau2 = 1e-4, ind_train_cal = c(1, 3, 4, 6))
#' out <- do.call(update_y_eta_one, c(args, ind_sam = 2))
#' dim(out$y)  # nT_cal x N_sp_train
#' dim(out$Sigma)
#'
#' @section Lifecycle:
#' Deprecated since spDBL 1.2.0 and will be removed in a future release.
#' The first call in a session warns; set
#' \code{options(spDBL.warn_deprecated = FALSE)} to silence it.
#' @export
update_y_eta_one <- function(eta, ind_sam, Nx, Ny, N_people, nT_ori, nT, nT_cal, N_sp, N_sp_train, AR_choice, bncol,
                      res_ffbs, n_b_cal, ind_block_cal, t_VinvY_FTheta_ls, Vinv, R_cal_train, input, input_cal, phi_para,
                      gp_sigma2, gp_tau2, ind_train_cal){
  .legacy_deprecated("update_y_eta_one")
  res <- list()
  nsam <- 1 # generate only one sample
  # 1. predict y_eta using PDE
  yt_eta_tt <- gen_pde(eta = eta, Nx = Nx, Ny = Ny, N = N_people,
                       nT_ori = nT_ori)
  yt_eta_tt_ls <- list()
  for (i in 1:dim(yt_eta_tt)[1]) {
    yt_eta_tt_ls[[i]] <- t(yt_eta_tt[i,])
  }

  # 2. calculate post mean and variance
  # 2.1 generate F_eta
  if(AR_choice == 1){
    F_eta_ls <- gen_F_ls_AR1_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  } else if(AR_choice == 2){
    F_eta_ls <- gen_F_ls_AR2_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  }

  F_eta_Theta_EP <- list()
  for (t in 1:nT) {
    F_eta_ls_t <- F_eta_ls[[t]]
    for (i in ind_sam) {
      F_eta_Theta_EP[[t]] <- F_eta_ls_t %*% res_ffbs[[t]][,,i]
    }
  }
  F_eta_Theta <- recover_from_EP_exact(dat_EP = F_eta_Theta_EP, nT_ori = nT_ori, nT = nT)

  # 2.2 generate J_eta
  J_eta <- gen_Jt(input = input, input_new = input_cal, phi_para = phi_para,
                  gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2)
  # 2.3 calculate mu
  mu_arry <- array(dim = c(nT_cal, N_sp))
  for (t in 1:nT_cal) {
    t_VinvY_FTheta_t <- t_VinvY_FTheta_ls[[t]]
    for (i in ind_sam) {
      mu_arry[t,] <- F_eta_Theta[[t]] + as.vector(t_VinvY_FTheta_t[,,i] %*% J_eta)
    }
  }

  mu_arry_train <- array(dim = c(nT_cal, N_sp_train))
  for (t in 1:nT_cal) {
    mu_arry_train[t,] <- mu_arry[t,ind_train_cal]
  }
  # 2.4 calculate Sigma constant
  # Conditional variance of the new run, as in FFBS_predict_exact().
  Sigma_eta_c <- gp_sigma2 + gp_tau2 - crossprod(J_eta, Vinv) %*% J_eta
  sigma2_cal <- as.numeric(res_ffbs$Sigma[ind_sam] * Sigma_eta_c)
  Sigma_eta_train <- sigma2_cal * R_cal_train
  # 3. sampling
  y_arry <- array(dim = c(nT_cal, N_sp_train))
  for (t in 1:nT_cal) {
    for (i in ind_sam) {
      y_arry[t,] <- mniw::rmNorm(n = 1, mu = mu_arry_train[t,],
                                  Sigma = Sigma_eta_train)
    }
  }

  res$y <- y_arry
  res$mu <- mu_arry_train
  res$Sigma <- Sigma_eta_train
  return(res)
}

#' Compute posterior predictive mean and covariance without sampling (single sample)
#'
#' Like \code{\link{update_y_eta_one}} but skips the final sampling step,
#' returning only the predictive mean and covariance. Useful when sampling is
#' done separately (e.g., in \code{sample_y_eta_one}).
#'
#' @param eta Numeric vector. PDE parameters (length 5).
#' @param ind_sam Integer. Index of the FFBS sample to use.
#' @param Nx Integer. Grid cells in x-direction.
#' @param Ny Integer. Grid cells in y-direction.
#' @param N_people Numeric. Total population size.
#' @param nT_ori Integer. Number of original time steps.
#' @param nT Integer. Number of episode time steps.
#' @param nT_cal Integer. Number of calibration time steps.
#' @param N_sp Integer. Total number of spatial locations.
#' @param N_sp_train Integer. Number of training spatial locations.
#' @param AR_choice Integer. AR order (1 or 2).
#' @param bncol Integer. Number of columns per block.
#' @param res_ffbs List. FFBS posterior samples.
#' @param n_b_cal Integer. Number of calibration blocks.
#' @param ind_block_cal Data frame. Block indices for calibration.
#' @param t_VinvY_FTheta_ls List. Precomputed \eqn{V^{-1}(Y - F\Theta)} terms.
#' @param Vinv Numeric matrix. Inverse of the observation covariance.
#' @param R_cal_train Numeric matrix. Correlation matrix for training locations.
#' @param input Numeric matrix. Observed location coordinates.
#' @param input_cal Numeric matrix. Calibration location coordinates.
#' @param phi_para Numeric scalar. GP range parameter.
#' @param gp_sigma2 Numeric scalar. GP marginal variance.
#' @param gp_tau2 Numeric scalar. GP nugget variance.
#' @param ind_train_cal Integer vector. Indices of training locations.
#'
#' @return A named list with \code{mu} and \code{Sigma}.
#'
#' @seealso \code{\link{update_y_eta_one}}, \code{\link{sample_y_eta_one}}
#' @examples
#' set.seed(1)
#' input <- matrix(runif(10), nrow = 5, ncol = 2)
#' res_ffbs <- lapply(1:12, function(e) array(rnorm(8), c(2, 2, 2)))
#' res_ffbs$Sigma <- c(0.2, 0.3)
#' args <- list(eta = c(0.3, 0.1, 0.01, 0.01, 0.01), Nx = 2, Ny = 3, N_people = 100,
#'   nT_ori = 4, nT = 12, nT_cal = 3, N_sp = 6, N_sp_train = 4, AR_choice = 1,
#'   bncol = 2, res_ffbs = res_ffbs, n_b_cal = 3,
#'   ind_block_cal = generate_grid(1, 6, 1, 2),
#'   t_VinvY_FTheta_ls = lapply(1:3, function(t) array(rnorm(60), c(6, 5, 2))),
#'   Vinv = solve(gen_exp_kernel(input, phi = 2, sigma2 = 1.1, tau2 = 1e-4)),
#'   R_cal_train = diag(4), input = input, input_cal = c(0.4, 0.6), phi_para = 2,
#'   gp_sigma2 = 1.1, gp_tau2 = 1e-4, ind_train_cal = c(1, 3, 4, 6))
#' out <- do.call(update_muSigma_eta_one, c(args, ind_sam = 1))
#' dim(out$mu)  # nT_cal x N_sp_train
#' y <- sample_y_eta_one(1, nT_cal = 3, N_sp_train = 4, out$mu, out$Sigma)
#' dim(y)
#'
#' @section Lifecycle:
#' Deprecated since spDBL 1.2.0 and will be removed in a future release.
#' The first call in a session warns; set
#' \code{options(spDBL.warn_deprecated = FALSE)} to silence it.
#' @export
update_muSigma_eta_one <- function(eta, ind_sam, Nx, Ny, N_people, nT_ori, nT, nT_cal, N_sp, N_sp_train, AR_choice, bncol,
                         res_ffbs, n_b_cal, ind_block_cal, t_VinvY_FTheta_ls, Vinv, R_cal_train, input, input_cal, phi_para,
                         gp_sigma2, gp_tau2, ind_train_cal){
  .legacy_deprecated("update_muSigma_eta_one")
  res <- list()
  nsam <- 1 # generate only one sample
  # 1. predict y_eta using PDE
  yt_eta_tt <- gen_pde(eta = eta, Nx = Nx, Ny = Ny, N = N_people,
                       nT_ori = nT_ori)
  yt_eta_tt_ls <- list()
  for (i in 1:dim(yt_eta_tt)[1]) {
    yt_eta_tt_ls[[i]] <- t(yt_eta_tt[i,])
  }

  # 2. calculate post mean and variance
  # 2.1 generate F_eta
  if(AR_choice == 1){
    F_eta_ls <- gen_F_ls_AR1_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  } else if(AR_choice == 2){
    F_eta_ls <- gen_F_ls_AR2_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  }

  F_eta_Theta_EP <- list()
  for (t in 1:nT) {
    F_eta_ls_t <- F_eta_ls[[t]]
    for (i in ind_sam) {
      F_eta_Theta_EP[[t]] <- F_eta_ls_t %*% res_ffbs[[t]][,,i]
    }
  }
  F_eta_Theta <- recover_from_EP_exact(dat_EP = F_eta_Theta_EP, nT_ori = nT_ori, nT = nT)

  # 2.2 generate J_eta
  J_eta <- gen_Jt(input = input, input_new = input_cal, phi_para = phi_para,
                  gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2)
  # 2.3 calculate mu
  mu_arry <- array(dim = c(nT_cal, N_sp))
  for (t in 1:nT_cal) {
    t_VinvY_FTheta_t <- t_VinvY_FTheta_ls[[t]]
    for (i in ind_sam) {
      mu_arry[t,] <- F_eta_Theta[[t]] + as.vector(t_VinvY_FTheta_t[,,i] %*% J_eta)
    }
  }

  mu_arry_train <- array(dim = c(nT_cal, N_sp_train))
  for (t in 1:nT_cal) {
    mu_arry_train[t,] <- mu_arry[t,ind_train_cal]
  }
  # 2.4 calculate Sigma constant
  # Conditional variance of the new run, as in FFBS_predict_exact().
  Sigma_eta_c <- gp_sigma2 + gp_tau2 - crossprod(J_eta, Vinv) %*% J_eta
  sigma2_cal <- as.numeric(res_ffbs$Sigma[ind_sam] * Sigma_eta_c)
  Sigma_eta_train <- sigma2_cal * R_cal_train

  res$mu <- mu_arry_train
  res$Sigma <- Sigma_eta_train
  return(res)
}

#' Draw predictive samples from a precomputed mean and covariance
#'
#' Draws one sample from the predictive distribution
#' \eqn{y_t \sim N(\mu_t, \Sigma)} at each calibration time step, given
#' precomputed mean and covariance (e.g., from
#' \code{\link{update_muSigma_eta_one}}).
#'
#' @param ind_sam Integer. Sample index (used only for loop iteration; single
#'   sample is drawn).
#' @param nT_cal Integer. Number of calibration time steps.
#' @param N_sp_train Integer. Number of training spatial locations.
#' @param mu_arry_train Numeric array of dimension \code{c(nT_cal, N_sp_train)}.
#'   Predictive means.
#' @param Sigma_eta_train Numeric matrix. Predictive covariance
#'   (\code{N_sp_train x N_sp_train}).
#'
#' @return Numeric matrix of dimension \code{c(nT_cal, N_sp_train)} with one
#'   draw per time step.
#'
#' @seealso \code{\link{update_muSigma_eta_one}}
#' @examples
#' set.seed(1)
#' mu <- matrix(1:6, nrow = 2, ncol = 3)
#' y <- sample_y_eta_one(ind_sam = 1, nT_cal = 2, N_sp_train = 3,
#'                       mu_arry_train = mu, Sigma_eta_train = 0.01 * diag(3))
#' dim(y)
#' round(y - mu, 2)
#'
#' @section Lifecycle:
#' Deprecated since spDBL 1.2.0 and will be removed in a future release.
#' The first call in a session warns; set
#' \code{options(spDBL.warn_deprecated = FALSE)} to silence it.
#' @export
sample_y_eta_one <- function(ind_sam, nT_cal, N_sp_train, mu_arry_train, Sigma_eta_train){
  .legacy_deprecated("sample_y_eta_one")
  res <- list()
  # 3. sampling
  y_arry <- array(dim = c(nT_cal, N_sp_train))
  for (t in 1:nT_cal) {
    for (i in ind_sam) {
      y_arry[t,] <- mniw::rmNorm(n = 1, mu = mu_arry_train[t,],
                                 Sigma = Sigma_eta_train)
    }
  }

  res = y_arry
  return(res)
}

#' Update the posterior mean and covariance of the discrepancy field
#'
#' Computes the posterior mean vector \eqn{B_t b_t} and covariance matrix
#' \eqn{B_t} of the discrepancy \eqn{u_t} given the current predictive mean
#' \eqn{\mu_t}, observations \eqn{z_t - u_t}, and variance \eqn{\tau^2_t}.
#'
#' @param mut Numeric vector. Predictive mean at calibration locations.
#' @param invSigma Numeric matrix. Inverse of the predictive spatial covariance.
#' @param tau2t Numeric scalar. Current variance \eqn{\tau^2_t}.
#' @param zt_ut Numeric vector. Observation residuals \eqn{z_t - u_t}.
#' @param I_Stilde Numeric matrix. Identity-like matrix for the discrepancy
#'   prior precision (\code{N_sp_train x N_sp_train}).
#'
#' @return A named list with:
#'   \describe{
#'     \item{Bt}{Numeric matrix. Posterior covariance of the discrepancy.}
#'     \item{Btbt}{Numeric vector. Posterior mean of the discrepancy.}
#'   }
#'
#' @examples
#' set.seed(1)
#' Sigma <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.3, 0.1, 0.3, 1), 3)
#' mu <- rnorm(3)
#' zu <- rnorm(3)
#' post <- cal_Bt_bt(mu, invSigma = solve(Sigma), tau2t = 0.4, zt_ut = zu, I_Stilde = diag(3))
#' K <- Sigma %*% solve(Sigma + 0.4 * diag(3))  # Kalman gain
#' all.equal(post$Bt, Sigma - K %*% Sigma)
#' all.equal(as.vector(post$Btbt), as.vector(mu + K %*% (zu - mu)))
#'
#' @section Lifecycle:
#' Deprecated since spDBL 1.2.0 and will be removed in a future release.
#' The first call in a session warns; set
#' \code{options(spDBL.warn_deprecated = FALSE)} to silence it.
#' @export
cal_Bt_bt <- function(mut, invSigma, tau2t, zt_ut, I_Stilde){
  .legacy_deprecated("cal_Bt_bt")
  res <- list()
  invBt <- invSigma + 1/tau2t * I_Stilde
  Bt <- inv_chol(invBt)
  bt <- invSigma %*% mut + 1/tau2t * zt_ut
  Btbt <- Bt %*% bt
  res$Bt <- Bt
  res$Btbt <- Btbt
  return(res)
}
