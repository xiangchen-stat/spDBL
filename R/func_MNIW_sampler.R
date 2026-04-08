#' Check and repair a matrix to be positive definite and symmetric
#'
#' Verifies that a matrix is approximately symmetric (within a relative and
#' absolute tolerance). If so, symmetrises it and, if still not positive
#' definite, adds a small diagonal ridge. Raises an error if the asymmetry
#' exceeds the tolerance.
#'
#' @param C Numeric matrix to check.
#' @param eps Numeric scalar. Absolute tolerance for the Frobenius norm of
#'   \code{C - t(C)}. Defaults to \code{1e-5}.
#' @param per Numeric scalar. Relative tolerance: the ratio of the asymmetry
#'   norm to the mean of \code{C}. Defaults to \code{0.05}.
#'
#' @return The (possibly corrected) positive definite symmetric matrix.
#'
#' @seealso \code{\link{make_pds}}
#' @export
check_pds <- function(C, eps = 10^(-5), per = 0.05) {
  norm_f <- norm(C - t(C), type="F")
  condition <- ((norm_f / mean(C)) < per) || (norm_f < eps)
  if (condition) {
    C <- (C + t(C))/2
    if (!matrixcalc::is.positive.definite(C)) {
      C <- C + 3 * 10^(-4) * diag(dim(C)[1])
      warning("Matrix is not positive definite.")
    }
    return(C)
  } else {
    print(norm_f)
    stop("Matrix symmetry is outside of tolerance: ",eps)
  }
}


#' Force a matrix to be positive definite and symmetric
#'
#' Symmetrises a matrix by averaging it with its transpose and, if the result
#' is not positive definite, adds a small diagonal ridge.
#'
#' @param C Numeric matrix.
#' @param eps Numeric scalar. Size of the diagonal ridge added when \code{C}
#'   is not positive definite. Defaults to \code{1e-4}.
#'
#' @return The corrected positive definite symmetric matrix.
#'
#' @seealso \code{\link{check_pds}}
#' @export
make_pds <- function(C, eps = 10^(-4)) {
  C <- (C + t(C))/2
  if (!matrixcalc::is.positive.definite(C)) {
    C <- C + eps * diag(dim(C)[1])
    warning("Matrix is not positive definite.")
  }
  return(C)
}

#' Sample from the Matrix Normal Inverse Wishart (MNIW) distribution
#'
#' Draws \code{nsam} posterior samples of the regression coefficient matrix
#' \eqn{B} and the covariance matrix \eqn{\Sigma} under the MNIW model:
#' \deqn{Y = X B + E, \quad E \sim MN(0, H, \Sigma)}
#' \deqn{\Sigma \sim IW(v, S), \quad B | \Sigma \sim MN(C, V_b, \Sigma).}
#'
#' @param nsam Integer. Number of samples to draw.
#' @param X Numeric matrix. Covariate (design) matrix (\eqn{n \times p}).
#' @param v Numeric scalar. Degrees of freedom of the inverse-Wishart prior on
#'   \eqn{\Sigma}.
#' @param S Numeric matrix. Scale matrix of the inverse-Wishart prior
#'   (\eqn{q \times q}, positive definite).
#' @param C Numeric matrix. Prior mean of \eqn{B | \Sigma} (\eqn{p \times q}).
#' @param Vb Numeric matrix. Prior row covariance of \eqn{B | \Sigma}
#'   (\eqn{p \times p}, positive definite).
#'
#' @return A named list with:
#'   \describe{
#'     \item{sigma}{Array of dimension \code{c(q, q, nsam)} with samples of
#'       \eqn{\Sigma}.}
#'     \item{B}{Array of dimension \code{c(p, q, nsam)} with samples of
#'       \eqn{B}.}
#'   }
#'
#' @seealso \code{\link{MNIG_sampler}}, \code{\link{FFBS_sampling}}
#' @export
MNIW_sampler <- function(nsam, X, v, S, C, Vb){
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(S)

  sigma <- matrixsampling::rinvwishart(nsam, nu = v, Omega = S, epsilon = 0, checkSymmetry = F)
  B <- array(dim = c(p, q, nsam))

  for (i in 1:nsam) {
    B[,,i] <- rmn_cpp(m = C, U = Vb, V = sigma[,,i])
  }

  out <- list(sigma = sigma, B = B)
  return(out)
}

#' Sample from the Matrix Normal Inverse Gamma (MNIG) distribution
#'
#' Draws \code{nsam} posterior samples of the regression coefficient matrix
#' \eqn{B} and the scalar variance \eqn{\sigma^2} under the MNIG model where
#' the right covariance of \eqn{B | \sigma^2} is \eqn{\sigma^2 R}:
#' \deqn{\sigma^2 \sim IG(v, S), \quad B | \sigma^2 \sim MN(C, V_b, \sigma^2 R).}
#'
#' @param nsam Integer. Number of samples to draw.
#' @param X Numeric matrix. Covariate (design) matrix (\eqn{n \times p}).
#' @param v Numeric scalar. Shape parameter of the inverse-gamma prior on
#'   \eqn{\sigma^2}.
#' @param S Numeric scalar. Rate parameter of the inverse-gamma prior on
#'   \eqn{\sigma^2}.
#' @param C Numeric matrix. Prior mean of \eqn{B | \sigma^2} (\eqn{p \times q}).
#' @param Vb Numeric matrix. Prior row covariance of \eqn{B | \sigma^2}
#'   (\eqn{p \times p}, positive definite).
#' @param R Numeric matrix. Fixed right correlation matrix (\eqn{q \times q},
#'   positive definite).
#'
#' @return A named list with:
#'   \describe{
#'     \item{sigma2}{Numeric vector of length \code{nsam} with samples of
#'       \eqn{\sigma^2}.}
#'     \item{B}{Array of dimension \code{c(p, q, nsam)} with samples of
#'       \eqn{B}.}
#'   }
#'
#' @seealso \code{\link{MNIW_sampler}}, \code{\link{FFBS_sampling_sigma2R}}
#' @export
MNIG_sampler <- function(nsam, X, v, S, C, Vb, R){
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(C)

  sigma2 <- invgamma::rinvgamma(n = nsam, shape = v, rate = S)
  B <- array(dim = c(p, q, nsam))

  for (i in 1:nsam) {
    B[,,i] <- rmn_cpp(m = C, U = Vb, V = sigma2[i] * R)
  }

  out <- list(sigma2 = sigma2, B = B)
  return(out)
}
