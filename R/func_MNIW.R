# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

#' Naive MNIW posterior update
#'
#' Computes the posterior parameters of the Matrix Normal Inverse Wishart (MNIW)
#' model using explicit matrix transposes. Equivalent to \code{\link{MNIW_R}} but
#' without \code{crossprod} optimisations. Intended for reference and testing.
#'
#' The model is \eqn{Y = X B + E}, where \eqn{E \sim MN(0, H, \Sigma)},
#' \eqn{\Sigma \sim IW(v, S)}, and \eqn{B | \Sigma \sim MN(C, V_b, \Sigma)}.
#'
#' @param X Numeric matrix. Covariate (design) matrix of dimension \eqn{n \times p}.
#' @param Y Numeric matrix. Response matrix of dimension \eqn{n \times q}.
#' @param newX Unused legacy parameter. Defaults to \code{NULL}.
#' @param newY Unused legacy parameter. Defaults to \code{NULL}.
#' @param H Numeric matrix. Row covariance matrix of the error term (\eqn{n \times n}).
#' @param v Numeric scalar. Prior degrees of freedom of the inverse-Wishart on \eqn{\Sigma}.
#' @param S Numeric matrix. Prior scale matrix of the inverse-Wishart on \eqn{\Sigma}.
#' @param C Numeric matrix. Prior mean matrix of \eqn{B | \Sigma} (\eqn{p \times q}).
#' @param Vb Numeric matrix. Prior row covariance matrix of \eqn{B | \Sigma} (\eqn{p \times p}).
#'
#' @return A named list with posterior parameters:
#'   \describe{
#'     \item{vnew}{Posterior degrees of freedom.}
#'     \item{Snew}{Posterior scale matrix of the inverse-Wishart.}
#'     \item{Cnew}{Posterior mean matrix of \eqn{B | \Sigma}.}
#'     \item{Vbnew}{Posterior row covariance matrix of \eqn{B | \Sigma}.}
#'   }
#'
#' @seealso \code{\link{MNIW_R}}
#' @examples
#' set.seed(1)
#' n <- 30
#' X <- cbind(1, rnorm(n))
#' Y <- X %*% matrix(c(1, -1, 0.5, 2), 2) + matrix(rnorm(n * 2, sd = 0.1), n, 2)
#' prior <- list(H = diag(n), v = 4, S = diag(2), C = matrix(0, 2, 2), Vb = diag(10, 2))
#' post <- do.call(MNIW_R_naiive, c(list(X = X, Y = Y), prior))
#' round(post$Cnew, 2)
#' # Same result as the crossprod version
#' all.equal(post, do.call(MNIW_R, c(list(X = X, Y = Y), prior)))
#'
#' @export
MNIW_R_naiive <- function(X, Y, newX = NULL, newY = NULL, H, v, S, C, Vb){
  out <- list()
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)

  Vbinv <- solve(Vb)
  Hinv <- solve(H)
  Xt <- t(X)
  Yt <- t(Y)
  Ct <- t(C)
  Oinv <- Vbinv + Xt %*% Hinv %*% X
  O <- solve(Oinv)
  M <- Vbinv %*% C + Xt %*% Hinv %*% Y
  Mt <- t(M)
  Snew <- S + Ct %*% Vbinv %*% C + Yt %*% Hinv %*% Y - Mt %*% O %*% M
  vnew <- v + n
  Cnew <- O %*% M
  Vbnew <- O

  out[["vnew"]] <- vnew
  out[["Snew"]] <- Snew
  out[["Cnew"]] <- Cnew
  out[["Vbnew"]] <- Vbnew
  return(out)
}

#' MNIW posterior update
#'
#' Computes the posterior parameters of the Matrix Normal Inverse Wishart (MNIW)
#' model using \code{crossprod} for improved numerical efficiency.
#'
#' The model is \eqn{Y = X B + E}, where \eqn{E \sim MN(0, H, \Sigma)},
#' \eqn{\Sigma \sim IW(v, S)}, and \eqn{B | \Sigma \sim MN(C, V_b, \Sigma)}.
#'
#' @param X Numeric matrix. Covariate (design) matrix of dimension \eqn{n \times p}.
#' @param Y Numeric matrix. Response matrix of dimension \eqn{n \times q}.
#' @param newX Unused legacy parameter. Defaults to \code{NULL}.
#' @param newY Unused legacy parameter. Defaults to \code{NULL}.
#' @param H Numeric matrix. Row covariance matrix of the error term (\eqn{n \times n}).
#' @param v Numeric scalar. Prior degrees of freedom of the inverse-Wishart on \eqn{\Sigma}.
#' @param S Numeric matrix. Prior scale matrix of the inverse-Wishart on \eqn{\Sigma}.
#' @param C Numeric matrix. Prior mean matrix of \eqn{B | \Sigma} (\eqn{p \times q}).
#' @param Vb Numeric matrix. Prior row covariance matrix of \eqn{B | \Sigma} (\eqn{p \times p}).
#'
#' @return A named list with posterior parameters:
#'   \describe{
#'     \item{vnew}{Posterior degrees of freedom.}
#'     \item{Snew}{Posterior scale matrix of the inverse-Wishart.}
#'     \item{Cnew}{Posterior mean matrix of \eqn{B | \Sigma}.}
#'     \item{Vbnew}{Posterior row covariance matrix of \eqn{B | \Sigma}.}
#'   }
#'
#' @seealso \code{\link{MNIW_R_naiive}}
#' @examples
#' set.seed(1)
#' n <- 30
#' X <- cbind(1, rnorm(n))
#' B <- matrix(c(1, -1, 0.5, 2, 0, 1), nrow = 2)
#' Y <- X %*% B + matrix(rnorm(n * 3, sd = 0.1), n, 3)
#' post <- MNIW_R(X, Y, H = diag(n), v = 5, S = diag(3),
#'                C = matrix(0, 2, 3), Vb = diag(10, 2))
#' post$vnew
#' round(post$Cnew, 2)  # close to B
#'
#' @export
MNIW_R <- function(X, Y, newX = NULL, newY = NULL, H, v, S, C, Vb){
  out <- list()
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)

  Vbinv <- solve(Vb)
  Hinv <- solve(H)
  XtHinv <- crossprod(X, Hinv)
  Oinv <- Vbinv + XtHinv %*% X
  O <- solve(Oinv)
  VbinvC <- Vbinv %*% C
  M <- VbinvC + XtHinv %*% Y
  Snew <- S + crossprod(C, VbinvC) + crossprod(Y, Hinv) %*% Y - crossprod(M, O) %*% M
  vnew <- v + n
  Cnew <- O %*% M
  Vbnew <- O

  out[["vnew"]] <- vnew
  out[["Snew"]] <- Snew
  out[["Cnew"]] <- Cnew
  out[["Vbnew"]] <- Vbnew
  return(out)
}
