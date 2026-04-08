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
#' @param newX Unused legacy parameter. Defaults to \code{none}.
#' @param newY Unused legacy parameter. Defaults to \code{none}.
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
#' @param newX Unused legacy parameter. Defaults to \code{none}.
#' @param newY Unused legacy parameter. Defaults to \code{none}.
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
