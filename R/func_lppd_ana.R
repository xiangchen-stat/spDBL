## lppd analytical functions

#' Log density of the matrix-T distribution with inverse-gamma right covariance
#'
#' Evaluates the log (or raw) density of the matrix-T distribution that arises
#' when the right covariance is \eqn{\sigma^2 R} and \eqn{\sigma^2 \sim IG(\nu, d)}.
#'
#' @param Y Numeric matrix. Observation matrix (\eqn{p \times S}).
#' @param m Numeric matrix. Mean matrix (\eqn{p \times S}).
#' @param M Numeric matrix. Left covariance matrix (\eqn{p \times p},
#'   positive definite).
#' @param nu Numeric scalar. Degrees of freedom (shape parameter of the
#'   inverse-gamma on \eqn{\sigma^2}).
#' @param d Numeric scalar. Rate parameter of the inverse-gamma on
#'   \eqn{\sigma^2}).
#' @param R Numeric matrix. Right correlation matrix (\eqn{S \times S},
#'   positive definite).
#' @param log Logical. If \code{TRUE} (default), returns the log density.
#'
#' @return Numeric scalar. Log density (or density if \code{log = FALSE}).
#'
#' @export
dMTig <- function(Y, m, M, nu, d, R, log=TRUE) {
  p <- nrow(Y)
  S <- ncol(Y)
  Rinv <- solve(R)
  Minv <- solve(M)
  q <- Y - m
  Q <- Rinv %*% t(q) %*% Minv %*% q
  temp <- 1 + 1/(2*d) * LaplacesDemon::tr(Q)
  out <- lgamma(p*S/2 + nu) - lgamma(nu) - p*S/2 * log(2 * pi * d) -
    S/2 * LaplacesDemon::logdet(M) - p/2 * LaplacesDemon::logdet(R) -
    (p*S/2 + nu) * log(temp)

  if (log == TRUE) {
    return(out)
  } else{
    return(exp(out))
  }
}


#' One-step log posterior predictive density (MNIW / inverse-Wishart model)
#'
#' Computes the log posterior predictive density of \eqn{Y_t} at a single time
#' step under the MNIW model, evaluated at the smoothed state parameters from
#' \code{\link{FFBS}}.
#'
#' @param Yt Numeric matrix. Observed data at time \eqn{t} (\eqn{N \times S}).
#' @param Ft Numeric matrix. Covariate matrix at time \eqn{t} (\eqn{N \times p}).
#' @param Vt Numeric matrix. Observation noise left-covariance (\eqn{N \times N}).
#' @param st Numeric matrix. Smoothed state mean at time \eqn{t}
#'   (\eqn{p \times S}).
#' @param St Numeric matrix. Smoothed state left-covariance at time \eqn{t}
#'   (\eqn{p \times p}).
#' @param nt Numeric scalar. Filtered degrees of freedom at time \eqn{t}.
#' @param Dt Numeric matrix. Filtered scale matrix at time \eqn{t}
#'   (\eqn{S \times S}).
#'
#' @return Numeric scalar. Log posterior predictive density.
#'
#' @seealso \code{\link{lppd_IG_1t}}, \code{\link{lppd_id_1t}}
#' @export
lppd_IW_1t <- function(Yt, Ft, Vt, st, St, nt, Dt) {
  # note nu = nt - S + 1; hyper T is not matrix t
  S <- dim(Yt)[2]
  lppd <- mniw::dMT(X = Yt, Lambda = Ft %*% st, SigmaR = Ft %*% St %*% t(Ft) + Vt,
                    nu = nt - S + 1, SigmaC = Dt, log = TRUE)
  return(lppd)
}



#' One-step log posterior predictive density (scalar sigma-squared-times-R model)
#'
#' Computes the log posterior predictive density of \eqn{Y_t} at a single time
#' step under the scalar-sigma model (\eqn{\sigma^2 \sim IG}), using
#' \code{\link{dMTig}}.
#'
#' @param Yt Numeric matrix. Observed data at time \eqn{t} (\eqn{N \times S}).
#' @param Ft Numeric matrix. Covariate matrix at time \eqn{t} (\eqn{N \times p}).
#' @param Vt Numeric matrix. Observation noise left-covariance (\eqn{N \times N}).
#' @param st Numeric matrix. Smoothed state mean at time \eqn{t}
#'   (\eqn{p \times S}).
#' @param St Numeric matrix. Smoothed state left-covariance at time \eqn{t}
#'   (\eqn{p \times p}).
#' @param nt Numeric scalar. Filtered shape parameter at time \eqn{t}.
#' @param Dt Numeric scalar. Filtered rate parameter at time \eqn{t}.
#' @param R Numeric matrix. Fixed right correlation matrix (\eqn{S \times S}).
#'
#' @return Numeric scalar. Log posterior predictive density.
#'
#' @seealso \code{\link{dMTig}}, \code{\link{lppd_IW_1t}}
#' @export
lppd_IG_1t <- function(Yt, Ft, Vt, st, St, nt, Dt, R) {
  # note nu = nt - S + 1; hyper T is not matrix t
  S <- dim(Yt)[2]
  lppd <- dMTig(Y = Yt, m = Ft %*% st, M = Ft %*% St %*% t(Ft) + Vt,
                nu = nt, d = Dt, R = R, log = TRUE)
  return(lppd)
}

#' One-step log posterior predictive density (identity right-covariance model)
#'
#' Computes the log posterior predictive density of \eqn{Y_t} at a single time
#' step under the model with an identity right-covariance matrix, using a
#' matrix-normal density.
#'
#' @param Yt Numeric matrix. Observed data at time \eqn{t} (\eqn{N \times S}).
#' @param Ft Numeric matrix. Covariate matrix at time \eqn{t} (\eqn{N \times p}).
#' @param Vt Numeric matrix. Observation noise left-covariance (\eqn{N \times N}).
#' @param st Numeric matrix. Smoothed state mean at time \eqn{t}
#'   (\eqn{p \times S}).
#' @param St Numeric matrix. Smoothed state left-covariance at time \eqn{t}
#'   (\eqn{p \times p}).
#'
#' @return Numeric scalar. Log posterior predictive density.
#'
#' @seealso \code{\link{lppd_IW_1t}}, \code{\link{lppd_IG_1t}}
#' @export
lppd_id_1t <- function (Yt, Ft, Vt, st, St) {
  lppd <- mniw::dMNorm(X = Yt, Lambda = Ft %*% st,
                       SigmaR = make_pds(Ft %*% St %*% t(Ft) + Vt), SigmaC = diag(ncol(Yt)), log=TRUE)
  return(lppd)
}
