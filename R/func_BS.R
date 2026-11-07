# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

## MNIW----

#' Backward sampler for the Forward Filter Backward Sampler (FFBS)
#'
#' Performs the backward sampling pass of the FFBS algorithm using the filtered
#' distributions produced by \code{\link{FF}}. Iterates from time \eqn{T} back
#' to \eqn{1}, drawing smoothed state samples at each step via
#' \code{BS_1step_cpp}.
#'
#' @param res_ff List of length \code{nT} (plus a \code{prior} element) returned
#'   by \code{\link{FF}}. Each element \code{res_ff[[t]]} must contain
#'   \code{mt}, \code{Mt}, \code{at}, and \code{At}.
#' @param G_ls Either a single numeric matrix (state transition matrix \eqn{G},
#'   constant over time) or a list of \code{nT} matrices (time-varying).
#' @param nT Integer. Number of time steps.
#' @param delta Discount factor. Not implemented: values other than 1 give a
#'   warning and are ignored. Defaults to \code{1}.
#' @param verbose Logical. If \code{TRUE}, emits progress messages. Defaults
#'   to \code{FALSE}.
#'
#' @return A named list with:
#'   \describe{
#'     \item{st}{Array of dimension \code{c(p, q, nT)} containing the smoothed
#'       state means.}
#'     \item{St}{Array of dimension \code{c(p, p, nT)} containing the smoothed
#'       state left-covariance matrices.}
#'   }
#'
#' @seealso \code{\link{FF}}, \code{\link{FFBS}}
#' @examples
#' set.seed(1)
#' dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
#' par <- dat$para
#' res_ff <- FF(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
#'              par$n0, par$D0, nT = 4)
#' res_bs <- BS(res_ff, par$G0, nT = 4)
#' dim(res_bs$st)
#' # At the last time the smoothed and filtered means agree
#' all.equal(res_bs$st[, , 4], res_ff$T4$mt)
#' @export
BS <- function(res_ff, G_ls, nT, delta = 1, verbose = FALSE){
  delta <- .check_delta(delta)
  out <- list()
  out$st <- array(dim = c(dim(res_ff[[nT]]$mt), nT))
  out$St <- array(dim = c(dim(res_ff[[nT]]$Mt), nT))

  # initialize matirces
  Gt1 <- G_ls

  # Firstly, calculate t = nT
  st1 <- res_ff[[nT]]$mt
  St1 <- res_ff[[nT]]$Mt
  out$st[,,nT] <- st1
  out$St[,,nT] <- St1

  if (nT == 1L) return(out)

  # then backward t
  for (i in (nT-1):1) {
    if(verbose && i %% max(1, round(0.5*nT)) == 0){
      message("BS: ", i, "/", nT)
    }
    # Check if G_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(G_ls)){
      Gt1 <- as.matrix(G_ls[[i+1]])
    }

    para_nt <- BS_1step_cpp(mt = res_ff[[i]]$mt, Mt = res_ff[[i]]$Mt,
                            st1 = st1, St1 = St1,
                            at1 = res_ff[[i+1]]$at, At1 = res_ff[[i+1]]$At,
                            Gt1 = Gt1, delta = delta)
    # update st1, St1
    st1 <- para_nt$st
    St1 <- para_nt$St
    out$st[,,i] <- st1
    out$St[,,i] <- St1
  }

  return(out)
}
