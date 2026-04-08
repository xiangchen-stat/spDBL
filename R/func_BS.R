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
#' @param delta Numeric scalar. Right-variance discount factor. Defaults to
#'   \code{1}.
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
#' @export
BS <- function(res_ff, G_ls, nT, delta = 1){
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

  # then backward t
  for (i in (nT-1):1) {
    if(i %% round(0.1*nT) == 0){
      print(paste("BS:", i, "/", nT))
      print(Sys.time())
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
