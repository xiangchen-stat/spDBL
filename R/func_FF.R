#' Read a rectangular block from a large CSV file
#'
#' Efficiently reads a contiguous rectangular block of rows and columns from a
#' CSV file using \code{readr::read_csv}. Useful for processing big data in
#' chunks without loading the entire file into memory.
#'
#' @param filename Character. Path to the CSV file.
#' @param rows Integer vector of length 2 giving the first and last row indices
#'   to read (1-based, excluding the header row if \code{header = TRUE}).
#'   Use \code{c(1, Inf)} to read all rows. Defaults to \code{c(1, Inf)}.
#' @param cols Integer vector of length 2 giving the first and last column
#'   indices to read. Use \code{c(1, Inf)} to read all columns. Defaults to
#'   \code{c(1, Inf)}.
#' @param header Logical. Whether the file has a header row. Defaults to
#'   \code{TRUE}.
#'
#' @return A \code{tibble} (from \code{readr}) containing the requested block.
#'
#' @export
read_big_csv_quick <- function(filename, rows = c(1,Inf), cols = c(1,Inf), header=TRUE) {

  if (is.finite(cols[2])) {
    #TODO: How to suppress progress bar here?
    out <- read_csv(file = filename, col_names = FALSE,
                    skip = rows[1] - 1 + as.integer(header),
                    n_max = rows[2] - rows[1] + 1,
                    col_select = cols[1]:cols[2],
                    show_col_types = FALSE)

  } else {
    #TODO: something is going on with this read_csv that is throwing an error, but when I run it in the R shell, it works fine. Oddly, neglecting this option does not throw an error.
    out <- read_csv(file = filename, col_names = FALSE,
                    skip = rows[1] - 1 + as.integer(header),
                    n_max = rows[2] - rows[1] + 1,
                    show_col_types = FALSE)
  }

  return(out)
}


## MNIW----

# TODO FF_cpp
#' Forward Filter for the MNIW dynamic linear model
#'
#' Runs the full forward filtering pass over \code{nT} time steps under the
#' Matrix Normal Inverse Wishart (MNIW) model. At each step, calls
#' \code{FF_1step_cpp} to update the state mean and covariance as well as the
#' inverse-Wishart parameters for the right-covariance matrix \eqn{\Sigma}.
#'
#' @param Y List of length \code{nT}. Each element \code{Y[[t]]} is the
#'   \eqn{N \times q} data matrix at time \eqn{t}.
#' @param F_ls Either a single \eqn{N \times p} covariate matrix (constant
#'   over time) or a list of \code{nT} such matrices (time-varying).
#' @param G_ls Either a single \eqn{p \times p} state transition matrix or a
#'   list of \code{nT} such matrices.
#' @param W_ls Either a single \eqn{p \times p} state noise left-covariance
#'   matrix or a list of \code{nT} such matrices.
#' @param V_ls Either a single \eqn{N \times N} observation noise left-covariance
#'   matrix or a list of \code{nT} such matrices.
#' @param m0 Numeric matrix. Prior mean of the state \eqn{\beta_0 | D_0}
#'   (\eqn{p \times q}).
#' @param M0 Numeric matrix. Prior left-covariance of \eqn{\beta_0 | D_0}
#'   (\eqn{p \times p}).
#' @param n0 Numeric scalar. Prior degrees of freedom of the inverse-Wishart on
#'   \eqn{\Sigma | D_0}.
#' @param D0 Numeric matrix. Prior scale matrix of the inverse-Wishart on
#'   \eqn{\Sigma | D_0} (\eqn{q \times q}).
#' @param nT Integer. Number of time steps.
#' @param delta Numeric scalar. Discount factor for the right-variance matrix.
#'   Defaults to \code{1.0}.
#'
#' @return A named list of length \code{nT + 1}. Elements \code{"T1"} through
#'   \code{"T<nT>"} each contain a list with filtered parameters
#'   \code{nt}, \code{Dt}, \code{at}, \code{At}, \code{mt}, \code{Mt}.
#'   The additional element \code{prior} stores \code{m0}, \code{M0},
#'   \code{n0}, \code{D0}.
#'
#' @seealso \code{\link{BS}}, \code{\link{FFBS}}
#' @export
FF <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, delta = 1.0){
  out <- list()
  out_name <- c()
  # initialize matirces
  Ft <- F_ls
  Gt <- G_ls
  Wt <- W_ls
  Vt <- V_ls

  for (i in 1:nT) {
    if(i %% round(0.1*nT) == 0){
      print(paste("FF:", i, "/", nT))
      print(Sys.time())
    }
    # Check if F_ls, G_ls, W_ls, V_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(F_ls)){
      Ft <- as.matrix(F_ls[[i]])
    }
    if(is.list(G_ls)){
      Gt <- as.matrix(G_ls[[i]])
    }
    if(is.list(W_ls)){
      Wt <- as.matrix(W_ls[[i]])
    }
    if(is.list(V_ls)){
      Vt <- as.matrix(V_ls[[i]])
    }

    Yt <- Y[[i]]
    # Set T-1 values
    if (i == 1) {
      mt_1 = m0
      Mt_1 = M0
      nt_1 = n0
      Dt_1 = D0
    } else{
      mt_1 = one_step$mt
      Mt_1 = one_step$Mt
      nt_1 = one_step$nt
      Dt_1 = one_step$Dt
    }

    # perform 1 step forward
    one_step <- FF_1step_cpp(Yt = Yt, Ft = Ft, Gt = Gt,
                             Wt = Wt, Vt = Vt,
                             mt_1 = mt_1, Mt_1 = Mt_1,
                             nt_1 = nt_1, Dt_1 = Dt_1, delta = delta)
    one_step <- one_step[c("nt", "Dt", "at", "At", "mt", "Mt")] # extract parameters
    out[[i]] <- one_step # save results
  }

  out_name <- paste("T", seq(1:nT), sep = "") # accumulate output list names
  names(out) <- out_name
  out$prior <- list(m0 = m0, M0 = M0, n0 = n0, D0 = D0)
  return(out)
}



#' Forward Filter for big data stored in CSV files (MNIW model)
#'
#' Runs the forward filtering pass for large datasets that are partitioned into
#' spatial blocks and stored as CSV files on disk. At each time step the spatial
#' domain is traversed block-by-block using a snake traversal order (see
#' \code{\link{generate.grid.rowsnake}}), and results are written to disk to
#' avoid exhausting memory.
#'
#' @param Y_ls Character vector of length \code{nT}. File paths to the response
#'   CSV files, one per time step.
#' @param F_ls Either a single file path (constant \eqn{F_t}) or a character
#'   vector of length \code{nT} (time-varying). Can also contain a second set
#'   of columns (AR2 lags) for the same file at offset \code{fncol}.
#' @param G_ls Either a single file path (constant \eqn{G_t}) or a character
#'   vector of length \code{nT} (time-varying).
#' @param W_ls Either a single file path (constant \eqn{W_t}) or a character
#'   vector of length \code{nT} (time-varying).
#' @param V_ls Either a single file path (constant \eqn{V_t}) or a character
#'   vector of length \code{nT} (time-varying). Expected to be a square block
#'   submatrix indexed by the row range of each spatial block.
#' @param m0 Numeric matrix. Prior state mean (\eqn{p \times q}).
#' @param M0 Numeric matrix. Prior left-covariance of the state (\eqn{p \times p}).
#' @param n0 Numeric scalar. Prior degrees of freedom of the inverse-Wishart.
#' @param D0 Numeric matrix. Prior scale matrix of the inverse-Wishart.
#' @param nT Integer. Number of time steps (files).
#' @param fnrow Integer. Total number of rows in each data file.
#' @param fncol Integer. Total number of columns in each data file.
#' @param bnrow Integer. Number of rows per spatial block.
#' @param bncol Integer. Number of columns per spatial block.
#' @param path_out Character. Directory path where output CSV files are written.
#' @param delta Numeric scalar. Right-variance discount factor. Defaults to
#'   \code{1.0}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages. Defaults to
#'   \code{FALSE}.
#'
#' @return A named list of character vectors with paths to the output CSV files:
#'   \describe{
#'     \item{nt_ls}{Paths to filtered degrees-of-freedom files.}
#'     \item{Dt_ls}{Paths to filtered scale-matrix files.}
#'     \item{mt_ls}{Paths to filtered state-mean files.}
#'     \item{Mt_ls}{Paths to filtered left-covariance files.}
#'     \item{at_ls}{Paths to predicted state-mean files.}
#'     \item{At_ls}{Paths to predicted left-covariance files.}
#'     \item{FT}{Path to the covariate matrix at the final time step.}
#'     \item{VT}{Path to the noise covariance matrix at the final time step.}
#'   }
#'
#' @seealso \code{\link{FF}}, \code{\link{generate.grid.rowsnake}}
#' @export
FF_bigdata_R <- function(Y_ls, F_ls, G_ls, W_ls, V_ls,
                                m0, M0, n0, D0, nT,
                                fnrow, fncol, bnrow, bncol, path_out,
                                delta = 1.0, verbose = FALSE){
  # generate index
  ind <- generate.grid.rowsnake(fnrow = fnrow, fncol = fncol, bnrow = bnrow, bncol = bncol)
  n_b <- dim(ind)[1] # number of blocks

  # initialize matirces if Gt and Wt never changes
  if(length(G_ls) == 1){
    Gt <- as.matrix(read_big_csv_quick(filename = G_ls))
  }

  if(length(W_ls) == 1){
    Wt <- as.matrix(read_big_csv_quick(filename = W_ls))
  }

  # loop through nT time/files
  for (f in 1:nT) {
    out <- list()

    # Check if G_ls, W_ls are changing over file
    # If they are array of filename, change accordingly, otherwise keep unchanged
    if(length(G_ls) != 1){
      Gt <- as.matrix(read_big_csv_quick(filename = G_ls[f]))
    }

    if(length(W_ls) != 1){
      Wt <- as.matrix(read_big_csv_quick(filename = W_ls[f]))
    }

    # indicator to avoid read Vt each time
    row_vt_old <- 0

    # loop through blocks
    for (i in 1:n_b) {
      # Extract Yt, Ft, Vt from files
      Yt <- as.matrix(read_big_csv_quick(filename = Y_ls[f],
                                       rows = c(ind[i,2], ind[i,3]),
                                       cols = c(ind[i,4], ind[i,5])))
      # extract Ft
      if(length(F_ls) != 1){
        # read ar1
        Ft <- as.matrix(read_big_csv_quick(filename = F_ls[f],
                                           rows = c(ind[i,2], ind[i,3]),
                                           cols = c(ind[i,4], ind[i,5])))
        # try read ar2
        try(expr = {
            Ft_ar2 <- as.matrix(read_big_csv_quick(filename = F_ls[f],
                                                   rows = c(ind[i,2], ind[i,3]),
                                                   cols = c((fncol + ind[i,4]), (fncol + ind[i,5]))));
            Ft <- cbind(Ft, Ft_ar2);
          },
          silent = T)
      } else{
        # read ar1
        Ft <- as.matrix(read_big_csv_quick(filename = F_ls,
                                           rows = c(ind[i,2], ind[i,3]),
                                           cols = c(ind[i,4], ind[i,5])))
        # try read ar2
        try(expr = {
          Ft_ar2 <- as.matrix(read_big_csv_quick(filename = F_ls,
                                                 rows = c(ind[i,2], ind[i,3]),
                                                 cols = c((fncol + ind[i,4]), (fncol + ind[i,5]))));
          Ft <- cbind(Ft, Ft_ar2);
        },
        silent = T)
      }

      # extract Vt
      row_vt_new <- ind[i,2]
      if(row_vt_new != row_vt_old){
        # extract Vt
        if(length(V_ls) != 1){
          Vt <- as.matrix(read_big_csv_quick(filename = V_ls[f],
                                             rows = c(ind[i,2], ind[i,3]),
                                             cols = c(ind[i,2], ind[i,3])))
        } else{
          Vt <- as.matrix(read_big_csv_quick(filename = V_ls,
                                             rows = c(ind[i,2], ind[i,3]),
                                             cols = c(ind[i,2], ind[i,3])))
        }
      }
      row_vt_old <- row_vt_new # update indicator for Ft and Vt 1st row


      # Set T-1 values
      if (i == 1 && f == 1) {
        mt_1 = m0
        Mt_1 = M0
        nt_1 = n0
        Dt_1 = D0
      } else{
        mt_1 = one_step$mt
        Mt_1 = one_step$Mt
        nt_1 = one_step$nt
        Dt_1 = one_step$Dt
      }

      # perform 1 step forward
      one_step <- FF_1step_cpp(Yt = Yt, Ft = Ft, Gt = Gt,
                               Wt = Wt, Vt = Vt,
                               mt_1 = mt_1, Mt_1 = Mt_1,
                               nt_1 = nt_1, Dt_1 = Dt_1, delta = delta)
      if(i == 1){
        out$nt <- one_step$nt
        out$Dt <- one_step$Dt
        out$at <- one_step$at
        out$At <- one_step$At
        out$mt <- one_step$mt
        out$Mt <- one_step$Mt
      } else{
        out$nt <- cbind(out$nt, one_step$nt)
        out$Dt <- cbind(out$Dt, one_step$Dt)
        out$at <- cbind(out$at, one_step$at)
        out$At <- cbind(out$At, one_step$At)
        out$mt <- cbind(out$mt, one_step$mt)
        out$Mt <- cbind(out$Mt, one_step$Mt)
      }
    }

    # after looping one time/file, save results
    write.csv(out$nt, file = paste0(path_out, "/output_FF_nt_", f, ".csv",
                                    sep = ""), row.names = F)
    write.csv(out$Dt, file = paste0(path_out, "/output_FF_Dt_", f, ".csv",
                                    sep = ""), row.names = F)
    write.csv(out$mt, file = paste0(path_out, "/output_FF_mt_", f, ".csv",
                                    sep = ""), row.names = F)
    write.csv(out$Mt, file = paste0(path_out, "/output_FF_MMt_", f, ".csv",
                                    sep = ""), row.names = F)
    write.csv(out$at, file = paste0(path_out, "/output_FF_at_", f, ".csv",
                                    sep = ""), row.names = F)
    write.csv(out$At, file = paste0(path_out, "/output_FF_AAt_", f, ".csv",
                                    sep = ""), row.names = F)
    if(verbose == TRUE){
      print(paste0("Forward filter", f, "/", nT))
      print(Sys.time())
    }
  }

  # Save last Ft and Vt for BS 1st sample
  write.csv(Ft, file = paste0(path_out, "/output_FF_FT.csv",
                              sep = ""), row.names = F)
  write.csv(Vt, file = paste0(path_out, "/output_FF_VT.csv",
                              sep = ""), row.names = F)

  file_ls <- list()
  file_ls$nt_ls <- paste0(path_out, "/output_FF_nt_", seq(1:nT), ".csv", sep = "")
  file_ls$Dt_ls <- paste0(path_out, "/output_FF_Dt_", seq(1:nT), ".csv", sep = "")
  file_ls$mt_ls <- paste0(path_out, "/output_FF_mt_", seq(1:nT), ".csv", sep = "")
  file_ls$Mt_ls <- paste0(path_out, "/output_FF_MMt_", seq(1:nT), ".csv", sep = "")
  file_ls$at_ls <- paste0(path_out, "/output_FF_at_", seq(1:nT), ".csv", sep = "")
  file_ls$At_ls <- paste0(path_out, "/output_FF_AAt_", seq(1:nT), ".csv", sep = "")
  file_ls$FT <- paste0(path_out, "/output_FF_FT.csv", sep = "")
  file_ls$VT <- paste0(path_out, "/output_FF_VT.csv", sep = "")

  return(file_ls)
}

##sigma2R----

#' Single forward filter step (scalar right-covariance, sigma-squared times R)
#'
#' Performs one step of the forward filter under the model where the right
#' covariance of \eqn{Y} is \eqn{\sigma^2 R} with \eqn{\sigma^2 \sim IG(n, d)}.
#' The inverse of \eqn{R} (\code{Rinv}) must be precomputed and passed in.
#'
#' @param Yt Numeric matrix. Data matrix at time \eqn{t} (\eqn{N \times S}).
#' @param Ft Numeric matrix. Covariate matrix at time \eqn{t} (\eqn{N \times p}).
#' @param Gt Numeric matrix. State transition matrix at time \eqn{t}
#'   (\eqn{p \times p}).
#' @param Wt Numeric matrix. State noise left-covariance at time \eqn{t}
#'   (\eqn{p \times p}).
#' @param Vt Numeric matrix. Observation noise left-covariance at time \eqn{t}
#'   (\eqn{N \times N}).
#' @param mt_1 Numeric matrix. Filtered state mean at \eqn{t-1} (\eqn{p \times S}).
#' @param Mt_1 Numeric matrix. Filtered state left-covariance at \eqn{t-1}
#'   (\eqn{p \times p}).
#' @param nt_1 Numeric scalar. Shape parameter of \eqn{\sigma^2 | D_{t-1}}.
#' @param Dt_1 Numeric scalar. Rate parameter of \eqn{\sigma^2 | D_{t-1}}.
#' @param Rinv Numeric matrix. Precomputed inverse of the spatial correlation
#'   matrix \eqn{R} (\eqn{S \times S}).
#' @param delta Numeric scalar. Right-variance discount factor. Defaults to
#'   \code{1.0}.
#'
#' @return A named list with updated filtering parameters:
#'   \describe{
#'     \item{nt}{Updated shape parameter.}
#'     \item{Dt}{Updated rate parameter.}
#'     \item{at}{One-step-ahead state mean (\eqn{p \times S}).}
#'     \item{At}{One-step-ahead state left-covariance (\eqn{p \times p}).}
#'     \item{mt}{Filtered state mean (\eqn{p \times S}).}
#'     \item{Mt}{Filtered state left-covariance (\eqn{p \times p}).}
#'     \item{delta}{The discount factor passed in.}
#'   }
#'
#' @seealso \code{\link{FF_sigma2R}}
#' @export
FF_1step_R_sigma2R <- function(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, nt_1, Dt_1, Rinv, delta = 1.0){
  N <- dim(Yt)[1]
  S <- dim(Yt)[2]
  at <- Gt %*% mt_1
  At <- tcrossprod(Gt %*% Mt_1, Gt) + Wt
  qt <- Ft %*% at
  FtAt <- Ft %*% At
  Qt <- tcrossprod(FtAt, Ft) + Vt
  Qtinv <- solve(Qt)
  AttFtQtinv <- tcrossprod(At, Ft)%*% Qtinv
  Yt_qt <- Yt - qt
  mt <- at + AttFtQtinv %*% Yt_qt
  Mt <- At - AttFtQtinv %*% FtAt
  nt <- nt_1 + N * S / 2
  Dt <- Dt_1 + 1/2 * matrixcalc::matrix.trace(crossprod(Yt_qt, Qtinv) %*% Yt_qt %*% Rinv)

  out <- list(nt = nt, Dt = Dt, at = at, At = At,
              mt = mt, Mt = Mt, delta = delta)
  return(out)
}



#' Forward Filter for the scalar-sigma-squared-times-R model
#'
#' Runs the full forward filtering pass over \code{nT} time steps under the
#' model where the right covariance of \eqn{Y} is \eqn{\sigma^2 R} with
#' \eqn{\sigma^2 \sim IG(n_0, d_0)}.
#'
#' @param Y List of length \code{nT}. Each element is the \eqn{N \times S}
#'   data matrix at the corresponding time step.
#' @param F_ls Covariate matrix or list of matrices (see \code{\link{FF}}).
#' @param G_ls State transition matrix or list (see \code{\link{FF}}).
#' @param W_ls State noise left-covariance matrix or list (see \code{\link{FF}}).
#' @param V_ls Observation noise left-covariance matrix or list
#'   (see \code{\link{FF}}).
#' @param m0 Numeric matrix. Prior state mean (\eqn{p \times S}).
#' @param M0 Numeric matrix. Prior state left-covariance (\eqn{p \times p}).
#' @param n0 Numeric scalar. Prior shape of \eqn{\sigma^2}.
#' @param D0 Numeric scalar. Prior rate of \eqn{\sigma^2}.
#' @param nT Integer. Number of time steps.
#' @param R Numeric matrix. Fixed spatial correlation matrix (\eqn{S \times S}).
#'   Its Cholesky inverse is used internally.
#' @param delta Numeric scalar. Right-variance discount factor. Defaults to
#'   \code{1.0}.
#'
#' @return A named list of length \code{nT + 1} (same structure as
#'   \code{\link{FF}}) with an additional scalar \code{Dt} in place of a matrix.
#'
#' @seealso \code{\link{FF_1step_R_sigma2R}}, \code{\link{FFBS_sigma2R}}
#' @export
FF_sigma2R <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, R, delta = 1.0){
  if (!is.null(dim(D0)) || length(D0) != 1) {
    stop("D0 is not a number")
  }

  out <- list()
  out_name <- c()
  # initialize matirces
  Ft <- F_ls
  Gt <- G_ls
  Wt <- W_ls
  Vt <- V_ls

  # Rinv = solve(R)
  RR <- chol(R)
  Rinv <- chol2inv(RR)

  for (i in 1:nT) {
    if(i %% round(0.1*nT) == 0){
      print(paste("FF:", i, "/", nT))
      print(Sys.time())
    }
    # Check if F_ls, G_ls, W_ls, V_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(F_ls)){
      Ft <- as.matrix(F_ls[[i]])
    }
    if(is.list(G_ls)){
      Gt <- as.matrix(G_ls[[i]])
    }
    if(is.list(W_ls)){
      Wt <- as.matrix(W_ls[[i]])
    }
    if(is.list(V_ls)){
      Vt <- as.matrix(V_ls[[i]])
    }

    Yt <- Y[[i]]
    # Set T+1 values
    if (i == 1) {
      mt_1 = m0
      Mt_1 = M0
      nt_1 = n0
      Dt_1 = D0
    } else{
      mt_1 = one_step$mt
      Mt_1 = one_step$Mt
      nt_1 = one_step$nt
      Dt_1 = one_step$Dt
    }

    # perform 1 step forward
    one_step <- FF_1step_R_sigma2R(Yt = Yt, Ft = Ft, Gt = Gt,
                                   Wt = Wt, Vt = Vt,
                                   mt_1 = mt_1, Mt_1 = Mt_1,
                                   nt_1 = nt_1, Dt_1 = Dt_1,
                                   Rinv = Rinv, delta = delta)
    out[[i]] <- one_step # save results
  }

  out_name <- paste("T", seq(1:nT), sep = "") # accumulate output list names
  names(out) <- out_name
  out$prior <- list(m0 = m0, M0 = M0, n0 = n0, D0 = D0)
  return(out)
}

# I----

#' Single forward filter step (identity right-covariance)
#'
#' Performs one step of the forward filter under the model where the right
#' covariance of the state and observation is the identity matrix (no
#' inverse-Wishart update for \eqn{\Sigma}).
#'
#' @param Yt Numeric matrix. Data matrix at time \eqn{t} (\eqn{N \times S}).
#' @param Ft Numeric matrix. Covariate matrix at time \eqn{t} (\eqn{N \times p}).
#' @param Gt Numeric matrix. State transition matrix (\eqn{p \times p}).
#' @param Wt Numeric matrix. State noise left-covariance (\eqn{p \times p}).
#' @param Vt Numeric matrix. Observation noise left-covariance (\eqn{N \times N}).
#' @param mt_1 Numeric matrix. Filtered state mean at \eqn{t-1} (\eqn{p \times S}).
#' @param Mt_1 Numeric matrix. Filtered state left-covariance at \eqn{t-1}
#'   (\eqn{p \times p}).
#' @param delta Numeric scalar. Discount factor. Defaults to \code{1.0}.
#'
#' @return A named list with:
#'   \describe{
#'     \item{at}{One-step-ahead state mean.}
#'     \item{At}{One-step-ahead state left-covariance.}
#'     \item{mt}{Filtered state mean.}
#'     \item{Mt}{Filtered state left-covariance.}
#'     \item{delta}{The discount factor passed in.}
#'   }
#'
#' @seealso \code{\link{FF_I}}
#' @export
FF_1step_R_I <- function(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, delta = 1.0){
  N <- dim(Yt)[1]
  S <- dim(Yt)[2]
  at <- Gt %*% mt_1
  At <- tcrossprod(Gt %*% Mt_1, Gt) + Wt
  qt <- Ft %*% at
  FtAt <- Ft %*% At
  Qt <- tcrossprod(FtAt, Ft) + Vt
  Qtinv <- solve(Qt)
  AttFtQtinv <- tcrossprod(At, Ft)%*% Qtinv
  Yt_qt <- Yt - qt
  mt <- at + AttFtQtinv %*% Yt_qt
  Mt <- At - AttFtQtinv %*% FtAt

  out <- list(at = at, At = At,
              mt = mt, Mt = Mt, delta = delta)
  return(out)
}


#' Forward Filter with identity right-covariance
#'
#' Runs the full forward filtering pass over \code{nT} time steps under the
#' model where the right covariance is fixed at the identity (no
#' inverse-Wishart update). Calls \code{\link{FF_1step_R_I}} at each step.
#'
#' @param Y List of length \code{nT}. Each element is the data matrix at time
#'   \eqn{t}.
#' @param F_ls Covariate matrix or list of matrices (see \code{\link{FF}}).
#' @param G_ls State transition matrix or list (see \code{\link{FF}}).
#' @param W_ls State noise left-covariance matrix or list (see \code{\link{FF}}).
#' @param V_ls Observation noise left-covariance matrix or list
#'   (see \code{\link{FF}}).
#' @param m0 Numeric matrix. Prior state mean.
#' @param M0 Numeric matrix. Prior state left-covariance.
#' @param nT Integer. Number of time steps.
#' @param delta Numeric scalar. Discount factor. Defaults to \code{1.0}.
#'
#' @return A named list of length \code{nT + 1}. Elements \code{"T1"} through
#'   \code{"T<nT>"} contain \code{at}, \code{At}, \code{mt}, \code{Mt}. The
#'   element \code{prior} stores \code{m0} and \code{M0}.
#'
#' @seealso \code{\link{FF_1step_R_I}}, \code{\link{FFBS_I}}
#' @export
FF_I <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, nT, delta = 1.0){
  out <- list()
  out_name <- c()
  # initialize matirces
  Ft <- F_ls
  Gt <- G_ls
  Wt <- W_ls
  Vt <- V_ls

  for (i in 1:nT) {
    if(i %% round(0.1*nT) == 0){
      print(paste("FF:", i, "/", nT))
      print(Sys.time())
    }
    # Check if F_ls, G_ls, W_ls, V_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(F_ls)){
      Ft <- as.matrix(F_ls[[i]])
    }
    if(is.list(G_ls)){
      Gt <- as.matrix(G_ls[[i]])
    }
    if(is.list(W_ls)){
      Wt <- as.matrix(W_ls[[i]])
    }
    if(is.list(V_ls)){
      Vt <- as.matrix(V_ls[[i]])
    }

    Yt <- Y[[i]]
    # Set T+1 values
    if (i == 1) {
      mt_1 = m0
      Mt_1 = M0
    } else{
      mt_1 = one_step$mt
      Mt_1 = one_step$Mt
    }

    # perform 1 step forward
    one_step <- FF_1step_R_I(Yt = Yt, Ft = Ft, Gt = Gt,
                                   Wt = Wt, Vt = Vt,
                                   mt_1 = mt_1, Mt_1 = Mt_1,
                                   delta = delta)
    out[[i]] <- one_step # save results
  }

  out_name <- paste("T", seq(1:nT), sep = "") # accumulate output list names
  names(out) <- out_name
  out$prior <- list(m0 = m0, M0 = M0)
  return(out)
}
