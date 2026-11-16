# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Read one design block: the lag-1 columns `cols` and, when the file holds a
# second set of fncol columns, the matching lag-2 columns at offset fncol.
.read_design_block <- function(file, rows, cols, fncol) {
  n_file <- length(strsplit(readLines(file, n = 1L), ",", fixed = TRUE)[[1]])
  out <- as.matrix(read_big_csv_quick(file, rows = rows, cols = cols))
  if (n_file >= fncol + cols[2])
    out <- cbind(out, as.matrix(read_big_csv_quick(file, rows = rows,
                                                   cols = fncol + cols)))
  unname(out)
}

# Split a matrix whose blocks are bound by column into a list of blocks.
.split_columns <- function(x, widths) {
  ends <- cumsum(widths)
  lapply(seq_along(widths), function(i)
    x[, seq.int(ends[i] - widths[i] + 1L, ends[i]), drop = FALSE])
}

.file_for_time <- function(files, f) if (length(files) == 1L) files else files[f]

#' Out-of-core backward smoothing
#'
#' Runs the backward smoothing pass over the filtered moments that
#' [FF_bigdata_R()] wrote to disk, one time step at a time, and writes the
#' smoothed moments to CSV files. The result equals [BS()] applied to the
#' in-memory filter over the same blocks, without holding all episodes in
#' memory.
#'
#' @param ff_files The list of file paths returned by [FF_bigdata_R()].
#' @param G_ls The state transition file path(s) given to [FF_bigdata_R()]:
#'   one path, or one per time step.
#' @param nT,fnrow,fncol,bnrow,bncol As given to [FF_bigdata_R()].
#' @param path_out Directory where the output CSV files are written.
#' @param delta Discount factor. Not implemented: values other than 1 give a
#'   warning and are ignored. Defaults to \code{1}.
#' @param verbose Logical. If \code{TRUE}, reports progress.
#'
#' @return A list of file paths with one file per time step:
#'   \describe{
#'     \item{st_ls}{Smoothed state means; the \eqn{p \times C_b} blocks are
#'       bound by column in traversal order, as in the filter's `mt` files.}
#'     \item{St_ls}{Smoothed left covariances, \eqn{p \times p} blocks bound
#'       by column.}
#'   }
#'
#' @seealso [FF_bigdata_R()], [FFBS_predict_bigdata_R()], [BS()]
#' @examples
#' set.seed(1)
#' dir <- file.path(tempdir(), "bs_big")
#' dir.create(dir, showWarnings = FALSE)
#' csv <- function(x, name) {
#'   path <- file.path(dir, name)
#'   utils::write.csv(x, path, row.names = FALSE)
#'   path
#' }
#' # 5 runs x 4 sites over 3 times, in two blocks of 2 sites (AR(1) designs)
#' Y <- lapply(1:3, function(t) matrix(rnorm(20, t), 5, 4))
#' F <- c(list(matrix(1, 5, 4)), Y[1:2])
#' ff <- FF_bigdata_R(Y_ls = mapply(csv, Y, paste0("Y", 1:3, ".csv")),
#'   F_ls = mapply(csv, F, paste0("F", 1:3, ".csv")),
#'   G_ls = csv(diag(2), "G.csv"), W_ls = csv(diag(2), "W.csv"),
#'   V_ls = csv(diag(5), "V.csv"), m0 = matrix(0, 2, 2), M0 = diag(2),
#'   n0 = 4, D0 = diag(2), nT = 3, fnrow = 5, fncol = 4, bnrow = 5,
#'   bncol = 2, path_out = dir)
#' bs <- BS_bigdata_R(ff, file.path(dir, "G.csv"), nT = 3, fnrow = 5,
#'   fncol = 4, bnrow = 5, bncol = 2, path_out = dir)
#' utils::read.csv(bs$st_ls[1])  # smoothed means at time 1, both blocks
#' unlink(dir, recursive = TRUE)
#' @export
BS_bigdata_R <- function(ff_files, G_ls, nT, fnrow, fncol, bnrow, bncol,
                         path_out, delta = 1, verbose = FALSE) {
  delta <- .check_delta(delta)
  ind <- generate.grid.rowsnake(fnrow, fncol, bnrow, bncol)
  K <- nrow(ind)
  widths <- ind$R - ind$L + 1
  read <- function(path) unname(as.matrix(read_big_csv_quick(path)))
  st_files <- file.path(path_out, paste0("output_BS_st_", seq_len(nT), ".csv"))
  St_files <- file.path(path_out, paste0("output_BS_SSt_", seq_len(nT), ".csv"))
  following <- NULL  # smoothed and one-step moments of the next episode
  for (f in seq.int(nT, 1L)) {
    mt <- .split_columns(read(ff_files$mt_ls[f]), widths)
    p <- nrow(mt[[1]])
    Mt <- .split_columns(read(ff_files$Mt_ls[f]), rep(p, K))
    at <- .split_columns(read(ff_files$at_ls[f]), widths)
    At <- .split_columns(read(ff_files$At_ls[f]), rep(p, K))
    G_now <- as.matrix(read_big_csv_quick(.file_for_time(G_ls, f)))
    st <- St <- vector("list", K)
    for (i in seq.int(K, 1L)) {
      # The next episode is block i + 1 at this time, entered with this time's
      # G, or block 1 at the next time, entered with the next time's G.
      nxt <- if (i < K) list(st = st[[i + 1L]], St = St[[i + 1L]],
                             at = at[[i + 1L]], At = At[[i + 1L]], G = G_now) else following
      if (is.null(nxt)) {
        st[[i]] <- mt[[i]]
        St[[i]] <- Mt[[i]]
      } else {
        one <- BS_1step_cpp(mt = mt[[i]], Mt = Mt[[i]], st1 = nxt$st, St1 = nxt$St,
                            at1 = nxt$at, At1 = nxt$At, Gt1 = nxt$G, delta = delta)
        st[[i]] <- one$st
        St[[i]] <- one$St
      }
    }
    utils::write.csv(do.call(cbind, st), st_files[f], row.names = FALSE)
    utils::write.csv(do.call(cbind, St), St_files[f], row.names = FALSE)
    following <- list(st = st[[1]], St = St[[1]], at = at[[1]], At = At[[1]], G = G_now)
    if (verbose) message("Backward smoother ", f, "/", nT)
  }
  list(st_ls = st_files, St_ls = St_files)
}

#' Out-of-core exact prediction
#'
#' Predicts simulator output at new inputs from the smoothed coefficients
#' written by [BS_bigdata_R()], streaming one block of data at a time. For
#' each block it computes the same mean as [FFBS_predict_exact()]:
#' \eqn{F_{new}\hat\Theta + J^\top V^{-1}(Y - F\hat\Theta)}, with the
#' exponential input kernel. Every block must contain all training runs
#' (`bnrow = fnrow` in [FF_bigdata_R()]), as in [emulator_learn()].
#'
#' @param Y_ls,F_ls The response and design file paths given to
#'   [FF_bigdata_R()].
#' @param F_new_ls Design file path(s) for the new inputs, laid out like
#'   `F_ls`: one row per new input, `fncol` lag-1 columns and, for AR(2)
#'   designs, `fncol` lag-2 columns after them. One path or one per time step.
#' @param input Matrix of training inputs, one row per run (`fnrow` rows).
#' @param input_new Matrix (or vector, for one input) of new inputs.
#' @param bs_files The list returned by [BS_bigdata_R()].
#' @param nT,fnrow,fncol,bncol As given to [FF_bigdata_R()].
#' @param phi_para,gp_sigma2,gp_tau2 Range, amplitude and nugget of the
#'   exponential input kernel, as in [FFBS_predict_exact()].
#' @param path_out Directory where the prediction files are written.
#' @param verbose Logical. If \code{TRUE}, reports progress.
#'
#' @return A character vector with one CSV path per time step; each file holds
#'   the predicted means, one row per new input and one column per site.
#'
#' @seealso [FF_bigdata_R()], [BS_bigdata_R()], [FFBS_predict_exact()]
#' @examples
#' set.seed(1)
#' dir <- file.path(tempdir(), "pred_big")
#' dir.create(dir, showWarnings = FALSE)
#' csv <- function(x, name) {
#'   path <- file.path(dir, name)
#'   utils::write.csv(x, path, row.names = FALSE)
#'   path
#' }
#' # 5 training runs and 2 new inputs; 4 sites in two blocks; 3 times
#' input <- matrix(runif(10), 5, 2)
#' input_new <- matrix(runif(4), 2, 2)
#' Y <- lapply(1:3, function(t) matrix(rnorm(20, t), 5, 4))
#' F <- c(list(matrix(1, 5, 4)), Y[1:2])
#' V <- 1.1 * exp(-2 * as.matrix(dist(input))) + diag(1e-4, 5)
#' ff <- FF_bigdata_R(Y_ls = mapply(csv, Y, paste0("Y", 1:3, ".csv")),
#'   F_ls = mapply(csv, F, paste0("F", 1:3, ".csv")),
#'   G_ls = csv(diag(2), "G.csv"), W_ls = csv(diag(2), "W.csv"),
#'   V_ls = csv(V, "V.csv"), m0 = matrix(0, 2, 2), M0 = diag(2), n0 = 4,
#'   D0 = diag(2), nT = 3, fnrow = 5, fncol = 4, bnrow = 5, bncol = 2,
#'   path_out = dir)
#' bs <- BS_bigdata_R(ff, file.path(dir, "G.csv"), nT = 3, fnrow = 5,
#'   fncol = 4, bnrow = 5, bncol = 2, path_out = dir)
#' # Designs for the new inputs, e.g. built from their own lagged outputs
#' F_new <- c(list(matrix(1, 2, 4)), lapply(1:2, function(t) matrix(t, 2, 4)))
#' pred <- FFBS_predict_bigdata_R(file.path(dir, paste0("Y", 1:3, ".csv")),
#'   file.path(dir, paste0("F", 1:3, ".csv")),
#'   mapply(csv, F_new, paste0("Fnew", 1:3, ".csv")), input, input_new, bs,
#'   nT = 3, fnrow = 5, fncol = 4, bncol = 2, phi_para = 2, path_out = dir)
#' utils::read.csv(pred[3])  # predicted means at time 3
#' unlink(dir, recursive = TRUE)
#' @export
FFBS_predict_bigdata_R <- function(Y_ls, F_ls, F_new_ls, input, input_new,
    bs_files, nT, fnrow, fncol, bncol, phi_para, gp_sigma2 = 1.1,
    gp_tau2 = 1e-4, path_out, verbose = FALSE) {
  input <- as.matrix(input)
  if (is.null(dim(input_new))) input_new <- matrix(input_new, nrow = 1)
  input_new <- as.matrix(input_new)
  N <- nrow(input)
  N_new <- nrow(input_new)
  if (N != fnrow || ncol(input_new) != ncol(input))
    stop("input must have one row per training run (fnrow) and match input_new")
  if (!is.finite(phi_para) || phi_para <= 0 || gp_sigma2 <= 0 || gp_tau2 < 0)
    stop("Invalid input kernel settings")
  ind <- generate.grid.rowsnake(fnrow, fncol, fnrow, bncol)
  widths <- ind$R - ind$L + 1
  distances <- as.matrix(stats::dist(rbind(input, input_new)))
  J <- gp_sigma2 * exp(-phi_para * distances[seq_len(N), N + seq_len(N_new), drop = FALSE])
  V_chol <- chol(gp_sigma2 * exp(-phi_para * distances[seq_len(N), seq_len(N)]) +
                 diag(gp_tau2, N))
  out_files <- file.path(path_out, paste0("output_pred_", seq_len(nT), ".csv"))
  for (f in seq_len(nT)) {
    st <- .split_columns(unname(as.matrix(read_big_csv_quick(bs_files$st_ls[f]))), widths)
    pred <- matrix(NA_real_, N_new, fncol)
    for (i in seq_len(nrow(ind))) {
      cols <- c(ind$L[i], ind$R[i])
      Y <- unname(as.matrix(read_big_csv_quick(.file_for_time(Y_ls, f),
                                               rows = c(1, N), cols = cols)))
      Ft <- .read_design_block(.file_for_time(F_ls, f), c(1, N), cols, fncol)
      F_new <- .read_design_block(.file_for_time(F_new_ls, f), c(1, N_new), cols, fncol)
      if (ncol(Ft) != nrow(st[[i]]) || ncol(F_new) != nrow(st[[i]]))
        stop("Design files do not match the smoothed coefficients at time ", f)
      weights <- backsolve(V_chol, forwardsolve(t(V_chol), Y - Ft %*% st[[i]]))
      pred[, cols[1]:cols[2]] <- F_new %*% st[[i]] + crossprod(J, weights)
    }
    utils::write.csv(pred, out_files[f], row.names = FALSE)
    if (verbose) message("Prediction ", f, "/", nT)
  }
  out_files
}
