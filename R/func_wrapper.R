#' Prepare PDE emulator training and testing data from CSV files
#'
#' Reads a matrix of PDE input parameters and a set of PDE output CSV files,
#' splits the inputs into training and testing sets, and reshapes the simulator
#' outputs into time-indexed lists for emulator fitting and evaluation.
#'
#' The function assumes that each file in \code{file_data} contains a matrix
#' whose rows correspond to time steps and whose columns correspond to spatial
#' locations. The number of time steps is inferred from the first PDE output
#' file.
#'
#' @param file_para Character scalar. Path to a CSV file containing the PDE
#'   input parameter matrix. Rows correspond to simulator runs.
#' @param file_data Character vector. Paths to CSV files containing PDE outputs,
#'   one file per simulator run. The length should match the number of rows in
#'   \code{file_para}.
#' @param prop_train Numeric scalar in \eqn{(0, 1)}. Proportion of simulator
#'   runs used for training. Defaults to \code{0.8}.
#'
#' @return A named list with:
#'   \describe{
#'     \item{pde_para_train}{Numeric matrix of training input parameters.}
#'     \item{pde_para_test}{Numeric matrix of testing input parameters.}
#'     \item{dt_pde_train}{Named list of training PDE outputs by time step.}
#'     \item{dt_pde_test}{Named list of testing PDE outputs by time step.}
#'     \item{n_train}{Integer. Number of training simulator runs.}
#'     \item{n_test}{Integer. Number of testing simulator runs.}
#'   }
#'
#' @seealso \code{\link{emulator_learn}}, \code{\link{emulator_predict}}
#' @importFrom utils read.csv
#' @export
prepare_data <- function(file_para, file_data, prop_train = 0.8){
  pde_para <- as.matrix(read.csv(file = file_para))
  n_input <- dim(pde_para)[1]
  n_train <- round(n_input * prop_train)
  n_test <- n_input - n_train
  pde_para_train <- pde_para[1:n_train,]
  pde_para_test <- pde_para[(n_train + 1):n_input,]
  
  # manipulate new data
  dat_pde <- list()
  for (i in 1:n_input) {
    dat_pde[[i]] <- as.matrix(read.csv(file = file_data[i]))
  }
  nT <- nrow(dat_pde[[1]])
  
  # transform the pde train data
  # read in transformed train data
  dt_pde_train <- list()
  for (i in 1:nT) {
    temp <- c()
    for (j in 1:n_train) {
      temp <- rbind(temp, dat_pde[[j]][i,])
    }
    dt_pde_train[[i]] <- temp
  }
  names(dt_pde_train) <- paste0("T", seq(1:nT))
  
  # transform the pde test data
  # read in transformed test data
  dt_pde_test <- list()
  for (i in 1:nT) {
    temp <- c()
    for (j in (n_train + 1):n_input) {
      temp <- rbind(temp, dat_pde[[j]][i,])
    }
    dt_pde_test[[i]] <- temp
  }
  names(dt_pde_test) <- paste0("T", seq(1:nT))
  
  return(list(
    pde_para_train = pde_para_train,
    pde_para_test = pde_para_test,
    dt_pde_train = dt_pde_train,
    dt_pde_test = dt_pde_test,
    n_train = n_train,
    n_test = n_test
  ))
}


#' Fit an FFBS-based dynamic emulator
#'
#' Fits an episode-based emulator to saved simulator runs. The default
#' inverse-Wishart model is retained for compatibility. Use
#' covariance = "sigma2R" for scalar-variance coefficient samples that can
#' be passed to [calibration_prepare()].
#'
#' @param pde_para_train Finite numeric matrix of training parameter inputs.
#' @param dt_pde_train Time-indexed list of training response matrices, with
#'   simulator runs in rows and original-grid locations in columns.
#' @param Nx,Ny Positive integer grid dimensions.
#' @param F_ls_train "default" or a list of episode design matrices.
#' @param F_ls_test Unused compatibility argument.
#' @param N_people Population scale recorded in the emulator metadata.
#' @param nsam Positive integer number of posterior coefficient draws.
#' @param AR_choice Autoregressive order, 1 or 2.
#' @param episode_window Number of spatial strips per block; must divide Nx.
#' @param gp_tune Positive range tuning factor.
#' @param gp_sigma2,gp_tau2 Positive input/spatial kernel amplitude and nugget.
#' @param covariance "iw" (default) or "sigma2R".
#' @param R_episode,R_field Optional spatial covariance matrices for the block
#'   and full grid, respectively. Defaults use the same exponential kernel.
#' @param seed Optional local RNG seed. NULL follows the session RNG.
#' @param verbose Whether to report progress.
#' @param sampling "marginal" (default) uses smoothed episode marginals
#'   conditional on a shared covariance draw, following Eq. (9) and Algorithm 4.
#'   "joint" uses backward conditionals to preserve dependence across episodes.
#' @return An object of class spDBL_emulator, retaining para_ffbs,
#'   res_ffbs, and a complete setup list for prediction and calibration.
#' @seealso [emulator_predict()], [calibration_prepare()]
#' @export
emulator_learn <- function(pde_para_train, dt_pde_train, Nx, Ny,
    F_ls_train = "default", F_ls_test = "default", N_people = 10000,
    nsam = 10, AR_choice = 2, episode_window = 1, gp_tune = 0.5,
    gp_sigma2 = 1.1, gp_tau2 = 1e-4, covariance = c("iw", "sigma2R"),
    R_episode = NULL, R_field = NULL, seed = NULL, verbose = FALSE,
    sampling = c("marginal", "joint")) {
  covariance <- match.arg(covariance)
  sampling <- match.arg(sampling)
  restore_seed <- .cal_seed(seed)
  on.exit(restore_seed(), add = TRUE)
  pde_para_train <- as.matrix(pde_para_train)
  ints <- c(Nx, Ny, nsam, episode_window)
  if (any(!is.finite(ints)) || any(ints < 1) || any(ints != as.integer(ints)) ||
      Nx %% episode_window != 0 || !AR_choice %in% c(1, 2) ||
      any(!is.finite(c(gp_tune, gp_sigma2, gp_tau2))) ||
      min(gp_tune, gp_sigma2, gp_tau2) <= 0)
    stop("Invalid grid, blocking, sampling, or kernel settings")
  n_train <- nrow(pde_para_train)
  nT_ori <- length(dt_pde_train)
  S <- as.integer(Nx * Ny)
  if (n_train < 2L || !ncol(pde_para_train) || any(!is.finite(pde_para_train)) ||
      nT_ori < 1L || !all(vapply(dt_pde_train, function(y)
        is.matrix(y) && identical(dim(y), c(n_train, S)) &&
          all(is.finite(y)), logical(1))))
    stop("Training inputs and time-indexed response matrices have incompatible dimensions")
  C <- as.integer(Ny * episode_window)
  B <- as.integer(S / C)
  ind_sp <- data.frame(row = rep(seq_len(Ny), times = Nx),
                       col = rep(seq_len(Nx), each = Ny))
  ind <- generate.grid.rowsnake(n_train, S, n_train, C)
  block_cols <- t(vapply(seq_len(B), function(b)
    seq.int(ind[b, 4], ind[b, 5]), integer(C)))
  input_distance <- max(stats::dist(pde_para_train))
  if (input_distance <= 0) stop("At least two distinct training inputs are required")
  phi_para <- 3 / (gp_tune * input_distance)
  block_distance <- if (C > 1L) max(stats::dist(ind_sp[block_cols[1, ], ])) else 0
  phi_sp <- if (block_distance > 0) 3 / (gp_tune * block_distance) else 1
  kernel <- function(x, phi) gp_sigma2 * exp(-phi * as.matrix(stats::dist(x))) +
    diag(gp_tau2, nrow(x))
  V <- kernel(pde_para_train, phi_para)
  if (is.null(R_episode)) R_episode <- kernel(ind_sp[block_cols[1, ], , drop = FALSE], phi_sp)
  if (is.null(R_field)) R_field <- kernel(ind_sp, phi_sp)
  if (!identical(dim(R_episode), c(C, C)) || !identical(dim(R_field), c(S, S)))
    stop("Spatial covariance dimensions do not match the block and grid")
  .cal_spd(R_episode, "R_episode")
  .cal_spd(R_field, "R_field")
  if (identical(F_ls_train, "default")) {
    if (AR_choice == 2 && nT_ori < 3L)
      stop("The default AR2 initialization requires at least three observations")
    F_ls_train <- if (AR_choice == 2)
      gen_F_ls_AR2_EP(dt_pde_train, nT_ori, B, ind) else
      gen_F_ls_AR1_EP(dt_pde_train, nT_ori, B, ind)
  }
  Y <- unlist(lapply(dt_pde_train, function(y) lapply(seq_len(B),
    function(b) y[, block_cols[b, ], drop = FALSE])), recursive = FALSE)
  E <- length(Y)
  if (!is.list(F_ls_train) || length(F_ls_train) != E ||
      !all(vapply(F_ls_train, function(f)
        is.matrix(f) && nrow(f) == n_train && all(is.finite(f)), logical(1))) ||
      length(unique(vapply(F_ls_train, ncol, integer(1)))) != 1L)
    stop("Supply one compatible design matrix per episode")
  P <- ncol(F_ls_train[[1]])
  G <- W <- diag(P)
  prior <- list(m0 = matrix(1, P, C), M0 = diag(P), n0 = P + 3,
                D0 = if (covariance == "iw") diag(C) else P + C)
  args <- c(list(Y = Y, F_ls = F_ls_train, G_ls = G, W_ls = W,
    V_ls = V, nT = E, verbose = verbose), prior)
  para_ffbs <- if (covariance == "iw") do.call(FFBS, args) else
    do.call(FFBS_sigma2R, c(args, list(R = R_episode)))
  res_ffbs <- if (covariance == "iw")
    FFBS_sampling(nsam, para_ffbs, F_ls_train, G, E, verbose = verbose,
                  sampling = sampling) else
    FFBS_sampling_sigma2R(nsam, para_ffbs, F_ls_train, G, E, R_episode,
                          verbose = verbose, sampling = sampling)
  structure(list(para_ffbs = para_ffbs, res_ffbs = res_ffbs,
    setup = list(version = 1L, covariance = covariance, sampling_method = sampling,
      Nx = Nx, Ny = Ny,
      pde_para_train = pde_para_train, N_people = N_people,
      AR_choice = AR_choice, episode_window = episode_window,
      n_train = n_train, nT = E, nT_ori = nT_ori, n_b = B, N_sp = S,
      bnrow = n_train, bncol = C, fncol = S, ind_sp = ind_sp, ind = ind,
      block_cols = block_cols, nsam = nsam, Y = Y, F_ls = F_ls_train,
      G_ls = G, W_ls = W, V = V, R_episode = R_episode, R_field = R_field,
      prior = prior, phi_para = phi_para, phi_sp = phi_sp,
      gp_tune = gp_tune, gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2)),
    class = "spDBL_emulator")
}

#' Predict simulator output from a fitted emulator
#'
#' With an ordinary [emulator_learn()] result, supply either new design
#' matrices or saved simulator outputs used to construct them. With a
#' [calibration_prepare()] result, the offline GP supplies those rows and
#' dt_pde_test is unnecessary.
#'
#' @param emulator An emulator or prepared calibration-input object.
#' @param input_new Matrix or vector of new mechanistic parameter inputs.
#' @param dt_pde_test Optional saved responses at new inputs, used only to
#'   construct default AR designs for an ordinary emulator.
#' @param F_new_ls "default" or explicit episode design matrices.
#' @param MC If TRUE, return coefficient-conditioned mean samples, not full
#'   draws of the emulator field. If FALSE, return smoothed means for an
#'   ordinary emulator, or draw-averaged means for a prepared emulator.
#' @return A time-indexed list of matrices, or arrays with a final posterior
#'   draw dimension when MC = TRUE.
#' @seealso [FFBS_predict_exact()], [calibration_predict()]
#' @export
emulator_predict <- function(emulator, input_new, dt_pde_test = NULL,
                             F_new_ls = "default", MC = FALSE) {
  if (is.null(dim(input_new))) input_new <- matrix(input_new, nrow = 1)
  input_new <- as.matrix(input_new)
  if (inherits(emulator, "spDBL_calibration_input")) {
    values <- lapply(seq_len(nrow(input_new)), function(i)
      lapply(seq_len(emulator$M), function(m)
        predict_calibration(input_new[i, ], m, emulator)$mu))
    return(lapply(seq_len(emulator$T), function(t) {
      out <- array(NA_real_, c(nrow(input_new), emulator$S, emulator$M))
      for (i in seq_len(nrow(input_new))) for (m in seq_len(emulator$M))
        out[i, , m] <- values[[i]][[m]][t, ]
      if (MC) out else matrix(apply(out, c(1, 2), mean), nrow(input_new), emulator$S)
    }))
  }
  s <- emulator$setup
  nT_ori <- s$nT_ori
  if (identical(F_new_ls, "default")) {
    if (is.null(dt_pde_test))
      stop("Supply F_new_ls or dt_pde_test, or use a calibration_prepare() result")
    n_test <- nrow(input_new)
    index <- generate.grid.rowsnake(n_test, s$fncol, n_test, s$bncol)
    F_new_ls <- if (s$AR_choice == 2)
      gen_F_ls_AR2_EP(dt_pde_test, nT_ori, s$n_b, index) else
      gen_F_ls_AR1_EP(dt_pde_test, nT_ori, s$n_b, index)
  }
  args <- list(Y = s$Y, F_ls = s$F_ls, F_new_ls = F_new_ls,
    input = s$pde_para_train, input_new = input_new, nT = s$nT,
    phi_para = s$phi_para, gp_sigma2 = s$gp_sigma2, gp_tau2 = s$gp_tau2)
  if (MC) {
    episodic <- do.call(FFBS_predict_MC,
      c(args, list(nsam = s$nsam, res_ffbs = emulator$res_ffbs)))
    recover_from_EP_MC(episodic, nT_ori, s$nT, s$nsam)
  } else {
    episodic <- do.call(FFBS_predict_exact, c(args, list(para_ffbs = emulator$para_ffbs)))
    recover_from_EP_exact(episodic, nT_ori, s$nT, block_cols = s$block_cols)
  }
}
