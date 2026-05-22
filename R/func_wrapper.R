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
#' Fits a dynamic linear emulator to training PDE output using Forward Filtering
#' Backward Sampling (FFBS). The wrapper constructs the input-parameter Gaussian
#' process covariance, partitions spatial output into episode blocks, builds
#' autoregressive design matrices, runs \code{\link{FFBS}}, and draws posterior
#' state samples with \code{\link{FFBS_sampling}}.
#'
#' @param pde_para_train Numeric matrix. Training input parameters, with rows
#'   corresponding to simulator runs and columns to parameters.
#' @param dt_pde_train List of length \code{nT}. Each element is a numeric
#'   matrix of training PDE outputs at one time step, with rows corresponding to
#'   simulator runs and columns to spatial locations.
#' @param Nx Integer. Number of spatial grid cells in the x-direction.
#' @param Ny Integer. Number of spatial grid cells in the y-direction.
#' @param F_ls_train Either \code{"default"} or a list of autoregressive design
#'   matrices for the training data. If \code{"default"}, the design matrices
#'   are generated from \code{dt_pde_train}.
#' @param F_ls_test Deprecated/unused argument retained for compatibility.
#' @param N_people Numeric scalar. Population size associated with the PDE
#'   simulator. Stored in the returned setup list. Defaults to \code{10000}.
#' @param nsam Integer. Number of posterior state samples drawn by
#'   \code{\link{FFBS_sampling}}. Defaults to \code{10}.
#' @param AR_choice Integer. Autoregressive design order. Use \code{1} for
#'   \code{\link{gen_F_ls_AR1_EP}} or \code{2} for
#'   \code{\link{gen_F_ls_AR2_EP}}. Defaults to \code{2}.
#' @param episode_window Integer. Number of y-direction columns per episode
#'   block. Defaults to \code{1}.
#' @param gp_tune Numeric scalar. Tuning factor used to set the exponential
#'   kernel range from the maximum pairwise input distance. Defaults to
#'   \code{0.5}.
#' @param gp_sigma2 Numeric scalar. Marginal variance parameter for the input
#'   Gaussian process kernel. Defaults to \code{1.1}.
#' @param gp_tau2 Numeric scalar. Nugget variance parameter for the input
#'   Gaussian process kernel. Defaults to \code{10^(-4)}.
#'
#' @return A named list with:
#'   \describe{
#'     \item{para_ffbs}{Output of \code{\link{FFBS}}.}
#'     \item{res_ffbs}{Posterior state samples returned by
#'       \code{\link{FFBS_sampling}}.}
#'     \item{setup}{A list of training data, blocking information, kernel
#'       parameters, generated design matrices, and other quantities needed by
#'       \code{\link{emulator_predict}}.}
#'   }
#'
#' @seealso \code{\link{emulator_predict}}, \code{\link{FFBS}},
#'   \code{\link{FFBS_sampling}}, \code{\link{gen_F_ls_AR1_EP}},
#'   \code{\link{gen_F_ls_AR2_EP}}
#' @export
emulator_learn <- function(pde_para_train,
                           # pde_para_test,
                           dt_pde_train,
                           # dt_pde_test, 
                           Nx,
                           Ny,
                           F_ls_train = "default",
                           F_ls_test = "default",
                           N_people = 10000,
                           nsam = 10,
                           AR_choice = 2,
                           episode_window = 1,
                           gp_tune = 0.5,
                           gp_sigma2 = 1.1,
                           gp_tau2 = 10^(-4)){
  n_train <- dim(pde_para_train)[1]
  # n_test <- dim(pde_para_test)[1]
  # prop_train <- n_train / (n_train + n_test)
  nT <- length(dt_pde_train)
  nT_ori <- nT
  N_sp = Nx * Ny
  
  bnrow <- n_train 
  # bnrow_test <- n_test 
  bncol <- Ny * episode_window 
  ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
  
  
  # Initialize para for FFBS
  # Get covariance matrix
  dist_para <- as.matrix(stats::dist(pde_para_train, method = "euclidean", diag = T, upper = T))
  phi_para <- 3 / (gp_tune * max(dist_para))
  V_para <- gen_exp_kernel(loc = pde_para_train, phi = phi_para, sigma2 = gp_sigma2, tau2 = gp_tau2) # exponential kernel
  
  ind_sp_EP <- ind_sp[1:bncol, ]
  dist_sp <- as.matrix(stats::dist(ind_sp_EP, method = "euclidean", diag = T, upper = T))
  phi_sp <- 3 / (gp_tune * max(dist_sp))
  
  # generate season-episode data
  V_para_full <- V_para 
  Y_full = dt_pde_train 
  # Y_test_full <- dt_pde_test 
  
  fnrow_train <- n_train 
  # fnrow_test <- n_test 
  fncol <- Nx * Ny
  N <- bnrow
  S <- bncol
  n_block <- fnrow_train / bnrow * fncol / bncol
  nT_block <- n_block * nT
  
  # generate index
  ind <- generate.grid.rowsnake(fnrow = fnrow_train, fncol = fncol, bnrow = bnrow, bncol = bncol)
  # ind_test <- generate.grid.rowsnake(fnrow = fnrow_test, fncol = fncol, bnrow = bnrow_test, bncol = bncol)
  n_b <- n_block
  # n_b_test <- dim(ind_test)[1]
  Y <- list()
  
  # set F
  if(F_ls_train == "default"){
    if(AR_choice == 1){
      F_ls_train <- gen_F_ls_AR1_EP(nT = nT_ori, n_b = n_b, ind = ind, Y = Y_full)
      # F_ls_test <- gen_F_ls_AR1_EP(nT = nT_ori, n_b = n_b, ind = ind_test, Y = Y_test_full)
    } else if(AR_choice == 2){
      F_ls_train <- gen_F_ls_AR2_EP(nT = nT_ori, n_b = n_b, ind = ind, Y = Y_full)
      # F_ls_test <- gen_F_ls_AR2_EP(nT = nT_ori, n_b = n_b, ind = ind_test, Y = Y_test_full)
    }
  }
  
  # set Y
  for (i in 1:nT) {
    temp_Y <- Y_full[[i]]
    for (j in 1:n_b) {
      Y[[(i-1)*n_b + j]] <- temp_Y[ind[j,2]:ind[j,3], ind[j,4]:ind[j,5]]
    }
  }
  
  F_ls <- F_ls_train
  V_ls <- V_para
  N <- bnrow
  S <- bncol
  nT_ori <- nT
  nT <- nT_block
  p <- dim(F_ls[[1]])[2]
  G_ls <- diag(p)
  W_ls <- diag(p)
  n0 <- p + 3
  m0 <- matrix(1, nrow = p, ncol = S)
  M0 <- diag(p)
  D0 <- diag(S)
  
  emulator <- list()
  ## Fit the emulator and predict at held-out inputs ----
  para_ffbs <- FFBS(Y = Y, F_ls = F_ls, G_ls = G_ls,
                    W_ls = W_ls, V_ls = V_ls,
                    m0 = m0, M0 = M0,
                    n0 = n0, D0 = D0,
                    nT = nT)
  
  # sampling
  res_ffbs <- FFBS_sampling(nsam = nsam, para_ffbs = para_ffbs, 
                            F_ls = F_ls, G_ls = G_ls,
                            nT = nT)
  emulator$para_ffbs <- para_ffbs
  emulator$res_ffbs <- res_ffbs
  emulator$setup <- list(
    pde_para_train = pde_para_train,
    # pde_para_test = pde_para_test,
    N_people = N_people,
    AR_choice = AR_choice,
    episode_window = episode_window,
    n_train = n_train,
    # n_test = n_test,
    # prop_train = prop_train,
    nT = nT,
    nT_ori = nT_ori,
    n_b = n_b,
    N_sp = N_sp,
    bnrow = bnrow,
    # bnrow_test = bnrow_test,
    bncol = bncol,
    fncol = fncol,
    ind_sp = ind_sp,
    nsam = nsam,
    Y = Y,
    F_ls = F_ls,
    # F_ls_test = F_ls_test,
    phi_para = phi_para,
    gp_tune = gp_tune,
    gp_sigma2 = gp_sigma2,
    gp_tau2 = gp_tau2
  )
  
  return(emulator)
}

#' Predict PDE output from a fitted FFBS emulator
#'
#' Uses a fitted emulator returned by \code{\link{emulator_learn}} to predict
#' PDE output at new input parameter values. Prediction can use either the exact
#' posterior predictive mean from \code{\link{FFBS_predict_exact}} or Monte
#' Carlo posterior predictive draws from \code{\link{FFBS_predict_MC}}.
#'
#' If \code{F_new_ls = "default"}, the function constructs the required
#' autoregressive design matrices for the new inputs from \code{dt_pde_test}
#' using the autoregressive order and episode blocking structure stored in
#' \code{emulator}.
#'
#' @param emulator Named list. Fitted emulator object returned by
#'   \code{\link{emulator_learn}}.
#' @param input_new Numeric matrix. New input parameter values at which to
#'   predict, with rows corresponding to new simulator settings.
#' @param dt_pde_test List of length \code{nT_ori}. PDE output matrices for the
#'   new inputs, used to construct default autoregressive design matrices when
#'   \code{F_new_ls = "default"}.
#' @param F_new_ls Either \code{"default"} or a list of autoregressive design
#'   matrices for the new inputs. Defaults to \code{"default"}.
#' @param MC Logical. If \code{FALSE}, return exact posterior predictive means.
#'   If \code{TRUE}, return Monte Carlo posterior predictive draws. Defaults to
#'   \code{FALSE}.
#'
#' @return If \code{MC = FALSE}, a named list of length \code{nT_ori}, where
#'   each element is a matrix of exact posterior predictive means for one
#'   original time step. If \code{MC = TRUE}, a named list of length
#'   \code{nT_ori}, where each element is an array of Monte Carlo predictive
#'   draws for one original time step.
#'
#' @seealso \code{\link{emulator_learn}}, \code{\link{FFBS_predict_exact}},
#'   \code{\link{FFBS_predict_MC}}, \code{\link{recover_from_EP_exact}},
#'   \code{\link{recover_from_EP_MC}}
#' @export
emulator_predict <- function(emulator,
                             input_new,
                             dt_pde_test,
                             F_new_ls = "default",
                             MC = FALSE){
  # Prediction----
  input_new <- input_new
  n_input_new <- dim(input_new)[1]
  if(F_new_ls == "default"){
    # set F
    AR_choice = emulator$setup$AR_choice
    nT_ori = emulator$setup$nT_ori
    n_b = emulator$setup$n_b
    n_test <- dim(input_new)[1]
    fnrow_test <- n_test 
    fncol <- emulator$setup$fncol
    bnrow_test <- n_test 
    bncol <- emulator$setup$bncol
    ind_test <- generate.grid.rowsnake(fnrow = fnrow_test, fncol = fncol, bnrow = bnrow_test, bncol = bncol)
    Y_test_full <- dt_pde_test
    if(AR_choice == 1){
      F_ls_test <- gen_F_ls_AR1_EP(nT = nT_ori, n_b = n_b, ind = ind_test, Y = Y_test_full)
      F_new_ls = F_ls_test
    } else if(AR_choice == 2){
      F_ls_test <- gen_F_ls_AR2_EP(nT = nT_ori, n_b = n_b, ind = ind_test, Y = Y_test_full)
      F_new_ls = F_ls_test
    }
  }
  
  nsam <- emulator$setup$nsam
  Y <- emulator$setup$Y
  input <- emulator$setup$pde_para_train
  F_ls <- emulator$setup$F_ls
  nT <- emulator$setup$nT
  phi_para <- emulator$setup$phi_para
  gp_sigma2 <- emulator$setup$gp_sigma2
  gp_tau2 <- emulator$setup$gp_tau2
  
  para_ffbs <- emulator$para_ffbs
  res_ffbs <- emulator$res_ffbs
  
  if(MC){
    res_pre_MC_EP <- FFBS_predict_MC(nsam = nsam, Y = Y, res_ffbs = res_ffbs,
                                     input = input, input_new = input_new,
                                     F_ls = F_ls, F_new_ls = F_new_ls,
                                     nT = nT, phi_para = phi_para, gp_sigma2 = gp_sigma2,
                                     gp_tau2 = gp_tau2)
    ### transfer back ----
    # from episode season to season only
    res_pre_MC <- recover_from_EP_MC(dat_EP = res_pre_MC_EP, nT_ori = nT_ori, nT = nT, nsam = nsam)
    return(res_pre_MC)
  } else{
    res_pre_exact_EP <- FFBS_predict_exact(Y = Y, para_ffbs = para_ffbs,
                                           input = input, input_new = input_new,
                                           F_ls = F_ls, F_new_ls = F_new_ls,
                                           nT = nT, phi_para = phi_para, gp_sigma2 = gp_sigma2,
                                           gp_tau2 = gp_tau2)
    ### transfer back ----
    # from episode season to season only
    res_pre_exact <- recover_from_EP_exact(dat_EP = res_pre_exact_EP, nT_ori = nT_ori, nT = nT)
    return(res_pre_exact)
  }
}
