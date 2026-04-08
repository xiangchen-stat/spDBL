#' Example emulation dataset
#'
#' A list of example training and testing objects used to illustrate the
#' emulation workflow in the package.
#'
#' @format A list with 4 elements:
#' \describe{
#'   \item{pde_para_train}{Training set of PDE parameter values.}
#'   \item{pde_para_test}{Test set of PDE parameter values.}
#'   \item{dt_pde_train}{Training set of PDE-based outputs corresponding to
#'   \code{pde_para_train}.}
#'   \item{dt_pde_test}{Test set of PDE-based outputs corresponding to
#'   \code{pde_para_test}.}
#' }
#' @source Simulated data / package example data
"dt_emulation"
