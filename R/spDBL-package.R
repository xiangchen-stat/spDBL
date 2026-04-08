#' spDBL: Dynamic Bayesian Learning for Spatiotemporal Mechanistic Models
#'
#' Provides tools for Bayesian learning of spatiotemporal dynamical mechanistic
#' models, including methods for parameter estimation, simulation, and inference
#' using hierarchical and state-space modelling approaches. Core computational
#' routines (forward filter, backward sampler, matrix-normal sampling) are
#' implemented in C++ via \pkg{RcppEigen}.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib spDBL, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dist median quantile rnorm
#' @importFrom utils write.csv
#' @importFrom readr read_csv
#' @importFrom ggplot2 aes element_text geom_raster ggplot ggsave labs
#'   scale_fill_gradientn theme unit
#' @importFrom ggpubr ggarrange
#' @importFrom matrixsampling rinvwishart
NULL

