# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

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
#' @importFrom matrixsampling rinvwishart
#' @importFrom rlang .data
NULL

# Version of spDBL that built an object, stored with fitted objects so saved
# results record the code that produced them.
.spdbl_version <- function() {
  tryCatch(as.character(utils::packageVersion("spDBL")),
           error = function(e) NA_character_)
}
