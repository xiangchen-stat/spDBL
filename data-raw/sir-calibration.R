# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Maintainer script: run from the package directory.
# Rscript data-raw/sir-calibration.R /path/to/saved/simulations
# Input files: pde_para.csv and pde_solution_1.csv,...,pde_solution_50.csv.
# They are the saved 6x6 SIR simulations from the development study.
# The CSV responses already contain log(I + 1); no transform is applied twice.
# This importer never runs a PDE solver. Original simulations were produced
# by gen_pde() using the SIR equations in JMLR (2025), Equation (26).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the directory containing saved SIR simulation CSV files")
source_dir <- normalizePath(args[1], mustWork = TRUE)
input_file <- file.path(source_dir, "pde_para.csv")
files <- file.path(source_dir, paste0("pde_solution_", seq_len(50L), ".csv"))
inputs <- as.matrix(read.csv(input_file))
simulations <- lapply(files, function(path) as.matrix(read.csv(path)))
stopifnot(identical(dim(inputs), c(50L, 5L)),
  all(vapply(simulations, function(y) identical(dim(y), c(26L, 36L)), logical(1))))
stack <- function(ids) lapply(seq_len(26L), function(t)
  do.call(rbind, lapply(simulations[ids], function(y) y[t, , drop = FALSE])))
grid <- data.frame(row = rep(seq_len(6L), times = 6L),
                   col = rep(seq_len(6L), each = 6L))
distance <- as.matrix(dist(grid))
U <- cov2cor(1.1 * exp(-3 / (0.5 * max(distance)) * distance) + diag(1e-4, 36L))
set.seed(1236L)
tau <- 1 / rgamma(27L, shape = 30, rate = 2)
u0 <- sqrt(tau[1]) * rnorm(36L)
previous <- u0
u <- z <- matrix(NA_real_, 26L, 36L)
for (t in seq_len(26L)) {
  previous <- as.numeric(previous + sqrt(tau[t + 1L]) * t(chol(U)) %*% rnorm(36L))
  u[t, ] <- previous
  z[t, ] <- simulations[[50L]][t, ] + previous + rnorm(36L, sd = sqrt(tau[t + 1L]))
}
dt_calibration <- list(
  training = list(input = inputs[1:40, , drop = FALSE], responses = stack(1:40)),
  validation = list(input = inputs[41:50, , drop = FALSE], responses = stack(41:50)),
  lower = setNames(c(2, 0.2, 0, 0, 0), colnames(inputs)),
  upper = setNames(c(4, 0.4, 0, 0.2, 0), colnames(inputs)),
  Nx = 6L, Ny = 6L, time = 0:25, grid = grid,
  field = list(z = z, locations = seq_len(36L), U = U,
    n0 = 30, d0 = 2, discount = 1, m0 = rep(0, 36L), M0 = diag(36L)),
  truth = list(input_id = 50L, eta = inputs[50L, ], y = simulations[[50L]],
    u = u, u0 = u0, tau2 = tau[-1], tau2_0 = tau[1]),
  provenance = list(response_scale = "log(I + 1)", field_seed = 1236L,
    train_ids = 1:40, validation_ids = 41:50, synthetic_input_id = 50L,
    source_md5 = setNames(unname(tools::md5sum(c(input_file, files))),
                         basename(c(input_file, files)))))
save(dt_calibration, file = "data/dt_calibration.rda", compress = "xz", version = 2)
