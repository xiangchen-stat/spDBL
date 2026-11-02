# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# SIR reaction-diffusion training/held-out designs (spDBL::gen_pde) and
# synthetic field data for calibration. Simulations are cached on disk so
# that every method and every worker process sees identical data.

# Random Latin hypercube on the active coordinates; fixed coordinates are set
# to their (equal) bounds. `inner` shrinks the box towards its centre.
sir_lhs <- function(n, lower, upper, seed, inner = 1) {
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv))
  set.seed(seed)
  active <- which(upper > lower)
  centre <- (lower + upper) / 2
  lo <- centre - inner * (upper - lower) / 2
  hi <- centre + inner * (upper - lower) / 2
  X <- matrix(rep(lower, each = n), n, length(lower),
              dimnames = list(NULL, names(lower)))
  for (d in active) {
    u <- (sample.int(n) - stats::runif(n)) / n
    X[, d] <- lo[d] + u * (hi[d] - lo[d])
  }
  X
}

sir_simulate <- function(X, grid, n_time, n_people) {
  lapply(seq_len(nrow(X)), function(i)
    spDBL::gen_pde(X[i, ], grid, grid, n_people, n_time))
}

# Runs (list of T x S matrices) -> time-indexed list of n x S matrices.
sir_stack <- function(runs) {
  lapply(seq_len(nrow(runs[[1]])), function(t)
    do.call(rbind, lapply(runs, function(y) y[t, , drop = FALSE])))
}

# Training and held-out SIR runs on a grid x grid lattice. Training and test
# inputs are independent Latin hypercubes in the same box; the design does
# not depend on the grid, so every grid size uses the same parameter values.
sir_dataset <- function(grid, n_train, n_test, n_time, seed, box, cache) {
  file <- file.path(cache, sprintf("sir_g%03d_n%d-%d_T%d_s%d.rds",
    as.integer(grid), as.integer(n_train), as.integer(n_test),
    as.integer(n_time), as.integer(seed)))
  if (file.exists(file)) return(file)
  X_train <- sir_lhs(n_train, box$lower, box$upper, seed)
  X_test <- sir_lhs(n_test, box$lower, box$upper, seed + 1L)
  t_sim <- system.time({
    runs_train <- sir_simulate(X_train, grid, n_time, box$n_people)
    runs_test <- sir_simulate(X_test, grid, n_time, box$n_people)
  })[["elapsed"]]
  d <- list(grid = grid, n_time = n_time, input_train = X_train,
    input_test = X_test, train = sir_stack(runs_train),
    test = sir_stack(runs_test), lower = box$lower, upper = box$upper,
    n_people = box$n_people, simulation_seconds = t_sim,
    generator = "spDBL::gen_pde (log(I + 1) of the 2-D SIR PDE)")
  saveRDS(d, file)
  message(sprintf("Simulated %d SIR runs on %d^2 grid in %.1f s -> %s",
    n_train + n_test, grid, t_sim, basename(file)))
  file
}

# Calibration replicates. Simulator design: the 40 training runs bundled in
# dt_calibration (6 x 6 grid, 26 times). Truths: a Latin hypercube in the
# inner part of the parameter box, simulated with gen_pde. Field data follow
# data-raw/sir-calibration.R: tau2_t ~ IG(n0, d0), random-walk discrepancy
# with increments N(0, tau2_t U), and N(0, tau2_t) observation noise.
sir_calibration_replicates <- function(n_rep, seed, box, cache, inner = 0.8,
                                       n0 = 30, d0 = 2) {
  file <- file.path(cache, sprintf("sircal_r%d_s%d_in%g_n%g_d%g.rds",
    as.integer(n_rep), as.integer(seed), inner, n0, d0))
  if (file.exists(file)) return(file)
  utils::data("dt_calibration", package = "spDBL", envir = environment())
  d <- get("dt_calibration", environment())
  truths <- sir_lhs(n_rep, d$lower, d$upper, seed + 7L, inner = inner)
  U <- d$field$U
  L <- t(chol(U))
  S <- d$Nx * d$Ny
  T <- length(d$time)
  reps <- lapply(seq_len(n_rep), function(r) {
    y <- spDBL::gen_pde(truths[r, ], d$Nx, d$Ny, box$n_people, T)
    set.seed(seed + 1000L + r)
    tau <- 1 / stats::rgamma(T + 1L, shape = n0, rate = d0)
    u0 <- sqrt(tau[1]) * stats::rnorm(S)
    previous <- u0
    u <- z <- matrix(NA_real_, T, S)
    for (t in seq_len(T)) {
      previous <- as.numeric(previous + sqrt(tau[t + 1L]) * L %*% stats::rnorm(S))
      u[t, ] <- previous
      z[t, ] <- y[t, ] + previous + stats::rnorm(S, sd = sqrt(tau[t + 1L]))
    }
    list(eta = truths[r, ], y = y, u = u, z = z, tau2 = tau[-1], tau2_0 = tau[1])
  })
  out <- list(training = d$training, lower = d$lower, upper = d$upper,
    Nx = d$Nx, Ny = d$Ny, time = d$time, grid = d$grid, U = U,
    prior = list(n0 = n0, d0 = d0), replicates = reps,
    generator = "gen_pde truth + random-walk discrepancy (data-raw/sir-calibration.R recipe)")
  saveRDS(out, file)
  file
}
