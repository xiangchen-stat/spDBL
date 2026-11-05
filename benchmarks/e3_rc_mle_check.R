# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Check of RobustCalibration's maximum-likelihood fit on one E3 replicate:
# where the estimate is stored and how long a multi-start fit takes. Used to
# start its MCMC chains at its own estimate in the mode-start E3 profile.
#   Rscript benchmarks/e3_rc_mle_check.R --rep=11 --starts=10

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))
args <- bench_parse_args()
rep_id <- if (is.null(args$rep)) 11L else as.integer(args$rep)
starts <- if (is.null(args$starts)) 10L else as.integer(args$starts)
bench_setup_lib(Sys.getenv("SPDBL_BENCH_LIB"))
cfg <- bench_config("e3", list(profile = "paper"))
suppressPackageStartupMessages({library(RobustCalibration); library(RobustGaSP)})
d <- readRDS(sir_calibration_replicates(cfg$n_rep, cfg$seed, cfg$box, cfg$cache,
  inner = cfg$inner, n0 = cfg$field_n0, d0 = cfg$field_d0))
rep <- d$replicates[[rep_id]]
active <- which(d$upper > d$lower)
S <- d$Nx * d$Ny; T <- length(d$time)
Y <- do.call(cbind, lapply(d$training$responses, as.matrix))
keep <- apply(Y, 2, function(v) max(v) - min(v) > 1e-10)
obs <- as.vector(t(rep$z))[keep]
x <- cbind(time = rep(d$time, each = S), row = rep(d$grid$row, T),
           col = rep(d$grid$col, T))[keep, , drop = FALSE]
X_sim <- d$training$input[, active, drop = FALSE]
theta_range <- cbind(d$lower[active], d$upper[active])
cat("truth:", round(rep$eta[active], 3), "\n")
for (disc in c("S-GaSP", "GaSP")) {
  set.seed(cfg$seed + rep_id)
  t0 <- proc.time()[["elapsed"]]
  fit <- rcalibration(design = x, observations = obs, p_theta = length(active),
    simul_type = 0, input_simul = X_sim, output_simul = Y[, keep],
    loc_index_emulator = seq_len(sum(keep)), theta_range = theta_range,
    discrepancy_type = disc, method = "mle", num_initial_values = starts)
  cat(sprintf("\n%s mle, %d starts: %.0f s\n", disc, starts, proc.time()[["elapsed"]] - t0))
  cat("slots:", paste(methods::slotNames(fit), collapse = ", "), "\n")
  for (sl in intersect(c("param_est", "opt_value", "initial_values", "post_sample"), methods::slotNames(fit))) {
    cat("--", sl, "\n"); v <- methods::slot(fit, sl); print(utils::head(v, 12))
  }
}
