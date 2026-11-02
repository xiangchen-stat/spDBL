# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Rscript data-raw/record-calibration.R /path/to/completed/oneiter500/profile
# Run from the package directory. Keep the full fit and field arrays in the
# reproduction output; distribute only samples and plot quantiles.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a completed oneiter500 output directory")
fit <- readRDS(file.path(args[1], "calibration.rds"))
predictions <- readRDS(file.path(args[1], "predictions.rds"))
stopifnot(inherits(fit, "spDBL_calibration"), isTRUE(fit$execution_only),
  fit$settings$chains == 1L, fit$settings$iter == 1L,
  fit$settings$warmup == 0L, fit$settings$draws_per_emulator == 1L,
  length(fit$settings$draws) == 500L,
  nrow(fit$mcmc_eta) == 500L, length(fit$field$locations) == 36L,
  identical(fit$settings$emulator_sampling, "marginal"),
  identical(fit$emulator$posterior$sampling_method, "marginal"),
  isFALSE(fit$emulator$posterior$joint_trajectories),
  all(table(fit$emulator_draw) == 1L),
  identical(sort(unique(fit$emulator_draw)), as.integer(fit$settings$draws)),
  all(fit$chain == 1L), all(fit$iteration == 1L),
  all(is.finite(fit$mcmc_eta)), all(is.finite(fit$mcmc_tau2)))
expected_dims <- c(500L, nrow(fit$field$z), ncol(fit$field$z))
stopifnot(all(vapply(predictions, function(x)
  identical(dim(x), expected_dims) && all(is.finite(x)), logical(1))))
diagnostics <- do.call(rbind, lapply(fit$diagnostics, function(x)
  data.frame(emulator_draw = x$emulator_draw, divergences = x$divergences,
    max_treedepth_hits = x$max_treedepth_hits, max_rhat = x$max_rhat,
    min_n_eff = x$min_n_eff)))
stopifnot(nrow(diagnostics) == 500L, all(is.na(diagnostics$max_rhat)),
          all(is.na(diagnostics$min_n_eff)))
recorded <- list(
  samples = fit[c("mcmc_eta", "mcmc_tau2", "tau2_0", "emulator_draw",
                  "chain", "iteration")],
  field_quantiles = lapply(predictions, function(x)
    apply(x, c(2, 3), stats::quantile, probs = c(0.025, 0.5, 0.975))),
  diagnostics = diagnostics, settings = fit$settings,
  locations = fit$field$locations,
  provenance = list(profile = "oneiter500", execution_only = TRUE,
    emulator_sampling_method = fit$settings$emulator_sampling,
    joint_trajectories = fit$emulator$posterior$joint_trajectories,
    design_uncertainty = isTRUE(fit$emulator$design_uncertainty),
    generated = as.character(Sys.Date()),
    package_version = as.character(utils::packageVersion("spDBL")),
    seed = fit$settings$seed, prediction_seed = 1234L,
    response_scale = "log(I + 1)", emulator_draws = 500L)
)
saveRDS(recorded, "inst/extdata/calibration-example.rds", compress = "xz", version = 2)
cat("Recorded 500 one-transition samples and field quantiles\n")
