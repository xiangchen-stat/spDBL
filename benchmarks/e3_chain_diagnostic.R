# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E3 chain diagnostic. For one replicate and one emulator draw of the paper
# E3 problem, run many NUTS chains from random starts and record where each
# chain starts and where it settles, to tell chains stuck near a start from
# separated posterior modes.
#   Rscript benchmarks/e3_chain_diagnostic.R --rep=11 --draw=1 --chains=24 --cores=8
# The emulator, design GP, field data, seeds and sampler settings are those of
# the paper profile (benchmarks/R/config.R). Output in
# results/e3_diagnostic/rep<r>_draw<m>_c<chains>/: chains.csv (one row per
# chain) and diagnostics.rds.

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))
args <- bench_parse_args()
own <- c("rep", "draw")
rep_id <- if (is.null(args$rep)) 11L else as.integer(args$rep)
draw <- if (is.null(args$draw)) 1L else as.integer(args$draw)
cfg_args <- args[setdiff(names(args), own)]
if (is.null(cfg_args$profile)) cfg_args$profile <- "paper"
bench_setup_lib(Sys.getenv("SPDBL_BENCH_LIB"))
cfg <- bench_config("e3", cfg_args)
suppressPackageStartupMessages(library(spDBL))
chains <- as.integer(cfg$chains); cores <- as.integer(cfg$cores)
out <- file.path(bench_dir(), "results", "e3_diagnostic",
                 sprintf("rep%03d_draw%02d_c%d", rep_id, draw, chains))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(sir_calibration_replicates(cfg$n_rep, cfg$seed, cfg$box, cfg$cache,
  inner = cfg$inner, n0 = cfg$field_n0, d0 = cfg$field_d0))
rep <- d$replicates[[rep_id]]
em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
  covariance = "sigma2R", episode_window = cfg$episode_window,
  nsam = cfg$n_draws, seed = cfg$seed)
cal <- calibration_prepare(em, d$lower, d$upper)
seed <- cfg$seed + rep_id
keep <- chains * (cfg$iter - cfg$warmup)
t0 <- proc.time()[["elapsed"]]
fit <- calibration_fit(cal, rep$z, d$U, draws = draw, chains = chains,
  iter = cfg$iter, warmup = cfg$warmup, draws_per_emulator = keep,
  n0 = d$prior$n0, d0 = d$prior$d0, seed = seed, cores = cores, refresh = 0L,
  init = cfg$init, keep_discrepancy = FALSE, cache_dir = cfg$stan_cache)
elapsed <- proc.time()[["elapsed"]] - t0
active <- which(d$upper > d$lower)
# The random starts calibration_fit() draws for its first conditional fit
# (unit scale, one row per chain), reproduced from the same seed.
set.seed(seed + 1L)
start_unit <- matrix(stats::runif(chains * length(active)), chains, length(active))
range <- (d$upper - d$lower)[active]
# Median log posterior density (R reference density of the marginal Stan
# program, on the unit scale of eta) over 50 of the chain's draws.
log_density <- function(idx) {
  idx <- idx[unique(round(seq(1, length(idx), length.out = 50)))]
  stats::median(vapply(idx, function(k) spDBL:::calibration_marginal_log_density(
    (fit$mcmc_eta[k, active] - d$lower[active]) / range,
    1 / c(fit$tau2_0[k], fit$mcmc_tau2[k, ]), cal, draw, rep$z, d$U,
    n0 = d$prior$n0, d0 = d$prior$d0), numeric(1)))
}
rows <- lapply(seq_len(chains), function(ch) {
  sel <- fit$chain == ch
  e <- fit$mcmc_eta[sel, active, drop = FALSE]
  data.frame(chain = ch,
    t(stats::setNames(start_unit[ch, ], paste0("start_", names(active)))),
    t(stats::setNames(colMeans(e), paste0("mean_", names(active)))),
    t(stats::setNames(apply(e, 2, stats::sd) / range, paste0("relsd_", names(active)))),
    mean_tau2 = mean(fit$mcmc_tau2[sel, , drop = FALSE]),
    min_start_distance_to_bound = min(pmin(start_unit[ch, ], 1 - start_unit[ch, ])),
    log_density = log_density(which(sel)))
})
tab <- do.call(rbind, rows)
utils::write.csv(tab, file.path(out, "chains.csv"), row.names = FALSE)
saveRDS(fit[c("mcmc_eta", "mcmc_tau2", "tau2_0", "chain", "iteration")], file.path(out, "draws.rds"))
saveRDS(list(diagnostics = fit$diagnostics, truth = rep$eta[active], elapsed_s = elapsed,
  settings = list(rep = rep_id, draw = draw, chains = chains, iter = cfg$iter,
    warmup = cfg$warmup, init = cfg$init, seed = seed)), file.path(out, "diagnostics.rds"))
cat(sprintf("rep %d, emulator draw %d, %d chains, %.0f s; truth %s\n", rep_id, draw, chains,
  elapsed, paste(round(rep$eta[active], 3), collapse = " ")))
print(tab, row.names = FALSE, digits = 3)
print(fit$diagnostics[[1]][c("divergences", "max_treedepth_hits", "max_rhat", "min_n_eff")])
