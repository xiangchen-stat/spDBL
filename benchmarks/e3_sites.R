# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# On-grid site prediction: in each E3 replicate, site_holdout grid sites are
# left out of the calibration and predicted with calibration_predict(sites =),
# then scored against the held-out observations and the noise-free field
# (e3_sites_worker.R). Uses the same 30-replicate data as E3.
#
#   Rscript benchmarks/e3_sites.R --profile=sites_check --reps=1
#   Rscript benchmarks/e3_sites.R --profile=sites --jobs=17
#
# Results: results/e3/<profile>/e3_sites_results.csv (one row per replicate,
# target and predictor) and e3_sites_summary.csv (pooled over replicates).

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))

args <- bench_parse_args()
bench_setup_lib(if (!is.null(args$lib)) args$lib else Sys.getenv("SPDBL_BENCH_LIB"))
cfg <- bench_config("e3", args)
if (!identical(cfg$methods, "spdbl_sites")) stop("Use a site-prediction profile (sites, sites_check)")
bench_set_threads(cfg$threads)
suppressPackageStartupMessages(library(spDBL))
bench_write_meta(cfg)

data_file <- sir_calibration_replicates(cfg$n_rep, cfg$seed, cfg$box, cfg$cache,
  inner = cfg$inner, n0 = cfg$field_n0, d0 = cfg$field_d0)

# Compile (or load) the Stan model once, outside the timed workers.
dir.create(cfg$stan_cache, recursive = TRUE, showWarnings = FALSE)
data("dt_calibration", package = "spDBL")
dc <- dt_calibration
invisible(suppressWarnings(utils::capture.output({
  em <- emulator_learn(dc$training$input, dc$training$responses, dc$Nx, dc$Ny,
    covariance = "sigma2R", episode_window = 2L, nsam = 1L, seed = 1L)
  calibration_fit(calibration_prepare(em, dc$lower, dc$upper), dc$field$z,
    dc$field$U, chains = 1L, iter = 2L, warmup = 1L, refresh = 0L,
    keep_discrepancy = FALSE, cache_dir = cfg$stan_cache)
})))

reps <- if (is.null(cfg$reps)) seq_len(cfg$n_rep) else as.integer(cfg$reps)
if (!length(reps) || anyNA(reps) || any(reps < 1L | reps > cfg$n_rep) ||
    anyDuplicated(reps)) stop("--reps must list distinct replicates in 1..n_rep")
tasks <- lapply(reps, function(r) list(id = sprintf("sites_r%03d", r), rep = r,
  method = "spdbl_sites", data_file = data_file, seed = cfg$seed + r))
results <- bench_run_tasks(tasks, "e3_sites_worker.R", cfg)
tab <- bench_write_results(results, cfg, "e3_sites_results")

ok <- tab[tab$status == "OK", ]
if (nrow(ok)) {
  # Pooled over replicates: every replicate has the same number of cells.
  key <- interaction(ok$target, ok$predictor, drop = TRUE, lex.order = TRUE)
  agg <- do.call(rbind, lapply(split(ok, key), function(x) data.frame(
    target = x$target[1], predictor = x$predictor[1], n_rep = nrow(x),
    rmse = sqrt(mean(x$rmse^2)), mae = mean(x$mae), crps = mean(x$crps),
    coverage95 = mean(x$coverage95), width95 = mean(x$width95),
    min_rep_coverage95 = suppressWarnings(min(x$coverage95)),
    mean_fit_s = mean(x$fit_s), mean_predict_s = mean(x$predict_s),
    mean_predict_sites_s = mean(x$predict_sites_s))))
  utils::write.csv(agg, file.path(cfg$out, "e3_sites_summary.csv"), row.names = FALSE)
  print(agg, row.names = FALSE, digits = 3)
}
