# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Off-grid prediction: the E3 truths simulated on a 3x finer grid of the same
# domain (sir_offgrid_replicates()), field data at the 36 coarse cell centres
# and at target points between them; spDBL calibrates on the coarse centres
# and predicts the targets with calibration_predict(coords =), scored against
# the held-out observations and the true field (e3_offgrid_worker.R).
#
#   Rscript benchmarks/e3_offgrid.R --profile=offgrid_check --reps=1
#   Rscript benchmarks/e3_offgrid.R --profile=offgrid --jobs=15
#
# Results: results/e3/<profile>/e3_offgrid_results.csv (one row per
# replicate, target and predictor) and e3_offgrid_summary.csv.

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))

args <- bench_parse_args()
bench_setup_lib(if (!is.null(args$lib)) args$lib else Sys.getenv("SPDBL_BENCH_LIB"))
cfg <- bench_config("e3", args)
if (!identical(cfg$methods, "spdbl_offgrid")) stop("Use an off-grid profile (offgrid, offgrid_check)")
bench_set_threads(cfg$threads)
suppressPackageStartupMessages(library(spDBL))
bench_write_meta(cfg)

data_file <- sir_offgrid_replicates(cfg$n_rep, cfg$seed, cfg$box, cfg$cache,
  inner = cfg$inner, n0 = cfg$field_n0, d0 = cfg$field_d0, refine = cfg$refine,
  n_targets = cfg$n_targets)

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
tasks <- lapply(reps, function(r) list(id = sprintf("offgrid_r%03d", r), rep = r,
  method = "spdbl_offgrid", data_file = data_file, seed = cfg$seed + r))
results <- bench_run_tasks(tasks, "e3_offgrid_worker.R", cfg)
tab <- bench_write_results(results, cfg, "e3_offgrid_results")

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
    mean_predict_cell_s = mean(x$predict_cell_s),
    mean_predict_bilinear_s = mean(x$predict_bilinear_s))))
  utils::write.csv(agg, file.path(cfg$out, "e3_offgrid_summary.csv"), row.names = FALSE)
  print(agg, row.names = FALSE, digits = 3)
}
