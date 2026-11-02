# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E3: calibration of the 6 x 6 SIR problem, spDBL (RStan) vs
# RobustCalibration, over replicate field datasets at held-out eta.
#
#   Rscript benchmarks/e3_calibration.R --profile=validation
#   Rscript benchmarks/e3_calibration.R --profile=paper --jobs=4
#   Rscript benchmarks/e3_calibration.R --n_rep=2 --methods=spdbl
#
# Environment: SPDBL_BENCH_LIB, SPDBL_STAN_CACHE_DIR (compiled Stan model),
# SPDBL_BENCH_CACHE. Each (replicate, method) runs in its own Rscript process.
# Results: results/e3/<profile>/e3_results.csv (one row per parameter).

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))

args <- bench_parse_args()
bench_setup_lib(if (!is.null(args$lib)) args$lib else Sys.getenv("SPDBL_BENCH_LIB"))
cfg <- bench_config("e3", args)
bench_set_threads(cfg$threads)
suppressPackageStartupMessages(library(spDBL))
bench_write_meta(cfg)

data_file <- sir_calibration_replicates(cfg$n_rep, cfg$seed, cfg$box, cfg$cache,
  inner = cfg$inner, n0 = cfg$field_n0, d0 = cfg$field_d0)

# Compile (or load) the Stan model once, outside the timed workers, through
# the public API: a two-iteration fit on the bundled example.
if ("spdbl" %in% cfg$methods) {
  dir.create(cfg$stan_cache, recursive = TRUE, showWarnings = FALSE)
  data("dt_calibration", package = "spDBL")
  dc <- dt_calibration
  warm <- system.time(suppressWarnings(utils::capture.output({
    em <- emulator_learn(dc$training$input, dc$training$responses, dc$Nx, dc$Ny,
      covariance = "sigma2R", episode_window = 2L, nsam = 1L, seed = 1L)
    calibration_fit(calibration_prepare(em, dc$lower, dc$upper), dc$field$z,
      dc$field$U, chains = 1L, iter = 2L, warmup = 1L, refresh = 0L,
      keep_discrepancy = FALSE, cache_dir = cfg$stan_cache)
  })))[["elapsed"]]
  message(sprintf("Stan model ready in %.1f s (cache: %s)", warm, cfg$stan_cache))
}

tasks <- list()
for (r in seq_len(cfg$n_rep)) for (m in cfg$methods)
  tasks[[length(tasks) + 1]] <- list(id = sprintf("e3_r%03d_%s", r, m),
    rep = r, method = m, data_file = data_file, seed = cfg$seed + r)
results <- bench_run_tasks(tasks, "e3_worker.R", cfg)
tab <- bench_write_results(results, cfg, "e3_results")

ok <- tab[tab$status == "OK", ]
if (nrow(ok)) {
  agg <- do.call(rbind, lapply(split(ok, ok$method), function(x) data.frame(
    method = x$method[1], n_rep = length(unique(x$rep)),
    rmse_rel_error = sqrt(mean(x$rel_error^2)), coverage95 = mean(x$covered),
    mean_rel_width = mean(x$rel_width), max_rhat = max(x$rhat_max),
    min_ess_bulk = min(x$ess_bulk_min),
    mean_wall_s = mean(tapply(x$wall_total_s, x$rep, `[`, 1)),
    mean_ess_per_s = mean(tapply(x$ess_per_s, x$rep, `[`, 1)))))
  utils::write.csv(agg, file.path(cfg$out, "e3_summary.csv"), row.names = FALSE)
  print(agg, row.names = FALSE, digits = 3)
}
