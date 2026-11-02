# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E1: spDBL emulator scaling with grid size and episode window.
#
#   Rscript benchmarks/e1_scaling.R --profile=validation
#   Rscript benchmarks/e1_scaling.R --profile=paper --grids=10,25,50,100
#
# Environment: SPDBL_BENCH_LIB (library holding spDBL), SPDBL_BENCH_CACHE
# (simulated-data cache). Every configuration runs in its own Rscript
# process under /usr/bin/time to record peak RSS. Results:
# results/e1/<profile>/e1_results.csv plus raw/, logs/ and meta.{rds,txt}.

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))

args <- bench_parse_args()
bench_setup_lib(if (!is.null(args$lib)) args$lib else Sys.getenv("SPDBL_BENCH_LIB"))
cfg <- bench_config("e1", args)
bench_set_threads(cfg$threads)
suppressPackageStartupMessages(library(spDBL))
bench_write_meta(cfg)

tasks <- list()
for (g in cfg$grids) {
  data_file <- sir_dataset(g, cfg$n_train, cfg$n_test, cfg$n_time, cfg$seed,
                           cfg$box, cfg$cache)
  for (w in cfg$windows) {
    if (g %% w != 0) {
      message("Skipping grid ", g, " with window ", w, " (window must divide Nx)")
      next
    }
    tasks[[length(tasks) + 1]] <- list(id = sprintf("e1_g%03d_w%02d", g, w),
      grid = g, window = w, data_file = data_file, seed = cfg$seed)
  }
}
results <- bench_run_tasks(tasks, "e1_worker.R", cfg)
tab <- bench_write_results(results, cfg, "e1_results")
print(tab[, intersect(c("grid", "episode_window", "status", "fit_s", "predict_s",
  "rfield_share", "object_mb", "peak_rss_mb"), names(tab))], row.names = FALSE)
