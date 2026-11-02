# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E2: emulator accuracy on held-out SIR runs, spDBL vs RobustGaSP::ppgasp
# (and an mlegp PCA-GP baseline). All methods receive the same training and
# held-out designs (cached simulations shared with E1).
#
#   Rscript benchmarks/e2_accuracy.R --profile=validation
#   Rscript benchmarks/e2_accuracy.R --profile=paper --grids=10,25,50
#   Rscript benchmarks/e2_accuracy.R --methods=ppgasp_default,spdbl_s2r_designgp
#
# Each (grid, method) pair runs in its own Rscript process under
# /usr/bin/time. Results: results/e2/<profile>/e2_results.csv.

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "sir_data.R"))

args <- bench_parse_args()
bench_setup_lib(if (!is.null(args$lib)) args$lib else Sys.getenv("SPDBL_BENCH_LIB"))
cfg <- bench_config("e2", args)
bench_set_threads(cfg$threads)
suppressPackageStartupMessages(library(spDBL))
bench_write_meta(cfg)

tasks <- list()
for (g in cfg$grids) {
  if (g %% cfg$episode_window != 0) stop("episode_window must divide every grid size")
  data_file <- sir_dataset(g, cfg$n_train, cfg$n_test, cfg$n_time, cfg$seed,
                           cfg$box, cfg$cache)
  for (m in cfg$methods)
    tasks[[length(tasks) + 1]] <- list(id = sprintf("e2_g%03d_%s", g, m),
      grid = g, method = m, data_file = data_file, seed = cfg$seed)
}
results <- bench_run_tasks(tasks, "e2_worker.R", cfg)
tab <- bench_write_results(results, cfg, "e2_results")
print(tab[, intersect(c("grid", "method", "status", "rmse", "coverage95",
  "width95", "crps", "fit_s", "predict_s", "peak_rss_mb"), names(tab))],
  row.names = FALSE, digits = 3)
