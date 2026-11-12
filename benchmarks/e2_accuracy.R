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

# Profiles with `designs` repeat E2 over independent training/test designs
# (design d uses seed + 1000 (d - 1)); --reps=1,2 runs only those designs.
# `domain` fixes the physical domain (see sir_pde_domain()).
replicated <- !is.null(cfg$designs)
designs <- if (replicated) seq_len(cfg$designs) else 1L
if (!is.null(cfg$reps)) {
  designs <- as.integer(cfg$reps)
  if (!replicated || anyNA(designs) || any(designs < 1 | designs > cfg$designs) ||
      anyDuplicated(designs)) stop("--reps must list distinct designs in 1..designs")
}
for (g in cfg$grids)
  if (g %% cfg$episode_window != 0) stop("episode_window must divide every grid size")
combos <- expand.grid(grid = cfg$grids, design = designs)
data_files <- parallel::mclapply(seq_len(nrow(combos)), function(k)
  sir_dataset(combos$grid[k], cfg$n_train, cfg$n_test, cfg$n_time,
    cfg$seed + 1000L * (combos$design[k] - 1L), cfg$box, cfg$cache, domain = cfg$domain),
  mc.cores = max(1L, min(as.integer(cfg$jobs), nrow(combos))))
tasks <- list()
for (k in seq_len(nrow(combos))) {
  g <- combos$grid[k]; dsg <- combos$design[k]
  for (m in cfg$methods)
    tasks[[length(tasks) + 1]] <- list(
      id = if (replicated) sprintf("e2_g%03d_d%02d_%s", g, dsg, m) else sprintf("e2_g%03d_%s", g, m),
      grid = g, design = dsg, method = m, data_file = data_files[[k]],
      seed = cfg$seed + 1000L * (dsg - 1L))
}
results <- bench_run_tasks(tasks, "e2_worker.R", cfg)
tab <- bench_write_results(results, cfg, "e2_results")
print(tab[, intersect(c("grid", "design", "method", "status", "rmse", "coverage95",
  "width95", "crps", "fit_s", "predict_s", "peak_rss_mb"), names(tab))],
  row.names = FALSE, digits = 3)
