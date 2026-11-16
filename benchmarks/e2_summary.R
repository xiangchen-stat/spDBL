# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Summary of a replicated E2 run (a profile with `designs`, e.g. "fixed"):
# for each grid and method, the mean over designs of RMSE, CRPS, 95% coverage
# and width, and paired differences from spDBL's simulator-free prediction
# with t intervals over designs.
#   Rscript benchmarks/e2_summary.R [results directory]
# Two time ranges: t >= eval_from (as in E2), and eval_from <= t <= T - 1, the
# common range that includes the design-GP ablation (no prediction at T).
# Writes e2_summary.csv and e2_paired.csv next to e2_results.csv.

args <- commandArgs(TRUE)
local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
dir <- if (length(args)) args[1] else
  file.path(getOption("spdbl_bench_dir"), "results", "e2", "fixed")
tab <- utils::read.csv(file.path(dir, "e2_results.csv"))
tab <- tab[tab$status == "OK", ]
eval_from <- unique(tab$eval_from)
reference <- "spdbl_s2r_designgp"

# Per-design scores over a time range, from the per-time scores in raw/.
per_design <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) {
  r <- tab[i, ]
  pt <- readRDS(file.path(dir, "raw", paste0(r$task_id, ".rds")))$per_time
  T <- max(pt$time)
  common <- pt[pt$time >= eval_from & pt$time <= T - 1, ]
  data.frame(grid = r$grid, design = r$design, method = r$method,
    rmse = r$rmse, crps = r$crps, coverage95 = r$coverage95, width95 = r$width95,
    rmse_common = sqrt(mean(common$rmse^2)), crps_common = mean(common$crps),
    coverage_common = mean(common$coverage95), width_common = mean(common$width95))
}))
metrics <- c("rmse", "crps", "coverage95", "width95", "rmse_common", "crps_common",
             "coverage_common", "width_common")

summary <- do.call(rbind, lapply(split(per_design, list(per_design$grid, per_design$method), drop = TRUE),
  function(x) cbind(data.frame(grid = x$grid[1], method = x$method[1], designs = nrow(x)),
    as.data.frame(lapply(x[metrics], function(v) mean(v, na.rm = TRUE))))))
summary <- summary[order(summary$grid, summary$method), ]

paired <- do.call(rbind, lapply(split(per_design, per_design$grid), function(x) {
  ref <- x[x$method == reference, ]
  do.call(rbind, lapply(setdiff(unique(x$method), reference), function(m) {
    o <- merge(ref, x[x$method == m, ], by = "design", suffixes = c("_ref", "_m"))
    do.call(rbind, lapply(metrics, function(k) {
      dlt <- o[[paste0(k, "_m")]] - o[[paste0(k, "_ref")]]
      dlt <- dlt[is.finite(dlt)]
      if (length(dlt) < 2) return(NULL)
      ci <- mean(dlt) + c(-1, 1) * stats::qt(0.975, length(dlt) - 1) * stats::sd(dlt) / sqrt(length(dlt))
      data.frame(grid = x$grid[1], method = m, metric = k, designs = length(dlt),
        difference = mean(dlt), lo95 = ci[1], hi95 = ci[2],
        method_better_in = sum(if (grepl("coverage", k)) abs(o[[paste0(k, "_m")]] - 0.95) <
          abs(o[[paste0(k, "_ref")]] - 0.95) else dlt < 0, na.rm = TRUE))
    }))
  }))
}))
utils::write.csv(per_design, file.path(dir, "e2_per_design.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(dir, "e2_summary.csv"), row.names = FALSE)
utils::write.csv(paired, file.path(dir, "e2_paired.csv"), row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
cat("\nPaired differences, method minus", reference, "(t intervals over designs):\n")
print(paired[paired$metric %in% c("rmse", "crps", "coverage95"), ], row.names = FALSE, digits = 3)
