# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Print compact Markdown tables of E1-E3 results for one profile.
#
#   Rscript benchmarks/summarise.R --profile=validation [--root=benchmarks/results]

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
args <- bench_parse_args()
profile <- if (is.null(args$profile)) "validation" else args$profile
root <- if (is.null(args$root)) file.path(bench_dir(), "results") else args$root

md <- function(df, digits = 3) {
  df[] <- lapply(df, function(x) if (is.numeric(x)) signif(x, digits) else x)
  cat("|", paste(names(df), collapse = " | "), "|\n")
  cat("|", paste(rep("---", ncol(df)), collapse = " | "), "|\n")
  for (i in seq_len(nrow(df)))
    cat("|", paste(vapply(df[i, ], as.character, ""), collapse = " | "), "|\n")
  cat("\n")
}
read <- function(exp) {
  f <- file.path(root, exp, profile, paste0(exp, "_results.csv"))
  if (file.exists(f)) utils::read.csv(f) else NULL
}

e1 <- read("e1")
if (!is.null(e1)) {
  cat("## E1 emulator scaling (", e1$label[1], ")\n\n", sep = "")
  md(e1[order(e1$grid, e1$episode_window), intersect(c("grid", "S",
    "episode_window", "status", "fit_s", "predict_s", "rfield_build_chol_s",
    "rfield_share", "object_mb", "rfield_mb", "heap_peak_fit_mb",
    "peak_rss_mb"), names(e1))])
}
e2 <- read("e2")
if (!is.null(e2)) {
  cat("## E2 emulator accuracy (", e2$label[1], ")\n\n", sep = "")
  md(e2[order(e2$grid, e2$method), intersect(c("grid", "method", "status",
    "rmse", "coverage95", "width95", "crps", "fit_s", "predict_s",
    "peak_rss_mb"), names(e2))])
}
e3 <- read("e3")
if (!is.null(e3)) {
  cat("## E3 calibration (", e3$label[1], ")\n\n", sep = "")
  ok <- e3[e3$status == "OK", ]
  md(ok[order(ok$method, ok$rep, ok$parameter), intersect(c("method", "rep",
    "parameter", "truth", "post_mean", "q025", "q975", "covered", "rhat_max",
    "ess_bulk_min", "divergences", "wall_total_s", "ess_per_s"), names(ok))])
  s <- file.path(root, "e3", profile, "e3_summary.csv")
  if (file.exists(s)) md(utils::read.csv(s))
  bad <- e3[e3$status != "OK", intersect(c("task_id", "status", "message"), names(e3))]
  if (nrow(bad)) md(bad)
}
