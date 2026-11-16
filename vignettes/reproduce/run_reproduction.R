#!/usr/bin/env Rscript
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Rscript spDBL/vignettes/reproduce/run_reproduction.R [emulation|calibration|all] [full|quick] [--dry-run]
args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args[startsWith(args, "--file=")][1])
directory <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(directory, "reproduction_helpers.R"))
root <- reproduction_setup(reproduction_root(directory))
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("emulation", "calibration", "all")) else "all"
profile <- if (length(args) >= 2L) match.arg(args[2], c("full", "quick")) else "full"
dry_run <- length(args) >= 3L && identical(args[3], "--dry-run")
if (length(args) > 3L || (length(args) == 3L && !dry_run)) stop("Unknown arguments")
if (dry_run) {
  if (mode != "emulation") {
    print(calibration_reproduction_settings(profile))
    cat("Calibration output:", file.path(directory, "_results",
      if (profile == "full") "calibration-full" else "calibration"), "\n")
    cat("Report:", if (profile == "full") "calibration-full.html" else "calibration.html", "\n")
  }
  if (mode != "calibration") cat("Emulation: 100 x 100 grid, 51 times, 10 draws; existing output directory.\n")
  quit(save = "no", status = 0L)
}
for (experiment in if (mode == "all") c("emulation", "calibration") else mode) {
  if (experiment == "emulation") run_emulation_reproduction(root) else run_calibration_reproduction(root, profile)
  input <- file.path(directory, paste0(experiment, ".Rmd"))
  if (experiment == "calibration") {
    rmarkdown::render(input, quiet = TRUE, params = list(profile = profile),
      output_file = if (profile == "full") "calibration-full.html" else "calibration.html",
      envir = new.env(parent = globalenv()))
  } else rmarkdown::render(input, quiet = TRUE, envir = new.env(parent = globalenv()))
}
