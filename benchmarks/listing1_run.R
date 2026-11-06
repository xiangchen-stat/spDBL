# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Run Listing 1 of the MLOSS manuscript (benchmarks/listing1.R) exactly as
# printed, echoing it as at the console, and record its output, run time and
# session information:
#   Rscript benchmarks/listing1_run.R [output directory]
# Default output: benchmarks/results/listing1/. Files:
#   listing1_output.txt  the echoed listing and everything it prints
#   listing1.rds         run time, summary(fit), per-draw diagnostics, session
#   session.txt          sessionInfo()

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
bench <- getOption("spdbl_bench_dir")
out <- commandArgs(TRUE)[1]
if (is.na(out)) out <- file.path(bench, "results", "listing1")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

env <- new.env()
con <- file(file.path(out, "listing1_output.txt"), open = "wt")
sink(con, split = TRUE)
started <- Sys.time()
timing <- system.time(source(file.path(bench, "listing1.R"), local = env,
  echo = TRUE, print.eval = TRUE, max.deparse.length = Inf, spaced = FALSE,
  skip.echo = 5, keep.source = TRUE))
sink()
close(con)

cpu_model <- tryCatch(suppressWarnings(sub(".*: ", "", grep("model name",
  readLines("/proc/cpuinfo"), value = TRUE)[1])), error = function(e) NA_character_)
info <- list(
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  elapsed_s = timing[["elapsed"]],
  cpu_s = sum(timing[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE),
  summary = summary(env$fit),
  diagnostics = env$fit$diagnostics,
  draws = nrow(env$fit$mcmc_eta),
  fields_dim = lapply(env$fields, dim),
  spdbl_version = as.character(utils::packageVersion("spDBL")),
  git_commit = tryCatch({
    sha <- suppressWarnings(system2("git", c("-C", shQuote(bench), "rev-parse", "HEAD"),
                                    stdout = TRUE, stderr = FALSE))
    dirty <- length(suppressWarnings(system2("git", c("-C", shQuote(bench), "status",
      "--porcelain", "--untracked-files=no"), stdout = TRUE, stderr = FALSE))) > 0
    paste0(sha[1], if (dirty) "-dirty")
  }, error = function(e) NA_character_),
  r_version = R.version.string, cpu = cpu_model,
  sys = Sys.info()[c("sysname", "release", "machine")],
  slurm = Sys.getenv(c("SLURM_JOB_ID", "SLURM_CPUS_ON_NODE")),
  session = utils::sessionInfo())
saveRDS(info, file.path(out, "listing1.rds"))
writeLines(c(paste("Commit:", info$git_commit), "",
  utils::capture.output(print(info$session))), file.path(out, "session.txt"))
message(sprintf("Listing 1 finished in %.1f min (CPU %.1f min); results in %s",
  info$elapsed_s / 60, info$cpu_s / 60, out))
