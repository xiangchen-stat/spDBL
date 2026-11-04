# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Shared helpers for the spDBL benchmark harness: command-line parsing,
# library and environment setup, provenance capture, peak-memory measurement
# of worker processes, and scoring rules. Base R only (ps is optional).

bench_dir <- function() {
  if (!is.null(getOption("spdbl_bench_dir"))) return(getOption("spdbl_bench_dir"))
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) stop("Run benchmark scripts with Rscript")
  dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  options(spdbl_bench_dir = dir)
  dir
}

# --key=value pairs; comma-separated values become vectors; numbers and
# TRUE/FALSE are converted. Positional arguments are returned in $positional.
bench_parse_args <- function(args = commandArgs(TRUE)) {
  out <- list(positional = character())
  for (a in args) {
    if (!startsWith(a, "--")) {
      out$positional <- c(out$positional, a)
      next
    }
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    value <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"
    parts <- strsplit(value, ",", fixed = TRUE)[[1]]
    num <- suppressWarnings(as.numeric(parts))
    out[[key]] <- if (all(parts %in% c("TRUE", "FALSE"))) as.logical(parts) else
      if (!anyNA(num)) num else parts
  }
  out
}

# Prepend the benchmark library (env SPDBL_BENCH_LIB or --lib=) so that the
# spDBL copy under test and the competitors are found before system packages.
bench_setup_lib <- function(lib = Sys.getenv("SPDBL_BENCH_LIB")) {
  if (nzchar(lib)) {
    if (!dir.exists(lib)) stop("Benchmark library not found: ", lib)
    Sys.setenv(SPDBL_BENCH_LIB = normalizePath(lib))
    .libPaths(c(normalizePath(lib), .libPaths()))
  }
  invisible(.libPaths())
}

bench_config <- function(experiment, args = bench_parse_args()) {
  source(file.path(bench_dir(), "R", "config.R"), local = TRUE)
  profile <- if (is.null(args$profile)) "validation" else args$profile
  cfg <- bench_profiles[[experiment]][[profile]]
  if (is.null(cfg)) stop("Unknown profile '", profile, "' for ", experiment)
  control <- c("profile", "out", "lib", "positional", "cache", "resume")
  unknown <- setdiff(names(args), c(names(cfg), control))
  if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "))
  for (k in intersect(names(args), names(cfg))) cfg[[k]] <- args[[k]]
  cfg$experiment <- experiment
  cfg$profile <- profile
  cfg$resume <- isTRUE(args$resume)
  cfg$box <- bench_sir_box
  cfg$out <- if (is.null(args$out))
    file.path(bench_dir(), "results", experiment, profile) else args$out
  cfg$cache <- if (!is.null(args$cache)) args$cache else {
    env <- Sys.getenv("SPDBL_BENCH_CACHE")
    if (nzchar(env)) env else file.path(bench_dir(), "cache")
  }
  cfg$stan_cache <- Sys.getenv("SPDBL_STAN_CACHE_DIR")
  if (!nzchar(cfg$stan_cache)) cfg$stan_cache <- file.path(cfg$cache, "stan")
  for (d in c(cfg$out, file.path(cfg$out, c("raw", "logs", "tasks")), cfg$cache))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  cfg$out <- normalizePath(cfg$out)
  cfg$cache <- normalizePath(cfg$cache)
  cfg
}

# Single-threaded linear algebra unless the profile says otherwise. The
# variables are inherited by worker processes.
bench_set_threads <- function(n) {
  n <- as.character(n)
  Sys.setenv(OMP_NUM_THREADS = n, OPENBLAS_NUM_THREADS = n,
             MKL_NUM_THREADS = n, VECLIB_MAXIMUM_THREADS = n,
             RCPP_PARALLEL_NUM_THREADS = n)
}

bench_pkg_versions <- function(pkgs = c("spDBL", "RobustGaSP",
    "RobustCalibration", "mlegp", "rstan", "StanHeaders", "deSolve",
    "ReacTran", "Rcpp", "RcppEigen", "mniw", "matrixsampling")) {
  vapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)),
    error = function(e) NA_character_), character(1))
}

bench_hardware <- function() {
  sys <- Sys.info()[["sysname"]]
  run <- function(cmd, args) tryCatch(paste(system2(cmd, args, stdout = TRUE,
    stderr = FALSE), collapse = " "), error = function(e) NA_character_,
    warning = function(w) NA_character_)
  # Batch-scheduler job (Grid Engine / Slurm), empty when run interactively.
  scheduler <- Sys.getenv(c("JOB_ID", "NSLOTS", "QUEUE", "SLURM_JOB_ID",
                            "SLURM_CPUS_ON_NODE"))
  if (sys == "Darwin") {
    list(cpu = run("sysctl", c("-n", "machdep.cpu.brand_string")),
         cores = run("sysctl", c("-n", "hw.ncpu")),
         memory_bytes = run("sysctl", c("-n", "hw.memsize")),
         load = run("uptime", character()), scheduler = scheduler)
  } else {
    list(cpu = run("sh", c("-c", shQuote("grep -m1 'model name' /proc/cpuinfo"))),
         cores = as.character(parallel::detectCores()),
         memory_bytes = run("sh", c("-c", shQuote("grep MemTotal /proc/meminfo"))),
         load = run("uptime", character()), scheduler = scheduler)
  }
}

# Provenance for one run: session, platform, BLAS/LAPACK, versions, config.
bench_meta <- function(cfg = NULL) {
  spdbl_path <- tryCatch(find.package("spDBL"), error = function(e) NA_character_)
  list(time = format(Sys.time(), tz = "UTC", usetz = TRUE),
    label = if (!is.null(cfg)) cfg$label else NA,
    config = cfg,
    session_info = utils::capture.output(print(utils::sessionInfo())),
    sys_info = Sys.info(),
    r_version = R.version.string,
    blas = unname(extSoftVersion()["BLAS"]),
    la_library = La_library(),
    la_version = La_version(),
    packages = bench_pkg_versions(),
    spdbl_path = spdbl_path,
    lib_paths = .libPaths(),
    threads_env = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
      "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")),
    hardware = bench_hardware())
}

bench_write_meta <- function(cfg, name = "meta") {
  meta <- bench_meta(cfg)
  saveRDS(meta, file.path(cfg$out, paste0(name, ".rds")))
  txt <- c(paste("Label:", meta$label), paste("Time:", meta$time),
    paste("Profile:", cfg$profile), "", "Config:",
    utils::capture.output(utils::str(cfg[setdiff(names(cfg), "box")])), "",
    "Hardware:", utils::capture.output(utils::str(meta$hardware)), "",
    paste("BLAS:", meta$blas), paste("LAPACK:", meta$la_library, meta$la_version),
    "", "Packages:", paste(names(meta$packages), meta$packages), "",
    paste("spDBL path:", meta$spdbl_path), "", meta$session_info)
  writeLines(txt, file.path(cfg$out, paste0(name, ".txt")))
  invisible(meta)
}

# Wall-clock and CPU seconds of an expression evaluated in the caller.
timed <- function(expr) {
  start <- proc.time()
  value <- expr
  d <- proc.time() - start
  list(value = value, elapsed = d[["elapsed"]],
       cpu = d[["user.self"]] + d[["sys.self"]] +
         sum(d[c("user.child", "sys.child")], na.rm = TRUE))
}

# Resident set size of this process in MB (NA without the ps package).
rss_now_mb <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  tryCatch(ps::ps_memory_info(ps::ps_handle())[["rss"]] / 2^20,
           error = function(e) NA_real_)
}

# R-heap high-water mark (MB) since the last gc(reset = TRUE). This counts
# R-allocated vectors only; Eigen/C++ temporaries appear only in peak RSS.
gc_peak_mb <- function() {
  g <- gc()
  sum(g[, ncol(g)])  # last column is "max used (Mb)"
}

# --- Worker processes -------------------------------------------------------

# GNU or BSD time binary: SPDBL_TIME if set (e.g. a GNU time built on a
# cluster without /usr/bin/time), else /usr/bin/time.
bench_time_command <- function() {
  bin <- Sys.getenv("SPDBL_TIME", "/usr/bin/time")
  if (!file.exists(bin)) return(NULL)
  if (Sys.info()[["sysname"]] == "Darwin") c(bin, "-l") else c(bin, "-v")
}

# Peak RSS (MB) reported by /usr/bin/time: bytes on macOS (-l), kB on Linux (-v).
bench_parse_time_log <- function(lines) {
  mac <- grep("maximum resident set size", lines, value = TRUE)
  if (length(mac)) {
    rss <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1", utils::tail(mac, 1))) / 2^20
    foot <- grep("peak memory footprint", lines, value = TRUE)
    footprint <- if (length(foot)) as.numeric(sub("^\\s*([0-9]+).*$", "\\1",
      utils::tail(foot, 1))) / 2^20 else NA_real_
    return(c(peak_rss_mb = rss, peak_footprint_mb = footprint))
  }
  gnu <- grep("Maximum resident set size", lines, value = TRUE)
  if (length(gnu))
    return(c(peak_rss_mb = as.numeric(sub(".*:\\s*", "", utils::tail(gnu, 1))) / 1024,
             peak_footprint_mb = NA_real_))
  c(peak_rss_mb = NA_real_, peak_footprint_mb = NA_real_)
}

# Settings that may change between a run and its resumption without
# changing any result.
bench_volatile <- c("jobs", "timeout_sec", "resume", "out", "cache", "stan_cache")

# With --resume: the earlier result of this task if it succeeded under the
# same task definition and the same settings (volatile ones aside), else NULL.
# Results written before resume support lack process_wall_s and are rerun.
bench_reusable <- function(task, cfg) {
  file <- file.path(cfg$out, "raw", paste0(task$id, ".rds"))
  old <- if (file.exists(file)) tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(old) || !identical(old$row$status[1], "OK") ||
      is.null(old$row$process_wall_s)) return(NULL)
  stable <- function(x) x[setdiff(names(x), bench_volatile)]
  same <- identical(old$task[names(task)], task) &&
    identical(stable(old$meta$config), stable(cfg))
  if (same) old else NULL
}

# Run every task in its own Rscript process under /usr/bin/time. Each task is
# a list with at least $id; the worker writes raw/<id>.rds, which is then
# rewritten with the process measurements added. Returns a list of worker
# results (status FAILED/TIMEOUT rows are kept, never dropped). With
# cfg$resume, successful earlier results of identical tasks are reused.
bench_run_tasks <- function(tasks, worker, cfg) {
  rscript <- file.path(R.home("bin"), "Rscript")
  time_cmd <- bench_time_command()
  worker <- file.path(bench_dir(), worker)
  run_one <- function(task) {
    if (isTRUE(cfg$resume)) {
      old <- bench_reusable(task, cfg)
      if (!is.null(old)) {
        message(sprintf("[%s] %s: reused earlier result", format(Sys.time(), "%H:%M:%S"),
          task$id))
        return(old)
      }
    }
    task$cfg <- cfg
    task$result_file <- file.path(cfg$out, "raw", paste0(task$id, ".rds"))
    task_file <- file.path(cfg$out, "tasks", paste0(task$id, ".rds"))
    log_file <- file.path(cfg$out, "logs", paste0(task$id, ".log"))
    saveRDS(task, task_file)
    unlink(task$result_file)
    cmd <- c(time_cmd, rscript, shQuote(worker), shQuote(task_file))
    started <- Sys.time()
    status <- suppressWarnings(system2(cmd[1], cmd[-1], stdout = log_file,
      stderr = log_file, timeout = cfg$timeout_sec))
    process_wall <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    log <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character()
    mem <- bench_parse_time_log(log)
    res <- if (file.exists(task$result_file)) readRDS(task$result_file) else
      list(row = data.frame(status = if (identical(status, 124L)) "TIMEOUT" else "FAILED",
        message = paste(utils::tail(log, 5), collapse = " | ")))
    res$row$task_id <- task$id
    res$row$exit_status <- status
    res$row$process_wall_s <- process_wall
    res$row$peak_rss_mb <- mem[["peak_rss_mb"]]
    res$row$peak_footprint_mb <- mem[["peak_footprint_mb"]]
    res$row$rss_method <- if (is.null(time_cmd)) "none" else paste(time_cmd, collapse = " ")
    res$row$label <- cfg$label
    if (file.exists(task$result_file)) saveRDS(res, task$result_file)
    message(sprintf("[%s] %s: %s (%.1f s, peak RSS %.0f MB)", format(Sys.time(), "%H:%M:%S"),
      task$id, res$row$status[1], process_wall, mem[["peak_rss_mb"]]))
    res
  }
  if (cfg$jobs > 1 && .Platform$OS.type == "unix")
    parallel::mclapply(tasks, run_one, mc.cores = cfg$jobs, mc.preschedule = FALSE) else
    lapply(tasks, run_one)
}

# Body of every worker script: load the task, run fun(task), and always write
# a result file with provenance, even on error.
bench_worker_main <- function(fun) {
  args <- commandArgs(TRUE)
  task <- readRDS(args[1])
  bench_setup_lib()
  set.seed(task$seed)
  out <- tryCatch({
    res <- fun(task)
    res$row$status <- if (is.null(res$row$status)) "OK" else res$row$status
    res
  }, error = function(e) list(row = data.frame(status = "FAILED",
      message = conditionMessage(e))))
  out$row$message <- if (is.null(out$row$message)) "" else out$row$message
  out$meta <- bench_meta(task$cfg)
  out$task <- task[setdiff(names(task), "cfg")]
  saveRDS(out, task$result_file)
  if (out$row$status[1] != "OK") {
    message("Worker failed: ", out$row$message[1])
    quit(status = 1, save = "no")
  }
  invisible(out)
}

bench_bind_rows <- function(rows) {
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    for (k in setdiff(cols, names(r))) r[[k]] <- NA
    r[cols]
  }))
}

bench_write_results <- function(results, cfg, name) {
  tab <- bench_bind_rows(lapply(results, `[[`, "row"))
  file <- file.path(cfg$out, paste0(name, ".csv"))
  utils::write.csv(tab, file, row.names = FALSE)
  message("Wrote ", file)
  invisible(tab)
}

# --- Scoring rules ---------------------------------------------------------

# Closed-form CRPS of a Gaussian predictive N(mu, sd^2) at y (Gneiting and
# Raftery 2007, Eq. 21). A degenerate predictive (sd = 0) gives |y - mu|.
crps_gaussian <- function(y, mu, sd) {
  out <- abs(y - mu)
  ok <- sd > 0
  z <- (y[ok] - mu[ok]) / sd[ok]
  out[ok] <- sd[ok] * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
  out
}

# Pointwise scores over the evaluation cells. Arrays share one shape.
score_predictions <- function(truth, mean, sd, lower, upper, mask = NULL) {
  if (is.null(mask)) mask <- array(TRUE, dim(truth))
  y <- truth[mask]; m <- mean[mask]; s <- sd[mask]
  lo <- lower[mask]; hi <- upper[mask]
  data.frame(n_cells = length(y), rmse = sqrt(mean((y - m)^2)),
    mae = mean(abs(y - m)), coverage95 = mean(y >= lo & y <= hi),
    width95 = mean(hi - lo), crps = mean(crps_gaussian(y, m, s)),
    mean_sd = mean(s))
}

# Time-indexed list of n x S matrices -> array [n, S, T].
list_to_array <- function(x) {
  array(unlist(lapply(x, as.matrix)), c(nrow(as.matrix(x[[1]])),
    ncol(as.matrix(x[[1]])), length(x)))
}
