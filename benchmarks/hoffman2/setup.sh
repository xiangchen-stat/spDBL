#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# One-time setup on Hoffman2: R packages into $SPDBL_BENCH_LIB, spDBL from
# this clone, GNU time if the node lacks it, and the compiled Stan model.
# Takes 30-60 min (rstan and StanHeaders compile from source). Runs inside
# the prepare job (bash benchmarks/hoffman2/submit.sh prepare), or by hand
# from the repository root in an interactive session, not on a login node:
#
#   srun -t 4:00:00 -c 8 --mem=32G --pty bash -l
#   cd ~/spDBL && bash benchmarks/hoffman2/setup.sh
#
# Safe to rerun: installed packages are kept and spDBL is reinstalled.
set -eo pipefail
. benchmarks/hoffman2/env.sh
export NCPU=$SPDBL_CPUS
TMP=$(mktemp -d "${TMPDIR:-/tmp}/spdbl-setup.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

echo "== R packages"
Rscript - <<'EOF'
lib <- Sys.getenv("SPDBL_BENCH_LIB")
d <- read.dcf("DESCRIPTION", c("Imports", "LinkingTo"))
deps <- trimws(sub("\\(.*", "", unlist(strsplit(d[!is.na(d)], ","))))
base <- rownames(installed.packages(priority = "base"))
need <- unique(c(setdiff(deps[nzchar(deps)], base), "rstan", "StanHeaders",
  "BH", "RobustGaSP", "RobustCalibration", "mlegp", "ps"))
todo <- setdiff(need, rownames(installed.packages(lib.loc = lib)))
if (length(todo)) install.packages(todo, lib = lib,
  repos = "https://cloud.r-project.org", Ncpus = as.integer(Sys.getenv("NCPU")))
missing <- setdiff(need, rownames(installed.packages()))
if (length(missing)) stop("Not installed: ", paste(missing, collapse = ", "))
stopifnot(packageVersion("rstan") >= "2.32.7")
EOF

echo "== spDBL from this clone"
rsync -a --exclude .git --exclude benchmarks --exclude paper ./ "$TMP/spDBL/"
R CMD INSTALL --no-docs -l "$SPDBL_BENCH_LIB" "$TMP/spDBL"

echo "== GNU time (peak memory per task)"
if ! "$SPDBL_TIME" -v true 2>&1 | grep -q "Maximum resident set size"; then
  echo "No GNU time at $SPDBL_TIME; building GNU time 1.9 into $SPDBL_WORK/bin"
  (cd "$TMP" && curl -fsSLO https://ftp.gnu.org/gnu/time/time-1.9.tar.gz &&
     tar xzf time-1.9.tar.gz && cd time-1.9 &&
     ./configure --prefix="$SPDBL_WORK" >/dev/null && make >/dev/null &&
     make install >/dev/null)
  export SPDBL_TIME=$SPDBL_WORK/bin/time
fi
"$SPDBL_TIME" -v true 2>&1 | grep "Maximum resident set size"

echo "== Stan model (compiled once, cached in $SPDBL_STAN_CACHE_DIR)"
Rscript - <<'EOF'
suppressPackageStartupMessages(library(spDBL))
stopifnot(packageVersion("spDBL") >= "1.2.0")
data("dt_calibration", package = "spDBL")
dc <- dt_calibration
t <- system.time(suppressWarnings(utils::capture.output({
  em <- emulator_learn(dc$training$input, dc$training$responses, dc$Nx, dc$Ny,
    covariance = "sigma2R", episode_window = 2L, nsam = 1L, seed = 1L)
  calibration_fit(calibration_prepare(em, dc$lower, dc$upper), dc$field$z,
    dc$field$U, chains = 1L, iter = 2L, warmup = 1L, refresh = 0L,
    keep_discrepancy = FALSE, cache_dir = Sys.getenv("SPDBL_STAN_CACHE_DIR"))
})))[["elapsed"]]
cat(sprintf("Stan model ready (%.0f s)\n", t))
cat("spDBL", format(packageVersion("spDBL")), "| rstan", format(packageVersion("rstan")),
    "| RobustCalibration", format(packageVersion("RobustCalibration")),
    "| RobustGaSP", format(packageVersion("RobustGaSP")), "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")
EOF
echo "== Setup done. Next: bash benchmarks/hoffman2/smoke.sh"
