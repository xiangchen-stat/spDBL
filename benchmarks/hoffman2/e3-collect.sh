#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Collect E3 after it was split over nodes (see e3.sh): rerun the driver for
# all replicates with --resume, which reuses every finished task and writes
# e3_results.csv and e3_summary.csv for the whole run. Submit it with a
# dependency on the split jobs:
#   SPDBL_SBATCH_ARGS="--dependency=afterok:<ids>" \
#     bash benchmarks/hoffman2/submit.sh e3-collect
# If a task is missing it runs here, so check the split jobs' logs first.
#SBATCH --job-name=spdbl-e3-collect
#SBATCH --partition=campus
#SBATCH --constraint=intel-gold-6240
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
. benchmarks/hoffman2/env.sh
# Install this checkout from a clean copy in job-local storage, ahead of the
# benchmark library, so the run uses exactly the checked-out package.
WORK_TMP=$(mktemp -d "${TMPDIR:-/tmp}/spdbl-e3collect.XXXXXX")
trap 'rm -rf "$WORK_TMP"' EXIT
mkdir -p "$WORK_TMP/src" "$WORK_TMP/lib"
git archive HEAD | tar -x -C "$WORK_TMP/src"
R CMD INSTALL --no-test-load --library="$WORK_TMP/lib" "$WORK_TMP/src"
export SPDBL_BENCH_LIB=$WORK_TMP/lib
Rscript benchmarks/e3_calibration.R --profile=${SPDBL_E3_PROFILE:-paper} --jobs=$SPDBL_CPUS --resume
