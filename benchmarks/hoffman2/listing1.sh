#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Listing 1 of the MLOSS manuscript, run exactly as printed (its
# calibration_fit() call uses parallel_draws = 20: 20 conditional fits at
# once, each running its chains in sequence): 20 emulator draws x 4 chains x 2000 iterations. This checkout is installed into a
# job-local library first, so the run uses exactly the checked-out package;
# dependencies come from the benchmark library.
# Submit with:   bash benchmarks/hoffman2/submit.sh listing1
#SBATCH --job-name=spdbl-listing1
#SBATCH --constraint=intel-gold-6240
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
. benchmarks/hoffman2/env.sh
# Install from a clean copy of the commit in job-local storage: jobs that
# start together from one clone must not compile in the same src/.
WORK_TMP=$(mktemp -d "${TMPDIR:-/tmp}/spdbl-listing1.XXXXXX")
trap 'rm -rf "$WORK_TMP"' EXIT
mkdir -p "$WORK_TMP/src" "$WORK_TMP/lib"
git archive HEAD | tar -x -C "$WORK_TMP/src"
LIB=$WORK_TMP/lib
R CMD INSTALL --no-test-load --library="$LIB" "$WORK_TMP/src"
# One BLAS thread: parallel_draws forks R, and forking after a multithreaded
# BLAS (oneMKL with OpenMP here) leaves the forked processes hung.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R_LIBS="$LIB" Rscript benchmarks/listing1_run.R
