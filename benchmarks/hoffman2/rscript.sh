#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Run one R script of this checkout as a Slurm job (diagnostics, pilots). The
# checkout is installed into a job-local library from a clean copy of the
# commit, ahead of the benchmark library that holds the dependencies.
#   SPDBL_RSCRIPT=benchmarks/e3_chain_diagnostic.R \
#   SPDBL_RARGS="--rep=11 --draw=1 --chains=24 --cores=8" \
#   SPDBL_SBATCH_ARGS="--cpus-per-task=8 --mem=32G --time=3:00:00" \
#     bash benchmarks/hoffman2/submit.sh rscript
#SBATCH --job-name=spdbl-rscript
#SBATCH --partition=campus
#SBATCH --constraint=intel-gold-6240
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
[ -n "$SPDBL_RSCRIPT" ] || { echo "Set SPDBL_RSCRIPT to the R script to run" >&2; exit 1; }
. benchmarks/hoffman2/env.sh
WORK_TMP=$(mktemp -d "${TMPDIR:-/tmp}/spdbl-rscript.XXXXXX")
trap 'rm -rf "$WORK_TMP"' EXIT
mkdir -p "$WORK_TMP/src" "$WORK_TMP/lib"
git archive HEAD | tar -x -C "$WORK_TMP/src"
R CMD INSTALL --no-test-load --library="$WORK_TMP/lib" "$WORK_TMP/src"
export SPDBL_BENCH_LIB=$WORK_TMP/lib
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
echo "Running $SPDBL_RSCRIPT $SPDBL_RARGS"
Rscript "$SPDBL_RSCRIPT" $SPDBL_RARGS
