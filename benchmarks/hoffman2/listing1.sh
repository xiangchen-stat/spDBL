#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Listing 1 of the MLOSS manuscript, run exactly as printed (its
# calibration_fit() call uses cores = 4): 20 emulator draws x 4 chains x 2000
# iterations, about an hour on 4 cores. This checkout is installed into a
# job-local library first, so the run uses exactly the checked-out package;
# dependencies come from the benchmark library.
# Submit with:   bash benchmarks/hoffman2/submit.sh listing1
#SBATCH --job-name=spdbl-listing1
#SBATCH --partition=campus
#SBATCH --constraint=intel-gold-6240
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
. benchmarks/hoffman2/env.sh
LIB=$(mktemp -d "${TMPDIR:-/tmp}/spdbl-listing1.XXXXXX")
trap 'rm -rf "$LIB"' EXIT
R CMD INSTALL --no-test-load --library="$LIB" .
R_LIBS="$LIB" Rscript benchmarks/listing1_run.R
