#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# On-grid site prediction (e3_sites.R), profile sites: 30 replicates, one
# task per replicate on two cores (its two chains), 15 tasks side by side on
# a whole 36-core node. Each task calibrates on 30 of the 36 sites (100
# training runs, 10 emulator draws, mode starts), then predicts the other 6:
# about 1 CPU-hour, so the job takes about 3 hours.
# Submit with:   bash benchmarks/hoffman2/submit.sh e3_sites
# SPDBL_SITES_REPS lists the replicates (default all), SPDBL_SITES_JOBS the
# tasks run side by side, SPDBL_SITES_PROFILE the profile (default sites).
# A job stopped by the time limit is continued by submitting it again.
#SBATCH --job-name=spdbl-e3-sites
#SBATCH --constraint=intel-gold-6240
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
. benchmarks/hoffman2/env.sh
# Install this checkout from a clean copy in job-local storage, ahead of the
# benchmark library, so the run uses exactly the checked-out package.
WORK_TMP=$(mktemp -d "${TMPDIR:-/tmp}/spdbl-e3-sites.XXXXXX")
trap 'rm -rf "$WORK_TMP"' EXIT
mkdir -p "$WORK_TMP/src" "$WORK_TMP/lib"
git archive HEAD | tar -x -C "$WORK_TMP/src"
R CMD INSTALL --no-test-load --library="$WORK_TMP/lib" "$WORK_TMP/src"
export SPDBL_BENCH_LIB=$WORK_TMP/lib
JOBS=${SPDBL_SITES_JOBS:-$(( SPDBL_CPUS / 2 > 3 ? SPDBL_CPUS / 2 - 3 : 1 ))}
Rscript benchmarks/e3_sites.R --profile=${SPDBL_SITES_PROFILE:-sites} --jobs=$JOBS --resume \
  ${SPDBL_SITES_REPS:+--reps=$SPDBL_SITES_REPS}
