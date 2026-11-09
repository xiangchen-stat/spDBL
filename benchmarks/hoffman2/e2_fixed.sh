#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E2 on the fixed physical domain with 10 training/test designs (profile
# "fixed"). This checkout is installed into a job-local library first, so
# the workers use exactly the checked-out package (and its
# predict(..., interval = "prediction")); dependencies come from the
# benchmark library. SPDBL_E2_REPS selects designs and SPDBL_E2_JOBS the
# number of tasks side by side (default 1). The paper run:
#   SPDBL_E2_REPS=1 bash benchmarks/hoffman2/submit.sh e2_fixed        # timing: one task at a time
#   SPDBL_E2_REPS=2,3,4,5,6,7,8,9,10 SPDBL_E2_JOBS=16 \
#     bash benchmarks/hoffman2/submit.sh e2_fixed                      # accuracy, side by side
#   SPDBL_SBATCH_ARGS="--dependency=afterok:<both job ids>" \
#     bash benchmarks/hoffman2/submit.sh e2_fixed                      # collect (--resume)
# Then: Rscript benchmarks/e2_summary.R
#SBATCH --job-name=spdbl-e2-fixed
#SBATCH --constraint=intel-gold-6240
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
. benchmarks/hoffman2/env.sh
# Install from a clean copy of the commit in job-local storage: jobs that
# start together from one clone must not compile in the same src/.
WORK_TMP=$(mktemp -d "${TMPDIR:-/tmp}/spdbl-e2fixed.XXXXXX")
trap 'rm -rf "$WORK_TMP"' EXIT
mkdir -p "$WORK_TMP/src" "$WORK_TMP/lib"
git archive HEAD | tar -x -C "$WORK_TMP/src"
LIB=$WORK_TMP/lib
R CMD INSTALL --no-test-load --library="$LIB" "$WORK_TMP/src"
export SPDBL_BENCH_LIB=$LIB
Rscript benchmarks/e2_accuracy.R --profile=${SPDBL_E2_PROFILE:-fixed} --resume --jobs=${SPDBL_E2_JOBS:-1} \
  ${SPDBL_E2_REPS:+--reps=$SPDBL_E2_REPS}
