#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E3 (calibration), paper profile: 30 replicates x 3 methods, every task on
# one core, (node CPUs - 2) tasks side by side on a whole node: 34 on a
# 36-core Xeon Gold 6240. About 180 CPU-hours at Apple M3 Max per-core
# speed, more on these cores; the longest single task (spDBL, 160,000 NUTS
# iterations) sets the floor. Each task needs about 1 GB (spDBL) or 2 GB
# (RobustCalibration).
# Submit with:   bash benchmarks/hoffman2/submit.sh e3
# A job stopped by the 24 h limit is continued by submitting it again
# (--resume keeps finished tasks and reruns only the unfinished ones).
#
# With 34 tasks per node, memory-bandwidth contention slows every task
# (RobustCalibration about 5x, spDBL about 3x against a workstation), which
# also distorts the timing comparison. The paper run therefore splits E3 over
# five nodes, six replicates (18 tasks) each, then collects the results:
#   for r in 1-6 7-12 13-18 19-24 25-30; do
#     SPDBL_E3_REPS=$(seq -s, ${r%-*} ${r#*-}) SPDBL_E3_JOBS=18 \
#       bash benchmarks/hoffman2/submit.sh e3
#   done
#   SPDBL_SBATCH_ARGS="--dependency=afterok:<the five job ids>" \
#     bash benchmarks/hoffman2/submit.sh e3-collect
# SPDBL_E3_REPS lists this job's replicates; SPDBL_E3_JOBS sets the number of
# tasks run side by side (default: CPUs - 2); SPDBL_E3_PROFILE selects the
# profile (default paper; "mode" starts every chain at its posterior mode).
#SBATCH --job-name=spdbl-e3
#SBATCH --partition=campus
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
WORK_TMP=$(mktemp -d "${TMPDIR:-/tmp}/spdbl-e3.XXXXXX")
trap 'rm -rf "$WORK_TMP"' EXIT
mkdir -p "$WORK_TMP/src" "$WORK_TMP/lib"
git archive HEAD | tar -x -C "$WORK_TMP/src"
R CMD INSTALL --no-test-load --library="$WORK_TMP/lib" "$WORK_TMP/src"
export SPDBL_BENCH_LIB=$WORK_TMP/lib
JOBS=${SPDBL_E3_JOBS:-$(( SPDBL_CPUS > 2 ? SPDBL_CPUS - 2 : 1 ))}
Rscript benchmarks/e3_calibration.R --profile=${SPDBL_E3_PROFILE:-paper} --jobs=$JOBS --resume \
  ${SPDBL_E3_REPS:+--reps=$SPDBL_E3_REPS}
