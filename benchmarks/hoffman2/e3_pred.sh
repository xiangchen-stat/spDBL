#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Site and off-grid prediction studies (e3_sites.R, e3_offgrid.R) as small
# jobs on shared nodes. These studies score accuracy, not time, so they need
# neither a whole node nor a particular CPU model, and small jobs start
# sooner. Each job runs the replicates in SPDBL_PRED_REPS side by side, two
# cores each (two chains). SPDBL_PRED_STUDY is sites or offgrid;
# SPDBL_PRED_PROFILE defaults to the study's name. For example, 30
# replicates in eight jobs of up to four:
#   for r in 1-4 5-8 9-12 13-16 17-20 21-24 25-28 29-30; do
#     SPDBL_PRED_STUDY=sites SPDBL_PRED_REPS=$(seq -s, ${r%-*} ${r#*-}) \
#       bash benchmarks/hoffman2/submit.sh e3_pred
#   done
# Then collect every replicate (each is reused, none rerun) with a one-core
# job after them:
#   SPDBL_PRED_STUDY=sites SPDBL_SBATCH_ARGS="--cpus-per-task=1 --mem=4G \
#     --dependency=afterok:<job ids>" bash benchmarks/hoffman2/submit.sh e3_pred
#SBATCH --job-name=spdbl-e3-pred
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
. benchmarks/hoffman2/env.sh
STUDY=${SPDBL_PRED_STUDY:?Set SPDBL_PRED_STUDY to sites or offgrid}
WORK_TMP=$(mktemp -d "${TMPDIR:-/tmp}/spdbl-e3-pred.XXXXXX")
trap 'rm -rf "$WORK_TMP"' EXIT
mkdir -p "$WORK_TMP/src" "$WORK_TMP/lib"
git archive HEAD | tar -x -C "$WORK_TMP/src"
R CMD INSTALL --no-test-load --library="$WORK_TMP/lib" "$WORK_TMP/src"
export SPDBL_BENCH_LIB=$WORK_TMP/lib
JOBS=${SPDBL_PRED_JOBS:-$(( SPDBL_CPUS / 2 > 0 ? SPDBL_CPUS / 2 : 1 ))}
Rscript "benchmarks/e3_$STUDY.R" --profile=${SPDBL_PRED_PROFILE:-$STUDY} --jobs=$JOBS --resume \
  ${SPDBL_PRED_REPS:+--reps=$SPDBL_PRED_REPS}
