#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E3 (calibration), paper profile: 30 replicates x 3 methods, every task on
# one core, NSLOTS - 2 tasks side by side. About 180 CPU-hours in total at
# Apple M3 Max per-core speed, so about 13 h with 16 slots; server cores may
# be slower. With more cores per node, raise the slot count, e.g.
# -pe shared 32 (about 6-7 h). Each task needs about 2 GB
# (RobustCalibration) or 1 GB (spDBL): 16 slots x 4 GB = 64 GB.
# Submit from the repository root:   qsub benchmarks/hoffman2/e3.sh
# If the job reaches h_rt, submit it again: --resume keeps finished tasks
# and reruns only the unfinished ones.
#$ -N spdbl-e3
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -o $HOME/spdbl-bench/logs/$JOB_NAME.$JOB_ID.log
#$ -pe shared 16
#$ -l h_data=4G,h_rt=24:00:00,exclusive
#$ -m ea
set -eo pipefail
. benchmarks/hoffman2/env.sh
JOBS=$(( ${NSLOTS:-16} > 2 ? ${NSLOTS:-16} - 2 : 1 ))
Rscript benchmarks/e3_calibration.R --profile=paper --jobs=$JOBS --resume
