#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E2 (accuracy), paper profile. One task at a time on a whole node (fit and
# prediction times are reported). 4 slots x 8 GB = 32 GB: spDBL at 50^2
# with 20 draws is projected at about 15 GB peak RSS.
# Submit from the repository root:   qsub benchmarks/hoffman2/e2.sh
# If the job reaches h_rt, submit it again: --resume keeps finished tasks.
#$ -N spdbl-e2
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -o $HOME/spdbl-bench/logs/$JOB_NAME.$JOB_ID.log
#$ -pe shared 4
#$ -l h_data=8G,h_rt=24:00:00,exclusive
#$ -m ea
set -eo pipefail
. benchmarks/hoffman2/env.sh
Rscript benchmarks/e2_accuracy.R --profile=paper --resume
