#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E1 (scaling), paper profile. One task at a time on a whole node, so the
# timings are not disturbed by other jobs. 8 slots x 16 GB = 128 GB: the
# 100^2 grid with episode window 5 is projected at about 68 GB peak RSS.
# Submit from the repository root:   qsub benchmarks/hoffman2/e1.sh
# If the job reaches h_rt, submit it again: --resume keeps finished tasks.
#$ -N spdbl-e1
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -o $HOME/spdbl-bench/logs/$JOB_NAME.$JOB_ID.log
#$ -pe shared 8
#$ -l h_data=16G,h_rt=24:00:00,exclusive
#$ -m ea
set -eo pipefail
. benchmarks/hoffman2/env.sh
Rscript benchmarks/e1_scaling.R --profile=paper --resume
