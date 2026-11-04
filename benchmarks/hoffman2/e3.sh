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
JOBS=$(( SPDBL_CPUS > 2 ? SPDBL_CPUS - 2 : 1 ))
Rscript benchmarks/e3_calibration.R --profile=paper --jobs=$JOBS --resume
