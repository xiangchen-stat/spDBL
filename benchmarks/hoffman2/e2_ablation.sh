#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E2 ablation (benchmarks/e2_ablation.R) on the paper-profile E2 data: the
# design GP alone and ppgasp per time step against spDBL's simulator-free
# prediction. One core, about 15 minutes; uses the E2 simulated data in
# $SPDBL_BENCH_CACHE and the E2 raw results committed in the repository.
# Submit with:   bash benchmarks/hoffman2/submit.sh e2_ablation
#SBATCH --job-name=spdbl-e2-ablation
#SBATCH --partition=campus
#SBATCH --constraint=intel-gold-6240
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
. benchmarks/hoffman2/env.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript benchmarks/e2_ablation.R
