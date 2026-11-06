#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E1 (scaling), paper profile: one task at a time on a whole node, so the
# timings are not disturbed by other jobs. The node's 192 GB cover the
# 100^2 grid with episode window 5 (projected at about 68 GB peak RSS).
# Submit with:   bash benchmarks/hoffman2/submit.sh e1
# A job stopped by the 24 h limit is continued by submitting it again
# (--resume keeps finished tasks).
#SBATCH --job-name=spdbl-e1
#SBATCH --constraint=intel-gold-6240
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
. benchmarks/hoffman2/env.sh
Rscript benchmarks/e1_scaling.R --profile=paper --resume
