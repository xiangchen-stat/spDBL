#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# E2 (accuracy), paper profile: one task at a time on a whole node (fit and
# prediction times are reported). spDBL at 50^2 with 20 draws is projected
# at about 15 GB peak RSS.
# Submit with:   bash benchmarks/hoffman2/submit.sh e2
# A job stopped by the 24 h limit is continued by submitting it again
# (--resume keeps finished tasks).
#SBATCH --job-name=spdbl-e2
#SBATCH --constraint=intel-gold-6240
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
. benchmarks/hoffman2/env.sh
Rscript benchmarks/e2_accuracy.R --profile=paper --resume
