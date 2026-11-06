#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# One-time preparation as a batch job: setup.sh (packages, spDBL, Stan
# model) followed by smoke.sh (tiny E1/E2/E3 runs). Same node type as the
# paper jobs, so the compiled Stan model and packages match them.
# Submit with:   bash benchmarks/hoffman2/submit.sh prepare
#SBATCH --job-name=spdbl-prepare
#SBATCH --constraint=intel-gold-6240
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=spdbl-%x.%j.log
set -eo pipefail
bash benchmarks/hoffman2/setup.sh
bash benchmarks/hoffman2/smoke.sh
