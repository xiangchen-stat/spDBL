#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Submit Hoffman2 jobs from the repository root, with the job log in
# ~/spdbl-bench/logs/ (the default SPDBL_WORK) and the clone as working
# directory:
#   bash benchmarks/hoffman2/submit.sh prepare
#   bash benchmarks/hoffman2/submit.sh e1 e2 e3
# Extra sbatch options can be given in SPDBL_SBATCH_ARGS, e.g.
#   SPDBL_SBATCH_ARGS="--constraint=intel-gold-6442Y" bash benchmarks/hoffman2/submit.sh e3
set -eo pipefail
if [ ! -f DESCRIPTION ] || ! grep -q "^Package: spDBL" DESCRIPTION; then
  echo "Run from the root of the spDBL repository" >&2
  exit 1
fi
[ $# -gt 0 ] || { echo "Usage: $0 prepare|e1|e2|e3 ..." >&2; exit 1; }
logs=${SPDBL_WORK:-$HOME/spdbl-bench}/logs
mkdir -p "$logs"
for job in "$@"; do
  script=benchmarks/hoffman2/$job.sh
  [ -f "$script" ] || { echo "No job script $script" >&2; exit 1; }
  sbatch --chdir="$PWD" --output="$logs/%x.%j.log" $SPDBL_SBATCH_ARGS "$script"
done
echo "Logs: $logs (squeue -u \$USER to follow the jobs)"
