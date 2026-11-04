# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Environment for the spDBL benchmarks on UCLA Hoffman2 (Slurm). Sourced
# from the repository root by setup.sh, smoke.sh and the job scripts. Any
# variable can be set beforehand, e.g.
#   export SPDBL_R_MODULE=R/4.6.1      # see: module avail R

if [ ! -f DESCRIPTION ] || ! grep -q "^Package: spDBL" DESCRIPTION; then
  echo "Run from the root of the spDBL repository (found no spDBL DESCRIPTION in $PWD)" >&2
  return 1 2>/dev/null || exit 1
fi
export SPDBL_ROOT=$PWD

if ! command -v module >/dev/null 2>&1 && [ -f /u/local/Modules/default/init/modules.sh ]; then
  . /u/local/Modules/default/init/modules.sh
fi
# R 4.4.1 matches the workstation runs; the module also sets the BLAS
# (Intel MKL) thread count to 1.
SPDBL_R_MODULE=${SPDBL_R_MODULE:-R/4.4.1}
if command -v module >/dev/null 2>&1; then
  # Not in a pipeline: module changes the environment of this shell.
  module load $SPDBL_EXTRA_MODULES "$SPDBL_R_MODULE" || true
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found after 'module load $SPDBL_R_MODULE'. Available R modules:" >&2
  module avail R 2>&1 | grep -i "R/" >&2 || true
  echo "Set SPDBL_R_MODULE to one of them." >&2
  return 1 2>/dev/null || exit 1
fi

# Everything the benchmarks write outside the repository. $HOME is kept
# ($SCRATCH is purged); the total is a few GB.
export SPDBL_WORK=${SPDBL_WORK:-$HOME/spdbl-bench}
export SPDBL_BENCH_LIB=${SPDBL_BENCH_LIB:-$SPDBL_WORK/rlib}
export SPDBL_BENCH_CACHE=${SPDBL_BENCH_CACHE:-$SPDBL_WORK/sim-cache}
export SPDBL_STAN_CACHE_DIR=${SPDBL_STAN_CACHE_DIR:-$SPDBL_WORK/stan-cache}
export R_USER_CACHE_DIR=${R_USER_CACHE_DIR:-$SPDBL_WORK/r-cache}
export R_LIBS_USER=$SPDBL_BENCH_LIB
# Peak memory per task comes from GNU time; setup.sh builds one if the node
# has no /usr/bin/time.
if [ -z "$SPDBL_TIME" ]; then
  if [ -x /usr/bin/time ]; then SPDBL_TIME=/usr/bin/time
  else SPDBL_TIME=$SPDBL_WORK/bin/time; fi
fi
export SPDBL_TIME
mkdir -p "$SPDBL_BENCH_LIB" "$SPDBL_BENCH_CACHE" "$SPDBL_STAN_CACHE_DIR" \
  "$R_USER_CACHE_DIR" "$SPDBL_WORK/logs"

# CPUs this job may use: the Slurm allocation, else the machine.
SPDBL_CPUS=${SLURM_CPUS_ON_NODE:-${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 1)}}
export SPDBL_CPUS

echo "spDBL benchmarks: host $(hostname), job ${SLURM_JOB_ID:-interactive}, CPUs $SPDBL_CPUS"
echo "  $(Rscript --version 2>&1 | head -1)"
echo "  CPU: $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | sed 's/.*: //')"
echo "  library: $SPDBL_BENCH_LIB"
