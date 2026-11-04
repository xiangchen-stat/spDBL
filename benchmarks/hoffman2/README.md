# Running the paper benchmarks on UCLA Hoffman2

Scripts for Hoffman2's Slurm scheduler. They run the same harness and
profiles as on a workstation (`../README.md`); only the environment and the
job requests are cluster-specific. Everything the runs write outside the
repository goes to `~/spdbl-bench/` (R library, simulation and Stan caches,
job logs); results go to `benchmarks/results/<exp>/paper/` in the clone, as
usual.

| File | Purpose |
|---|---|
| `env.sh` | Loads R (`SPDBL_R_MODULE`, default `R/4.4.1`, the version used on the workstation) and sets the library, cache and GNU time paths. Sourced by the others. |
| `setup.sh` | One-time install: dependencies, rstan and the competitors into `~/spdbl-bench/rlib`, spDBL from the clone, GNU time if missing, and the compiled Stan model. |
| `smoke.sh` | Tiny E1/E2/E3 runs (a few minutes) that check the whole pipeline, including peak-memory measurement, before the long jobs. |
| `prepare.sh` | Batch job running `setup.sh` then `smoke.sh` (8 CPUs, 32 GB, 4 h). |
| `e1.sh`, `e2.sh`, `e3.sh` | Batch jobs for the paper profiles. |
| `submit.sh` | Submits any of the above with the clone as working directory and the log in `~/spdbl-bench/logs/`. |

## Steps

1. Log in and clone the branch under test:

   ```sh
   ssh <user>@hoffman2.idre.ucla.edu
   git clone -b mloss-readiness https://github.com/xiangchen-stat/spDBL.git ~/spDBL
   cd ~/spDBL
   ```

2. Prepare (30-60 min, mostly compiling rstan), and check the log ends with
   "Smoke test passed":

   ```sh
   bash benchmarks/hoffman2/submit.sh prepare
   tail -f ~/spdbl-bench/logs/spdbl-prepare.*.log
   ```

3. Submit the paper jobs. Each takes a whole node, so they can run at the
   same time:

   ```sh
   bash benchmarks/hoffman2/submit.sh e1 e2 e3
   ```

4. Monitor: `squeue -u $USER`; the driver log
   `~/spdbl-bench/logs/spdbl-e<k>.<job id>.log` has one line per finished
   task; per-task logs are in `benchmarks/results/e<k>/paper/logs/`.

5. A job stopped by the 24 h limit (or anything else) is continued by
   submitting the same job again: every job passes `--resume`, so
   successful tasks with unchanged settings are kept and only the rest run.

6. Copy the results back, from the workstation clone:

   ```sh
   for e in e1 e2 e3; do
     rsync -av <user>@hoffman2.idre.ucla.edu:spDBL/benchmarks/results/$e/paper/ \
       benchmarks/results/$e/paper/
   done
   ```

   `meta.txt` and `raw/*.rds` record the host name, user and paths; review
   them before committing.

## Resource requests

All jobs use the `campus` partition (24 h limit) and Intel Xeon Gold 6240
nodes (36 cores, 192 GB), the most common node type there, so the three
experiments are timed on the same hardware.

| Job | Request | Why |
|---|---|---|
| E1 | whole node | One task at a time for clean timings; 100² with window 5 is projected at about 68 GB. |
| E2 | whole node | One task at a time; spDBL at 50² is projected at about 15 GB. |
| E3 | whole node | 34 single-core tasks side by side. |

To use another node type, pass its feature (listed by
`sinfo -p campus -o "%f"`), e.g.
`SPDBL_SBATCH_ARGS="--constraint=intel-gold-6442Y" bash benchmarks/hoffman2/submit.sh e3`.
The E3 task count follows the node's CPU count.

## Notes

* Threads: the R module and the harness set `OMP_NUM_THREADS`,
  `MKL_NUM_THREADS` and related variables to 1, so the module's Intel MKL
  BLAS runs single-threaded. `meta.txt` records the BLAS, the CPU model, the
  node's memory and the Slurm job ID and CPU count.
* Peak memory per task needs GNU time (`time -v`), which the nodes have in
  `/usr/bin/time`; `setup.sh` builds GNU time 1.9 into `~/spdbl-bench/bin`
  otherwise, and `smoke.sh` fails if peak memory is missing.
* The compiled Stan model is cached per R, rstan and platform version, so
  changing the R module triggers one recompilation.
