# Running the paper benchmarks on UCLA Hoffman2

Scripts for the Grid Engine (UGE) scheduler on Hoffman2. They run the same
harness and profiles as on a workstation (`../README.md`); only the
environment and the job requests are cluster-specific. Everything the runs
write outside the repository goes to `~/spdbl-bench/` (R library, simulation
and Stan caches, job logs); results go to `benchmarks/results/<exp>/paper/`
in the clone, as usual.

| File | Purpose |
|---|---|
| `env.sh` | Loads R (`SPDBL_R_MODULE`, default `R/4.4.0`) and sets the library, cache and GNU time paths. Sourced by the others. |
| `setup.sh` | One-time install: dependencies, rstan and the competitors into `~/spdbl-bench/rlib`, spDBL from the clone, GNU time if missing, and the compiled Stan model. |
| `smoke.sh` | Tiny E1/E2/E3 runs (10-15 min) that check the whole pipeline, including peak-memory measurement, before the long jobs. |
| `e1.sh`, `e2.sh`, `e3.sh` | Job scripts for the paper profiles. |

## Steps

1. Log in and clone the branch under test:

   ```sh
   ssh <user>@hoffman2.idre.ucla.edu
   git clone -b mloss-readiness https://github.com/xiangchen-stat/spDBL.git ~/spDBL
   cd ~/spDBL
   module avail R            # pick an R >= 4.4 module if R/4.4.0 is absent
   ```

2. Set up in an interactive compute session (30-60 min, mostly compiling
   rstan). If the R module name differs, export it first, for example
   `export SPDBL_R_MODULE=R/4.4.1`.

   ```sh
   qrsh -l h_rt=4:00:00,h_data=4G -pe shared 4
   cd ~/spDBL
   bash benchmarks/hoffman2/setup.sh
   bash benchmarks/hoffman2/smoke.sh
   exit
   ```

   If package downloads fail on the compute node, run `setup.sh` on the
   login node instead, then `smoke.sh` in the `qrsh` session.

3. Submit the paper jobs from the repository root. They use separate whole
   nodes (`exclusive`) and can run at the same time:

   ```sh
   cd ~/spDBL
   qsub benchmarks/hoffman2/e1.sh
   qsub benchmarks/hoffman2/e2.sh
   qsub benchmarks/hoffman2/e3.sh
   ```

   If `SPDBL_R_MODULE` is not the default, pass it on: `qsub -v
   SPDBL_R_MODULE=R/4.4.1 benchmarks/hoffman2/e1.sh`.

4. Monitor: `qstat -u $USER`; the driver log is
   `~/spdbl-bench/logs/spdbl-e<k>.<job id>.log` (one line per finished
   task); per-task logs are in `benchmarks/results/e<k>/paper/logs/`.

5. A job stopped by its 24 h limit (or anything else) is resumed by
   submitting the same script again: every script passes `--resume`, so
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

| Job | Request | Why |
|---|---|---|
| E1 | 8 slots x 16 GB, whole node, 24 h | One task at a time for clean timings; 100² with window 5 is projected at about 68 GB. |
| E2 | 4 slots x 8 GB, whole node, 24 h | One task at a time; spDBL at 50² is projected at about 15 GB. |
| E3 | 16 slots x 4 GB, whole node, 24 h | 14 single-core tasks side by side; about 13 h at M3 Max per-core speed. Use more slots on larger nodes (`qhost` lists them); the task count follows `NSLOTS`. |

Jobs without `highp` are limited to 24 h; `--resume` makes resubmission
cheap. Edit the `#$` lines to change the requests.

## Notes

* Threads: the harness sets `OMP_NUM_THREADS` and related variables to 1, so
  an optimised BLAS in the R module runs single-threaded, as on the
  workstation. `meta.txt` records the BLAS, the CPU model, the node's
  memory and the job's `JOB_ID` and `NSLOTS`.
* Peak memory per task needs GNU time (`time -v`). `setup.sh` builds GNU
  time 1.9 into `~/spdbl-bench/bin` when the node has no `/usr/bin/time`,
  and `smoke.sh` fails if peak memory is missing.
* The compiled Stan model is cached per R, rstan and platform version, so
  changing the R module triggers one recompilation.
