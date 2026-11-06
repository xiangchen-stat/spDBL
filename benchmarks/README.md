# spDBL benchmark harness

Reproducible benchmarks for the JMLR MLOSS submission of spDBL. They compare
spDBL with related implementations on run time, memory, accuracy and
calibration quality. The scripts live outside the package build
(`^benchmarks$` is in `.Rbuildignore`) and use only R.

| Experiment | Script | What it measures |
|---|---|---|
| E1 emulator scaling | `e1_scaling.R` (+ `e1_worker.R`) | spDBL `emulator_learn()` and `emulator_predict()` wall and CPU time, peak RSS, R-heap peak, `object.size` of the fit, and the share of fit time spent building and factorising the dense S x S `R_field`, across grid sizes and `episode_window`. |
| E2 emulator accuracy | `e2_accuracy.R` (+ `e2_worker.R`) | Held-out RMSE, MAE, 95% interval coverage and mean width, Gaussian CRPS, fit and predict time, and peak RSS for spDBL (3 settings), RobustGaSP `ppgasp` (default and tuned), and an mlegp PCA-GP. |
| E3 calibration | `e3_calibration.R` (+ `e3_worker.R`) | Point error, 95% credible-interval coverage and width for eta, R-hat, bulk/tail ESS, ESS/s, divergences, and wall time for spDBL (`calibration_fit`, RStan) vs RobustCalibration (S-GaSP and GaSP discrepancy), over replicate field datasets with known truth. |
| Check | `check_predictive.R` | Verifies the analytic spDBL predictive variance used in E2 against Monte Carlo built from the package's own samplers. |
| E2 summary | `e2_summary.R` | Means over designs and design-paired t intervals against spDBL (`e2_summary.csv`, `e2_paired.csv`). |
| E2 ablation | `e2_ablation.R` | (First E2 profile; profile `fixed` runs the same ablations as methods.) Where spDBL's simulator-free accuracy comes from: the design GP alone and `ppgasp` fitted per time step, against spDBL, on the same runs. |
| E3 chain diagnostic | `e3_chain_diagnostic.R` | 24 NUTS chains from random starts on one replicate and emulator draw: where each chain starts and settles, and its log density (the local modes in the manuscript). |
| E3 MLE check | `e3_rc_mle_check.R` | Where RobustCalibration stores its maximum-likelihood estimate and how long a multi-start fit takes (used to start its chains in profile `mode`). |
| Listing 1 | `listing1.R`, `listing1_run.R` | The manuscript's Listing 1, run exactly as printed, with its output, run time, git commit and session. |

Shared code is in `R/`: `config.R` (profiles), `common.R` (CLI parsing,
provenance, per-process peak memory, scoring rules), `sir_data.R` (designs
and simulations), and `spdbl_predictive.R` (spDBL predictive distributions).

## Running

Install spDBL and the competitors into a dedicated library rather than the
system library, and point the scripts at it:

```sh
LIB=/path/to/rlib-bench
mkdir -p "$LIB"
Rscript -e 'install.packages(c("RobustGaSP", "RobustCalibration", "mlegp"),
  lib = "'"$LIB"'", repos = "https://cloud.r-project.org")'
rsync -a --exclude .git spDBL/ /tmp/spDBL-bench-src/     # build from a copy
R CMD INSTALL -l "$LIB" /tmp/spDBL-bench-src
export SPDBL_BENCH_LIB="$LIB"                 # prepended to .libPaths()
export SPDBL_STAN_CACHE_DIR=/path/to/stan-cache   # compiled Stan model (E3)
export SPDBL_BENCH_CACHE=/path/to/sim-cache        # cached SIR simulations
```

rstan (>= 2.32.7) and a C++ toolchain are needed for E3. `ps` is optional and
used only for in-process RSS snapshots. Peak memory per task comes from
`/usr/bin/time`, or from the binary named in `SPDBL_TIME`.

`hoffman2/` has setup, smoke-test and job scripts for running the paper
profiles on a Slurm cluster (UCLA Hoffman2); see `hoffman2/README.md`.

One command per experiment (run from the package root):

```sh
Rscript benchmarks/e1_scaling.R     --profile=validation
Rscript benchmarks/e2_accuracy.R    --profile=validation
Rscript benchmarks/e3_calibration.R --profile=validation
Rscript benchmarks/check_predictive.R
```

Use `--profile=paper` for the manuscript settings. Any profile field can be
overridden, for example `--grids=10,20,50,100`, `--windows=1,2,5`,
`--methods=spdbl,rc_gasp`, `--n_rep=20`, `--jobs=4`, `--out=/some/dir`,
`--lib=/path/to/lib`. Unknown options are rejected. Grid sizes of 100 and
above are not in the validation profile; request them explicitly or use
`--profile=paper`.

Add `--resume` to continue an interrupted run: a task is skipped, and its
earlier result reused, when its result file holds a successful run of the
same task under the same settings (`jobs`, `timeout_sec` and the output and
cache paths may differ). Anything else is rerun.

The runs behind the manuscript (E1 needs an otherwise idle machine; on
Hoffman2 use the job scripts in `hoffman2/`):

```sh
Rscript benchmarks/e1_scaling.R     --profile=paper
Rscript benchmarks/e2_accuracy.R    --profile=fixed
Rscript benchmarks/e2_summary.R     --profile=fixed
Rscript benchmarks/e2_ablation.R
Rscript benchmarks/e3_calibration.R --profile=mode --jobs=14
Rscript benchmarks/listing1_run.R
```

Outputs go to `results/<experiment>/<profile>/` (override with `--out`):

* `<experiment>_results.csv`: one row per task (E3: one row per task and
  parameter), including `status`, `exit_status`, timings, `peak_rss_mb` and
  the run `label`.
* `e3_summary.csv`: per-method aggregates for E3.
* `meta.rds` / `meta.txt`: the configuration, the git commit of the
  benchmarks (with `-dirty` if tracked files differed; recorded from
  November 5, 2026), `sessionInfo()`, `Sys.info()`, BLAS/LAPACK, package
  versions, library paths, thread settings, CPU, memory and load.
* `raw/<task>.rds`: the worker result with its own provenance (the same
  fields, recorded inside the worker process) and the task definition.
* `logs/<task>.log`: worker output and the `/usr/bin/time` report.
* `tasks/<task>.rds`: the exact task given to the worker.

### Which results the manuscript uses

| Folder | Profile | Used for |
|---|---|---|
| `results/e1/paper/` | `paper` | E1, Figure 2(a), Table A1 |
| `results/e2/fixed/` | `fixed` | E2, Figure 2(b, c), Tables A2 and A3 |
| `results/e3/mode/` | `mode` | E3, Figure 2(d), Table A4 |
| `results/e3_diagnostic/rep011_draw01_c24/` | `paper` | the local modes in E3 |
| `results/listing1/` | -- | the run time and diagnostics of Listing 1 |

Superseded runs, kept for the record and cited in the manuscript's appendix
where they bear on a result:

* `results/e2/paper/`: the first E2 profile, with one design per grid and
  unit grid cells, so that, as in E1, larger grids were larger domains with
  different epidemics, which confounds resolution with the problem. There
  spDBL's intervals also covered 95.0-95.6% on every grid, but at 50² tuned
  `ppgasp` had lower RMSE than spDBL (0.66 against 0.70) and similar coverage
  (95.4% against 95.0%). The fixed-domain profile replaced it and was set
  before any of its results were seen.
* `results/e3/paper/`: E3 with every chain of both packages started at a
  uniform draw in the box (`init = "random"`). Chains settled in local modes,
  so the pooled draws were unreliable (spDBL Rhat above 1.1 in 26 of 30
  replicates); profile `mode` replaced it.
* `results/e2_ablation/paper/`: `e2_ablation.R` on the first E2 profile's
  data; profile `fixed` includes the same ablations (`designgp_only`,
  `ppgasp_per_time`) as methods.
* `results/e2/paper-sqexp-shared/`, `results/e3/paper-sqexp-shared/`: the
  paper-profile runs with the design GP used before the Matern 5/2 change.

`results/e2/fixed/` (profile `fixed`): the 10^2, 20^2 and 50^2 grids discretise the same 10 x 10
physical domain (`sir_pde_domain()` in `R/sir_data.R`, identical to
`gen_pde()` at 10^2), 10 independent training/test designs per grid, with the
design-GP and per-time `ppgasp` ablations. `e2_summary.R` writes the means
over designs and paired t intervals against spDBL (`e2_summary.csv`,
`e2_paired.csv`); run times are from design 1, one task per node.

All runs used Hoffman2 Xeon Gold 6240 nodes with R 4.6.1 and oneMKL: E1 and E2
one task per node, E3 split over five nodes with 18 one-core tasks each. The
spDBL version is 1.2.0 throughout, but the runs come from development
snapshots of November 3-5, 2026, before the release; from November 5 the
metadata records the commit. The sqexp-shared runs used
`design_kernel = "sqexp", design_groups = "shared"`; the E3 one ran 34 tasks
on one node, so its timings are not comparable, and its RobustCalibration
rows are identical to those in `results/e3/paper/`. Before committing, the
user name and home and scratch directories were replaced by `<user>`, `~`
and `$TMPDIR`.

Simulated data are cached (`SPDBL_BENCH_CACHE`, default `benchmarks/cache/`)
and keyed on grid, design size, number of times and seed, so all methods and
experiments see identical runs.

## Profiles

Defined in `R/config.R`.

| | validation | paper | fixed (E2) / mode (E3) |
|---|---|---|---|
| E1 grids | 10², 20², 30², 40², 50² | 10², 20², 30², 40², 50², 70², 100² | -- |
| E1 `episode_window` | 1, 2, 5 | 1, 2, 5 | -- |
| E1 design | 40 training + 10 held-out runs, T = 26 | same | -- |
| E2 grids | 10², 20² | 10², 20², 50², unit cells | 10², 20², 50² on one 10 x 10 domain |
| E2 design | 40 training + 10 held-out runs, T = 26 | 40 training + 50 held-out runs, T = 26, one design | same, 10 designs per grid |
| E2 spDBL draws (`nsam`) | 20 | 20 | 20 |
| E2 ablations | -- | -- | `designgp_only`, `ppgasp_per_time` |
| E3 replicates | 3 | 30 | 30 |
| E3 spDBL | 4 emulator draws x 2 chains x 500 iter (250 warmup) | 20 draws x 4 chains x 2000 iter (1000 warmup) | same |
| E3 starts | centre of the box | uniform random (`init = "random"`) | spDBL `init = "mode"` (best of 20 L-BFGS optima); RobustCalibration its 10-start MLE, jittered |
| E3 RobustCalibration | S = 2000, S_0 = 500, 2 chains | S = 10000, S_0 = 2000 (package defaults), 4 chains | S = 20000, S_0 = 5000, 4 chains |
| E3 cores per task | 1 | 1 | 1 |
| parallel workers (`jobs`) | E1: 1, E2: 2, E3: 3 | E1, E2: 1; E3: `--jobs=14` | E2: 1; E3: 18 per node |

E2 stops at 50² in the paper profile because spDBL keeps every coefficient
draw: with 20 draws, peak memory at 100² is projected at about 120 GB. E3
runs every task on one core, as RobustCalibration has no parallel chains, so
wall and CPU time agree for both packages; throughput comes from `jobs`.

All validation numbers are labelled `validation only` in the outputs: they come
from a laptop shared with other jobs and from reduced MCMC lengths.

## Data

* **E1/E2.** SIR reaction-diffusion fields from `spDBL::gen_pde()` (response
  `log(I + 1)`, population 10,000, T = 26 time points including t = 0); in
  profile `fixed`, from `sir_pde_domain()` (`R/sir_data.R`), which
  discretises one 10 x 10 domain at every grid size and equals `gen_pde()`
  at 10².
  Inputs: (eta1, eta2, alpha2) in [2, 4] x [0.2, 0.4] x [0, 0.2] with
  alpha1 = alpha3 = 0 (the `dt_calibration` box). Training and held-out inputs
  are independent random Latin hypercubes (seeds `seed` and `seed + 1`), the
  same for every grid size.
* **E3.** Simulator design: the 40 training runs bundled in `dt_calibration`
  (6 x 6 grid, T = 26). Truths: a Latin hypercube in the central 80% of the
  box (`inner`). Field data follow `data-raw/sir-calibration.R`:
  tau2_t ~ IG(30, 2), random-walk discrepancy with increments
  N(0, tau2_t U) (U from `dt_calibration`), and N(0, tau2_t) noise. This is
  spDBL's calibration model, so spDBL is correctly specified and the
  competitors are not (see caveats).

## Metrics

* **Peak RSS.** Every task runs in a fresh `Rscript` process under
  `/usr/bin/time -l` (macOS, bytes) or `/usr/bin/time -v` (Linux, kB); the
  maximum resident set size of that process is reported as `peak_rss_mb`.
  It includes the R baseline (about 100-200 MB with packages and data
  loaded; see `rss_start_mb` in E1). `heap_peak_*_mb` is R's own high-water
  mark (`gc(reset = TRUE)` then `gc()` "max used") for a phase; it counts
  R-allocated vectors only, not Eigen/C++ temporaries. With `cores > 1` in
  E3, rstan runs chains in child processes; `/usr/bin/time` then reports the
  largest single process, not the sum. Without `/usr/bin/time`, the column is
  NA and only the in-process measures remain.
* **Time.** Wall-clock (`proc.time()` elapsed) and CPU (user + system)
  seconds, fit and predict separately. Data loading, simulation and Stan
  compilation are excluded (the Stan model is compiled once by the driver
  through the public API before the workers start).
* **E2 scoring.** Computed over every held-out run, grid cell and time
  t >= `eval_from` = 3 (physical time >= 2). In spDBL 1.1.1 the default AR(2)
  design at t = 1, 2 used the response itself and, at t = 2, the future value
  Y_3, which leaked the target when AR rows were built from held-out
  trajectories. Since 1.2.0 the designs are causal; scoring still starts at
  t = 3 so results compare across versions. CRPS uses the Gaussian closed
  form with each method's predictive mean and standard deviation. Coverage
  and width use each method's own central 95% interval: ppgasp's
  `lower95`/`upper95`; spDBL known-AR: Student-t quantiles (analytic, see
  `R/spdbl_predictive.R`); spDBL design-GP: moment-matched Gaussian of the
  equal-weight mixture over emulator draws, using the package's own
  per-draw mean and variance (including the design-GP variance, which
  `calibration_prepare()` propagates by default); mlegp: mean +/- 1.96 se.
* **E3.** For each active parameter: posterior mean minus truth (also
  divided by the prior range), central 95% interval coverage and width (also
  relative to the range), rank-normalised split R-hat and bulk/tail ESS
  (`rstan::Rhat`, `rstan::ess_bulk`, `rstan::ess_tail`, applied identically to
  both packages). spDBL: per emulator-conditional fit; the table reports the
  worst R-hat and smallest ESS across conditionals, and ESS/s is the summed
  bulk ESS of the slowest parameter across conditionals divided by
  `calibration_fit()` wall time (`ess_per_s`) and CPU time
  (`ess_per_cpu_s`, which counts rstan's chain processes). RobustCalibration:
  chains are separate `rcalibration()` calls with different seeds, run one
  after another; ESS/s is the slowest parameter's bulk ESS divided by the
  summed run time, which is also its CPU time. `e3_summary.csv` adds, per
  method, the share of replicates with R-hat above 1.01 and the total number
  of divergences. spDBL's own
  diagnostics (`divergences`, `treedepth_hits`, `pkg_max_rhat_all`,
  `pkg_min_neff_all`, which include tau2) are also recorded.

## Fairness protocol

* Same simulated designs, held-out runs, field data, parameter bounds and
  hardware for every method; one method per process; single-threaded BLAS
  (`OMP_NUM_THREADS` etc. set to `threads`); R's reference BLAS unless the
  meta file says otherwise.
* Competitors run with their documented defaults plus one documented tuning
  pass, fixed before looking at results:
  * `ppgasp_tuned`: `nugget.est = TRUE`, `optimization = "nelder-mead"`
    (ppgasp itself warns that the derivative-based optimiser may be
    inaccurate with many outputs).
  * spDBL tuning pass `spdbl_s2r_designgp_unit`: inputs rescaled to [0, 1]
    (the default isotropic exponential input kernel on raw inputs is
    dominated by eta1, whose range is ten times larger).
  * RobustCalibration: S-GaSP discrepancy (default) and GaSP
    (Kennedy-O'Hagan) as the alternative.
* Failures, time-outs (`timeout_sec`) and out-of-memory exits are written as
  rows with `status` FAILED/TIMEOUT and the tail of the log; they are never
  dropped.
* Seeds are fixed (`seed` in the profile; task seeds are `seed + replicate`).

## Comparability caveats

* **spDBL known-AR rows (`spdbl_*_known`).** The AR design rows for held-out
  runs are built from the true held-out trajectories, so these are one-step-
  ahead predictions given the true past. They are a conditional/oracle
  reference, not comparable to ppgasp or mlegp, which predict the whole
  trajectory from eta alone. The comparable spDBL setting is
  `spdbl_s2r_designgp` (PDE-free design-GP rows). Its intervals include the
  design-GP variance; spDBL 1.1.1 omitted it, and results from that version
  are not comparable.
* **ppgasp.** Emulates the whole T x S trajectory as k = T S conditionally
  independent outputs with a shared kernel; no temporal structure. It
  rejects constant output columns, so the t = 0 field (identical for all
  runs) is removed before fitting and predicted exactly; this is excluded
  from scoring anyway.
* **mlegp PCA-GP.** Outputs are centred before the SVD (mlegp does not
  centre). The predictive variance propagates the PC-weight GP standard
  errors only; truncation error is not modelled.
* **RobustCalibration.** (i) The discrepancy is a stationary GaSP/S-GaSP over
  (time, row, column) with one homoscedastic noise variance, whereas the data
  are generated (and spDBL models them) as a random walk in time with
  time-varying variances; RobustCalibration has no dynamic discrepancy.
  (ii) It plugs in the ppgasp predictive mean, ignoring emulator
  uncertainty; spDBL integrates the emulator covariance and pools over
  emulator draws. (iii) t = 0 observations are dropped for RobustCalibration
  (900 vs 936 observations) because ppgasp rejects constant columns; t = 0
  carries no information about eta. (iv) `loc_index_emulator` must be passed
  explicitly: with the default `NULL`, `rcalibration()` stores
  `1:k` on the object but passes `NULL` to the sampler, and
  `predict.ppgasp()` then fails ("argument is of length zero").
  (v) Random-walk Metropolis vs NUTS: compare ESS/s, not iterations; ESS/s
  for RobustCalibration comes from chains that had not converged
  (Rhat above 1.01 in every replicate of profile `mode`).
  (vi) Starts follow the profile's `init`. `"midpoint"` (validation): every
  chain at the centre of the box. `"random"` (profile `paper`): every chain
  at a uniform draw in the box, through `calibration_fit(init = "random")`
  and `rcalibration(initial_values = )`; chains then settled in local modes
  (`e3_chain_diagnostic.R`: six modes in one replicate, 36 to 586
  log-density units below the main one, three on the boundary of the box).
  `"mode"` (profile `mode`, the manuscript): spDBL starts every chain of a
  conditional fit, slightly jittered, at the best of 20 L-BFGS optima
  (`calibration_fit()`'s default); RobustCalibration starts at its own
  10-start maximum-likelihood estimate, jittered. Because a fit's chains
  start near one optimum, Rhat then checks mixing within that mode; the
  optima of all starts are recorded (`diagnostics[[j]]$mode_search`).
  (vii) ESS for spDBL is summed over the 20 conditional fits, each keeping
  all 4,000 post-warmup draws; `calibration_fit()` keeps 200 per emulator
  draw by default.
* **Validation runs share the machine** with other jobs and use up to three
  parallel workers; timings are indicative only.

## Resource planning for the paper profile

Written before the paper runs, and kept as a record. Rough log-log
extrapolations from the validation runs (Apple M3 Max,
reference BLAS, one thread; spDBL 1.1.1 before the `R_field` fix):

* **E1, 100², N = 40, T = 26, nsam = 10.** `emulator_learn()` about 4 min
  (w = 1), 7 min (w = 2) and 25 min (w = 5). The dense `R_field` build and
  Cholesky is about 3 min of each, i.e. most of the w = 1 time. Fitted object
  about 8, 16 and 38 GB, and peak RSS about 17, 37 and 55-70 GB: w = 5 is
  likely to fail on a 64 GB machine (memory grows as S^1.5 w, dominated by
  the stored coefficient draws and the per-episode filter/smoother moments).
* **E2, 100².** ppgasp handles k = 260,000 outputs with modest memory. spDBL
  design-GP prediction cost grows with grid size x draws x held-out runs
  (17 s at 20² with 20 draws and 10 runs), so expect tens of minutes at 100²;
  memory follows E1.
* **E3 (6 x 6, 936 observations).** Measured on three paper-setting
  replicates (spDBL 1.2.0, five tasks sharing the machine): spDBL about
  3.1-4.0 CPU-hours per replicate (160,000 NUTS iterations);
  RobustCalibration about 1.8 h (S-GaSP) and 0.6 h (GaSP) for 4 x 10,000
  iterations. That is about 6 CPU-hours per replicate for the three
  methods, so 30 replicates take about 13 h with 14 single-core tasks on
  16 cores. Peak RSS is about 0.8 GB per spDBL process and about 2 GB per
  RobustCalibration process.

## Workarounds for spDBL (no package code is modified)

* Up to spDBL 1.1.1, `emulator_predict()` returned means (or
  coefficient-conditioned mean draws) only. `R/spdbl_predictive.R` derived
  the marginal predictive from the stored filter/smoother output and
  `calibration_prepare()` objects; `check_predictive.R` confirmed it against
  Monte Carlo from the package's own samplers (variance ratio 1.000, 5-95%
  range 0.987-1.016 with 400 draws). Since 1.2.0 the package's
  `predict(..., interval = "prediction")` gives the same distributions, and
  the harness uses it.
* In spDBL 1.1.1 the default AR design for held-out runs drew random numbers
  from the global RNG and used Y_1 (t = 1) and Y_3 (t = 2); the harness builds
  the rows once with a fixed seed. Since 1.2.0 the designs are deterministic
  and causal, and the seed has no effect.
* `calibration_fit()` keeps a random subset of draws by default; the
  harness keeps all post-warmup draws
  (`draws_per_emulator = chains * (iter - warmup)`) so chains can be
  reassembled from `$chain`/`$iteration` for rank-normalised diagnostics.
* The Stan model cache is passed explicitly as `cache_dir`
  (`SPDBL_STAN_CACHE_DIR`), which the package itself reads only in tests.

## TODO (not yet implemented)

* Previous implementation: the research code of Banerjee et al. (2025)
  (Gibbs + Metropolis with a PDE solve per iteration) is cited, not rerun;
  that paper reports 59 min for the 6 x 6 SIR on an Apple M3 Max.
* SEPIA (Python) calibration baseline.
* duqling `pollutant` space-time function as a neutral E2 problem.
* E4 ablations of E3's under-coverage: marginal vs joint emulator sampling,
  design uncertainty on and off, number of emulator draws M, and design rows
  from the simulator (oracle); also marginal vs explicit discrepancy Stan
  program.
* 200² grid via out-of-core filtering (`FF_bigdata_R`): no wrapper or
  backward pass exists yet, so E1 stops at what fits in memory.
* An E3 variant with a stationary-GP discrepancy generator, where
  RobustCalibration rather than spDBL is correctly specified.
* The E1 validation numbers predate the 1.2.0 `R_field` change
  (`emulator_learn()` no longer forms the dense S x S matrix); the paper
  profile is run on 1.2.0.
