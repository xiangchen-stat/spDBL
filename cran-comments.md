## Submission

Update of spDBL from 1.1.1 to 1.2.0. Main changes (details in NEWS.md):
calibration without simulator runs with propagated design-GP uncertainty,
posterior-mode chain starts and rank-normalised convergence diagnostics,
pointwise prediction intervals, out-of-core smoothing and prediction,
batched simulator-free prediction (about 40 times faster), coefficient draws
that can be regenerated on demand instead of stored (`keep_draws = FALSE`,
about half the memory on large grids), faster calibration gradients,
conditional fits that can run in parallel (`parallel_draws`), and fewer hard
dependencies (ggplot2, ggpubr, scales, deSolve and ReacTran moved to
Suggests; magrittr, invgamma, LaplacesDemon and matrixcalc removed). The
CITATION file now gives the package's CRAN DOI.
Some defaults change results; NEWS.md says how to restore the earlier
behaviour.

## Test environments

* local macOS (aarch64), R 4.4.1: `R CMD check --as-cran` (2026-11-15)
* GitHub Actions: Ubuntu (R-devel, release, oldrel-1), macOS (release),
  Windows (release)
* win-builder (R-devel, 2026-11-14 r90890 ucrt): <result of the rerun>

## R CMD check results

0 errors | 0 warnings | 1 note

* Possibly misspelled word in DESCRIPTION: "gridded". The word is spelled
  correctly (output on a spatial grid).

The first win-builder run also timed three examples at 16-17 s; they now fit
on 6 of the bundled data's 26 time steps.

(Locally, `R CMD check --as-cran` also reports "unable to verify current
time" and skips math rendering in the HTML manual because V8 is not
installed; both are local to the check machine.)

## Reverse dependencies

There are no reverse dependencies on CRAN.

## Notes for the CRAN team

* Calibration uses rstan, which is in Suggests; examples that compile Stan
  run only interactively, and the tests that compile Stan run only when
  `SPDBL_RUN_STAN_TESTS=true`.
