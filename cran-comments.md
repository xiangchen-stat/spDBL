## Submission

Update of spDBL from 1.1.1 to 1.2.0. Main changes (details in NEWS.md):
calibration without simulator runs with propagated design-GP uncertainty,
posterior-mode chain starts and rank-normalised convergence diagnostics,
pointwise prediction intervals, out-of-core smoothing and prediction, and
fewer hard dependencies (ggplot2, ggpubr, scales, deSolve and ReacTran moved
to Suggests; magrittr, invgamma, LaplacesDemon and matrixcalc removed).
Some defaults change results; NEWS.md says how to restore the earlier
behaviour.

## Test environments

* local macOS (aarch64), R 4.4.1: `R CMD check --as-cran`
* GitHub Actions: Ubuntu (R-devel, release, oldrel-1), macOS (release),
  Windows (release)
* win-builder (R-devel): <to run before submission>

## R CMD check results

0 errors | 0 warnings | 0 notes

(Locally, `R CMD check --as-cran` also reports "unable to verify current
time" and skips math rendering in the HTML manual because V8 is not
installed; both are local to the check machine.)

## Reverse dependencies

There are no reverse dependencies on CRAN.

## Notes for the CRAN team

* Calibration uses rstan, which is in Suggests; examples that compile Stan
  run only interactively, and the tests that compile Stan run only when
  `SPDBL_RUN_STAN_TESTS=true`.
