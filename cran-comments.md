## Resubmission

This is a resubmission of 1.2.0. The incoming pretest of the first
submission (2026-11-16) reported "Overall checktime 14 min > 10 min" on
Windows R-devel (Debian: OK). In this version:

* `vignettes/calibration.Rmd` skips its two slowest chunks (the emulator fit
  and the held-out validation) while `R CMD check` runs; `R CMD build` runs
  every chunk, so the vignette shipped in the tarball is unchanged.
* The six slowest tests (Monte Carlo cross-checks and two plotting tests)
  use `skip_on_cran()`; our GitHub Actions CI runs them.

Rebuilding the vignettes went from 279 s to <win-builder> and the tests from
194 s to <win-builder> on win-builder (Windows R-devel); the whole check now
takes <win-builder> there.

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

This update comes less than three weeks after 1.1.1 (October 29, 2026). It
fixes intervals that were too narrow in 1.1.1, which ignored the predictive
variance of the design Gaussian process (95% prediction intervals for
held-out simulator runs covered 54% before and 94% after), and it is the
version described in a manuscript for the JMLR Machine Learning Open Source
Software track, which archives the reviewed version. We do not plan
another update for at least two months.

## Test environments

* local macOS (aarch64), R 4.4.1: `R CMD check --as-cran` (2026-11-16)
* GitHub Actions: Ubuntu (release, oldrel-1; R-devel weekly), macOS
  (release), Windows (release)
* win-builder (R-devel, 2026-11-14 r90890 ucrt): Status OK

## R CMD check results

0 errors | 0 warnings | 0 notes

A first win-builder run gave two NOTEs, both fixed: "gridded" flagged as
possibly misspelled in DESCRIPTION (reworded), and three examples at 16-17 s
(they now fit on 6 of the bundled data's 26 time steps).

(Locally, `R CMD check --as-cran` also reports "unable to verify current
time" and skips math rendering in the HTML manual because V8 is not
installed; both are local to the check machine.)

## Reverse dependencies

There are no reverse dependencies on CRAN.

## Notes for the CRAN team

* Calibration uses rstan, which is in Suggests; examples that compile Stan
  run only interactively, and the tests that compile Stan run only when
  `SPDBL_RUN_STAN_TESTS=true`.
