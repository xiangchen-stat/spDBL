These research vignettes rerun both scripts in the repository's
`vignettes/reproduce/` using installed spDBL and the original CSV data.

From the repository root:

```sh
mkdir -p spDBL/.r-library
R CMD INSTALL --library=spDBL/.r-library spDBL
# Inspect the full configuration without starting computations:
Rscript spDBL/vignettes/reproduce/run_reproduction.R calibration full --dry-run
# Full calibration: 30,000 emulator draws, one transition per draw:
Rscript spDBL/vignettes/reproduce/run_reproduction.R calibration full
# Completed quick calibration: 300 emulator draws:
Rscript spDBL/vignettes/reproduce/run_reproduction.R calibration quick
```

Use `emulation` to run the emulation experiment, or `all full` / `all quick`
to run both experiments with the selected calibration profile.
The command-line runner defaults to `all full`.
The script works from any directory. Each Rmd can also be rendered directly
with `rmarkdown::render()`; it selects the local package installation.
Direct rendering of `calibration.Rmd` defaults to the quick profile; pass
`params = list(profile = "full")` to render the full profile.

Outputs are `emulation.html`, `calibration.html` (quick),
`calibration-full.html` (full), and `_results/` (fits, CSV
tables, checksums, figures, and session information). Completed models and
conditional fits are checkpointed. Cache fingerprints include inputs, source
scripts, settings, helper implementation, and package version. Move a result
directory aside before intentionally changing its setup.

Figure functions live separately in `reproduction_plots.R`, so plot edits do
not invalidate numerical fits. The emulation figures match the original
palette, colour limits, layouts, themes, and export dimensions/resolutions.
Calibration also uses the original histogram fills and variance/discrepancy/
field colours. Its `plot_heat_y_panel.png` follows the original 3×3 layout,
using GP means at the componentwise parameter quantile vectors instead of
PDE evaluations. Re-rendering uses the saved fit and does not repeat sampling.

A resumed emulation run reuses completed models. For the full uninterrupted
seed-1234 random stream, start with an empty emulation result directory;
unfinished models in a resumed run can use a different Monte Carlo stream.

Emulation retains the 100×100 grid, 51 times, 40/10 split, ten coefficient
draws, and all three covariance models. It is a substantial computation.
Calibration retains the 6×6 grid, 26 times, and 40/10 split, with
**30,000 emulator draws in full** (`nsam <- 3 * 10^4`) or 300 in quick,
three estimated parameters (eta1, eta2, alpha2), alpha1 and alpha3 fixed at zero,
and the unnormalised discrepancy covariance. It uses the
package's one-transition-per-emulator default, so its quantiles are execution
summaries rather than evidence of posterior convergence.

Full outputs use `_results/calibration-full/`; quick outputs remain in
`_results/calibration/`. The full run was launched in the background on
2026-10-25. Its current status is in `_results/calibration-full/job.json`,
with progress in `run.log` and per-draw files in `conditionals/`.
The detached wrapper records completion or failure and keeps the Mac awake
while the computation runs. The measured quick run suggests roughly 12–14
hours for 30,000 conditional fits, plus preparation,
validation, recovery, and saving. Increasing the number of emulator draws does
not add warmup or additional transitions within each conditional fit.

The original script generates the synthetic field after drawing emulator
coefficients. The profiles retain that RNG sequence, so changing the draw
count also changes the synthetic field realization and site permutation.
Both use the same saved simulator truth and generative settings.

The previous five-parameter calibration and reports are preserved under
`_archive/five-parameter-run/`.

The research folder is excluded by `.Rbuildignore`: routine installation does
not rerun experiments or bundle their external data. Ordinary package vignettes
remain in the parent `vignettes/` directory.
