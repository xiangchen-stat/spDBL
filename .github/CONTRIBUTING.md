# Contributing to spDBL

Thank you for considering a contribution. This guide covers reporting
problems, setting up a development copy, and the checks a change should pass.
Please follow the [Code of Conduct](CODE_OF_CONDUCT.md) in all project spaces.

## Reporting bugs and asking questions

Ask usage questions, and share how you use spDBL, in
[GitHub Discussions](https://github.com/xiangchen-stat/spDBL/discussions).
Report bugs at <https://github.com/xiangchen-stat/spDBL/issues>. A minimal
reproducible example (`reprex::reprex()`) and your `sessionInfo()` make it
much easier to help. For calibration problems, include the rstan and
StanHeaders versions and your compiler.

## Before a larger change

For a new feature, model component, or change to an exported function, please
open an issue first so we can agree on the design. Small fixes (typos,
documentation, clear bugs) can go straight to a pull request.

The [architecture article](https://xiangchen-stat.github.io/spDBL/articles/architecture.html)
(`vignettes/articles/architecture.Rmd`) maps each workflow function to the
code that implements it and lists the main extension points.

## Development setup

You need R (≥ 4.0), a C++ compiler for the Rcpp/RcppEigen code, and, for
calibration work, rstan.

```r
# install.packages(c("pak", "devtools"))
pak::local_install_dev_deps()   # from the package root; installs Imports and Suggests
devtools::load_all()            # compile and load the development version
```

## Making a change

1. Create a branch from `main`.
2. Edit the code. Documentation lives in roxygen comments next to each
   function; run `devtools::document()` to regenerate `man/` and `NAMESPACE`.
   Never edit `man/*.Rd` or `NAMESPACE` by hand.
3. Add or update tests (see below).
4. Add a bullet to `NEWS.md` under the development version.
5. Run `devtools::check()`. It should finish with no errors or warnings.
6. Open a pull request; the checklist in the template repeats these steps.

### Conventions

* **License header.** Every new R, C++, Stan, or shell file starts with

  ```
  # Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
  # SPDX-License-Identifier: MIT
  ```

  (use `//` in C++ and Stan files). `bash .github/scripts/check-license-headers.sh`
  checks this, and so does CI.
* **Names.** New functions use `snake_case`. Kernel functions follow the
  covariance-model suffixes: none for inverse-Wishart, `_sigma2R`, and `_I`.
  Some older exported names (`generate.grid.*`) are kept for compatibility.
* **Randomness.** Functions that take a `seed` set it locally and restore the
  caller's random-number state on exit; reuse the internal `.cal_seed()`
  helper for this.
* **Files.** Examples, tests, and vignettes must not write outside
  `tempdir()`. Plot helpers that can save files (`savei = TRUE`) need a
  temporary path in examples and tests.
* **C++.** Edit `src/*.cpp`. `devtools::document()` runs
  `Rcpp::compileAttributes()`, which regenerates `R/RcppExports.R` and
  `src/RcppExports.cpp`; do not edit those two files.
* **Stan.** The programs are in `inst/stan/`. The R functions
  `calibration_log_density()` and `calibration_marginal_log_density()` (in
  `R/calibration_recover.R`) are reference implementations of the same
  densities, and the Stan tests compare them. Change both together.

## Tests

Tests are in `tests/testthat/` (testthat, 3rd edition).

```r
devtools::test()                        # unit tests; no Stan compilation
Sys.setenv(SPDBL_RUN_STAN_TESTS = "true")
devtools::test()                        # also compile and run the Stan programs
```

The Stan tests are skipped by default because they need a C++ compiler and
are much slower than the unit tests. Set `SPDBL_STAN_CACHE_DIR` to a directory to keep compiled
models between runs.

Prefer tests that compare against an independent calculation (a closed-form
posterior, direct Gaussian conditioning, an in-memory version of an
out-of-core computation) over tests that only check that code runs. Keep test
dimensions small: the default suite runs on CRAN. Mark slow tests with
`skip_on_cran()`.

Coverage:

```r
Sys.setenv(SPDBL_RUN_STAN_TESTS = "true")
covr::report(covr::package_coverage())   # opens an HTML report
```

GitHub Actions runs `R CMD check` on macOS, Windows, and Ubuntu (R devel,
release, and the previous release), measures coverage with the Stan tests on,
checks license headers, and rebuilds the website.

## Vignettes

The `.Rmd` files in `vignettes/` are the sources. The Quarto copies (`.qmd`)
are generated from them; never edit them directly:

```sh
Rscript data-raw/sync-vignettes.R
```

To render a vignette against the current source rather than an older
installed version, run from the package directory:

```sh
bash vignettes/renderEmulation.sh quarto calibration
```

The helper installs the current source into the ignored `.r-library`
directory, then renders `calibration.qmd`. Omit `calibration` to render all
three vignettes, or use `r` instead of `quarto` for R Markdown.

Package builds of the calibration vignette never compile Stan; they show a
recorded result. To reproduce it, run the installed script
`system.file("examples", "run_calibration.R", package = "spDBL")`, then
record the output with `Rscript data-raw/record-calibration.R <output_dir>`.

## Example data

Scripts in `data-raw/` rebuild the shipped data:

* `data-raw/sir-calibration.R` builds `dt_calibration` from saved SIR
  simulations (pass their directory as the argument).
* `data-raw/record-calibration.R` records the calibration vignette's result.

## Website

The website is built with pkgdown from `_pkgdown.yml`, the README, the help
pages, and the vignettes. To preview it locally, run `pkgdown::build_site()` in a copy
of the repository with `vignettes/*.qmd` (generated duplicates) and
`vignettes/reproduce/` (research notebooks that need the original simulation
data) removed; the GitHub Actions workflow does the same before deploying the
site from `main`.
