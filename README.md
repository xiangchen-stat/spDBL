# spDBL <a href="https://github.com/xiangchen-stat/spDBL"><img src="man/figures/logo.png" align="right" height="138" alt="spDBL hex logo" /></a>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/spDBL)](https://CRAN.R-project.org/package=spDBL)
[![R-CMD-check](https://github.com/xiangchen-stat/spDBL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xiangchen-stat/spDBL/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/xiangchen-stat/spDBL/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/xiangchen-stat/spDBL/actions/workflows/test-coverage.yaml)
[![GitHub](https://img.shields.io/badge/GitHub-xiangchen--stat%2FspDBL-black?logo=github)](https://github.com/xiangchen-stat/spDBL)
[![Paper](https://img.shields.io/badge/JMLR-2025-blue)](https://jmlr.org/papers/v26/22-0896.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/xiangchen-stat/spDBL/blob/main/LICENSE.md)
<!-- badges: end -->

## Overview

**spDBL** is an R package for **Dynamic Bayesian Learning for Spatiotemporal Mechanistic Models**. It implements a scalable framework for **Gaussian process emulation** and **Bayesian inverse problems** in complex spatiotemporal dynamical systems, combining **mechanistic modeling** with **hierarchical state-space models**.

Mechanistic models arising from ordinary differential equations, partial differential equations, and computer simulators are often computationally expensive. Bayesian inverse problems seek to infer unknown mechanistic parameters from noisy observations and quantify their uncertainty through posterior distributions. `spDBL` makes this inference more tractable by emulating the mechanistic system and then melding emulator output with observations in a hierarchical Bayesian model.

The framework supports:

- **Gaussian process emulation** of mechanistic systems
- **Hierarchical state-space modeling** for spatiotemporal dynamics
- **Closed-form Bayesian updating** in matrix-variate Normal and Wishart models during emulator construction
- **Post-emulation MCMC** for Bayesian inverse problems and calibration of mechanistic parameters
- **Dynamic Bayesian transfer learning** for large-scale emulation

The methodological details are presented in [Banerjee, Chen, Frankenburg, and Zhou (2025)](https://jmlr.org/papers/v26/22-0896.html).

## Installation

### Release version from CRAN

```r
install.packages("spDBL")
```

CRAN provides binaries for Windows and macOS. On Linux the package is built from
source, which needs a C++ compiler.

### Development version from GitHub

```r
# install.packages("pak")
pak::pak("xiangchen-stat/spDBL")
```

You can also use `devtools::install_github("xiangchen-stat/spDBL")`.

### Optional: RStan for calibration

Emulation needs no extra software. `calibration_fit()` uses the suggested
package [rstan](https://mc-stan.org/rstan/), which compiles a Stan program on
first use (cached afterwards) and therefore needs a C++ toolchain:

```r
install.packages("rstan")
```

On macOS install the Xcode Command Line Tools, on Windows install Rtools, and
on Linux install the system compilers recommended for R package compilation.
See the [RStan Getting Started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
if compilation fails.

## Quick example

Fit a dynamic emulator to 40 runs of a spatial SIR simulator on a 10 × 10 grid
over 26 time steps, then predict 10 held-out runs:

```r
library(spDBL)
data("dt_emulation")
d <- dt_emulation

# Fit the emulator to the training simulator runs (10 x 10 grid)
em <- emulator_learn(d$pde_para_train, d$dt_pde_train, Nx = 10, Ny = 10,
                     seed = 1)

# Predict the held-out runs; their lagged outputs supply the AR design rows
pred <- emulator_predict(em, d$pde_para_test, dt_pde_test = d$dt_pde_test)

# Root mean squared error over held-out runs, locations, and times
sqrt(mean(unlist(Map(function(p, y) (p - y)^2, pred, d$dt_pde_test))))
#> [1] 0.1241679
```

The held-out outputs have standard deviation 3.17, so the error is small
relative to the signal. `vignette("emulation")` walks through this example with
plots and uncertainty intervals. The next section shows calibration, which
predicts without any new simulator runs.

## Package scope

`spDBL` is designed for Bayesian learning in settings where a mechanistic model governs spatiotemporal dynamics but repeated direct simulation is too costly for routine inference. Typical applications include:

- inverse problems for nonlinear ODE and PDE systems
- spatiotemporal epidemiologic models
- environmental and ecological dynamical systems
- black-box computer models generating spatial or spatiotemporal output
- transfer learning across related dynamical systems

## Methodological summary

The `spDBL` workflow has four main components:

1. **Mechanistic model runs**  
   Generate training runs from the underlying dynamical system over a finite collection of inputs.

2. **Gaussian process emulation**  
   Learn the relationship between model inputs and system outputs using Gaussian process regression.

3. **Hierarchical spatiotemporal modeling**  
   Embed the emulator in a state-space model to combine mechanistic knowledge with noisy observed data.

4. **Inference on model parameters**  
   Use analytically tractable posterior distributions for emulator construction, followed by MCMC-based inference for mechanistic parameters when needed.

A central advantage of this framework is that the emulator stage avoids expensive iterative Bayesian computation such as full MCMC or variational inference, while still preserving a coherent probabilistic learning framework.

## Related paper

**Dynamic Bayesian Learning for Spatiotemporal Mechanistic Models**  
Sudipto Banerjee, Xiang Chen, Ian Frankenburg, Daniel Zhou  
*Journal of Machine Learning Research*, 26(146):1-43, 2025

### Abstract

We develop an approach for Bayesian learning of spatiotemporal dynamical mechanistic models. Such learning consists of statistical emulation of the mechanistic system that can efficiently interpolate the output of the system from arbitrary inputs. The emulated learner can then be used to train the system from noisy data achieved by melding information from observed data with the emulated mechanistic system. This joint melding of mechanistic systems employ hierarchical state-space models with Gaussian process regression. Assuming the dynamical system is controlled by a finite collection of inputs, Gaussian process regression learns the effect of these parameters through a number of training runs, driving the stochastic innovations of the spatiotemporal state-space component. This enables efficient modeling of the dynamics over space and time.

This article details exact inference with analytically accessible posterior distributions in hierarchical matrix-variate Normal and Wishart models in designing the emulator. This step obviates expensive iterative algorithms such as Markov chain Monte Carlo or variational approximations. We also show how emulation is applicable to large-scale emulation by designing a dynamic Bayesian transfer learning framework. Inference on mechanistic model parameters proceeds using Markov chain Monte Carlo as a post-emulation step using the emulator as a regression component. We demonstrate this framework through solving inverse problems arising in the analysis of ordinary and partial nonlinear differential equations and, in addition, to a black-box computer model generating spatiotemporal dynamics across a graphical model.

Paper link: <https://jmlr.org/papers/v26/22-0896.html>

## Calibration from saved simulations

The calibration interface supports the scalar covariance model and uses an
offline GP to supply autoregressive design rows at new mechanistic inputs.
Prediction and calibration need no new PDE solves. The design GP's predictive
covariance is propagated into the field covariance, along with emulator
coefficient and predictive uncertainty.

    library(spDBL)
    data("dt_calibration")
    d <- dt_calibration
    em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
      covariance = "sigma2R", episode_window = 2, nsam = 20, seed = 1234)
    cal <- calibration_prepare(em, d$lower, d$upper)
    validation <- calibration_validate(cal, d$validation$input,
                                       d$validation$responses)
    # Optional RStan backend; requires a compiler on first use.
    # iter must be set; 4 chains x 2000 iterations per emulator draw for inference.
    fit <- calibration_fit(cal, d$field$z, d$field$U, chains = 4, iter = 2000)
    summary(fit)
    fields <- calibration_predict(fit)

The inverse-Wishart default of emulator_learn() is preserved. Calibration
requires an explicit sigma2R fit and does not convert covariance posteriors.
Both covariance models use sampling = "marginal" by default, drawing smoothed
episode marginals conditional on a shared covariance draw. Choose
sampling = "joint" for coefficient paths that preserve episode dependence.
calibration_prepare() inherits the fitted method; an explicit sampling option
changes methods using saved filtering results:

    cal_joint <- calibration_prepare(em, d$lower, d$upper, sampling = "joint")
    fit_long <- calibration_fit(cal_joint, d$field$z, d$field$U,
      chains = 4, iter = 2000, warmup = 1000, draws_per_emulator = 200)

`calibration_fit()` has no default for `iter`. Warmup defaults to half the
iterations, chains to one, and up to 200 available samples are retained per
emulator draw. The default sampler controls disable adaptation for zero warmup
and enable it when warmup is positive. `iter = 1` gives one transition per
emulator draw: a quick execution check that does not establish convergence.

emulator_predict(cal, input_new) also works without new PDE trajectories;
its mean samples are distinct from full conditional field draws.

See vignette("calibration", package = "spDBL") for the end-to-end example and
vignette("calibration-details", package = "spDBL") for the probability model.
Package builds run preparation and validation, then display a compact recorded
calibration result; they never compile Stan or launch MCMC. The calibration
vignette uses 500 emulator draws with one NUTS transition per draw, and shows
eta histograms and a variance/discrepancy/field comparison panel. Its recorded
sample quantiles describe an execution experiment, not converged inference.
Full reproduction is available through the installed
system.file("examples", "run_calibration.R", package = "spDBL") script.
Its default oneiter500 profile uses marginal emulator sampling. The small,
full, and smoke profiles remain available; an optional third argument selects
joint sampling, for example run_calibration.R full output_dir joint.
One-iteration runs are execution experiments and do not establish convergence.

## MNIW model comparison

For multivariate regression with independent observation rows,
`compute_analytical_waic()` evaluates WAIC directly from an MNIW posterior.
`compute_sample_waic()` evaluates the same row likelihood from paired posterior
coefficient and covariance draws. Both return `WAIC`, `lppd`, `p_waic`,
`elpd_waic`, and pointwise scores. The sampled function can also return a log
likelihood matrix for the optional `loo` package.

```r
X <- cbind(1, c(-1, 0, 1, 2))
Y <- cbind(c(0.1, 0.4, 1.2, 1.7), c(1.2, 0.8, 0.5, 0.1))
post <- spDBL::MNIW_R(X, Y, H = diag(nrow(Y)), v = 4, S = diag(2),
                     C = matrix(0, 2, 2), Vb = diag(2))
spDBL::compute_analytical_waic(Y, X, post$Vbnew, post$Cnew,
                              post$Snew, post$vnew)
```

`mn` is the posterior coefficient mean itself. Each multivariate response row
is scored jointly; these functions do not directly evaluate the correlated-row
GP likelihood used by the spatiotemporal emulator. See
`vignette("model-comparison", package = "spDBL")` for posterior sampling,
comparison with `loo`, and sequential covariate comparisons. To render it from
this checkout, use `bash vignettes/renderEmulation.sh quarto model-comparison`.

## Documentation

* Website with the function reference and all vignettes:
  <https://xiangchen-stat.github.io/spDBL/>
* `vignette("emulation")`: fit an emulator and check its predictions.
* `vignette("calibration")`: end-to-end calibration from saved simulations.
* `vignette("calibration-details")`: the probability model behind calibration.
* `vignette("model-comparison")`: WAIC for MNIW multivariate regression.
* [Package architecture](https://xiangchen-stat.github.io/spDBL/articles/architecture.html):
  how the code is organised, for contributors.
* [IMPROVEMENTS.md](https://github.com/xiangchen-stat/spDBL/blob/main/IMPROVEMENTS.md):
  what the package adds to the published method and its research code.

## Getting help and contributing

Ask questions and share how you use spDBL in
[GitHub Discussions](https://github.com/xiangchen-stat/spDBL/discussions), and
report bugs on the [issue tracker](https://github.com/xiangchen-stat/spDBL/issues),
with a minimal reproducible example. Contributions are welcome; the
[contributing guide](https://github.com/xiangchen-stat/spDBL/blob/main/.github/CONTRIBUTING.md) explains how to set up a
development copy, run the tests (including the optional Stan tests), and
render the vignettes. Please note that this project is released with a
[Contributor Code of Conduct](https://github.com/xiangchen-stat/spDBL/blob/main/.github/CODE_OF_CONDUCT.md). By participating you
agree to abide by its terms.

## Repository structure

```text
R/            R functions
src/          C++ kernels (Rcpp/RcppEigen)
inst/stan/    Stan programs for calibration
inst/examples/  reproduction script for the calibration vignette
tests/        testthat unit and integration tests
vignettes/    vignettes (R Markdown sources; Quarto copies are generated)
data/         example datasets (dt_emulation, dt_calibration)
data-raw/     scripts that build the example data and recorded results
man/          help pages generated by roxygen2
.github/      CI workflows, contributing guide, issue templates
```

## License

This package is released under the **MIT License** (see `LICENSE` and
`LICENSE.md`). Every hand-written source file states this with an
`SPDX-License-Identifier: MIT` header; the Rcpp-generated `RcppExports` files
are covered by the same license. All dependencies use OSI-approved open-source
licenses (GPL, MIT, or BSD), and the only optional dependency for
calibration, rstan, is GPL (>= 3).

## Acknowledgement

The hex logo for spDBL was generated with the assistance of ChatGPT (OpenAI) using AI-based image generation tools.
