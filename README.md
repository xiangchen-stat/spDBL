# spDBL <a href="https://github.com/xiangchen-stat/spDBL"><img src="man/figures/logo.png" align="right" height="138" alt="spDBL hex logo" /></a>

<!-- badges: start -->
[![CRAN](https://img.shields.io/badge/CRAN-spDBL-blue)](https://cran.r-project.org/web/packages/spDBL/index.html)
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

### Install from GitHub

```r
# install.packages("pak")
pak::pak("xiangchen-stat/spDBL")
```

You can also install with `devtools`:

```r
# install.packages("devtools")
devtools::install_github("xiangchen-stat/spDBL")
```

### Install from source

To install from a source tarball, download the package source file and run:

```r
install.packages("spDBL_X.X.X.tar.gz", repos = NULL, type = "source")
```

Since `spDBL` contains compiled code, you will need a working C/C++ toolchain on your system. On macOS, this typically means installing Xcode Command Line Tools. On Windows, use Rtools. On Linux, install the corresponding system compilers and development libraries recommended for R package compilation.

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

<!-- ## Quick example

```r
library(spDBL)

# Example workflow
# sim_data <- simulate_mechanistic_model(n_runs = 50)
# gp_model <- fit_gp_emulator(sim_data)
# fit <- run_inference(gp_model, observations)

# summary(fit)
```

Please replace the example above with package-specific functions as your API stabilizes. -->

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
Prediction and calibration need no new PDE solves. The GP design interpolation
uncertainty remains excluded; emulator coefficient and predictive uncertainty
are retained.

    library(spDBL)
    data("dt_calibration")
    d <- dt_calibration
    em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
      covariance = "sigma2R", episode_window = 2, nsam = 20, seed = 1234)
    cal <- calibration_prepare(em, d$lower, d$upper)
    validation <- calibration_validate(cal, d$validation$input,
                                       d$validation$responses)
    # Optional RStan backend; requires a compiler on first use.
    # Defaults: marginal emulator sampling, one calibration iteration per draw.
    fit <- calibration_fit(cal, d$field$z, d$field$U)
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

The calibration defaults are chains = 1, iter = 1, warmup = 0, and one retained
sample per emulator draw. Larger iter values default to half the iterations as
warmup and retain up to 200 available samples per emulator draw. The default
sampler controls disable adaptation for zero warmup and enable it when warmup
is positive. One-transition results do not establish convergence.

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

Maintain .Rmd vignette sources; Rscript data-raw/sync-vignettes.R generates
their Quarto counterparts. To render the calibration vignette from this source
checkout, run this command from the package directory:

```sh
bash vignettes/renderEmulation.sh quarto calibration
```

The helper installs the current source into the ignored .r-library directory,
then renders calibration.qmd. This avoids using an older system installation
that lacks dt_calibration and the calibration functions. Subsequent direct
Quarto renders use the same local library. Omit calibration to render all three
vignettes, or select r instead of quarto for R Markdown.

Set SPDBL_RUN_STAN_TESTS=true with the package's
testthat suite to enable compiler-dependent integration checks. The ordinary
unit suite does not require RStan.

## Repository structure

```text
R/          R functions
src/        C++ source code
man/        documentation files
tests/      tests
```

## License

This package is released under the **MIT License**.

## Acknowledgement

The hex logo for spDBL was generated with the assistance of ChatGPT (OpenAI) using AI-based image generation tools.
