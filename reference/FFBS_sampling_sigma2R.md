# Draw posterior samples from FFBS output (scalar sigma-squared-times-R model)

Draws `nsam` posterior samples of the state \\\Theta_t\\ and the scalar
variance \\\sigma^2\\ using the MNIG sampler, given the FFBS output
under the scalar-sigma model.

## Usage

``` r
FFBS_sampling_sigma2R(
  nsam,
  para_ffbs,
  F_ls,
  G_ls,
  nT,
  R,
  delta = 1,
  verbose = FALSE,
  seed = NULL,
  sampling = c("marginal", "joint"),
  keep_draws = TRUE
)
```

## Arguments

- nsam:

  Integer. Number of posterior samples.

- para_ffbs:

  List. Output of
  [`FFBS_sigma2R`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sigma2R.md).

- F_ls:

  Covariate matrix or list (see
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)).

- G_ls:

  State transition matrix or list (see
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)).

- nT:

  Integer. Number of time steps.

- R:

  Numeric matrix. Fixed spatial correlation matrix.

- delta:

  Discount factor. Must be 1: discounting is not implemented. Defaults
  to `1`.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

- seed:

  Optional local seed. NULL uses the session RNG.

- sampling:

  "marginal" (default) draws smoothed episode marginals conditional on a
  shared variance. "joint" preserves cross-episode dependence by
  sampling backward conditionals.

- keep_draws:

  If TRUE (default), store every episode's coefficient draws. If FALSE
  (marginal sampling only), store the random-number state at the start
  of each episode's draws instead; spDBL's prediction and calibration
  functions then regenerate identical draws on demand, which cuts the
  memory of large emulators. Indexing the returned draws by episode
  (`[[e]]`) then fails.

## Value

A named list of length `nT + 1`. Elements `"T1"` through `"T<nT>"` are
arrays of dimension `c(p, q, nsam)` with state samples. The element
`Sigma` is a vector of length `nsam` with samples of \\\sigma^2\\.
Attributes sampling_method and joint_trajectories identify the chosen
law.

## See also

[`FFBS_sigma2R`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sigma2R.md),
[`MNIG_sampler`](https://xiangchen-stat.github.io/spDBL/reference/MNIG_sampler.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
R <- 0.5^abs(outer(1:3, 1:3, "-"))
res <- FFBS_sigma2R(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
                    n0 = 3, D0 = 1, nT = 4, R = R)
sam <- FFBS_sampling_sigma2R(500, res, par$F0, par$G0, nT = 4, R = R,
                             seed = 1)
dim(sam$T2)
#> [1]   2   3 500
summary(sam$Sigma)  # draws of sigma^2
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1039  0.1638  0.1848  0.1880  0.2075  0.3099 
round(apply(sam$T2, 1:2, mean) - res$bs$st[, , 2], 2)
#>      [,1] [,2] [,3]
#> [1,] 0.01    0    0
#> [2,] 0.00    0    0
```
