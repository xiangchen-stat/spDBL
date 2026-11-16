# Exact posterior predictive mean using FFBS smoothed states (MNIW model)

Computes the analytical posterior predictive mean at new spatial
locations using the smoothed state means \\s_t\\ from
[`FFBS`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md) and
an exponential GP kernel for the cross-covariance.

## Usage

``` r
FFBS_predict_exact(
  Y,
  para_ffbs = NULL,
  F_ls,
  F_new_ls,
  input,
  input_new,
  nT,
  phi_para,
  gp_sigma2 = 1.1,
  gp_tau2 = 10^(-4),
  delta = 1,
  theta = NULL,
  V_chol = NULL,
  residual_weights = NULL,
  return_covariance = FALSE
)
```

## Arguments

- Y:

  List of length `nT`. Observed data matrices.

- para_ffbs:

  List. Output of
  [`FFBS`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md),
  containing `ff` and `bs`.

- F_ls:

  Covariate matrix or list for observed locations.

- F_new_ls:

  Covariate matrix or list for new locations.

- input:

  Numeric matrix. Coordinates of observed locations.

- input_new:

  Numeric matrix or vector. Coordinates of new locations.

- nT:

  Integer. Number of time steps.

- phi_para:

  Numeric scalar. Range parameter of the exponential kernel.

- gp_sigma2:

  Numeric scalar. Variance parameter of the GP kernel. Defaults to
  `1.1`.

- gp_tau2:

  Numeric scalar. Nugget variance of the GP kernel. Defaults to `1e-4`.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

- theta:

  Optional list of fixed coefficient matrices, one per episode.

- V_chol:

  Optional cached upper Cholesky factor of the input covariance.

- residual_weights:

  Optional cached inverse-covariance-weighted residuals for the supplied
  coefficient trajectory.

- return_covariance:

  If TRUE, return a list with mean and input_covariance before scaling
  by the emulator variance and spatial covariance.

## Value

A named list of length `nT`. Each element `"T<t>"` is a matrix of
dimension `c(N_new, q)` with the posterior predictive mean at the new
locations.

## See also

[`FFBS_predict_MC`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_MC.md),
[`gen_exp_kernel`](https://xiangchen-stat.github.io/spDBL/reference/gen_exp_kernel.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res <- FFBS(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
            par$n0, par$D0, nT = 4)
x_new <- matrix(c(-0.5, 0.5))  # two new inputs
pred <- FFBS_predict_exact(dat$Y, res, F_ls = par$F0,
                           F_new_ls = cbind(1, x_new), input = par$F0[, 2],
                           input_new = x_new, nT = 4, phi_para = 1)
pred$T4  # predictive mean at the new inputs (2 x S)
#>         [,1]       [,2]       [,3]
#> 6 -0.4612179 -1.5019804  0.2489372
#> 7 -0.4220251 -0.3272489 -0.1435978
```
