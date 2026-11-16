# WAIC from paired coefficient and covariance draws

Evaluates row-wise multivariate normal log likelihoods for paired
posterior draws, using log-sum-exp aggregation for the predictive
density and sample variances for the WAIC penalty.

## Usage

``` r
compute_sample_waic(Y, X, B_samples, Sigma_samples, return_log_lik = FALSE)
```

## Arguments

- Y:

  Finite numeric response matrix with dimensions n by q.

- X:

  Finite numeric design matrix with dimensions n by p. Include an
  intercept column explicitly when required.

- B_samples:

  Finite numeric array with dimensions p by q by S, containing posterior
  coefficient draws. At least two draws are required.

- Sigma_samples:

  Finite numeric array with dimensions q by q by S, containing symmetric
  positive-definite posterior response covariance draws. Slice s must be
  paired with slice s of `B_samples`.

- return_log_lik:

  Logical scalar. If `TRUE`, also return the S by n `log_lik` matrix,
  suitable for
  [`loo::waic()`](https://mc-stan.org/loo/reference/waic.html) or
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html). The loo
  package is optional and is not called by this function.

## Value

The same summaries and pointwise matrix as
[`compute_analytical_waic()`](https://xiangchen-stat.github.io/spDBL/reference/compute_analytical_waic.md),
with an additional `log_lik` matrix when `return_log_lik = TRUE`.

## Details

Each response row is one observation, with likelihood \\Y_i\mid
B,\Sigma\sim N_q(x_i B,\Sigma)\\ independently across rows. The draws
can come from any posterior with this likelihood; they need not follow
an MNIW law. The analytical comparison requires an MNIW posterior.
Arrays retain three dimensions even when p, q, or n equals one. The
penalty uses sample variance with denominator S - 1, matching
[`loo::waic()`](https://mc-stan.org/loo/reference/waic.html) for a
matrix of log likelihoods. Correlated draws can increase Monte Carlo
error; this function does not diagnose convergence.

## See also

[`compute_analytical_waic()`](https://xiangchen-stat.github.io/spDBL/reference/compute_analytical_waic.md)

## Author

Sudipto Banerjee (original development implementation).

## Examples

``` r
X <- matrix(1, 2, 1)
Y <- matrix(c(0.1, 0.3), 2, 1)
B <- array(c(0, 0.2, 0.4), c(1, 1, 3))
Sigma <- array(c(0.5, 1, 1.5), c(1, 1, 3))
score <- compute_sample_waic(Y, X, B, Sigma, return_log_lik = TRUE)
score$WAIC
#> [1] 3.759926
dim(score$log_lik)
#> [1] 3 2
```
