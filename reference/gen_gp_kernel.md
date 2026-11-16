# Compute a Gaussian Process covariance kernel matrix

Evaluates the isotropic GP exponential kernel \\k(s, s') = \sigma^2
\exp(-\phi \\s - s'\\) + \tau^2 \mathbf{1}\_{s=s'}\\ for all pairs of
rows in `loc`.

## Usage

``` r
gen_gp_kernel(loc, phi, sigma2, tau2)
```

## Arguments

- loc:

  Numeric matrix. Spatial locations (\\n \times d\\).

- phi:

  Numeric scalar. Range (decay) parameter of the exponential kernel.

- sigma2:

  Numeric scalar. Marginal variance parameter.

- tau2:

  Numeric scalar. Nugget (white-noise) variance.

## Value

A symmetric positive definite matrix of dimension `c(n, n)`.

## See also

[`gen_exp_kernel`](https://xiangchen-stat.github.io/spDBL/reference/gen_exp_kernel.md),
[`gen_expsq_kernel`](https://xiangchen-stat.github.io/spDBL/reference/gen_expsq_kernel.md)

## Examples

``` r
loc <- matrix(c(0, 1, 3), ncol = 1)
K <- gen_gp_kernel(loc, phi = 1, sigma2 = 2, tau2 = 0.1)
K
#>            1         2          3
#> 1 2.10000000 0.7357589 0.09957414
#> 2 0.73575888 2.1000000 0.27067057
#> 3 0.09957414 0.2706706 2.10000000
all.equal(K[1, 2], 2 * exp(-1))
#> [1] TRUE
```
