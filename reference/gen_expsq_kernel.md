# Compute a squared-exponential (Gaussian) GP kernel matrix

Evaluates the squared-exponential kernel \\k(s, s') = \sigma^2
\exp(-\phi \\s - s'\\^2) + \tau^2 \mathbf{1}\_{s=s'}\\ for all pairs of
rows in `loc`.

## Usage

``` r
gen_expsq_kernel(loc, phi, sigma2 = 1, tau2 = 0)
```

## Arguments

- loc:

  Numeric matrix. Spatial locations (\\n \times d\\).

- phi:

  Numeric scalar. Inverse length-scale parameter.

- sigma2:

  Numeric scalar. Marginal variance. Defaults to `1`.

- tau2:

  Numeric scalar. Nugget variance. Defaults to `0`.

## Value

A symmetric positive definite matrix of dimension `c(n, n)`.

## See also

[`gen_gp_kernel`](https://xiangchen-stat.github.io/spDBL/reference/gen_gp_kernel.md),
[`gen_exp_kernel`](https://xiangchen-stat.github.io/spDBL/reference/gen_exp_kernel.md)

## Examples

``` r
loc <- cbind(c(0, 1, 0), c(0, 0, 2))
K <- gen_expsq_kernel(loc, phi = 0.5, tau2 = 0.01)
round(K, 3)
#>       1     2     3
#> 1 1.010 0.607 0.135
#> 2 0.607 1.010 0.082
#> 3 0.135 0.082 1.010
all.equal(K[1, 3], exp(-0.5 * 2^2))
#> [1] TRUE
```
