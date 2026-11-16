# Compute an exponential GP kernel matrix

Evaluates the exponential kernel \\k(s, s') = \sigma^2 \exp(-\phi \\s -
s'\\) + \tau^2 \mathbf{1}\_{s=s'}\\ for all pairs of rows in `loc`.

## Usage

``` r
gen_exp_kernel(loc, phi, sigma2 = 1, tau2 = 0)
```

## Arguments

- loc:

  Numeric matrix. Spatial locations (\\n \times d\\).

- phi:

  Numeric scalar. Range (decay) parameter.

- sigma2:

  Numeric scalar. Marginal variance. Defaults to `1`.

- tau2:

  Numeric scalar. Nugget variance. Defaults to `0`.

## Value

A symmetric positive definite matrix of dimension `c(n, n)`.

## See also

[`gen_gp_kernel`](https://xiangchen-stat.github.io/spDBL/reference/gen_gp_kernel.md),
[`gen_expsq_kernel`](https://xiangchen-stat.github.io/spDBL/reference/gen_expsq_kernel.md)

## Examples

``` r
loc <- cbind(c(0, 1, 0), c(0, 0, 2))
K <- gen_exp_kernel(loc, phi = 0.5)
round(K, 3)
#>       1     2     3
#> 1 1.000 0.607 0.368
#> 2 0.607 1.000 0.327
#> 3 0.368 0.327 1.000
all.equal(K[1, 3], exp(-0.5 * 2))
#> [1] TRUE
```
