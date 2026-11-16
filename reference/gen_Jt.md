# Compute the cross-covariance matrix between observed and new locations

Evaluates the exponential GP kernel between observed and new locations,
giving the \\N \times N\_{\text{new}}\\ cross-covariance block \\J_t\\.
The nugget does not enter a cross-covariance, even when a new location
coincides with an observed one; this matches
[`FFBS_predict_exact`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md).

## Usage

``` r
gen_Jt(input, input_new, phi_para, gp_sigma2 = 1.1, gp_tau2 = 10^(-4))
```

## Arguments

- input:

  Numeric matrix. Coordinates of observed locations (\\N \times d\\).

- input_new:

  Numeric matrix or vector. Coordinates of new locations
  (\\N\_{\text{new}} \times d\\).

- phi_para:

  Numeric scalar. Range parameter of the exponential kernel.

- gp_sigma2:

  Numeric scalar. Marginal variance of the kernel. Defaults to `1.1`.

- gp_tau2:

  Numeric scalar. Nugget variance. Not used; kept for compatibility.
  Defaults to `1e-4`.

## Value

Numeric matrix of dimension `c(N, N_new)`.

## See also

[`gen_exp_kernel`](https://xiangchen-stat.github.io/spDBL/reference/gen_exp_kernel.md)

## Examples

``` r
set.seed(1)
input <- matrix(runif(10), nrow = 5, ncol = 2)
J <- gen_Jt(input, input_new = c(0.5, 0.5), phi_para = 2)
length(J)
#> [1] 5
# At an observed input the cross-covariance equals gp_sigma2
gen_Jt(input, input[1, ], phi_para = 2)[1]
#>   1 
#> 1.1 
```
