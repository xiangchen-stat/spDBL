# Compute mean and 95% credible interval across rows

For each row of the matrix `X`, computes the mean and the 2.5th and
97.5th percentiles across columns (samples).

## Usage

``` r
cal_errorbar_mean(X)
```

## Arguments

- X:

  Numeric matrix. Rows correspond to spatial locations or parameters;
  columns correspond to posterior samples.

## Value

A data frame with columns:

- med:

  Row-wise mean.

- lower:

  Row-wise 2.5th percentile.

- upper:

  Row-wise 97.5th percentile.

## See also

[`cal_errorbar`](https://xiangchen-stat.github.io/spDBL/reference/cal_errorbar.md)

## Examples

``` r
set.seed(1)
X <- matrix(rexp(3 * 200, rate = 1:3), nrow = 3)
cal_errorbar_mean(X)
#>         med      lower    upper
#> 1 1.0183048 0.01460505 4.384868
#> 2 0.4991048 0.01811867 1.699184
#> 3 0.3217754 0.01661676 1.004608
```
