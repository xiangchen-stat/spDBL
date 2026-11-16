# Compute median and 95% credible interval across rows

For each row of the matrix `X`, computes the median and the 2.5th and
97.5th percentiles across columns (samples).

## Usage

``` r
cal_errorbar(X)
```

## Arguments

- X:

  Numeric matrix. Rows correspond to spatial locations or parameters;
  columns correspond to posterior samples.

## Value

A data frame with columns:

- med:

  Row-wise median.

- lower:

  Row-wise 2.5th percentile.

- upper:

  Row-wise 97.5th percentile.

## See also

[`cal_errorbar_mean`](https://xiangchen-stat.github.io/spDBL/reference/cal_errorbar_mean.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(3 * 200, mean = 1:3), nrow = 3)
cal_errorbar(X)
#>         med       lower    upper
#> 1 0.9506307 -1.11860766 2.874439
#> 2 2.0673627 -0.03490053 3.766650
#> 3 2.8427249  1.01037797 5.089122
```
