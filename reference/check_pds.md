# Check and repair a matrix to be positive definite and symmetric

Verifies that a matrix is approximately symmetric (within a relative and
absolute tolerance). If so, symmetrises it and, if still not positive
definite, adds a small diagonal ridge. Raises an error if the asymmetry
exceeds the tolerance.

## Usage

``` r
check_pds(C, eps = 10^(-5), per = 0.05)
```

## Arguments

- C:

  Numeric matrix to check.

- eps:

  Numeric scalar. Absolute tolerance for the Frobenius norm of
  `C - t(C)`. Defaults to `1e-5`.

- per:

  Numeric scalar. Relative tolerance: the ratio of the asymmetry norm to
  the mean of `C`. Defaults to `0.05`.

## Value

The (possibly corrected) positive definite symmetric matrix.

## See also

[`make_pds`](https://xiangchen-stat.github.io/spDBL/reference/make_pds.md)

## Examples

``` r
# A covariance matrix with a tiny rounding asymmetry
A <- matrix(c(2, 1, 1 + 1e-8, 2), 2)
isSymmetric(A)
#> [1] FALSE
A_fixed <- check_pds(A)
isSymmetric(A_fixed)
#> [1] TRUE
A_fixed
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    1    2
```
