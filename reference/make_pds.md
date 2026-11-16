# Force a matrix to be positive definite and symmetric

Symmetrises a matrix by averaging it with its transpose and, if the
result is not positive definite, adds a small diagonal ridge.

## Usage

``` r
make_pds(C, eps = 10^(-4))
```

## Arguments

- C:

  Numeric matrix.

- eps:

  Numeric scalar. Size of the diagonal ridge added when `C` is not
  positive definite. Defaults to `1e-4`.

## Value

The corrected positive definite symmetric matrix.

## See also

[`check_pds`](https://xiangchen-stat.github.io/spDBL/reference/check_pds.md)

## Examples

``` r
A <- matrix(c(2, 0.9, 1.1, 2), 2)
make_pds(A)  # averaged with its transpose
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    1    2
isSymmetric(make_pds(A))
#> [1] TRUE
```
