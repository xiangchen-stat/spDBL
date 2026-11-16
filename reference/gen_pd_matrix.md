# Generate a random positive definite matrix

Constructs a positive definite matrix of dimension `dim` by forming
\\A^\top A\\ from a random normal matrix \\A\\ and then rescaling so the
maximum entry is 1.

## Usage

``` r
gen_pd_matrix(dim)
```

## Arguments

- dim:

  Integer. Number of rows (and columns) of the output matrix.

## Value

A symmetric positive definite numeric matrix of dimension `c(dim, dim)`.

## Examples

``` r
set.seed(1)
S <- gen_pd_matrix(4)
isSymmetric(S)
#> [1] TRUE
max(S)
#> [1] 1
all(eigen(S, symmetric = TRUE)$values > 0)
#> [1] TRUE
```
