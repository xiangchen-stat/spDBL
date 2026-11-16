# Generate a random matrix with entries scaled to \\\[-1, 1\]\\

Draws a matrix of independent standard normal entries and divides by the
largest absolute entry, so all entries lie in \\\[-1, 1\]\\. (Versions
up to 1.1.1 divided by the absolute value of the largest entry, so
entries could fall below -1.)

## Usage

``` r
gen_ran_matrix(nrow, ncol)
```

## Arguments

- nrow:

  Integer. Number of rows.

- ncol:

  Integer. Number of columns.

## Value

A numeric matrix of dimension `c(nrow, ncol)`.

## Examples

``` r
set.seed(1)
M <- gen_ran_matrix(3, 4)
dim(M)
#> [1] 3 4
max(abs(M))  # the largest absolute entry is scaled to 1
#> [1] 1
```
