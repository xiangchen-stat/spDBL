# Generate a flexible block grid with snake traversal

Iteratively computes block indices using a snake (boustrophedon)
traversal order. Handles cases where block dimensions do not evenly
divide file dimensions by allowing the last block in each row or column
to be smaller.

## Usage

``` r
generate.grid.rowsnake(fnrow, fncol, bnrow, bncol)
```

## Arguments

- fnrow:

  Integer. Number of rows in the file or matrix.

- fncol:

  Integer. Number of columns in the file or matrix.

- bnrow:

  Integer. Number of rows per block.

- bncol:

  Integer. Number of columns per block.

## Value

A `data.frame` with columns `block_no`, `U`, `D`, `L`, `R`. See
[`generate_grid`](https://xiangchen-stat.github.io/spDBL/reference/generate_grid.md)
for column descriptions.

## See also

[`generate_grid`](https://xiangchen-stat.github.io/spDBL/reference/generate_grid.md),
[`generate.grid.exact`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.exact.md),
[`generate.grid.lr`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.lr.md)

## Examples

``` r
# Last block row and column are smaller; even block rows run right to left
generate.grid.rowsnake(5, 7, 2, 3)
#>   block_no U D L R
#> 1        1 1 2 1 3
#> 2        2 1 2 4 6
#> 3        3 1 2 7 7
#> 4        4 3 4 7 7
#> 5        5 3 4 4 6
#> 6        6 3 4 1 3
#> 7        7 5 5 1 3
#> 8        8 5 5 4 6
#> 9        9 5 5 7 7
```
