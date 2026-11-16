# Generate an exact block grid analytically

Computes block indices when block dimensions divide file dimensions
exactly. Uses vectorised operations rather than a loop, making it fast
for large grids.

## Usage

``` r
generate.grid.exact(fnrow, fncol, bnrow, bncol, traversal.mode)
```

## Arguments

- fnrow:

  Integer. Number of rows in the file or matrix.

- fncol:

  Integer. Number of columns in the file or matrix.

- bnrow:

  Integer. Number of rows per block. Must divide `fnrow` exactly.

- bncol:

  Integer. Number of columns per block. Must divide `fncol` exactly.

- traversal.mode:

  Character. `"rowsnake"` or `"lr"`. See
  [`generate_grid`](https://xiangchen-stat.github.io/spDBL/reference/generate_grid.md)
  for details.

## Value

A `data.frame` with columns `block_no`, `U`, `D`, `L`, `R`. See
[`generate_grid`](https://xiangchen-stat.github.io/spDBL/reference/generate_grid.md)
for column descriptions.

## See also

[`generate_grid`](https://xiangchen-stat.github.io/spDBL/reference/generate_grid.md),
[`generate.grid.rowsnake`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.rowsnake.md),
[`generate.grid.lr`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.lr.md)

## Examples

``` r
# Block sizes divide the grid exactly
generate.grid.exact(4, 6, 2, 3, traversal.mode = "rowsnake")
#>   block_no U D L R
#> 1        1 1 2 1 3
#> 2        2 1 2 4 6
#> 3        3 3 4 4 6
#> 4        4 3 4 1 3
generate.grid.exact(4, 6, 2, 3, traversal.mode = "lr")
#>   block_no U D L R
#> 1        1 1 2 1 3
#> 2        2 1 2 4 6
#> 3        3 3 4 1 3
#> 4        4 3 4 4 6
```
