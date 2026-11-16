# Build AR(2) covariate list for the episode-block model

Generates a flattened covariate list for the episode-partition (EP)
model with AR(2) lags. Both lag-1 and lag-2 column blocks are extracted
and column-bound for each spatial block.

## Usage

``` r
gen_F_ls_AR2_EP(Y, nT, n_b, ind)
```

## Arguments

- Y:

  List of length `nT`. Each element is the full spatial data matrix at
  time \\t\\.

- nT:

  Integer. Number of time steps.

- n_b:

  Integer. Number of spatial blocks.

- ind:

  Data frame of block indices as returned by
  [`generate_grid`](https://xiangchen-stat.github.io/spDBL/reference/generate_grid.md).

## Value

A list of length `nT * n_b` with block-wise AR(2) covariate matrices
(each of width \\2 \times \text{bncol}\\).

## See also

[`gen_F_ls_AR2`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR2.md),
[`gen_F_ls_AR1_EP`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR1_EP.md)

## Examples

``` r
set.seed(1)
Y <- lapply(1:3, function(t) matrix(10 * t + 1:6, nrow = 1))
ind <- generate_grid(1, 6, 1, 2)
F_ep <- gen_F_ls_AR2_EP(Y, nT = 3, n_b = 3, ind = ind)
length(F_ep)
#> [1] 9
F_ep[[8]]  # time 3, block 2: sites 3-4 at times 1 and 2
#>      [,1] [,2] [,3] [,4]
#> [1,]   13   14   23   24
```
