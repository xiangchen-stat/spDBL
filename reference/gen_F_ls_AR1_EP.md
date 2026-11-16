# Build AR(1) covariate list for the episode-block model

Generates a flattened covariate list for the episode-partition (EP)
model with AR(1) lags. The output list is indexed by \\(t-1) \times
n_b + j\\, where \\j\\ is the block index and \\t\\ is the time index.

## Usage

``` r
gen_F_ls_AR1_EP(Y, nT, n_b, ind)
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

A list of length `nT * n_b` with block-wise AR(1) covariate matrices.

## See also

[`gen_F_ls_AR1`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR1.md),
[`gen_F_ls_AR2_EP`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR2_EP.md)

## Examples

``` r
# Three time steps on six sites, split into three blocks of two sites
Y <- lapply(1:3, function(t) matrix(10 * t + 1:6, nrow = 1))
ind <- generate_grid(1, 6, 1, 2)
F_ep <- gen_F_ls_AR1_EP(Y, nT = 3, n_b = 3, ind = ind)
length(F_ep)  # nT * n_b
#> [1] 9
F_ep[[5]]     # time 2, block 2: sites 3-4 at time 1
#>      [,1] [,2]
#> [1,]   13   14
```
