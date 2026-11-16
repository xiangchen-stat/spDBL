# Recover episode-partitioned data to original time dimension (exact)

Reassembles a flattened episode-partitioned list (as produced by the EP
model) back into a list indexed by the original `nT_ori` time steps.
Assumes the number of blocks per season `n_block = nT / nT_ori` is an
integer.

## Usage

``` r
recover_from_EP_exact(dat_EP, nT_ori, nT, block_cols = NULL)
```

## Arguments

- dat_EP:

  List of length `nT`. Flattened episode-partitioned data.

- nT_ori:

  Integer. Number of original (non-partitioned) time steps.

- nT:

  Integer. Total number of episode time steps (`nT = nT_ori * n_block`).

- block_cols:

  Optional block-by-column mapping to the original grid.

## Value

A named list of length `nT_ori` (`"T1"`, ..., `"T<nT_ori>"`), each
element a matrix covering one original time step.

## See also

[`recover_from_EP_MC`](https://xiangchen-stat.github.io/spDBL/reference/recover_from_EP_MC.md)

## Examples

``` r
# Two time steps on six sites, split into three episodes of two sites each
Y <- lapply(1:2, function(t) matrix(10 * t + 1:6, nrow = 1))
ind <- generate_grid(1, 6, 1, 2)
Y_ep <- unlist(lapply(Y, function(y) lapply(1:3, function(j)
  y[, ind$L[j]:ind$R[j], drop = FALSE])), recursive = FALSE)
length(Y_ep)
#> [1] 6
rec <- recover_from_EP_exact(Y_ep, nT_ori = 2, nT = 6)
all.equal(unname(rec$T2), Y[[2]])
#> [1] TRUE
```
