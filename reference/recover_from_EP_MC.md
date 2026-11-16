# Recover episode-partitioned posterior samples to original time dimension

Reassembles a flattened episode-partitioned list of posterior sample
arrays back into a list indexed by the original `nT_ori` time steps.

## Usage

``` r
recover_from_EP_MC(dat_EP, nT_ori, nT, nsam)
```

## Arguments

- dat_EP:

  List of length `nT`. Each element is an array of dimension
  `c(p, bncol, nsam)`.

- nT_ori:

  Integer. Number of original time steps.

- nT:

  Integer. Total number of episode time steps.

- nsam:

  Integer. Number of posterior samples.

## Value

A named list of length `nT_ori`. Each element is an array of dimension
`c(p, bncol * n_block, nsam)`.

## See also

[`recover_from_EP_exact`](https://xiangchen-stat.github.io/spDBL/reference/recover_from_EP_exact.md)

## Examples

``` r
set.seed(1)
# 2 original times x 3 episodes; each episode holds 4 draws of a 2 x 2 block
ep <- lapply(1:6, function(e) array(rnorm(2 * 2 * 4), c(2, 2, 4)))
rec <- recover_from_EP_MC(ep, nT_ori = 2, nT = 6, nsam = 4)
names(rec)
#> [1] "T1" "T2"
dim(rec$T1)
#> [1] 2 6 4
all.equal(rec$T2[, 3:4, 1], ep[[5]][, , 1])
#> [1] TRUE
```
