# Out-of-core backward smoothing

Runs the backward smoothing pass over the filtered moments that
[`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md)
wrote to disk, one time step at a time, and writes the smoothed moments
to CSV files. The result equals
[`BS()`](https://xiangchen-stat.github.io/spDBL/reference/BS.md) applied
to the in-memory filter over the same blocks, without holding all
episodes in memory.

## Usage

``` r
BS_bigdata_R(
  ff_files,
  G_ls,
  nT,
  fnrow,
  fncol,
  bnrow,
  bncol,
  path_out,
  delta = 1,
  verbose = FALSE
)
```

## Arguments

- ff_files:

  The list of file paths returned by
  [`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md).

- G_ls:

  The state transition file path(s) given to
  [`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md):
  one path, or one per time step.

- nT, fnrow, fncol, bnrow, bncol:

  As given to
  [`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md).

- path_out:

  Directory where the output CSV files are written.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

- verbose:

  Logical. If `TRUE`, reports progress.

## Value

A list of file paths with one file per time step:

- st_ls:

  Smoothed state means; the \\p \times C_b\\ blocks are bound by column
  in traversal order, as in the filter's `mt` files.

- St_ls:

  Smoothed left covariances, \\p \times p\\ blocks bound by column.

## See also

[`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md),
[`FFBS_predict_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_bigdata_R.md),
[`BS()`](https://xiangchen-stat.github.io/spDBL/reference/BS.md)

## Examples

``` r
set.seed(1)
dir <- file.path(tempdir(), "bs_big")
dir.create(dir, showWarnings = FALSE)
csv <- function(x, name) {
  path <- file.path(dir, name)
  utils::write.csv(x, path, row.names = FALSE)
  path
}
# 5 runs x 4 sites over 3 times, in two blocks of 2 sites (AR(1) designs)
Y <- lapply(1:3, function(t) matrix(rnorm(20, t), 5, 4))
F <- c(list(matrix(1, 5, 4)), Y[1:2])
ff <- FF_bigdata_R(Y_ls = mapply(csv, Y, paste0("Y", 1:3, ".csv")),
  F_ls = mapply(csv, F, paste0("F", 1:3, ".csv")),
  G_ls = csv(diag(2), "G.csv"), W_ls = csv(diag(2), "W.csv"),
  V_ls = csv(diag(5), "V.csv"), m0 = matrix(0, 2, 2), M0 = diag(2),
  n0 = 4, D0 = diag(2), nT = 3, fnrow = 5, fncol = 4, bnrow = 5,
  bncol = 2, path_out = dir)
bs <- BS_bigdata_R(ff, file.path(dir, "G.csv"), nT = 3, fnrow = 5,
  fncol = 4, bnrow = 5, bncol = 2, path_out = dir)
utils::read.csv(bs$st_ls[1])  # smoothed means at time 1, both blocks
#>          V1        V2         V3        V4
#> 1 0.2511867 0.5446719 0.09863948 0.6876736
#> 2 0.8224550 0.5645969 0.95554192 0.7175611
unlink(dir, recursive = TRUE)
```
