# Forward Filter for big data stored in CSV files (MNIW model)

Runs the forward filtering pass for large datasets that are partitioned
into spatial blocks and stored as CSV files on disk. At each time step
the spatial domain is traversed block-by-block using a snake traversal
order (see
[`generate.grid.rowsnake`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.rowsnake.md)),
and results are written to disk to avoid exhausting memory.

## Usage

``` r
FF_bigdata_R(
  Y_ls,
  F_ls,
  G_ls,
  W_ls,
  V_ls,
  m0,
  M0,
  n0,
  D0,
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

- Y_ls:

  Character vector of length `nT`. File paths to the response CSV files,
  one per time step.

- F_ls:

  Either a single file path (constant \\F_t\\) or a character vector of
  length `nT` (time-varying). Can also contain a second set of columns
  (AR2 lags) for the same file at offset `fncol`. Each design block is
  read with the same rows and columns as the response block, so the
  design must have as many columns as a block (`bncol`), as with the
  autoregressive designs.

- G_ls:

  Either a single file path (constant \\G_t\\) or a character vector of
  length `nT` (time-varying).

- W_ls:

  Either a single file path (constant \\W_t\\) or a character vector of
  length `nT` (time-varying).

- V_ls:

  Either a single file path (constant \\V_t\\) or a character vector of
  length `nT` (time-varying). Expected to be a square block submatrix
  indexed by the row range of each spatial block.

- m0:

  Numeric matrix. Prior state mean (\\p \times q\\).

- M0:

  Numeric matrix. Prior left-covariance of the state (\\p \times p\\).

- n0:

  Numeric scalar. Prior degrees of freedom of the inverse-Wishart.

- D0:

  Numeric matrix. Prior scale matrix of the inverse-Wishart.

- nT:

  Integer. Number of time steps (files).

- fnrow:

  Integer. Total number of rows in each data file.

- fncol:

  Integer. Total number of columns in each data file.

- bnrow:

  Integer. Number of rows per spatial block.

- bncol:

  Integer. Number of columns per spatial block.

- path_out:

  Character. Directory path where output CSV files are written.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

## Value

A named list of character vectors with paths to the output CSV files:

- nt_ls:

  Paths to filtered degrees-of-freedom files.

- Dt_ls:

  Paths to filtered scale-matrix files.

- mt_ls:

  Paths to filtered state-mean files.

- Mt_ls:

  Paths to filtered left-covariance files.

- at_ls:

  Paths to predicted state-mean files.

- At_ls:

  Paths to predicted left-covariance files.

- FT:

  Path to the covariate matrix at the final time step.

- VT:

  Path to the noise covariance matrix at the final time step.

## See also

[`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md),
[`generate.grid.rowsnake`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.rowsnake.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 2, p = 2, nT = 3)
par <- dat$para
dir <- file.path(tempdir(), "ff_big")
dir.create(dir, showWarnings = FALSE)
csv <- function(x, name) {
  f <- file.path(dir, name)
  write.csv(x, f, row.names = FALSE)
  f
}
Y_files <- mapply(csv, dat$Y, paste0("Y", 1:3, ".csv"))
files <- FF_bigdata_R(Y_files, csv(par$F0, "F.csv"), csv(par$G0, "G.csv"),
                      csv(par$W0, "W.csv"), csv(par$V0, "V.csv"), par$m0,
                      par$M0, par$n0, par$D0, nT = 3, fnrow = 5, fncol = 2,
                      bnrow = 5, bncol = 2, path_out = dir)
# With a single block this matches the in-memory filter
res_ff <- FF(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
             par$n0, par$D0, nT = 3)
all.equal(unname(as.matrix(read.csv(files$mt_ls[3]))), unname(res_ff$T3$mt))
#> [1] TRUE
unlink(dir, recursive = TRUE)
```
