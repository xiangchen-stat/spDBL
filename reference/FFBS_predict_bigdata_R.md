# Out-of-core exact prediction

Predicts simulator output at new inputs from the smoothed coefficients
written by
[`BS_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/BS_bigdata_R.md),
streaming one block of data at a time. For each block it computes the
same mean as
[`FFBS_predict_exact()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md):
\\F\_{new}\hat\Theta + J^\top V^{-1}(Y - F\hat\Theta)\\, with the
exponential input kernel. Every block must contain all training runs
(`bnrow = fnrow` in
[`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md)),
as in
[`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md).

## Usage

``` r
FFBS_predict_bigdata_R(
  Y_ls,
  F_ls,
  F_new_ls,
  input,
  input_new,
  bs_files,
  nT,
  fnrow,
  fncol,
  bncol,
  phi_para,
  gp_sigma2 = 1.1,
  gp_tau2 = 1e-04,
  path_out,
  verbose = FALSE
)
```

## Arguments

- Y_ls, F_ls:

  The response and design file paths given to
  [`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md).

- F_new_ls:

  Design file path(s) for the new inputs, laid out like `F_ls`: one row
  per new input, `fncol` lag-1 columns and, for AR(2) designs, `fncol`
  lag-2 columns after them. One path or one per time step.

- input:

  Matrix of training inputs, one row per run (`fnrow` rows).

- input_new:

  Matrix (or vector, for one input) of new inputs.

- bs_files:

  The list returned by
  [`BS_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/BS_bigdata_R.md).

- nT, fnrow, fncol, bncol:

  As given to
  [`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md).

- phi_para, gp_sigma2, gp_tau2:

  Range, amplitude and nugget of the exponential input kernel, as in
  [`FFBS_predict_exact()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md).

- path_out:

  Directory where the prediction files are written.

- verbose:

  Logical. If `TRUE`, reports progress.

## Value

A character vector with one CSV path per time step; each file holds the
predicted means, one row per new input and one column per site.

## See also

[`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md),
[`BS_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/BS_bigdata_R.md),
[`FFBS_predict_exact()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md)

## Examples

``` r
set.seed(1)
dir <- file.path(tempdir(), "pred_big")
dir.create(dir, showWarnings = FALSE)
csv <- function(x, name) {
  path <- file.path(dir, name)
  utils::write.csv(x, path, row.names = FALSE)
  path
}
# 5 training runs and 2 new inputs; 4 sites in two blocks; 3 times
input <- matrix(runif(10), 5, 2)
input_new <- matrix(runif(4), 2, 2)
Y <- lapply(1:3, function(t) matrix(rnorm(20, t), 5, 4))
F <- c(list(matrix(1, 5, 4)), Y[1:2])
V <- 1.1 * exp(-2 * as.matrix(dist(input))) + diag(1e-4, 5)
ff <- FF_bigdata_R(Y_ls = mapply(csv, Y, paste0("Y", 1:3, ".csv")),
  F_ls = mapply(csv, F, paste0("F", 1:3, ".csv")),
  G_ls = csv(diag(2), "G.csv"), W_ls = csv(diag(2), "W.csv"),
  V_ls = csv(V, "V.csv"), m0 = matrix(0, 2, 2), M0 = diag(2), n0 = 4,
  D0 = diag(2), nT = 3, fnrow = 5, fncol = 4, bnrow = 5, bncol = 2,
  path_out = dir)
bs <- BS_bigdata_R(ff, file.path(dir, "G.csv"), nT = 3, fnrow = 5,
  fncol = 4, bnrow = 5, bncol = 2, path_out = dir)
# Designs for the new inputs, e.g. built from their own lagged outputs
F_new <- c(list(matrix(1, 2, 4)), lapply(1:2, function(t) matrix(t, 2, 4)))
pred <- FFBS_predict_bigdata_R(file.path(dir, paste0("Y", 1:3, ".csv")),
  file.path(dir, paste0("F", 1:3, ".csv")),
  mapply(csv, F_new, paste0("Fnew", 1:3, ".csv")), input, input_new, bs,
  nT = 3, fnrow = 5, fncol = 4, bncol = 2, phi_para = 2, path_out = dir)
utils::read.csv(pred[3])  # predicted means at time 3
#>         V1       V2       V3       V4
#> 1 3.992888 4.078902 2.134635 2.556459
#> 2 3.393615 3.673938 2.652602 1.888386
unlink(dir, recursive = TRUE)
```
