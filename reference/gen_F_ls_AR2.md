# Build AR(2) covariate list from a list of response matrices

Constructs the time-varying covariate list \\F_t\\ for a second-order
autoregressive model by column-binding the two most recent lags: \\F_t =
\[Y\_{t-2}, Y\_{t-1}\]\\. A lag that does not exist yet is replaced by a
matrix of ones (a per-site intercept): \\F_1 = \[1, 1\]\\ and \\F_2 =
\[1, Y_1\]\\. The design is causal and deterministic. (Versions up to
1.1.1 used \\Y_1\\ plus random noise at \\t = 1\\, and the future value
\\Y_3\\ at \\t = 2\\.)

## Usage

``` r
gen_F_ls_AR2(Y, nT)
```

## Arguments

- Y:

  List of length `nT`. Each element is the \\N \times S\\ data matrix at
  time \\t\\.

- nT:

  Integer. Number of time steps.

## Value

A list of length `nT` where element \\t\\ is the \\N \times 2S\\
covariate matrix \\F_t\\.

## See also

[`gen_F_ls_AR1`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR1.md),
[`gen_F_ls_AR2_EP`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR2_EP.md)

## Examples

``` r
Y <- lapply(1:4, function(t) matrix(t, 2, 3))
F_ls <- gen_F_ls_AR2(Y, nT = 4)
dim(F_ls[[4]])
#> [1] 2 6
F_ls[[2]][1, ]  # ones (no lag 2 yet), then lag 1
#> [1] 1 1 1 1 1 1
F_ls[[4]][1, ]  # lags 2 and 3
#> [1] 2 2 2 3 3 3
```
