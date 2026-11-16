# Build AR(1) covariate list from a list of response matrices

Constructs the time-varying covariate list \\F_t\\ for a first-order
autoregressive model: \\F_t = Y\_{t-1}\\. No earlier output exists at
\\t = 1\\, so \\F_1\\ is a matrix of ones (a per-site intercept). The
design is causal and deterministic: \\F_t\\ uses only outputs before
time \\t\\. (Versions up to 1.1.1 used \\F_1 = Y_1\\, the response
itself.)

## Usage

``` r
gen_F_ls_AR1(Y, nT)
```

## Arguments

- Y:

  List of length `nT`. Each element is the \\N \times S\\ data matrix at
  time \\t\\.

- nT:

  Integer. Number of time steps.

## Value

A list of length `nT` where element \\t\\ is the covariate matrix
\\F_t\\.

## See also

[`gen_F_ls_AR2`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR2.md),
[`gen_F_ls_AR1_EP`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR1_EP.md)

## Examples

``` r
Y <- lapply(1:4, function(t) matrix(t, 2, 3))
F_ls <- gen_F_ls_AR1(Y, nT = 4)
length(F_ls)
#> [1] 4
sapply(F_ls, function(Ft) Ft[1, 1])  # ones at t = 1, then lags: 1 1 2 3
#> [1] 1 1 2 3
```
