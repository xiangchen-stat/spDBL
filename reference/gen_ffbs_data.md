# Generate synthetic FFBS data in memory

Simulates data from the MNIW dynamic linear model for `nT` time steps
and returns the generated observations and model parameters as in-memory
objects.

## Usage

``` r
gen_ffbs_data(N, S, p, nT)
```

## Arguments

- N:

  Integer. Number of spatial locations (rows of \\Y_t\\).

- S:

  Integer. Number of response variables (columns of \\Y_t\\).

- p:

  Integer. Dimension of the state vector (must be \\\geq 2\\).

- nT:

  Integer. Number of time steps.

## Value

A named list with:

- Y:

  Named list of length `nT` (`"Y1"`, ..., `"Y<nT>"`), each element an
  \\N \times S\\ data matrix.

- para:

  Named list of model parameters: `loc`, `dist`, `m0`, `M0`, `G0`, `F0`,
  `Sigma`, `RSigma`, `V0`, `RV0`, `W0`, `n0`, `D0`.

## See also

[`gen_ffbs_csv`](https://xiangchen-stat.github.io/spDBL/reference/gen_ffbs_csv.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 10, S = 3, p = 2, nT = 4)
names(dat$Y)
#> [1] "Y1" "Y2" "Y3" "Y4"
dim(dat$Y$Y1)
#> [1] 10  3
dim(dat$para$V0)
#> [1] 10 10
```
