# Generate synthetic FFBS data and write to CSV files

Simulates data from the MNIW dynamic linear model for `nT` time steps
and writes the observations \\Y_1, \ldots, Y_T\\ and all model
parameters to CSV files in `path`.

## Usage

``` r
gen_ffbs_csv(N, S, p, nT, path)
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

- path:

  Character. Directory path (with trailing `/`) where CSV files are
  written.

## Value

Invisibly returns `NULL`. Files written: `m0.csv`, `MM0.csv`, `G0.csv`,
`F0.csv`, `Sigma.csv`, `RSigma.csv`, `V0.csv`, `RV0.csv`, `W0.csv`, and
`Y1.csv` through `Y<nT>.csv`.

## See also

[`gen_ffbs_data`](https://xiangchen-stat.github.io/spDBL/reference/gen_ffbs_data.md)

## Examples

``` r
set.seed(1)
out_dir <- file.path(tempdir(), "ffbs_csv")
dir.create(out_dir)
gen_ffbs_csv(N = 10, S = 3, p = 2, nT = 3, path = paste0(out_dir, "/"))
list.files(out_dir)
#>  [1] "F0.csv"     "G0.csv"     "MM0.csv"    "RSigma.csv" "RV0.csv"   
#>  [6] "Sigma.csv"  "V0.csv"     "W0.csv"     "Y1.csv"     "Y2.csv"    
#> [11] "Y3.csv"     "m0.csv"    
dim(read.csv(file.path(out_dir, "Y1.csv")))
#> [1] 10  3
unlink(out_dir, recursive = TRUE)
```
