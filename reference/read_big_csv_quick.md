# Read a rectangular block from a large CSV file

Efficiently reads a contiguous rectangular block of rows and columns
from a CSV file using
[`readr::read_csv`](https://readr.tidyverse.org/reference/read_delim.html).
Useful for processing big data in chunks without loading the entire file
into memory.

## Usage

``` r
read_big_csv_quick(filename, rows = c(1, Inf), cols = c(1, Inf), header = TRUE)
```

## Arguments

- filename:

  Character. Path to the CSV file.

- rows:

  Integer vector of length 2 giving the first and last row indices to
  read (1-based, excluding the header row if `header = TRUE`). Use
  `c(1, Inf)` to read all rows. Defaults to `c(1, Inf)`.

- cols:

  Integer vector of length 2 giving the first and last column indices to
  read. Use `c(1, Inf)` to read all columns. Defaults to `c(1, Inf)`.

- header:

  Logical. Whether the file has a header row. Defaults to `TRUE`.

## Value

A `tibble` (from `readr`) containing the requested block.

## Examples

``` r
f <- file.path(tempdir(), "big.csv")
write.csv(matrix(1:20, 5, 4), f, row.names = FALSE)
read_big_csv_quick(f, rows = c(2, 4), cols = c(2, 3))
#> # A tibble: 3 × 2
#>      X2    X3
#>   <dbl> <dbl>
#> 1     7    12
#> 2     8    13
#> 3     9    14
unlink(f)
```
