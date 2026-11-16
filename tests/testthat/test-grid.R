# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

cell_counts <- function(ind, fnrow, fncol) {
  counts <- matrix(0L, fnrow, fncol)
  for (k in seq_len(nrow(ind))) {
    rows <- ind$U[k]:ind$D[k]
    cols <- ind$L[k]:ind$R[k]
    counts[rows, cols] <- counts[rows, cols] + 1L
  }
  counts
}

expect_partition <- function(ind, fnrow, fncol, bnrow, bncol) {
  expect_named(ind, c("block_no", "U", "D", "L", "R"))
  expect_equal(ind$block_no, seq_len(nrow(ind)))
  expect_equal(nrow(ind), ceiling(fnrow / bnrow) * ceiling(fncol / bncol))
  expect_true(all(ind$U <= ind$D & ind$L <= ind$R))
  expect_true(all(ind$D - ind$U + 1 <= bnrow & ind$R - ind$L + 1 <= bncol))
  expect_true(all(cell_counts(ind, fnrow, fncol) == 1L))
}

# Consecutive blocks of a snake traversal share an edge.
expect_snake <- function(ind) {
  for (k in seq_len(nrow(ind) - 1L)) {
    a <- ind[k, ]
    b <- ind[k + 1L, ]
    same_rows <- a$U == b$U && a$D == b$D
    same_cols <- a$L == b$L && a$R == b$R
    expect_true((same_rows && (b$L == a$R + 1 || b$R == a$L - 1)) ||
                (same_cols && b$U == a$D + 1))
  }
}

dims <- list(c(4, 6, 2, 2), c(6, 4, 2, 2), c(6, 6, 3, 2), c(2, 8, 2, 2),
             c(10, 9, 5, 3))

test_that("exact grids partition every cell once in snake or row order", {
  for (d in dims) {
    snake <- generate_grid(d[1], d[2], d[3], d[4])
    expect_partition(snake, d[1], d[2], d[3], d[4])
    expect_snake(snake)
    expect_equal(snake, generate.grid.rowsnake(d[1], d[2], d[3], d[4]), ignore_attr = TRUE)
    lr <- generate_grid(d[1], d[2], d[3], d[4], traversal.mode = "lr")
    expect_partition(lr, d[1], d[2], d[3], d[4])
    expect_equal(lr$L, rep(seq(1, d[2], d[4]), d[1] / d[3]))
    expect_equal(lr, generate.grid.lr(d[1], d[2], d[3], d[4]), ignore_attr = TRUE)
  }
  # A single column of blocks (exact path).
  expect_partition(generate_grid(8, 2, 2, 2), 8, 2, 2, 2)
  expect_partition(generate.grid.lr(8, 2, 2, 2), 8, 2, 2, 2)
})

test_that("flexible grids cover ragged edges exactly once", {
  for (d in list(c(5, 5, 2, 2), c(7, 10, 3, 4), c(105, 103, 10, 10), c(3, 7, 2, 3))) {
    snake <- generate_grid(d[1], d[2], d[3], d[4], is.flexible = TRUE)
    expect_partition(snake, d[1], d[2], d[3], d[4])
    expect_snake(snake)
    lr <- generate_grid(d[1], d[2], d[3], d[4], traversal.mode = "lr", is.flexible = TRUE)
    expect_partition(lr, d[1], d[2], d[3], d[4])
  }
  rows_only <- generate_grid(5, 4, 2, 2, is.flexible.row = TRUE)
  expect_partition(rows_only, 5, 4, 2, 2)
  cols_only <- generate_grid(4, 5, 2, 2, is.flexible.col = TRUE, traversal.mode = "lr")
  expect_partition(cols_only, 4, 5, 2, 2)
})

test_that("inexact block sizes are rejected unless flexible", {
  expect_error(generate_grid(5, 4, 2, 2), "rows per block 2 must divide number of rows in file 5")
  expect_error(generate_grid(4, 5, 2, 2), "cols per block 2 must divide number of cols in file 5")
  expect_error(generate_grid(5, 5, 2, 2, is.flexible.row = TRUE), "cols per block")
})

test_that("grids with one block or one block column are partitions", {
  expect_partition(generate.grid.rowsnake(2, 3, 2, 3), 2, 3, 2, 3)
  expect_partition(generate.grid.lr(2, 3, 2, 3), 2, 3, 2, 3)
  expect_partition(generate.grid.rowsnake(6, 2, 2, 2), 6, 2, 2, 2)
  expect_partition(generate_grid(5, 2, 2, 2, is.flexible = TRUE), 5, 2, 2, 2)
})

test_that("grid helpers reject non-positive or fractional sizes", {
  expect_error(generate.grid.rowsnake(0, 3, 1, 1), "positive whole numbers")
  expect_error(generate.grid.lr(4, 3, 1.5, 1), "positive whole numbers")
  expect_error(generate.grid.lr(4, NA, 1, 1), "positive whole numbers")
})
