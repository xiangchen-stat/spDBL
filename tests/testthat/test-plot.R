# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

skip_if_no_png <- function() {
  skip_if_not(isTRUE(capabilities("png")), "No PNG device available")
}

# The heatmap helpers need the suggested plotting packages.
skip_if_no_plot_packages <- function() {
  for (pkg in c("ggplot2", "ggpubr", "scales")) skip_if_not_installed(pkg)
}

test_that("error-bar summaries use row-wise quantiles", {
  set.seed(81)
  X <- matrix(rnorm(3 * 200), 3, 200)
  med <- cal_errorbar(X)
  avg <- cal_errorbar_mean(X)
  expect_named(med, c("med", "lower", "upper"))
  expect_equal(med$med, apply(X, 1, stats::median))
  expect_equal(avg$med, rowMeans(X))
  for (i in 1:3) {
    q <- stats::quantile(X[i, ], c(0.025, 0.975))
    expect_equal(c(med$lower[i], med$upper[i]), unname(q))
    expect_equal(c(avg$lower[i], avg$upper[i]), unname(q))
  }
})

test_that("quick_heat maps data to a raster and quick_save writes into path_fig", {
  skip_if_no_plot_packages()
  dt <- data.frame(x = rep(1:3, 2), y = rep(1:2, each = 3), fill = c(0, 1, 2, 3, 4, 9))
  p <- quick_heat(dt, max_y = 5)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(p$layers[[1]]$geom, "GeomRaster")
  # Values at 0 and above max_y take the end colours of the palette.
  palette <- toupper(spDBL:::spdbl_default_col_bgr())
  expect_identical(toupper(built$data[[1]]$fill[c(1, 6)]), palette[c(1, 14)])
  skip_if_no_png()
  dir <- tempfile("quick")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  quick_save("heat", dir, p)
  expect_length(list.files(dir, pattern = "^heat_.*\\.png$"), 1L)
})

nine_panel_data <- function() {
  set.seed(82)
  lapply(1:9, function(t) matrix(runif(2 * 6), 2, 6))
}

test_that("emulation panel helper returns one heatmap per time stamp", {
  skip_if_no_plot_packages()
  dat <- nine_panel_data()
  plots <- plot_panel_heatmap_9(dat, input_num = 2, tstamp = 1:9, max_y = 1,
                                Nx = 3, Ny = 2, nT = 9, savei = FALSE)
  expect_length(plots, 9L)
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_equal(plots[[4]]$data$sol, dat[[4]][2, ])
  expect_equal(plots[[4]]$data$row, rep(1:2, times = 3))
  expect_equal(plots[[4]]$data$col, rep(1:3, each = 2))
  skip_if_no_png()
  dir <- tempfile("panel")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  plot_panel_heatmap_9(dat, 1, 1:9, 1, 3, 2, 9, filename = "emu", path_fig = dir)
  expect_true(file.exists(file.path(dir, "emu.png")))
})

test_that("calibration panel helpers plot rows of a time-by-location matrix", {
  skip_if_no_plot_packages()
  set.seed(83)
  dat <- matrix(runif(9 * 6), 9, 6)
  loc <- data.frame(row = rep(1:2, times = 3), col = rep(1:3, each = 2))
  cal <- plot_panel_heatmap_9_cal(dat, 1:9, 1, loc, 3, 2, savei = FALSE)
  nolab <- plot_panel_heatmap_9_cal_nolab(dat, 1:9, min_y = 0.2, max_y = 0.8,
                                          loc_cal = loc, Nx = 3, Ny = 2, savei = FALSE)
  for (plots in list(cal, nolab)) {
    expect_length(plots, 9L)
    expect_equal(plots[[7]]$data$sol, dat[7, ])
  }
  expect_equal(nolab[[1]]$scales$get_scales("fill")$limits, c(0.2, 0.8))
  skip_if_no_png()
  dir <- tempfile("calpanel")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  plot_panel_heatmap_9_cal(dat, 1:9, 1, loc, 3, 2, filename = "cal", path_fig = dir)
  plot_panel_heatmap_9_cal_nolab(dat, 1:9, max_y = 1, loc_cal = loc, Nx = 3, Ny = 2,
                                 filename = "nolab", path_fig = dir)
  expect_setequal(list.files(dir), c("cal.png", "nolab.png"))
})
