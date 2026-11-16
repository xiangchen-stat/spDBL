# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

interval_data <- function() {
  input <- cbind(a = seq(0.05, 0.95, length.out = 8),
                 b = c(0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4, 0.8))
  responses <- lapply(1:3, function(t)
    outer(input[, 1] + t * input[, 2], 1:4, function(x, y) sin(x) + y / 10))
  list(input = input, responses = responses)
}

test_that("interval = 'none' keeps the previous predict() results", {
  d <- interval_data()
  em <- emulator_learn(d$input, d$responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
    covariance = "sigma2R", nsam = 3L, seed = 6L)
  new <- d$input[c(2, 5), ] + 0.01
  held_out <- lapply(d$responses, function(y) y[c(2, 5), ] + 0.01)
  expect_identical(predict(em, new, dt_pde_test = held_out, interval = "none"),
                   emulator_predict(em, new, dt_pde_test = held_out))
  cal <- calibration_prepare(em, c(0, 0), c(1, 1))
  expect_identical(predict(cal, new, interval = "none"), emulator_predict(cal, new))
})

test_that("known-row prediction intervals are Student-t around the emulator mean", {
  d <- interval_data()
  new <- d$input[c(2, 5), ] + 0.01
  held_out <- lapply(d$responses, function(y) y[c(2, 5), ] + 0.01)
  for (covariance in c("iw", "sigma2R")) {
    em <- emulator_learn(d$input, d$responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
      covariance = covariance, nsam = 3L, seed = 6L)
    band <- predict(em, new, dt_pde_test = held_out, interval = "prediction",
                    level = 0.9)
    mean <- emulator_predict(em, new, dt_pde_test = held_out)
    expect_equal(unname(band$mean), unname(mean))
    expect_identical(band$distribution, "student_t")
    expect_gt(band$df, 2)
    for (t in seq_along(mean)) {
      expect_true(all(band$sd[[t]] > 0))
      half <- (band$upper[[t]] - band$lower[[t]]) / 2
      expect_equal(half, stats::qt(0.95, band$df) * band$sd[[t]] *
                     sqrt((band$df - 2) / band$df))
      expect_equal((band$upper[[t]] + band$lower[[t]]) / 2, unname(band$mean[[t]]))
    }
    wide <- predict(em, new, dt_pde_test = held_out, interval = "prediction",
                    level = 0.99)
    expect_true(all(wide$upper[[3]] - wide$lower[[3]] > band$upper[[3]] - band$lower[[3]]))
  }
})

test_that("known-row predictive variance matches Monte Carlo from the samplers", {
  # Monte Carlo reference built only from the package's own draws: the
  # coefficient-conditioned means of emulator_predict(MC = TRUE) plus the
  # input-GP conditional variance q * Sigma_m[c, c] of each draw.
  data("dt_emulation", package = "spDBL", envir = environment())
  d <- dt_emulation
  X_new <- d$pde_para_test[1:3, , drop = FALSE]
  test <- lapply(d$dt_pde_test, function(y) y[1:3, , drop = FALSE])
  for (covariance in c("iw", "sigma2R")) {
    em <- emulator_learn(d$pde_para_train, d$dt_pde_train, 10, 10, nsam = 300L,
      episode_window = 2, covariance = covariance, seed = 1L)
    band <- predict(em, X_new, dt_pde_test = test, interval = "prediction")
    mc <- predict(em, X_new, dt_pde_test = test, MC = TRUE)
    s <- em$setup
    input <- as.matrix(s$pde_para_train); N <- nrow(input)
    dd <- as.matrix(stats::dist(rbind(input, X_new)))
    J <- s$gp_sigma2 * exp(-s$phi_para * dd[seq_len(N), N + seq_len(nrow(X_new)), drop = FALSE])
    q <- s$gp_sigma2 + s$gp_tau2 - colSums(J * solve(s$V, J))
    sig <- em$res_ffbs$Sigma
    diag_draw <- if (covariance == "iw") apply(sig, 3, diag) else
      outer(diag(s$R_episode), as.numeric(sig))
    diag_draw <- diag_draw[rep(seq_len(s$bncol), s$n_b), , drop = FALSE]
    diag_draw <- diag_draw[order(as.vector(t(s$block_cols))), , drop = FALSE]
    ratio <- numeric()
    for (t in 3:length(mc)) {
      a <- mc[[t]]
      mc_var <- outer(q, rowMeans(diag_draw)) +
        apply(a, c(1, 2), function(v) mean((v - mean(v))^2))
      expect_lt(max(abs(band$mean[[t]] - apply(a, c(1, 2), mean))), 0.05)
      ratio <- c(ratio, as.vector(band$sd[[t]]^2 / mc_var))
    }
    expect_gt(stats::median(ratio), 0.95)
    expect_lt(stats::median(ratio), 1.05)
  }
})

test_that("design-GP prediction intervals follow the law of total variance", {
  d <- interval_data()
  em <- emulator_learn(d$input, d$responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
    covariance = "sigma2R", nsam = 3L, seed = 6L)
  cal <- calibration_prepare(em, c(0, 0), c(1, 1))
  new <- d$input[c(2, 5), ] + 0.01
  band <- predict(cal, new, interval = "prediction")
  expect_identical(band$distribution, "gaussian_mixture_moments")
  expect_identical(band$n_draws, cal$M)
  expect_true(band$design_uncertainty)
  expect_equal(unname(band$mean), unname(emulator_predict(cal, new)))
  # Independent law-of-total-variance calculation from the per-draw predictions.
  for (i in seq_len(nrow(new))) {
    per <- lapply(seq_len(cal$M), function(m)
      predict_calibration(new[i, ], m, cal, covariance = "marginal"))
    mu <- sapply(per, function(p) p$mu, simplify = "array")
    v <- sapply(per, function(p) p$variance, simplify = "array")
    total <- apply(v, c(1, 2), mean) + apply(mu, c(1, 2), function(x) mean((x - mean(x))^2))
    for (t in seq_len(cal$T)) {
      expect_equal(band$sd[[t]][i, ], sqrt(total[t, ]))
      expect_equal(band$upper[[t]][i, ] - band$mean[[t]][i, ],
                   stats::qnorm(0.975) * sqrt(total[t, ]))
    }
  }
  # Without design uncertainty the intervals can only narrow.
  plain <- calibration_prepare(em, c(0, 0), c(1, 1), design_uncertainty = FALSE)
  narrow <- predict(plain, new, interval = "prediction")
  expect_false(narrow$design_uncertainty)
  expect_true(all(unlist(narrow$sd) <= unlist(band$sd) + 1e-12))
})

test_that("prediction intervals reject bad arguments", {
  d <- interval_data()
  em <- emulator_learn(d$input, d$responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
    covariance = "sigma2R", nsam = 3L, seed = 6L)
  held_out <- lapply(d$responses, function(y) y[1:2, ])
  expect_error(predict(em, d$input[1:2, ], dt_pde_test = held_out,
                       interval = "prediction", level = 1), "level")
  expect_error(predict(em, d$input[1:2, ], dt_pde_test = held_out,
                       interval = "prediction", MC = TRUE), "MC = FALSE")
  expect_error(predict(em, d$input[1:2, ], interval = "prediction"), "dt_pde_test")
  cal <- calibration_prepare(em, c(0, 0), c(1, 1))
  expect_error(predict(cal, d$input[1:2, ], interval = "prediction", level = 0), "level")
  expect_error(predict(cal, d$input[1:2, ], interval = "prediction", MC = TRUE), "MC = FALSE")
  # A single new input given as a vector.
  one <- predict(cal, d$input[2, ] + 0.01, interval = "prediction")
  expect_identical(dim(one$mean[[1]]), c(1L, 4L))
})
