# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# keep_draws = FALSE stores where each episode's draws start in the random-
# number stream instead of the draws; everything downstream must be identical.

small_problem <- function() {
  set.seed(4)
  input <- cbind(a = seq(0.05, 0.95, length.out = 10), b = runif(10))
  responses <- lapply(1:3, function(t)
    outer(input[, 1] + t * input[, 2], 1:4, function(x, y) sin(x) + y / 10))
  list(train = input[1:8, ], responses = lapply(responses, function(y) y[1:8, ]),
       new = input[9:10, ], new_responses = lapply(responses, function(y) y[9:10, ]))
}

learn <- function(d, covariance, keep_draws)
  emulator_learn(d$train, d$responses, Nx = 2L, Ny = 2L, AR_choice = 1L,
    covariance = covariance, nsam = 4L, seed = 17L, keep_draws = keep_draws)

test_that("regenerated draws equal stored draws for both covariance models", {
  d <- small_problem()
  for (covariance in c("iw", "sigma2R")) {
    stored <- learn(d, covariance, TRUE)
    lazy <- learn(d, covariance, FALSE)
    expect_s3_class(lazy$res_ffbs, "spDBL_lazy_draws")
    expect_identical(lazy$res_ffbs$Sigma, stored$res_ffbs$Sigma)
    moments <- spDBL:::.draw_moments(lazy$para_ffbs)
    before <- .Random.seed
    for (e in seq_len(stored$setup$nT))
      expect_identical(spDBL:::.regenerate_episode(lazy$res_ffbs, moments, e),
                       stored$res_ffbs[[e]])
    # Regenerating leaves the session's random numbers where they were.
    expect_identical(.Random.seed, before)
    expect_identical(emulator_predict(lazy, d$new, dt_pde_test = d$new_responses, MC = TRUE),
                     emulator_predict(stored, d$new, dt_pde_test = d$new_responses, MC = TRUE))
    expect_identical(spDBL:::.emulator_settings(lazy), spDBL:::.emulator_settings(stored))
  }
})

test_that("calibration on regenerated draws matches stored draws", {
  d <- small_problem()
  stored <- calibration_prepare(learn(d, "sigma2R", TRUE), c(0, 0), c(1, 1))
  lazy <- calibration_prepare(learn(d, "sigma2R", FALSE), c(0, 0), c(1, 1))
  expect_s3_class(lazy$posterior$theta, "spDBL_lazy_draws")
  expect_identical(predict(lazy, d$new, interval = "prediction"),
                   predict(stored, d$new, interval = "prediction"))
  expect_identical(emulator_predict(lazy, d$new), emulator_predict(stored, d$new))
  expect_identical(spDBL:::predict_calibration(d$new[1, ], 2, lazy),
                   spDBL:::predict_calibration(d$new[1, ], 2, stored))
  expect_identical(calibration_validate(lazy, d$new, d$new_responses)$field,
                   calibration_validate(stored, d$new, d$new_responses)$field)
  z <- matrix(rnorm(12), 3, 4)
  expect_identical(spDBL:::make_stan_calibration_data(lazy, 3, z, diag(4)),
                   spDBL:::make_stan_calibration_data(stored, 3, z, diag(4)))
  # Fresh marginal draws for a calibration are regenerated the same way.
  fresh_lazy <- calibration_prepare(learn(d, "sigma2R", FALSE), c(0, 0), c(1, 1),
                                    n_draws = 3L, seed = 5L)
  fresh_stored <- calibration_prepare(learn(d, "sigma2R", TRUE), c(0, 0), c(1, 1),
                                      n_draws = 3L, seed = 5L)
  expect_s3_class(fresh_lazy$posterior$theta, "spDBL_lazy_draws")
  for (e in seq_len(fresh_stored$E))
    expect_identical(spDBL:::.cal_theta(fresh_lazy, e), fresh_stored$posterior$theta[[e]])
  # Joint draws for a calibration are stored even when the emulator's are not.
  joint <- calibration_prepare(learn(d, "sigma2R", FALSE), c(0, 0), c(1, 1),
                               sampling = "joint")
  expect_false(inherits(joint$posterior$theta, "spDBL_lazy_draws"))
})

test_that("keep_draws is validated and unstored draws cannot be indexed", {
  d <- small_problem()
  expect_error(emulator_learn(d$train, d$responses, 2L, 2L, AR_choice = 1L,
    nsam = 2L, keep_draws = NA), "keep_draws must be TRUE or FALSE")
  expect_error(emulator_learn(d$train, d$responses, 2L, 2L, AR_choice = 1L,
    nsam = 2L, sampling = "joint", keep_draws = FALSE), "marginal")
  lazy <- learn(d, "iw", FALSE)
  expect_error(lazy$res_ffbs[[1]], "not stored")
  expect_identical(length(lazy$res_ffbs), lazy$setup$nT)
  expect_output(print(lazy$res_ffbs), "regenerated on demand")
  expect_error(FFBS_sampling(2L, lazy$para_ffbs, lazy$setup$F_ls, lazy$setup$G_ls,
    lazy$setup$nT, sampling = "joint", keep_draws = FALSE), "marginal")
})
