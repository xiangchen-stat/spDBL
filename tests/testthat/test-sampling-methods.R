test_that("scalar-variance samplers have correct marginals and dependence", {
  ff <- list(
    list(mt = matrix(0), Mt = matrix(2/3), at = matrix(0), At = matrix(2),
         nt = 1000, Dt = 1000),
    list(mt = matrix(0), Mt = matrix(5/8), at = matrix(0), At = matrix(5/3),
         nt = 1000, Dt = 1000))
  bs <- list(st = array(0, c(1,1,2)),
             St = array(c(0.5, 5/8), c(1,1,2)))
  object <- list(ff = ff, bs = bs)
  set.seed(38)
  before <- .Random.seed
  marginal <- FFBS_sampling_sigma2R(8000L, object, matrix(1), matrix(1),
    2L, matrix(1), seed = 200L)
  joint <- FFBS_sampling_sigma2R(8000L, object, matrix(1), matrix(1),
    2L, matrix(1), seed = 200L, sampling = "joint")
  expect_identical(before, .Random.seed)
  expect_identical(marginal$Sigma, joint$Sigma)
  expect_identical(marginal[[2]], joint[[2]])
  standardized <- function(x, e) as.numeric(x[[e]]) / sqrt(x$Sigma)
  expect_equal(var(standardized(marginal, 1)), 0.5, tolerance = 0.025)
  expect_equal(var(standardized(marginal, 2)), 5/8, tolerance = 0.025)
  expect_equal(cov(standardized(marginal, 1), standardized(marginal, 2)),
               0, tolerance = 0.025)
  expect_equal(cov(standardized(joint, 1), standardized(joint, 2)),
               0.25, tolerance = 0.025)
  expect_false(attr(marginal, "joint_trajectories"))
  expect_true(attr(joint, "joint_trajectories"))
})

test_that("IW sampling supports both posterior laws", {
  ff <- list(
    list(mt = matrix(0,2,2), Mt = diag(2/3,2), at = matrix(0,2,2),
         At = diag(2,2), nt = 1001, Dt = diag(1000,2)),
    list(mt = matrix(0,2,2), Mt = diag(5/8,2), at = matrix(0,2,2),
         At = diag(5/3,2), nt = 1003, Dt = diag(1000,2)))
  object <- list(ff = ff, bs = list(st = array(0,c(2,2,2)),
    St = array(c(diag(0.5,2), diag(5/8,2)), c(2,2,2))))
  marginal <- FFBS_sampling(8000L, object, diag(2), diag(2), 2L, seed = 201L)
  joint <- FFBS_sampling(8000L, object, diag(2), diag(2), 2L,
                         seed = 201L, sampling = "joint")
  standardized <- function(x, e) x[[e]][1,1,] / sqrt(x$Sigma[1,1,])
  expect_identical(marginal$Sigma, joint$Sigma)
  expect_identical(marginal[[2]], joint[[2]])
  expect_equal(var(standardized(marginal, 1)), 0.5, tolerance = 0.025)
  expect_equal(cov(standardized(marginal, 1), standardized(marginal, 2)),
               0, tolerance = 0.025)
  expect_equal(cov(standardized(joint, 1), standardized(joint, 2)),
               0.25, tolerance = 0.025)
  expect_identical(attr(marginal, "sampling_method"), "marginal")
  expect_identical(attr(joint, "sampling_method"), "joint")
})

test_that("both sampling methods support singleton state and output dimensions", {
  for (model in c("iw", "sigma2R")) for (method in c("marginal", "joint")) {
    em <- emulator_learn(matrix(c(0,1),2), list(matrix(c(0.1,0.3),2)),
      1L, 1L, AR_choice = 1L, covariance = model, nsam = 2L,
      sampling = method, seed = 2L)
    expect_identical(dim(em$res_ffbs[[1]]), c(1L,1L,2L))
    expect_identical(attr(em$res_ffbs, "sampling_method"), method)
    expect_identical(attr(em$res_ffbs, "joint_trajectories"), method == "joint")
    prediction <- emulator_predict(em, 0.4, F_new_ls = list(matrix(0.2)))
    expect_identical(dim(prediction[[1]]), c(1L,1L))
    expect_true(all(is.finite(prediction[[1]])))
  }
})
