# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

test_that("FF matches exact joint Gaussian conditioning of the MNIW DLM", {
  d <- dlm_problem(11)
  ff <- FF(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$n0, d$D0, d$nT)
  ref <- batch_dlm(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0)
  expect_named(ff, c("T1", "T2", "T3", "prior"))
  for (t in seq_len(d$nT)) {
    expect_named(ff[[t]], c("nt", "Dt", "at", "At", "mt", "Mt"))
    for (k in c("at", "At", "mt", "Mt"))
      expect_equal(ff[[t]][[k]], ref[[t]][[k]], tolerance = 1e-9, ignore_attr = TRUE)
    expect_equal(ff[[t]]$nt, d$n0 + sum(d$N[seq_len(t)]))
    expect_equal(ff[[t]]$Dt, d$D0 + ref[[t]]$quad, tolerance = 1e-9, ignore_attr = TRUE)
  }
  expect_identical(ff$prior, list(m0 = d$m0, M0 = d$M0, n0 = d$n0, D0 = d$D0))
})

test_that("constant system matrices equal their replicated time-varying lists", {
  d <- dlm_problem(12, N = c(3L, 3L, 3L), time_varying = FALSE)
  rep3 <- function(x) rep(list(x), d$nT)
  constant <- FF(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$n0, d$D0, d$nT)
  listed <- FF(d$Y, rep3(d$F_ls), rep3(d$G_ls), rep3(d$W_ls), rep3(d$V_ls),
               d$m0, d$M0, d$n0, d$D0, d$nT)
  expect_equal(constant, listed)
  ref <- batch_dlm(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0)
  expect_equal(constant$T3$mt, ref[[3]]$mt, tolerance = 1e-9)
  msgs <- capture_messages(
    FF(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$n0, d$D0, d$nT,
       verbose = TRUE))
  expect_true(any(grepl("FF: 2/3", msgs)))
})

test_that("the C++ one-step update equals the closed-form MNIW posterior", {
  set.seed(13)
  N <- 5; p <- 3; q <- 2
  X <- matrix(rnorm(N * p), N, p)
  Y <- matrix(rnorm(N * q), N, q)
  H <- rand_spd(N)
  C <- matrix(rnorm(p * q), p, q)
  Vb <- rand_spd(p)
  S <- rand_spd(q)
  v <- 6
  # Precision form of the conjugate update.
  O <- solve(solve(Vb) + t(X) %*% solve(H) %*% X)
  Cnew <- O %*% (solve(Vb) %*% C + t(X) %*% solve(H) %*% Y)
  r <- Y - X %*% C
  Snew <- S + t(r) %*% solve(X %*% Vb %*% t(X) + H) %*% r
  # With G = I and W = 0 the filter step is exactly this update.
  one <- spDBL:::FF_1step_cpp(Yt = Y, Ft = X, Gt = diag(p), Wt = matrix(0, p, p),
                              Vt = H, mt_1 = C, Mt_1 = Vb, nt_1 = v, Dt_1 = S)
  expect_equal(one$mt, Cnew, tolerance = 1e-10)
  expect_equal(one$Mt, O, tolerance = 1e-10)
  expect_equal(one$Dt, Snew, tolerance = 1e-10)
  expect_equal(one$nt, v + N)
  expect_equal(one$at, C)
  expect_equal(one$At, Vb)
  for (fun in list(MNIW_R, MNIW_R_naiive)) {
    post <- fun(X = X, Y = Y, H = H, v = v, S = S, C = C, Vb = Vb)
    expect_named(post, c("vnew", "Snew", "Cnew", "Vbnew"))
    expect_equal(post$Cnew, Cnew, tolerance = 1e-10)
    expect_equal(post$Vbnew, O, tolerance = 1e-10)
    expect_equal(post$Snew, Snew, tolerance = 1e-10)
    expect_equal(post$vnew, v + N)
  }
})

test_that("scalar-variance and identity one-step updates share the state recursion", {
  set.seed(14)
  N <- 4; p <- 2; S <- 3
  Yt <- matrix(rnorm(N * S), N, S)
  Ft <- matrix(rnorm(N * p), N, p)
  Gt <- diag(p) + 0.1
  Wt <- rand_spd(p, 0.2)
  Vt <- rand_spd(N)
  mt_1 <- matrix(rnorm(p * S), p, S)
  Mt_1 <- rand_spd(p)
  R <- rand_spd(S)
  iw <- spDBL:::FF_1step_cpp(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, 4, diag(S))
  ig <- FF_1step_R_sigma2R(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, nt_1 = 2, Dt_1 = 0.5,
                           Rinv = solve(R))
  id <- FF_1step_R_I(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1)
  for (k in c("at", "At", "mt", "Mt")) {
    expect_equal(ig[[k]], iw[[k]], tolerance = 1e-10)
    expect_equal(id[[k]], iw[[k]], tolerance = 1e-10)
  }
  # sigma^2 | Y ~ IG(n + NS/2, d + tr(R^{-1} r' Q^{-1} r) / 2).
  r <- Yt - Ft %*% Gt %*% mt_1
  Q <- Ft %*% (Gt %*% Mt_1 %*% t(Gt) + Wt) %*% t(Ft) + Vt
  expect_equal(ig$nt, 2 + N * S / 2)
  expect_equal(ig$Dt, 0.5 + sum(diag(solve(R) %*% t(r) %*% solve(Q) %*% r)) / 2,
               tolerance = 1e-10)
  expect_false(any(c("nt", "Dt") %in% names(id)))
})

test_that("FF_sigma2R and FF_I match exact conditioning", {
  d <- dlm_problem(15)
  set.seed(16)
  R <- rand_spd(d$q)
  ref <- batch_dlm(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0)
  ig <- FF_sigma2R(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0,
                   n0 = 3, D0 = 0.7, nT = d$nT, R = R)
  id <- FF_I(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$nT)
  expect_named(ig, c("T1", "T2", "T3", "prior"))
  expect_named(id, c("T1", "T2", "T3", "prior"))
  for (t in seq_len(d$nT)) {
    for (k in c("at", "At", "mt", "Mt")) {
      expect_equal(ig[[t]][[k]], ref[[t]][[k]], tolerance = 1e-9, ignore_attr = TRUE)
      expect_equal(id[[t]][[k]], ref[[t]][[k]], tolerance = 1e-9, ignore_attr = TRUE)
    }
    expect_equal(ig[[t]]$nt, 3 + sum(d$N[seq_len(t)]) * d$q / 2)
    expect_equal(ig[[t]]$Dt, 0.7 + sum(diag(ref[[t]]$quad %*% solve(R))) / 2,
                 tolerance = 1e-9)
  }
  expect_identical(id$prior, list(m0 = d$m0, M0 = d$M0))
  expect_identical(ig$prior$D0, 0.7)
  expect_error(FF_sigma2R(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0,
                          3, diag(2), d$nT, R), "D0 is not a number")
})

test_that("verbose filters report progress", {
  d <- dlm_problem(17, N = rep(2L, 10), time_varying = FALSE)
  R <- diag(d$q)
  msgs <- capture_messages(FF_sigma2R(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls,
    d$m0, d$M0, 3, 1, d$nT, R, verbose = TRUE))
  expect_length(msgs, 10L)
  expect_match(msgs[10], "FF: 10/10")
  msgs <- capture_messages(FF_I(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls,
    d$m0, d$M0, d$nT, verbose = TRUE))
  expect_length(msgs, 10L)
})

test_that("verbose filters and smoothers handle one- and two-step series", {
  # The progress step used to round to zero for short series (i %% 0 is NA).
  for (n in 1:2) {
    d <- dlm_problem(40 + n, N = rep(2L, n), time_varying = FALSE)
    R <- diag(d$q)
    expect_length(capture_messages(FF(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls,
      d$m0, d$M0, d$n0, d$D0, d$nT, verbose = TRUE)), n)
    expect_length(capture_messages(FF_sigma2R(d$Y, d$F_ls, d$G_ls, d$W_ls,
      d$V_ls, d$m0, d$M0, 3, 1, d$nT, R, verbose = TRUE)), n)
    expect_length(capture_messages(FF_I(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls,
      d$m0, d$M0, d$nT, verbose = TRUE)), n)
    fit <- FFBS_I(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$nT)
    # The smoother reports the n - 1 backward steps.
    expect_length(capture_messages(BS(fit$ff, d$G_ls, d$nT, verbose = TRUE)), n - 1L)
    expect_length(capture_messages(FFBS_sampling_I(2L, fit, d$F_ls, d$G_ls,
      d$nT, verbose = TRUE)), n)
  }
})

test_that("read_big_csv_quick reads rectangular blocks", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  x <- matrix(round(rnorm(20), 6), 5, 4)
  utils::write.csv(x, path, row.names = FALSE)
  block <- read_big_csv_quick(path, rows = c(2, 4), cols = c(2, 3))
  expect_equal(unname(as.matrix(block)), x[2:4, 2:3])
  expect_equal(unname(as.matrix(read_big_csv_quick(path))), x)
  expect_equal(unname(as.matrix(read_big_csv_quick(path, rows = c(5, 5)))),
               x[5, , drop = FALSE])
  raw <- tempfile(fileext = ".csv")
  on.exit(unlink(raw), add = TRUE)
  utils::write.table(x, raw, sep = ",", row.names = FALSE, col.names = FALSE)
  expect_equal(unname(as.matrix(read_big_csv_quick(raw, rows = c(1, 2),
    cols = c(1, 4), header = FALSE))), x[1:2, ])
})

# Split full-grid matrices into the snake-ordered blocks used by FF_bigdata_R.
snake_blocks <- function(full, ind, ar2_offset = NULL) {
  lapply(seq_len(nrow(ind)), function(i) {
    rows <- ind$U[i]:ind$D[i]
    cols <- ind$L[i]:ind$R[i]
    if (!is.null(ar2_offset)) cols <- c(cols, ar2_offset + cols)
    full[rows, cols, drop = FALSE]
  })
}

write_csv_matrix <- function(x, dir, name) {
  path <- file.path(dir, name)
  utils::write.csv(x, path, row.names = FALSE)
  path
}

check_bigdata_output <- function(files, ff, nT, K) {
  read <- function(path) unname(as.matrix(utils::read.csv(path)))
  for (t in seq_len(nT)) {
    steps <- ff[(t - 1L) * K + seq_len(K)]
    bound <- function(k) unname(do.call(cbind, lapply(steps, `[[`, k)))
    expect_equal(read(files$mt_ls[t]), bound("mt"), tolerance = 1e-8)
    expect_equal(read(files$Mt_ls[t]), bound("Mt"), tolerance = 1e-8)
    expect_equal(read(files$at_ls[t]), bound("at"), tolerance = 1e-8)
    expect_equal(read(files$At_ls[t]), bound("At"), tolerance = 1e-8)
    expect_equal(read(files$Dt_ls[t]), bound("Dt"), tolerance = 1e-8)
    expect_equal(as.numeric(read(files$nt_ls[t])),
                 unname(vapply(steps, `[[`, numeric(1), "nt")))
  }
}

test_that("out-of-core FF_bigdata_R reproduces in-memory FF over snake blocks", {
  fnrow <- 4; fncol <- 4; bnrow <- 2; bncol <- 2; nT <- 2
  ind <- generate.grid.rowsnake(fnrow, fncol, bnrow, bncol)
  K <- nrow(ind)
  root <- tempfile("ffbig")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  # Time-varying design files with AR1 (fncol) or AR1 + AR2 (2 fncol) columns.
  for (lags in 1:2) {
    set.seed(18 + lags)
    dir <- file.path(root, paste0("lags", lags))
    dir.create(dir, recursive = TRUE)
    p <- lags * bncol
    Y <- lapply(seq_len(nT), function(t) matrix(rnorm(fnrow * fncol), fnrow, fncol))
    Fm <- lapply(seq_len(nT), function(t) matrix(rnorm(fnrow * lags * fncol), fnrow))
    V <- rand_spd(fnrow, 0.5)
    G <- 0.9 * diag(p)
    W <- rand_spd(p, 0.2)
    m0 <- matrix(0, p, bncol)
    files <- FF_bigdata_R(
      Y_ls = vapply(seq_len(nT), function(t) write_csv_matrix(Y[[t]], dir, paste0("Y", t, ".csv")), ""),
      F_ls = vapply(seq_len(nT), function(t) write_csv_matrix(Fm[[t]], dir, paste0("F", t, ".csv")), ""),
      G_ls = write_csv_matrix(G, dir, "G.csv"), W_ls = write_csv_matrix(W, dir, "W.csv"),
      V_ls = write_csv_matrix(V, dir, "V.csv"), m0 = m0, M0 = diag(p), n0 = 3,
      D0 = diag(bncol), nT = nT, fnrow = fnrow, fncol = fncol, bnrow = bnrow,
      bncol = bncol, path_out = dir)
    offset <- if (lags == 2) fncol
    V_blocks <- lapply(seq_len(K), function(i) V[ind$U[i]:ind$D[i], ind$U[i]:ind$D[i]])
    F_blocks <- unlist(lapply(Fm, snake_blocks, ind, ar2_offset = offset), recursive = FALSE)
    ff <- FF(Y = unlist(lapply(Y, snake_blocks, ind), recursive = FALSE),
             F_ls = F_blocks, G_ls = G, W_ls = W, V_ls = rep(V_blocks, nT), m0 = m0,
             M0 = diag(p), n0 = 3, D0 = diag(bncol), nT = nT * K)
    expect_true(all(file.exists(unlist(files))))
    check_bigdata_output(files, ff, nT, K)
    read <- function(path) unname(as.matrix(utils::read.csv(path)))
    expect_equal(read(files$FT), unname(F_blocks[[nT * K]]), tolerance = 1e-12)
    expect_equal(read(files$VT), unname(V_blocks[[K]]), tolerance = 1e-12)
  }
})

test_that("FF_bigdata_R reads AR2 lags and time-varying system files", {
  set.seed(19)
  dir <- tempfile("ffbig2")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  fnrow <- 4; fncol <- 4; bnrow <- 2; bncol <- 2; nT <- 2; p <- 2 * bncol
  ind <- generate.grid.rowsnake(fnrow, fncol, bnrow, bncol)
  K <- nrow(ind)
  Y <- lapply(seq_len(nT), function(t) matrix(rnorm(fnrow * fncol), fnrow, fncol))
  Fm <- matrix(rnorm(fnrow * 2 * fncol), fnrow, 2 * fncol)
  V <- lapply(seq_len(nT), function(t) rand_spd(fnrow, 0.5))
  G <- lapply(seq_len(nT), function(t) diag(p) * (0.8 + 0.1 * t))
  W <- lapply(seq_len(nT), function(t) rand_spd(p, 0.2))
  m0 <- matrix(0.1, p, bncol)
  wr <- function(x, prefix) vapply(seq_len(nT), function(t)
    write_csv_matrix(x[[t]], dir, paste0(prefix, t, ".csv")), "")
  msgs <- capture_messages(files <- FF_bigdata_R(wr(Y, "Y"),
    write_csv_matrix(Fm, dir, "F.csv"), wr(G, "G"), wr(W, "W"), wr(V, "V"),
    m0 = m0, M0 = diag(p), n0 = 4, D0 = diag(bncol), nT = nT, fnrow = fnrow,
    fncol = fncol, bnrow = bnrow, bncol = bncol, path_out = dir, verbose = TRUE))
  expect_length(msgs, nT)
  expect_match(msgs[1], "Forward filter 1/2")
  rep_t <- function(x) unlist(lapply(seq_len(nT), function(t) rep(list(x[[t]]), K)),
                              recursive = FALSE)
  ff <- FF(Y = unlist(lapply(Y, snake_blocks, ind), recursive = FALSE),
           F_ls = rep(snake_blocks(Fm, ind, ar2_offset = fncol), nT),
           G_ls = rep_t(G), W_ls = rep_t(W),
           V_ls = unlist(lapply(V, function(v) lapply(seq_len(K), function(i)
             v[ind$U[i]:ind$D[i], ind$U[i]:ind$D[i]])), recursive = FALSE),
           m0 = m0, M0 = diag(p), n0 = 4, D0 = diag(bncol), nT = nT * K)
  check_bigdata_output(files, ff, nT, K)
  expect_equal(dim(utils::read.csv(files$FT)), c(bnrow, p))
})

test_that("an unsupported delta warns and is ignored", {
  d <- dlm_problem(44, N = rep(2L, 3), time_varying = FALSE)
  plain <- FFBS(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0, d$M0, d$n0, d$D0, d$nT)
  expect_warning(discounted <- FFBS(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls, d$m0,
    d$M0, d$n0, d$D0, d$nT, delta = 0.9), "not implemented")
  expect_identical(discounted, plain)
  # One warning per call, not one per time step.
  expect_length(capture_warnings(FF_I(d$Y, d$F_ls, d$G_ls, d$W_ls, d$V_ls,
    d$m0, d$M0, d$nT, delta = 0.5)), 1L)
  expect_warning(spDBL:::.check_delta(NA), "not implemented")
  expect_identical(spDBL:::.check_delta(1L), 1)
})

test_that("out-of-core smoothing and prediction reproduce the in-memory results", {
  N <- 5; N_new <- 2; fncol <- 6; bncol <- 2; nT <- 3
  ind <- generate.grid.rowsnake(N, fncol, N, bncol)
  ind_new <- generate.grid.rowsnake(N_new, fncol, N_new, bncol)
  K <- nrow(ind)
  block_cols <- t(vapply(seq_len(K), function(b) ind$L[b]:ind$R[b], numeric(bncol)))
  read <- function(path) unname(as.matrix(utils::read.csv(path)))
  root <- tempfile("bigpred")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  for (lags in 1:2) {
    set.seed(30 + lags)
    dir <- file.path(root, paste0("lags", lags))
    dir.create(dir, recursive = TRUE)
    p <- lags * bncol
    input <- matrix(runif(N * 2), N, 2)
    input_new <- matrix(runif(N_new * 2), N_new, 2)
    Y <- lapply(seq_len(nT), function(t) matrix(rnorm(N * fncol), N, fncol))
    Fm <- lapply(seq_len(nT), function(t) matrix(rnorm(N * lags * fncol), N))
    F_new <- lapply(seq_len(nT), function(t) matrix(rnorm(N_new * lags * fncol), N_new))
    V <- 1.1 * exp(-2 * as.matrix(dist(input))) + diag(1e-4, N)
    G <- 0.9 * diag(p)
    W <- rand_spd(p, 0.2)
    m0 <- matrix(0, p, bncol)
    files <- function(x, stem) vapply(seq_along(x), function(t)
      write_csv_matrix(x[[t]], dir, paste0(stem, t, ".csv")), "")
    Y_files <- files(Y, "Y")
    F_files <- files(Fm, "F")
    F_new_files <- files(F_new, "Fnew")
    G_file <- write_csv_matrix(G, dir, "G.csv")
    ff_files <- FF_bigdata_R(Y_ls = Y_files, F_ls = F_files, G_ls = G_file,
      W_ls = write_csv_matrix(W, dir, "W.csv"), V_ls = write_csv_matrix(V, dir, "V.csv"),
      m0 = m0, M0 = diag(p), n0 = 3, D0 = diag(bncol), nT = nT, fnrow = N,
      fncol = fncol, bnrow = N, bncol = bncol, path_out = dir)
    bs_files <- BS_bigdata_R(ff_files, G_file, nT, N, fncol, N, bncol, dir)
    # In memory: the same blocks as one episode sequence.
    offset <- if (lags == 2) fncol
    Y_ep <- unlist(lapply(Y, snake_blocks, ind), recursive = FALSE)
    F_ep <- unlist(lapply(Fm, snake_blocks, ind, ar2_offset = offset), recursive = FALSE)
    F_new_ep <- unlist(lapply(F_new, snake_blocks, ind_new, ar2_offset = offset),
                       recursive = FALSE)
    fit <- FFBS(Y_ep, F_ep, G, W, V, m0, diag(p), 3, diag(bncol), nT * K)
    for (t in seq_len(nT)) {
      episodes <- (t - 1L) * K + seq_len(K)
      expect_equal(read(bs_files$st_ls[t]),
                   do.call(cbind, lapply(episodes, function(e) fit$bs$st[, , e])),
                   tolerance = 1e-8)
      expect_equal(read(bs_files$St_ls[t]),
                   do.call(cbind, lapply(episodes, function(e) fit$bs$St[, , e])),
                   tolerance = 1e-8)
    }
    pred_files <- FFBS_predict_bigdata_R(Y_files, F_files, F_new_files, input,
      input_new, bs_files, nT, N, fncol, bncol, phi_para = 2, path_out = dir)
    exact <- FFBS_predict_exact(Y = Y_ep, para_ffbs = fit, F_ls = F_ep,
      F_new_ls = F_new_ep, input = input, input_new = input_new, nT = nT * K,
      phi_para = 2)
    full <- recover_from_EP_exact(exact, nT_ori = nT, nT = nT * K, block_cols = block_cols)
    for (t in seq_len(nT))
      expect_equal(read(pred_files[t]), unname(full[[t]]), tolerance = 1e-8)
  }
  # A single new input and constant design files also work.
  msgs <- capture_messages(one <- FFBS_predict_bigdata_R(Y_files, F_files,
    F_new_files[1], input, input_new[1, ], bs_files, nT, N, fncol, bncol,
    phi_para = 2, path_out = dir, verbose = TRUE))
  expect_identical(msgs, paste0("Prediction ", seq_len(nT), "/", nT, "\n"))
  expect_identical(dim(read(one[1])), c(1L, as.integer(fncol)))
  expect_error(FFBS_predict_bigdata_R(Y_files, F_files, F_new_files, input[1:3, ],
    input_new, bs_files, nT, N, fncol, bncol, phi_para = 2, path_out = dir),
    "one row per training run")
  expect_error(FFBS_predict_bigdata_R(Y_files, F_files, F_new_files, input,
    input_new, bs_files, nT, N, fncol, bncol, phi_para = -1, path_out = dir),
    "Invalid input kernel")
  wide <- FFBS_predict_bigdata_R  # designs narrower than the coefficients fail
  expect_error(wide(Y_files, files(lapply(Fm, function(x) x[, 1:fncol]), "Fshort"),
    F_new_files, input, input_new, bs_files, nT, N, fncol, bncol, phi_para = 2,
    path_out = dir), "do not match the smoothed coefficients")
})

test_that("out-of-core smoothing handles blocks that split the runs", {
  fnrow <- 4; fncol <- 4; bnrow <- 2; bncol <- 2; nT <- 2
  ind <- generate.grid.rowsnake(fnrow, fncol, bnrow, bncol)
  K <- nrow(ind)
  dir <- tempfile("bsrows")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  set.seed(35)
  Y <- lapply(seq_len(nT), function(t) matrix(rnorm(fnrow * fncol), fnrow, fncol))
  Fm <- lapply(seq_len(nT), function(t) matrix(rnorm(fnrow * fncol), fnrow))
  V <- rand_spd(fnrow, 0.5)
  # A different transition matrix at each time.
  G <- lapply(seq_len(nT), function(t) (0.8 + 0.1 * t) * diag(bncol))
  G_files <- vapply(seq_len(nT), function(t) write_csv_matrix(G[[t]], dir, paste0("G", t, ".csv")), "")
  W <- rand_spd(bncol, 0.2)
  ff_files <- FF_bigdata_R(
    Y_ls = vapply(seq_len(nT), function(t) write_csv_matrix(Y[[t]], dir, paste0("Y", t, ".csv")), ""),
    F_ls = vapply(seq_len(nT), function(t) write_csv_matrix(Fm[[t]], dir, paste0("F", t, ".csv")), ""),
    G_ls = G_files, W_ls = write_csv_matrix(W, dir, "W.csv"),
    V_ls = write_csv_matrix(V, dir, "V.csv"), m0 = matrix(0, bncol, bncol),
    M0 = diag(bncol), n0 = 3, D0 = diag(bncol), nT = nT, fnrow = fnrow,
    fncol = fncol, bnrow = bnrow, bncol = bncol, path_out = dir)
  msgs <- capture_messages(bs_files <- BS_bigdata_R(ff_files, G_files, nT, fnrow,
    fncol, bnrow, bncol, dir, verbose = TRUE))
  expect_identical(msgs, c("Backward smoother 2/2\n", "Backward smoother 1/2\n"))
  V_blocks <- lapply(seq_len(K), function(i) V[ind$U[i]:ind$D[i], ind$U[i]:ind$D[i]])
  G_ep <- rep(G, each = K)
  ff <- FF(Y = unlist(lapply(Y, snake_blocks, ind), recursive = FALSE),
           F_ls = unlist(lapply(Fm, snake_blocks, ind), recursive = FALSE),
           G_ls = G_ep, W_ls = W, V_ls = rep(V_blocks, nT), m0 = matrix(0, bncol, bncol),
           M0 = diag(bncol), n0 = 3, D0 = diag(bncol), nT = nT * K)
  bs <- BS(ff, G_ep, nT * K)
  read <- function(path) unname(as.matrix(utils::read.csv(path)))
  for (t in seq_len(nT))
    expect_equal(read(bs_files$st_ls[t]),
                 do.call(cbind, lapply((t - 1L) * K + seq_len(K), function(e) bs$st[, , e])),
                 tolerance = 1e-8)
})
