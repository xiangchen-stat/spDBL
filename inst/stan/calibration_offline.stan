// Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
// SPDX-License-Identifier: MIT

// Equations (27)-(29), conditional on one fixed emulator posterior draw.
// y is integrated out. F(eta) is the mean of a GP fitted offline to AR2 designs.
// Its predictive variance enters as design_variance[t] * design_covariance[t]
// (zero if disabled), with design_variance[t] the GP variance factor of time
// t's group.
// Prediction matches FFBS_predict_exact(theta = one posterior trajectory),
// followed by recover_from_EP_exact(block_cols = ...). F(eta) comes from the
// offline AR2 GP; the original smoothed-mean predictor is not substituted here.
data {
  int<lower=1> N;
  int<lower=1> D_full;
  int<lower=1> D;
  int<lower=1> T;
  int<lower=1> S_full;
  int<lower=1> S;
  int<lower=1> B;
  int<lower=1> C;
  int<lower=1> P;
  int<lower=1> E;
  array[D] int<lower=1, upper=D_full> active;
  vector[D_full] eta_lower;
  vector[D_full] eta_upper;
  matrix[N, D_full] input;
  matrix[N, D] predictor_x;
  // Design GP: one kernel per group of design columns; the episodes of one
  // time step share a group. design_kernel 1 = squared exponential,
  // 2 = Matern 5/2, both on the unit-scaled free parameters.
  int<lower=1> G;
  array[T] int<lower=1, upper=G> time_group;
  int<lower=1, upper=2> design_kernel;
  array[G] vector<lower=0>[D] lengthscale;
  vector[E * P] predictor_center;
  matrix[N, E * P] predictor_alpha;
  array[G] matrix[N, N] predictor_L;
  array[T] matrix[S, S] design_covariance;
  array[E] matrix[P, C] theta;
  array[E] matrix[N, C] residual_weights;
  array[B, C] int<lower=1, upper=S_full> block_cols;
  array[S] int<lower=1, upper=S_full> locations;
  real<lower=0> phi_para;
  real<lower=0> gp_sigma2;
  real<lower=0> gp_tau2;
  real<lower=0> sigma2;
  matrix[N, N] L_V;
  matrix[S, S] R_field;
  matrix[S, S] L_U;
  matrix[S, S] L_M0;
  vector[S] m0;
  array[T] vector[S] z;
  real<lower=0> n0;
  real<lower=0> d0;
  real<lower=0, upper=1> discount;
}
parameters {
  vector<lower=0, upper=1>[D] eta_unit;
  vector<lower=0>[T + 1] precision;
  array[T + 1] vector[S] u_std;
}
transformed parameters {
  vector[D_full] eta = eta_lower;
  vector[T + 1] tau2 = inv(precision);
  array[T + 1] vector[S] u;
  array[T] vector[S] mu;
  matrix[S, S] Sigma;
  vector[T] design_variance;
  {
    vector[N] J;
    matrix[N, G] k;
    vector[G] group_variance;
    array[T] vector[S_full] mu_full;
    real conditional_variance;
    for (d in 1:D)
      eta[active[d]] = eta_lower[active[d]] +
        (eta_upper[active[d]] - eta_lower[active[d]]) * eta_unit[d];
    for (i in 1:N) {
      real distance2 = 0;
      for (d in 1:D_full) distance2 += square(input[i, d] - eta[d]);
      J[i] = gp_sigma2 * exp(-phi_para * sqrt(distance2 + 1e-24));
    }
    for (g in 1:G) {
      for (i in 1:N) {
        real r2 = 0;
        for (d in 1:D)
          r2 += square((predictor_x[i, d] - eta_unit[d]) / lengthscale[g][d]);
        if (design_kernel == 1) {
          k[i, g] = exp(-0.5 * r2);
        } else {
          // The offset keeps the gradient finite at a training input.
          real a = sqrt(5 * r2 + 1e-24);
          k[i, g] = (1 + a + square(a) / 3) * exp(-a);
        }
      }
      group_variance[g] = fmax(0, 1 - dot_self(mdivide_left_tri_low(predictor_L[g], col(k, g))));
    }
    for (t in 1:T) design_variance[t] = group_variance[time_group[t]];
    conditional_variance = gp_sigma2 + gp_tau2 -
      dot_self(mdivide_left_tri_low(L_V, J));
    Sigma = sigma2 * conditional_variance * R_field;
    for (t in 1:T) mu_full[t] = rep_vector(0, S_full);
    for (e in 1:E) {
      int t = (e - 1) / B + 1;
      int b = (e - 1) % B + 1;
      vector[P] design_row = segment(predictor_center, (e - 1) * P + 1, P)
        + block(predictor_alpha, 1, (e - 1) * P + 1, N, P)' * col(k, time_group[t]);
      vector[C] block_mu = theta[e]' * design_row + residual_weights[e]' * J;
      for (c in 1:C) mu_full[t][block_cols[b, c]] = block_mu[c];
    }
    for (t in 1:T) for (s in 1:S) mu[t][s] = mu_full[t][locations[s]];
    u[1] = m0 + sqrt(tau2[1]) * L_M0 * u_std[1];
    for (t in 1:T) u[t + 1] = u[t] + sqrt(tau2[t + 1]) * L_U * u_std[t + 1];
  }
}
model {
  // eta_unit has a uniform prior; Stan supplies its constraint Jacobian.
  precision[1] ~ gamma(n0, d0);
  for (t in 1:T) precision[t + 1] ~ gamma(pow(discount, t) * n0, pow(discount, t) * d0);
  for (t in 1:(T + 1)) u_std[t] ~ std_normal();
  for (t in 1:T) {
    z[t] ~ multi_normal_cholesky(mu[t] + u[t + 1], cholesky_decompose(
      add_diag(Sigma + design_variance[t] * design_covariance[t], tau2[t + 1])));
  }
}
