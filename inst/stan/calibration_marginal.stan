// Equations (27)-(29), conditional on one fixed emulator posterior draw.
// y and u are integrated out exactly with a Kalman filter.
// F(eta) is the mean of a GP fitted offline to AR2 designs.
// Prediction matches FFBS_predict_exact(theta = one posterior trajectory) and
// recover_from_EP_exact(block_cols = ...), retaining sigma2 and R_field.
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
  vector<lower=0>[D] lengthscale;
  vector[E * P] predictor_center;
  matrix[N, E * P] predictor_alpha;
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
transformed data {
  matrix[S, S] U = L_U * L_U';
  matrix[S, S] M0 = L_M0 * L_M0';
}
parameters {
  vector<lower=0, upper=1>[D] eta_unit;
  vector<lower=0>[T + 1] precision;
}
transformed parameters {
  vector[D_full] eta = eta_lower;
  vector[T + 1] tau2 = inv(precision);
  array[T] vector[S] mu;
  matrix[S, S] Sigma;
  {
    vector[N] J;
    vector[N] k;
    vector[E * P] design;
    array[T] vector[S_full] mu_full;
    real conditional_variance;
    for (d in 1:D)
      eta[active[d]] = eta_lower[active[d]] +
        (eta_upper[active[d]] - eta_lower[active[d]]) * eta_unit[d];
    for (i in 1:N) {
      real distance2 = 0;
      real predictor_distance2 = 0;
      for (d in 1:D_full) distance2 += square(input[i, d] - eta[d]);
      for (d in 1:D)
        predictor_distance2 += square((predictor_x[i, d] - eta_unit[d]) / lengthscale[d]);
      J[i] = gp_sigma2 * exp(-phi_para * sqrt(distance2 + 1e-24));
      k[i] = exp(-0.5 * predictor_distance2);
    }
    design = predictor_center + predictor_alpha' * k;
    conditional_variance = gp_sigma2 + gp_tau2 -
      dot_self(mdivide_left_tri_low(L_V, J));
    Sigma = sigma2 * conditional_variance * R_field;
    for (t in 1:T) mu_full[t] = rep_vector(0, S_full);
    for (e in 1:E) {
      int t = (e - 1) / B + 1;
      int b = (e - 1) % B + 1;
      vector[C] block_mu = theta[e]' * segment(design, (e - 1) * P + 1, P)
        + residual_weights[e]' * J;
      for (c in 1:C) mu_full[t][block_cols[b, c]] = block_mu[c];
    }
    for (t in 1:T) for (s in 1:S) mu[t][s] = mu_full[t][locations[s]];
  }
}
model {
  vector[S] filtered_mean = m0;
  matrix[S, S] filtered_covariance = tau2[1] * M0;
  precision[1] ~ gamma(n0, d0);
  for (t in 1:T) precision[t + 1] ~ gamma(pow(discount, t) * n0, pow(discount, t) * d0);
  for (t in 1:T) {
    matrix[S, S] prior_covariance = filtered_covariance + tau2[t + 1] * U;
    matrix[S, S] observation_covariance = prior_covariance + Sigma;
    matrix[S, S] L;
    matrix[S, S] whitened_prior;
    vector[S] residual = z[t] - mu[t] - filtered_mean;
    for (s in 1:S) observation_covariance[s, s] += tau2[t + 1];
    L = cholesky_decompose(observation_covariance);
    z[t] ~ multi_normal_cholesky(mu[t] + filtered_mean, L);
    whitened_prior = mdivide_left_tri_low(L, prior_covariance);
    filtered_mean += whitened_prior' * mdivide_left_tri_low(L, residual);
    filtered_covariance = prior_covariance - whitened_prior' * whitened_prior;
    filtered_covariance = 0.5 * (filtered_covariance + filtered_covariance');
  }
}
