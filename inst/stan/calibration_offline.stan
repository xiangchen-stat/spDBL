// Equations (27)-(29), conditional on one fixed emulator posterior draw.
// y is integrated out. F(eta) is the mean of a GP fitted offline to AR2 designs.
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
    matrix[S, S] observation_covariance = Sigma;
    for (s in 1:S) observation_covariance[s, s] += tau2[t + 1];
    z[t] ~ multi_normal_cholesky(mu[t] + u[t + 1], cholesky_decompose(observation_covariance));
  }
}
