functions {
  vector rep_each(vector x, int K) {
    int N = rows(x);
    vector[N * K] y;
    int pos = 1;
    for (n in 1:N) {
      for (k in 1:K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
  }
  real lambda_gp(real L, int m) {
    real lam;
    lam = ((m*pi())/(2*L))^2;
    return lam;
  }
  vector phi_SE(real L, int m, vector x) {
    vector[rows(x)] fi;
    fi = 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));
    return fi;
  }
  real spd_SE(real alpha, real rho, real w) {
    real S;
    S = (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));
    return S;
  }
}
data {
  int<lower=0> total_obs; // total number of observations
  int<lower=0> n; // number of timepoints per series
  int<lower=0> n_sp; // number of smoothing parameters
  int<lower=0> n_series; // number of series
  int<lower=0> num_basis; // total number of basis coefficients
  vector[num_basis] zero; // prior locations for basis coefficients
  matrix[total_obs, num_basis] X; // mgcv GAM design matrix
  int<lower=0> ytimes[n, n_series]; // time-ordered matrix (which col in X belongs to each [time, series] observation?)
  matrix[11,33] S3; // mgcv smooth penalty matrix S3
  matrix[12,36] S4; // mgcv smooth penalty matrix S4
  matrix[12,36] S5; // mgcv smooth penalty matrix S5
  matrix[12,36] S6; // mgcv smooth penalty matrix S6
  matrix[12,36] S7; // mgcv smooth penalty matrix S7
  int<lower=0> n_nonmissing; // number of nonmissing observations
  int<lower=0> flat_ys[n_nonmissing]; // flattened nonmissing observations
  matrix[n_nonmissing, num_basis] flat_xs; // X values for nonmissing observations
  int<lower=0> obs_ind[n_nonmissing]; // indices of nonmissing observations
}
transformed data {
  // timepoints for latent gp; need 12 lags in total
  vector<lower=1>[n + 12] times;
  vector[n + 12] times_cent;
  real mean_times;
  real<lower=0> boundary;
  // number of basis functions for Hilbert space approx. gp
  int<lower=1> num_gp_basis;
  num_gp_basis = 20;
  matrix[n + 12, num_gp_basis] gp_phi;
  for (t in 1:(n + 12)){
    times[t] = t;
  }
  // zero-centre timepoints
  mean_times = mean(times);
  times_cent = times - mean_times;
  boundary = (5.0/4) * (max(times_cent) - min(times_cent));
  for (m in 1:num_gp_basis){
    gp_phi[,m] = phi_SE(boundary, m, times_cent);
  }
  int<lower=1> M;
  M = n_series * 4;
}
parameters {
  // raw basis coefficients
  vector[num_basis] b_raw;
  // random effect variances
  vector<lower=0>[2] sigma_raw;
  // random effect means
  vector[2] mu_raw;
  // negative binomial overdispersion (inverse)
  vector<lower=0>[n_series] r_inv;
  // gp length scale (alpha fixed at 0.25)
  real<lower=0> rho_gp;
  // gp basis coefficient weights
  matrix[num_gp_basis, 1] b_gp;
  // smoothing parameters
  vector<lower=0>[n_sp] lambda;
  // dynamic factor loadings
  vector[M] L;
}
transformed parameters {
  // latent gp spectral densities
  matrix[n + 12, 1] LV_current;
  matrix[n, 4] LV_lags;
  matrix[num_gp_basis, 1] diag_SPD;
  matrix[num_gp_basis, 1] SPD_beta;
  // dynamic factor loading matrix
  matrix[n_series, 4] lv_coefs_raw;
  // total trend component
  matrix[n, n_series] trend;
  // basis coefficients
  vector[num_basis] b;
  b[1:9] = mu_raw[1] + b_raw[1:9] * sigma_raw[1];
  b[10:18] = mu_raw[2] + b_raw[10:18] * sigma_raw[2];
  b[19:num_basis] = b_raw[19:num_basis];
  {
  int index;
  index = 0;
  for (j in 1:4) {
  for (i in 1:n_series) {
  index = index + 1;
  lv_coefs_raw[i, j] = L[index];
  }
  }
  }
  // gp estimates
  for (m in 1:num_gp_basis){
    diag_SPD[m, 1] = sqrt(spd_SE(0.25, rho_gp, sqrt(lambda_gp(boundary, m))));
  }
  SPD_beta = diag_SPD .* b_gp;
  // gp latent variable at the current timepoint
  LV_current = gp_phi * SPD_beta;
  // gp at lagged timepoints (3, 6 and 12 months)
  LV_lags[, 1] = LV_current[1:n, 1];
    for (i in 4:(n + 3)){
    LV_lags[i - 3, 2] = LV_current[i, 1];
  }
  for (i in 7:(n + 6)){
    LV_lags[i - 6, 3] = LV_current[i, 1];
  }
  for (i in 13:(n + 12)){
    LV_lags[i - 12, 4] = LV_current[i, 1];
  }
  // derived latent trends
  for (i in 1:n){;
    for (s in 1:n_series){
      trend[i, s] = dot_product(lv_coefs_raw[s, 1:4], LV_lags[i, 1:4]);
    }
  }
}
model {
  // prior for random effect population variances
  sigma_raw ~ exponential(0.5);
  // prior for random effect population means
  mu_raw ~ std_normal();
  // prior (non-centred) for s(series)...
  b_raw[1:9] ~ std_normal();
  // prior (non-centred) for s(ndvi_ma12,series)...
  b_raw[10:18] ~ std_normal();
  // prior for te(mintemp,lag)...
  b_raw[19:29] ~ multi_normal_prec(zero[19:29],S3[1:11,1:11] * lambda[3] + S3[1:11,12:22] * lambda[4] + S3[1:11,23:33] * lambda[5]);
  // prior for te(mintemp,lag):weights_dm...
  b_raw[30:41] ~ multi_normal_prec(zero[30:41],S4[1:12,1:12] * lambda[6] + S4[1:12,13:24] * lambda[7] + S4[1:12,25:36] * lambda[8]);
  // prior for te(mintemp,lag):weights_do...
  b_raw[42:53] ~ multi_normal_prec(zero[42:53],S5[1:12,1:12] * lambda[9] + S5[1:12,13:24] * lambda[10] + S5[1:12,25:36] * lambda[11]);
  // prior for te(mintemp,lag):weights_ot...
  b_raw[54:65] ~ multi_normal_prec(zero[54:65],S6[1:12,1:12] * lambda[12] + S6[1:12,13:24] * lambda[13] + S6[1:12,25:36] * lambda[14]);
  // prior for te(mintemp,lag):weights_pp...
  b_raw[66:77] ~ multi_normal_prec(zero[66:77],S7[1:12,1:12] * lambda[15] + S7[1:12,13:24] * lambda[16] + S7[1:12,25:36] * lambda[17]);
  // priors for smoothing parameters
  lambda ~ normal(30, 25);
  // priors for overdispersion parameters
  r_inv ~ exponential(5);
  // priors for gp parameters
  b_gp[1:num_gp_basis, 1] ~ std_normal();
  rho_gp ~ normal(6, 2.5);
  // priors for dynamic factor loading coefficients
  L ~ student_t(5, 0, 1);
  // likelihood functions
  vector[n_nonmissing] flat_trends;
  real flat_rs[n_nonmissing];
  flat_trends = (to_vector(trend))[obs_ind];
  flat_rs = to_array_1d(rep_each(r_inv, n)[obs_ind]);
  flat_ys ~ neg_binomial_2(
    exp(append_col(flat_xs, flat_trends) * append_row(b, 1.0)),
    inv(flat_rs));
}
generated quantities {
  vector[total_obs] eta;
  matrix[n, n_series] mus;
  matrix[n, 4] LV;
  matrix[n_series, 4] lv_coefs;
  vector[n_sp] rho;
  real alpha_gp;
  array[n, n_series] int ypred;
  rho = log(lambda);
  alpha_gp = 0.5;
  // Sign correct factor loadings and the latent gp;
  // picking a point where we expect strongly nonzero 
  // estimates (timepoint 180)
  for(j in 1:4){
   if(LV_lags[180, j] > 0){
    lv_coefs[, j] = -1 * lv_coefs_raw[, j];
    LV[, j] = -1 * LV_lags[, j];
  } else {
    lv_coefs[, j] = lv_coefs_raw[, j];
    LV[, j] = LV_lags[, j];
  }
  }

  matrix[n, n_series] r_vec;
  vector[n_series] r;
  r = inv(r_inv);
  for (s in 1:n_series) {
    r_vec[1:n,s] = rep_vector(r[s], n);
  }
  // posterior predictions
  eta = X * b;
  for(s in 1:n_series){ 
    mus[1:n, s] = eta[ytimes[1:n, s]] + trend[1:n, s];
    ypred[1:n, s] = neg_binomial_2_rng(exp(mus[1:n, s]), r_vec[1:n, s]);
  }
}
