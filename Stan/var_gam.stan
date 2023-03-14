// Stan model code for a Negative Binomial GAM with latent VAR(1) dynamic process
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
  matrix[11,33] S2; // mgcv smooth penalty matrix S2
  matrix[12,36] S3; // mgcv smooth penalty matrix S3
  matrix[12,36] S4; // mgcv smooth penalty matrix S4
  matrix[12,36] S5; // mgcv smooth penalty matrix S5
  matrix[12,36] S6; // mgcv smooth penalty matrix S6
  int<lower=0> n_nonmissing; // number of nonmissing observations
  int<lower=0> flat_ys[n_nonmissing]; // flattened nonmissing observations
  matrix[n_nonmissing, num_basis] flat_xs; // X values for nonmissing observations
  int<lower=0> obs_ind[n_nonmissing]; // indices of nonmissing observations
}

parameters {
  // trend component parameters //
  // VAR coefficients for latent trends
  matrix<lower=-1,upper=1>[n_series, n_series] beta_var;
  // trend estimates
  vector[n_series] trend_raw[n];
  // trend SDs
  vector<lower=0,upper=1>[n_series] sigma;
  
  // GAM component parameters //
  // raw basis coefficients
  vector[num_basis] b_raw;
  // random effect variance
  vector<lower=0>[1] sigma_raw;
  // raw random effect mean
  vector[1] mu_raw;
  // smoothing parameters
  vector<lower=0>[n_sp] lambda;
  
  // negative binomial overdispersion //
  vector<lower=0>[n_series] r;
}

transformed parameters {
  vector[num_basis] b;
  matrix[n, n_series] trend;
  
  // trend estimates in matrix-form
  for(i in 1:n){
    trend[i, 1:n_series] = to_row_vector(trend_raw[i]);
  }
  
  // rescaled z-score random NDVI slope estimates
  b[1:9] = mu_raw[1] + b_raw[1:9] * sigma_raw[1];
  // remaining GAM beta coefficients are unchanged from the model block
  b[10:num_basis] = b_raw[10:num_basis];
}

model {
  // latent trend component priors //
  // VAR coefficients
  to_vector(beta_var) ~ normal(0, 0.25);
  // trend SDs
  sigma ~ beta(8, 12);
  // trend means
  vector[n_series] mu[n - 1];
  for(i in 2:n){
    mu[i - 1] = beta_var * trend_raw[i - 1];
  }
  // latent trends (contemporaneously uncorrelated)
  trend_raw[1] ~ normal(0, sigma);
  for(i in 2:n){
    trend_raw[i] ~ normal(mu[i - 1], sigma);
  }

  // GAM component priors //
  // containment prior for random NDVI slope population SD
  sigma_raw ~ inv_gamma(2.3693353, 0.7311319);
  // prior for random NDVI slope mean
  mu_raw ~ std_normal();
  // prior (non-centred) for s(ndvi_ma12,series)...
  b_raw[1:9] ~ std_normal();
  // prior for te(mintemp,lag)...
  b_raw[10:20] ~ multi_normal_prec(zero[10:20],S2[1:11,1:11] * lambda[2] + S2[1:11,12:22] * lambda[3] + S2[1:11,23:33] * lambda[4]);
  // prior for te(mintemp,lag):weights_dm...
  b_raw[21:32] ~ multi_normal_prec(zero[21:32],S3[1:12,1:12] * lambda[5] + S3[1:12,13:24] * lambda[6] + S3[1:12,25:36] * lambda[7]);
  // prior for te(mintemp,lag):weights_do...
  b_raw[33:44] ~ multi_normal_prec(zero[33:44],S4[1:12,1:12] * lambda[8] + S4[1:12,13:24] * lambda[9] + S4[1:12,25:36] * lambda[10]);
  // prior for te(mintemp,lag):weights_ot...
  b_raw[45:56] ~ multi_normal_prec(zero[45:56],S5[1:12,1:12] * lambda[11] + S5[1:12,13:24] * lambda[12] + S5[1:12,25:36] * lambda[13]);
  // prior for te(mintemp,lag):weights_pp...
  b_raw[57:68] ~ multi_normal_prec(zero[57:68],S6[1:12,1:12] * lambda[14] + S6[1:12,13:24] * lambda[15] + S6[1:12,25:36] * lambda[16]);
  // priors for smoothing parameters
  lambda ~ normal(30, 25);
  
  // containment prior for overdispersion parameters
  r ~ inv_gamma(0.5408871, 6.8770244);
  
  // likelihood functions
  vector[n_nonmissing] flat_trends;
  real flat_rs[n_nonmissing];
  flat_trends = (to_vector(trend))[obs_ind];
  flat_rs = to_array_1d(rep_each(r, n)[obs_ind]);
  flat_ys ~ neg_binomial_2(exp(append_col(flat_xs, flat_trends) * append_row(b, 1.0)),
                           flat_rs);
}

generated quantities {
  vector[total_obs] eta;
  matrix[n, n_series] mus;
  array[n, n_series] int ypred;
  vector[n_sp] rho;
  matrix[n, n_series] r_vec;
  
  // vector form of overdispersion parameters
  for (s in 1:n_series) {
   r_vec[1:n, s] = rep_vector(r[s], n);
  }
  
  // logged smoothing parameters
  rho = log(lambda);
  
  // posterior predictions
  eta = X * b;
  for(s in 1:n_series){ 
    mus[1:n, s] = eta[ytimes[1:n, s]] + trend[1:n, s];
    ypred[1:n, s] = neg_binomial_2_rng(exp(mus[1:n, s]), r_vec[1:n, s]);
  }
}
