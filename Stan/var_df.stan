// Stan model code for a Poisson VAR(1)
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
  vector[n_series] sigma = rep_vector(0.1, n_series);
}
parameters {
  // trend component parameters //
  // VAR coefficients for latent trends
  matrix[n_series, n_series] beta_var;
  // trend estimates
  vector[n_series] trend_raw[n];
  
  // GAM component parameters //
  // raw basis coefficients
  vector[num_basis] b_raw;
  // random effect variances
  vector<lower=0>[1] sigma_raw;
  // random effect means
  vector[1] mu_raw;
  // smoothing parameters
  vector<lower=0>[n_sp] lambda;

}
transformed parameters {
  // basis coefficients
  vector[num_basis] b;
  
  // trends in matrix form
  matrix[n, n_series] trend;
  
  // GAM beta coefficients
  b[1:9] = mu_raw[1] + b_raw[1:9] * sigma_raw[1];
  b[10:18] = mu_raw[2] + b_raw[10:18] * sigma_raw[2];
  b[19:num_basis] = b_raw[19:num_basis];
  
  // trend matrix-form estimates
  for(i in 1:n){
    trend[i, 1:n_series] = to_row_vector(trend_raw[i]);
  }
}
model {
  // latent trend component priors //
  // VAR coefficients
  to_vector(beta_var) ~ std_normal();
  // latent trend means
  vector[n_series] mu[n - 1];
  for(i in 2:n){
    mu[i - 1] = beta_var * trend_raw[i - 1];
  }
  // latent trends (contemporaneously uncorrelated errors)
  trend_raw[1] ~ normal(0, 2.5);
  for(i in 1:n){
    trend_raw[i] ~ normal(mu[i], sigma);
  }

  // GAM component priors //
  // prior for random effect population variances
  sigma_raw ~ exponential(1.5);
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
  
  // likelihood functions // 
  vector[n_nonmissing] flat_trends;
  flat_trends = (to_vector(trend))[obs_ind];
  flat_ys ~ poisson_log_glm(append_col(flat_xs, flat_trends),
                            0.0, append_row(b, 1.0));
}
generated quantities {
  vector[total_obs] eta;
  matrix[n, n_series] mus;
  array[n, n_series] int ypred;
  vector[n_sp] rho;
  rho = log(lambda);
  
  // posterior predictions
  eta = X * b;
  for(s in 1:n_series){ 
  mus[1:n, s] = eta[ytimes[1:n, s]] + trend[1:n, s];
  ypred[1:n, s] = poisson_log_rng(mus[1:n, s]);
  }
}