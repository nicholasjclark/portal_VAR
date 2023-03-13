// Stan code for stationary VAR(p) model
// This code is modified from Heaps, S
// Enforcing stationarity through the prior in vector autoregressions
// Journal of Computational and Graphical Statistics (2022): 1-24.
functions {
  /* Function to compute the matrix square root */
    matrix sqrtm(matrix A) {
      int m = rows(A);
      vector[m] root_root_evals = sqrt(sqrt(eigenvalues_sym(A)));
      matrix[m, m] evecs = eigenvectors_sym(A);
      matrix[m, m] eprod = diag_post_multiply(evecs, root_root_evals);
      return tcrossprod(eprod);
    }
  /* Function to transform A to P (inverse of part 2 of reparameterisation) */
    matrix AtoP(matrix A) {
      int m = rows(A);
      matrix[m, m] B = tcrossprod(A);
      for(i in 1:m) B[i, i] += 1.0;
      return mdivide_left_spd(sqrtm(B), A);
    }
  /* Function to perform reverse mapping. The details of
  how to perform Step 1 are in Section S1.3 of Heaps 2022 Supplementary Materials.
  Returned: a (2 x p) array of (m x m) matrices; the (1, i)-th component
  of the array is phi_i and the (2, i)-th component of the array
  is Gamma_{i-1}*/
    matrix[,] rev_mapping(matrix[] P, matrix Sigma) {
      int p = size(P);
      int m = rows(Sigma);
      matrix[m, m] phi_for[p, p];   matrix[m, m] phi_rev[p, p];
      matrix[m, m] Sigma_for[p+1];  matrix[m, m] Sigma_rev[p+1];
      matrix[m, m] S_for;           matrix[m, m] S_rev;
      matrix[m, m] S_for_list[p+1];
      matrix[m, m] Gamma_trans[p+1];
      matrix[m, m] phiGamma[2, p];
      // Step 1:
        Sigma_for[p+1] = Sigma;
        S_for_list[p+1] = sqrtm(Sigma);
        for(s in 1:p) {
          // In this block of code S_rev is B^{-1} and S_for is a working matrix
          S_for = - tcrossprod(P[p-s+1]);
          for(i in 1:m) S_for[i, i] += 1.0;
          S_rev = sqrtm(S_for);
          S_for_list[p-s+1] = mdivide_right_spd(mdivide_left_spd(S_rev, 
                                                                 sqrtm(quad_form_sym(Sigma_for[p-s+2], S_rev))), S_rev);
          Sigma_for[p-s+1] = tcrossprod(S_for_list[p-s+1]);
        }
        // Step 2:
          Sigma_rev[1] = Sigma_for[1];
          Gamma_trans[1] = Sigma_for[1];
          for(s in 0:(p-1)) {
            S_for = S_for_list[s+1];
            S_rev = sqrtm(Sigma_rev[s+1]);
            phi_for[s+1, s+1] = mdivide_right_spd(S_for * P[s+1], S_rev);
            phi_rev[s+1, s+1] = mdivide_right_spd(S_rev * P[s+1]', S_for);
      Gamma_trans[s+2] = phi_for[s+1, s+1] * Sigma_rev[s+1];
      if(s>=1) {
        for(k in 1:s) {
          phi_for[s+1, k] = phi_for[s, k] - phi_for[s+1, s+1] * phi_rev[s, s-k+1];
          phi_rev[s+1, k] = phi_rev[s, k] - phi_rev[s+1, s+1] * phi_for[s, s-k+1];
        }
        for(k in 1:s) Gamma_trans[s+2] = Gamma_trans[s+2] + phi_for[s, k] * 
                                                               Gamma_trans[s+2-k];
      }
      Sigma_rev[s+2] = Sigma_rev[s+1] - quad_form_sym(Sigma_for[s+1], 
                                                      phi_rev[s+1, s+1]');
          }
          for(i in 1:p) phiGamma[1, i] = phi_for[p, i];
          for(i in 1:p) phiGamma[2, i] = Gamma_trans[i]';
    return phiGamma;
  }
}
data {
  int<lower=0> total_obs; // total number of observations
  int<lower=0> n; // number of timepoints per series
  int<lower=0> n_series; // number of series
  int<lower=0> num_basis; // total number of basis coefficients
  matrix[total_obs, num_basis] X; // mgcv GAM design matrix
  int<lower=0> ytimes[n, n_series]; // time-ordered matrix (which col in X belongs to each [time, series] observation?)
  int<lower=0> n_nonmissing; // number of nonmissing observations
  int<lower=0> flat_ys[n_nonmissing]; // flattened nonmissing observations
  matrix[n_nonmissing, num_basis] flat_xs; // X values for nonmissing observations
  int<lower=0> obs_ind[n_nonmissing]; // indices of nonmissing observations
  int<lower=1> p; // order of latent trend VAR model
}
transformed data {
  // trend covariance  degrees of freedom to ensure finite variance
  real df = n_series + 4.0;
  // (zero)-mean of trend VAR process
  vector[n_series] mu = rep_vector(0.0, n_series);
  // scale-matrix in prior for Sigma
  matrix[n_series, n_series] scale_mat;            
  for(i in 1:n_series) {
    for(j in 1:n_series) {
      if(i==j) scale_mat[i, j] = 1.0;
      else scale_mat[i, j] = 0.0;
    }
  }
}
parameters {
  // raw basis coefficients
  vector[num_basis] b_raw;
  // random effect variances
  vector<lower=0>[1] sigma_raw;
  // random effect means
  vector[1] mu_raw;
  // square matrices map trend partial autocorrelations to unconstrained space as
  // a (p-length array of n_series * n_series matrices);
  // allows construction of prior distributions for the A's which encourage shrinkage 
  // toward meaningful parametric structures. This mapping and prior combination also
  // allows borrowing of strength between the diagonal elements and between 
  // the off-diagonal elements of each A matrix
  matrix[n_series, n_series] A[p];
  // trend covariance, Sigma
  cov_matrix[n_series] Sigma; 
  // means and precisions in hyperpriors for unconstrained A's
  vector[p] Amu;
  vector<lower=0>[p] Aomega;
  // raw trends for initial times 1:p
  vector[p*n_series] trend_raw_init;
  // raw trends from times (p+1):n
  vector[n_series] trend_raw[n-p];
}
transformed parameters {
  // basis coefficients
  vector[num_basis] b;
  // VAR coeficients (p-length array of n_series * n_series matrices)
  matrix[n_series, n_series] phi[p];  
  // reconstructed trends in matrix form
  matrix[n, n_series] trend;
  // initial trend stationary covariance;
  // a positive definite block Toeplitz matrix 
  // available as a by-product of the reverse mapping
  cov_matrix[p*n_series] Gamma; 
  {
    matrix[n_series, n_series] P[p];
    matrix[n_series, n_series] phiGamma[2, p];
    for(i in 1:p) P[i] = AtoP(A[i]);
    phiGamma = rev_mapping(P, Sigma);
    phi = phiGamma[1];
    for(i in 1:p) {
      for(j in 1:p) {
        if(i<=j) Gamma[((i-1)*n_series+1):(i*n_series), ((j-1)*n_series+1):(j*n_series)] = phiGamma[2, j-i+1];
        else Gamma[((i-1)*n_series+1):(i*n_series), ((j-1)*n_series+1):(j*n_series)] = phiGamma[2, i-j+1]';
    }
   }
  }
  
  // reconstruct trends in matrix form
  for (t in 1:p){
    trend[t, 1:n_series] = to_row_vector(trend_raw_init[((t-1)*n_series+1):(t*n_series)]);
  }
  for (t in (p+1):n){
    trend[t, 1:n_series] = to_row_vector(trend_raw[t - p]);
  }
  
  // GAM component betas
  b[1:9] = mu_raw[1] + b_raw[1:9] * sigma_raw[1];
}
model {
  // mean trend estimates for times 1:p
  vector[p*n_series] mu_t_init;
  
  // mean trend estimates for times (p+1):n
  vector[n_series] mu_t[n-p];
  
  // shrinkage hyperpriors for unconstrained A's induce:
  // prior marginal means = 0;
  // prior marginal SDs = 0.5;
  // prior diagonal correlations = 0.7;
  // prior off-diagonal correlations = 0.7;
  Amu ~ normal(0, sqrt(0.35));
  Aomega ~ gamma(1.05, 0.0075); 
  
  // unconstrained A matrices
  for(s in 1:p) {
    diagonal(A[s]) ~ normal(Amu[s], 1 / sqrt(Aomega[s]));
    for(i in 1:n_series) {
      for(j in 1:n_series) {
        if(i != j) A[s, i, j] ~ normal(Amu[s], 1 / sqrt(Aomega[s]));
      }
    }
  }
  
  // trend covariance
  Sigma ~ inv_wishart(df, scale_mat);
  
  // initial trend means
  for(t in 1:p){
    mu_t_init[((t-1)*n_series+1):(t*n_series)] = mu;
  } 
  
  // remaining trend means
  for(t in (p+1):n) {
    mu_t[t-p] = mu;
    for(i in 1:p) {
      mu_t[t-p] += phi[i] * (trend_raw[t-i] - mu);
    }
  }
  
  // raw trends (in vector form)
  trend_raw_init ~ multi_normal(mu_t_init, Gamma);
  trend_raw ~ multi_normal(mu_t, Sigma);
  
  // prior for random effect population variances
  sigma_raw ~ exponential(0.5);
  // prior for random effect population means
  mu_raw ~ std_normal();
  // prior (non-centred) for s(ndvi_ma12,series)...
  b_raw[1:9] ~ std_normal();
  // likelihood functions
  vector[n_nonmissing] flat_trends;
  flat_trends = (to_vector(trend))[obs_ind];
  flat_ys ~ poisson_log_glm(append_col(flat_xs, flat_trends),
  0.0, append_row(b, 1.0));
}
generated quantities {
  vector[total_obs] eta;
  matrix[n, n_series] mus;
  array[n, n_series] int ypred;
  
  // posterior predictions
  eta = X * b;
  for(s in 1:n_series){ 
  mus[1:n, s] = eta[ytimes[1:n, s]] + trend[1:n, s];
  ypred[1:n, s] = poisson_log_rng(mus[1:n, s]);
  }
}
