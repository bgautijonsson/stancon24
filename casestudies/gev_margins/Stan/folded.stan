functions {
  #include gev.stanfunctions
  #include folded.stanfunctions
  #include circulant.stanfunctions
  #include gev_margins.stanfunctions
}

data {
  int<lower = 1> dim1;
  int<lower = 1> dim2;
  int<lower=0> nu;  
  int<lower=1> n_obs;
  matrix[dim1 * dim2, n_obs] y;
  matrix[dim1 * dim2, n_obs] y_test;
}

transformed data {
  int D = dim1 * dim2;
  real min_y = min(y);
}

parameters {
  vector<lower = 0, upper = 1>[2] rho;
  real<lower = 0> sigma;
  real<lower = 0> xi;
  real<lower = 0, upper = min_y + sigma/xi> mu;
}

transformed parameters {
  complex_matrix[2 * dim2, 2 * dim1] eigs = create_base_matrix_and_rescale_eigenvalues(2 * dim1, 2 * dim2, rho[1], rho[2], nu);
}

model {
  matrix[D, n_obs] u;

  for (i in 1:D) {
    for (j in 1:n_obs) {
      u[i, j] = gev_cdf(y[i, j] | mu, sigma, xi);
      target+= gev_lpdf(y[i, j] | mu, sigma, xi);
    }
  }
  
  // Likelihood
  matrix[D, n_obs] Z = inv_Phi(u);
  for (i in 1:n_obs) {
    target += matern_circulant_copula_lpdf(fold_data(Z[ , i], dim1, dim2) | eigs) / 4.0;
  }

  // Priors
  target += priors(mu, sigma, xi, rho);
}

generated quantities {
  matrix[D, n_obs] u;
  real log_lik = 0;

  for (i in 1:D) {
    for (j in 1:n_obs) {
      u[i, j] = gev_cdf(y_test[i, j] | mu, sigma, xi);
      log_lik += gev_lpdf(y_test[i, j] | mu, sigma, xi);
    }
  }
  
  // Likelihood
  matrix[D, n_obs] Z = inv_Phi(u);
  for (i in 1:n_obs) {
    log_lik += matern_circulant_copula_lpdf(fold_data(Z[ , i], dim1, dim2) | eigs) / 4.0;
  }
}
