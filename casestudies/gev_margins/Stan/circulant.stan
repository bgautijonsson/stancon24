functions {
  #include gev.stanfunctions
  #include circulant.stanfunctions
}

data {
  int<lower=1> dim1;
  int<lower=1> dim2;
  int<lower=0> nu;  
  int<lower=1> n_obs;
  matrix[dim1*dim2, n_obs] y;
  matrix[dim1*dim2, n_obs] y_test;
}

transformed data {
  int D = dim1 * dim2;
  real min_y = min(y);
}

parameters {
  vector<lower=0, upper=1>[2] rho;
  real<lower = 0> sigma;
  real<lower = 0> xi;
  real<lower = 0, upper = min_y + sigma/xi> mu;
}

transformed parameters {
  matrix[dim2, dim1] c = create_base_matrix(dim1, dim2, rho[1], rho[2]);
  complex_matrix[dim2, dim1] eigs = compute_and_rescale_eigenvalues(c, nu);
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
    target += matern_circulant_copula_lpdf(to_vector(Z[ , i]) | eigs);
  }

  // Priors
  target += beta_lpdf(rho | 1, 1);
  target += std_normal_lpdf(xi);
  target += exponential_lpdf(sigma | 1);
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
  matrix[D, n_obs] Z = inv_Phi(u);
  for (i in 1:n_obs) {
    log_lik += matern_circulant_copula_lpdf(to_vector(Z[ , i]) | eigs);
  }
}
