functions {
#include exact.stanfunctions
#include gev.stanfunctions
}

data {
  int dim1;
  int dim2;
  int nu;
  int n_obs;
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
  matrix[dim1, dim1] Q1 = ar1_precision(dim1, rho[1]);
  matrix[dim2, dim2] Q2 = ar1_precision(dim2, rho[2]);
  tuple(matrix[dim1, dim1], vector[dim1]) E1 = eigendecompose_sym(Q1);
  tuple(matrix[dim2, dim2], vector[dim2]) E2 = eigendecompose_sym(Q2);
}

model {
  matrix[D, n_obs] u;

  for (i in 1:D) {
    for (j in 1:n_obs) {
      u[i, j] = gev_cdf(y[i, j] | mu, sigma, xi);
      target+= gev_lpdf(y[i, j] | mu, sigma, xi);
    }
  }

  matrix[D, n_obs] Z = inv_Phi(u);

  target += matern_copula_exact_lpdf(Z | E1, E2, nu);

  // Priors
  target += beta_lpdf(rho | 1, 1);
  target += std_normal_lpdf(xi);
  target += exponential_lpdf(sigma | 1);
}

generated quantities {
  matrix[D, n_obs] u;
  matrix[D, n_obs] Z;
  real log_lik = 0;

  for (i in 1:D) {
    for (j in 1:n_obs) {
      u[i, j] = gev_cdf(y_test[i, j] | mu, sigma, xi);
      log_lik += gev_lpdf(y_test[i, j] | mu, sigma, xi);
    }
  }

  Z = inv_Phi(u);
  log_lik += matern_copula_exact_lpdf(Z | E1, E2, nu);
}
