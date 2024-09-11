functions {
#include gev.stanfunctions
#include gev_margins.stanfunctions
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
  real<lower = 0> sigma;
  real<lower = 0> xi;
  real<lower = 0, upper = min_y + sigma/xi> mu;
}

model {

  target += gev_lpdf(to_vector(y) | mu, sigma, xi);

  // Priors
  target += exponential_lpdf(xi | 1);
}


generated quantities {
  real log_lik = 0;
  log_lik += gev_lpdf(to_vector(y_test) | mu, sigma, xi);
}
