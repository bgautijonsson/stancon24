functions {
#include exact.stanfunctions
#include gev.stanfunctions
#include gev_margins.stanfunctions
}

data {
  int dim1;
  int dim2;
  int nu;
  int n_obs;
  matrix[dim1 * dim2, n_obs] y;
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
  tuple(matrix[dim1, dim1], vector[dim1]) E1 = ar1_precision_eigen(dim1, rho[1]);
  tuple(matrix[dim2, dim2], vector[dim2]) E2 = ar1_precision_eigen(dim2, rho[2]);
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
  target += priors(mu, sigma, xi, rho);
}
