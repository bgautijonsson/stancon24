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

model {
  matrix[D, n_obs] Z;

  target += gev_lpdf(to_vector(y) | mu, sigma, xi);
  for (i in 1:D) {
    for (j in 1:n_obs) {
      Z[i, j] = inv_Phi(gev_cdf(y[i, j] | mu, sigma, xi));
    }
  }

  target += matern_copula_exact_lpdf(Z | dim1, rho[1], dim2, rho[2], nu);

  // Priors
  target += exponential_lpdf(xi | 1);
  target += std_normal_lpdf(rho);
}


generated quantities {
  real log_lik = 0;
  {
    matrix[D, n_obs] Z;
    
    log_lik += gev_lpdf(to_vector(y_test) | mu, sigma, xi);
    for (i in 1:D) {
      for (j in 1:n_obs) {
        Z[i, j] = inv_Phi(gev_cdf(y_test[i, j] | mu, sigma, xi));
      }
    }

    log_lik += matern_copula_exact_lpdf(Z | dim1, rho[1], dim2, rho[2], nu);
  }
}
