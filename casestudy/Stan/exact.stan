functions {
#include exact.stanfunctions
#include gev.stanfunctions
#include gev_margins.stanfunctions
}

data {
  int<lower = 1> dim1;
  int<lower = 1> dim2;
  int<lower = 0> nu;
  int<lower = 1> n_obs;
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

model {
  matrix[D, n_obs] Z;

  profile("gev") {
    target += gev_lpdf(to_vector(y) | mu, sigma, xi);
    for (i in 1:D) {
      for (j in 1:n_obs) {
        Z[i, j] = inv_Phi(gev_cdf(y[i, j] | mu, sigma, xi));
      }
    }
  }

  profile("matern_copula_exact_lpdf") {
    target += matern_copula_exact_lpdf(Z | dim1, rho[1], dim2, rho[2], nu);
  }

  // Priors
  target += exponential_lpdf(xi | 1);
  target += beta_lpdf(rho | 1, 1);
}
