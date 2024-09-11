functions {
  #include gev.stanfunctions
  #include ar1_copula.stanfunctions
}

data {
  int n_obs;
  vector[n_obs] y;
}

transformed data {
  real min_y = min(y);
}

parameters {
  real<lower = -1, upper = 1> rho;
  real<lower = 0> sigma;
  real<lower = 0> xi;
  real<lower = 0, upper = min_y + sigma / xi> mu;
}

model {
  vector[n_obs] U;
  target += gev_lpdf(y | mu, sigma, xi);
  for (i in 1:n_obs) {
    U[i] = gev_cdf(y[i] | mu, sigma, xi);
  }
  target += normal_copula_ar1_lpdf(U | rho);

  // Priors
  target += std_normal_lpdf(rho);
  target += exponential_lpdf(sigma | 1);
  target += exponential_lpdf(xi | 1);
}
