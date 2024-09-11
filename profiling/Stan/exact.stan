functions {
  /*
    Computes the Kronecker product of two vectors a and b.

    @param a The first vector
    @param b The second vector

    @return A vector representing the Kronecker product of a and b
  */
  vector kronecker(vector a, vector b) {
    int n_a = rows(a);
    int n_b = rows(b);
    vector[n_a * n_b] result;
    
    for (i in 1:n_a) {
      for (j in 1:n_b) {
        result[(i-1)*n_b + j] = a[i] * b[j];
      }
    }
    
    return result;
  }

  /*
    Calculates the marginal standard deviations of the matrix Q, defined as the
    Kronecker sum of Q1 and Q2 with smoothness parameter nu.

    @param E1 Tuple containing the eigendecomposition of Q1
    @param E2 Tuple containing the eigendecomposition of Q2
    @param nu The smoothness parameter

    @return A vector of marginal standard deviations
  */
  vector marginal_sd(tuple(matrix, vector) E1, tuple(matrix, vector) E2, real nu) {
    int dim1 = cols(E1.1);
    int dim2 = cols(E2.1);

    matrix[dim1, dim1] V1 = E1.1;
    vector[dim1] A1 = E1.2;
    matrix[dim2, dim2] V2 = E2.1;
    vector[dim2] A2 = E2.2;

    vector[dim1 * dim2] marginal_sds = rep_vector(0.0, dim1 * dim2);
    for (i in 1:dim1) {
      for (j in 1:dim2) {
        vector[dim1 * dim2] v = kronecker(V1[, i], V2[, j]);
        real lambda = pow(A1[i] + A2[j], nu + 1);
        marginal_sds += square(v) / lambda;
      }
    }
    return sqrt(marginal_sds);
  }

  /*
    Computes the eigendecomposition of the precision matrix for an AR(1) process.

    @param n The size of the AR(1) process
    @param rho The correlation parameter

    @return A tuple containing the eigenvectors and eigenvalues of the precision matrix
  */
  tuple(matrix, vector) ar1_precision_eigen(int n, real rho) {
    matrix[n, n] Q;
    real scaling = 1.0 / (1.0 - rho * rho);
    real off_diag = -rho * scaling;

    Q = rep_matrix(0, n, n);
    for (i in 1:n) {
      Q[i, i] = (i == 1 || i == n) ? scaling : (1.0 + rho * rho) * scaling;
      if (i > 1) Q[i, i-1] = off_diag;
      if (i < n) Q[i, i+1] = off_diag;
    }
    
    return eigendecompose_sym(Q);
  }
  /*
    Computes the log probability density of a Matérn copula using the connection between 
    the eigendecomposition of the precision matrix and the smaller AR(1) precision matrices.

    @param Z The matrix of standard normal variates
    @param E1 Tuple containing the eigendecomposition of Q1
    @param E2 Tuple containing the eigendecomposition of Q2
    @param nu The smoothness parameter

    @return The log probability density
  */
  real matern_copula_exact_lpdf(matrix Z, int dim1, real rho1, int dim2, real rho2, int nu) {
    int n_obs = cols(Z);
    int D = dim1 * dim2;
    tuple(matrix[dim1, dim1], vector[dim1]) E1 = ar1_precision_eigen(dim1, rho1);
    tuple(matrix[dim2, dim2], vector[dim2]) E2 = ar1_precision_eigen(dim2, rho2);


    real log_det = 0;
    real quadform_sum = 0;

    vector[D] marginal_sds = marginal_sd(E1, E2, nu);

    for (i in 1:dim1) {
      for (j in 1:dim2) {
        vector[D] v = kronecker(E1.1[, i], E2.1[, j]);
        v = v .* marginal_sds;  
        real norm_v = sqrt(sum(square(v)));
        v /= norm_v;  
        
        real lambda = pow(E1.2[i] + E2.2[j], nu + 1) * square(norm_v);
        log_det += log(lambda);
        
        row_vector[n_obs] q = v' * Z;  
        quadform_sum += dot_self(q) * lambda;
      }
    }

    real z_squared = sum(columns_dot_self(Z));

    return -0.5 * (quadform_sum - n_obs * log_det - z_squared);
  }

  /*
    This function calculates the conditional log probability of each observation given
    the rest of the observations, p(y_i | y_{-i}). We do this by calculating the conditional mean and
    conditional variance of the observation given the rest of the observations using the
    precision matrix.
    Conditional mean, mu_i|(-i) = y_i - {Qy}_{i}
    Conditional variance, sigma_i|(-i)^2 = 1 / {Q_{ii}}
    First we calculate g, where g_i = Qy = V diag(lambda) V' y,
    then we calculate sigma_tilde where sigma_tilde_ii = 1 / Q_{ii} = Hadamard(V, V) * lambda
  */
  vector matern_cond_loglik(vector y, int dim1, real rho1, int dim2, real rho2, int nu) {
    int D = dim1 * dim2;
    tuple(matrix[dim1, dim1], vector[dim1]) E1 = ar1_precision_eigen(dim1, rho1);
    tuple(matrix[dim2, dim2], vector[dim2]) E2 = ar1_precision_eigen(dim2, rho2);

    vector[D] marginal_sds = marginal_sd(E1, E2, nu);
    vector[D] g = rep_vector(0, D);
    vector[D] tau_tilde = rep_vector(0, D);

    for (i in 1:dim1) {
      for (j in 1:dim2) {
        vector[D] v = kronecker(E1.1[, i], E2.1[, j]);
        v = v .* marginal_sds;
        real norm_v = sqrt(sum(square(v)));
        v /= norm_v;

        real lambda = pow(E1.2[i] + E2.2[j], nu + 1) * square(norm_v);
        
        g += v * lambda * v' * y;
        tau_tilde += square(v) * lambda;
      }
    }
    
    vector[D] loglik = - 0.5 * (tau_tilde .* pow(g, 2) - log(tau_tilde) - pow(y, 2));

    return loglik;
  }

  /*
  This function computes the CDF of the Generalized Extreme Value distribution.
  */
  real gev_cdf(real y, real mu, real sigma, real xi) {
    if (abs(xi) < 1e-10) {
      return exp(-exp((y - mu) / sigma));
    } else {
      return exp(-pow(1 + xi * (y - mu) / sigma, -1/xi));
    }
  }

  /*
  This function computes the log-pdf of the Generalized Extreme Value distribution.
  */
  real gev_lpdf(vector y, real mu, real sigma, real xi) {
    vector[rows(y)] z;
    if (abs(xi) < 1e-8) {
      z = (y - mu) / sigma;
      return sum(-log(sigma) - z - exp(-z));
    } else {
      z = 1 + xi * (y - mu) / sigma;
      return sum(-log(sigma) - (1 + 1/xi) * log(z) - pow(z, -1/xi));
    }
  }

  /*
  This function computes the log-pdf of the Generalized Extreme Value distribution.
  */
  vector gev_loglik(vector y, real mu, real sigma, real xi) {
    vector[rows(y)] z;
    if (abs(xi) < 1e-8) {
      z = (y - mu) / sigma;
      return -log(sigma) - z - exp(-z);
    } else {
      z = 1 + xi * (y - mu) / sigma;
      return -log(sigma) - (1 + 1/xi) * log(z) - pow(z, -1/xi);
    }
  }
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
  target += std_normal_lpdf(rho);
}

generated quantities {
  matrix[D, n_obs] log_lik = rep_matrix(0, D, n_obs);
  matrix[D, n_obs] Z;
  profile("loglik") {
    for (j in 1:n_obs) {
      log_lik[ , j] += gev_loglik(y[ , j], mu, sigma, xi);
    for (i in 1:D) {
        Z[i, j] = inv_Phi(gev_cdf(y[i, j] | mu, sigma, xi));
      }
      log_lik[ , j] += matern_cond_loglik(Z[ , j], dim1, rho[1], dim2, rho[2], nu);
    }
  }
}
