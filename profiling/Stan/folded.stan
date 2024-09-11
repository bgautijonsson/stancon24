functions {
  /*
  Folds a vector representing a 2D grid into a larger vector with symmetric repetitions.
  This function is used to create a folded version of the data for the circulant approximation
  of a Matérn covariance.

  @param x The input vector representing the original 2D grid data
  @param n1 The first dimension of the original grid
  @param n2 The second dimension of the original grid

  @return A vector of length 4 * n1 * n2 containing the folded data
  */
  vector fold_data(vector x, int n1, int n2) {
    vector[4 * n1 * n2] folded;
    for (i in 1:n1) {
      for (j in 1:n2) {
        int idx = (i - 1) * n2 + j;
        folded[(i - 1) * 2 * n2 + j] = x[idx];
        folded[(i - 1) * 2 * n2 + (2 * n2 - j + 1)] = x[idx];
        folded[(2 * n1 - i) * 2 * n2 + j] = x[idx];
        folded[(2 * n1 - i) * 2 * n2 + (2 * n2 - j + 1)] = x[idx];
      }
    }
    return folded;
  }
  /*
  Creates a base matrix for the circulant approximation of a Matérn covariance
  defined as a Kronecker sum of two AR(1) processes approximated with 
  circulant matrices.

  @param dim1 The first dimension of the grid
  @param dim2 The second dimension of the grid
  @param rho1 The correlation parameter for the first dimension
  @param rho2 The correlation parameter for the second dimension

  @return A matrix representing the base for the circulant approximation
  */
  matrix make_base_matrix(int dim1, int dim2, real rho1, real rho2) {
    matrix[dim2, dim1] c = rep_matrix(0, dim2, dim1);
    vector[dim1] c1 = rep_vector(0, dim1);
    vector[dim2] c2 = rep_vector(0, dim2);

    real scale1 = 1.0 / (1.0 - square(rho1));  
    real scale2 = 1.0 / (1.0 - square(rho2));

    c1[1] = 1 + square(rho1);
    c1[2] = -rho1;
    c1[dim1] = -rho1;
    c1 *= scale1;

    c2[1] = 1 + square(rho2);
    c2[2] = -rho2;
    c2[dim2] = -rho2;
    c2 *= scale2;

    // Set the first row
    c[1, 1] = c1[1] + c2[1];  
    c[1, 2] = c1[2];
    c[1, dim1] = c1[dim1];

    // Set the rest of the first column
    c[2, 1] = c2[2];
    c[dim2, 1] = c2[dim2];
    return c;
  }

  /*
  Creates a base matrix and rescales its eigenvalues for the circulant approximation.
  The eigenvalues are rescaled so that the inverse of the precision matrix is
  a correlation matrix.

  The function returns the square root of the eigenvalues for use in quadratic forms.

  @param dim1 The first dimension of the grid
  @param dim2 The second dimension of the grid
  @param rho1 The correlation parameter for the first dimension
  @param rho2 The correlation parameter for the second dimension
  @param nu The smoothness parameter of the Matérn covariance

  @return A complex matrix of square roots of rescaled eigenvalues
  */
  complex_matrix create_base_matrix_and_rescale_eigenvalues(int dim1, int dim2, real rho1, real rho2, int nu) {
    
    matrix[dim2, dim1] c = make_base_matrix(dim1, dim2, rho1, rho2);

    // Compute the eigenvalues and marginal standard deviation
    complex_matrix[dim2, dim1] eigs = pow(fft2(c), (nu + 1.0));
    complex_matrix[dim2, dim1] inv_eigs = pow(eigs, -1);
    real mvar = get_real(inv_fft2(inv_eigs)[1, 1]);
    eigs *= mvar;
    
    return eigs;
  }

  /*
  Performs matrix-vector multiplication in the Fourier domain.

  @param sqrt_eigs The square root of the eigenvalues
  @param v The vector to multiply

  @return The result of the matrix-vector product
  */
  vector matvec_prod(complex_matrix eigs, vector v) {
    int dim2 = rows(eigs);
    int dim1 = cols(eigs);
    complex_matrix[dim2, dim1] v_mat = to_matrix(v, dim2, dim1);
    
    complex_matrix[dim2, dim1] fft_v = fft2(v_mat);
    complex_matrix[dim2, dim1] prod = eigs .* fft_v;
    complex_matrix[dim2, dim1] result_complex = inv_fft2(prod);
    
    return to_vector(get_real(result_complex));
  }

  /*
  Computes the log probability density of a Matérn copula that's approximated
  with a block circulant matrix with circulant blocks.

  @param Z The vector of standard normal variates
  @param sqrt_eigs The square root of the eigenvalues of the precision matrix

  @return The log probability density
  */
  real matern_folded_copula_lpdf(matrix Z, int dim1, int dim2, real rho1, real rho2, int nu) {
    int n_obs = cols(Z);
    complex_matrix[2 * dim2, 2 * dim1] eigs = create_base_matrix_and_rescale_eigenvalues(2 * dim1, 2 * dim2, rho1, rho2, nu);
    real quad_forms = 0;
    real log_det = sum(log(get_real(eigs)));
    for (i in 1:n_obs) {
      vector[4 * dim1 * dim2] Z_fold = fold_data(Z[, i], dim1, dim2);
      vector[4 * dim1 * dim2] Qz = matvec_prod(eigs, Z_fold);
      quad_forms += dot_product(Z_fold, Qz) - dot_self(Z_fold);
    } 
    return - 0.5 * (quad_forms - n_obs * log_det);
  }

  /*
    This function calculates the conditional log probability of each observation given
    the rest of the observations, p(y_i | y_{-i}). We do this by calculating the conditional mean and
    conditional variance of the observation given the rest of the observations using the folded 
    circulant approximation.
    Conditional mean, mu_i|(-i) = y_i - {Qy}_{i}
    Conditional variance, sigma_i|(-i)^2 = 1 / {Q_{ii}}
    First we calculate g, where g_i = Qy = V diag(lambda) V' y,
    then we calculate sigma_tilde where sigma_tilde_ii = 1 / Q_{ii} = Hadamard(V, V) * lambda
  */
  vector matern_cond_loglik(vector y, int dim1, real rho1, int dim2, real rho2, int nu) {
    
    
  }
  
  /*
  This function computes the CDF of the Generalized Extreme Value distribution.
  */
  real gev_cdf(real y, real mu, real sigma, real xi) {
    if (abs(xi) < 1e-8) {
      real z = (y - mu) / sigma;
      return exp(-exp(z));
    } else {
      real z = 1 + xi * (y - mu) / sigma;
      return exp(-pow(z, -1/xi));
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
}

data {
  int<lower = 1> dim1;
  int<lower = 1> dim2;
  int<lower=0> nu;  
  int<lower=1> n_obs;
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

  profile("gef_pdf") {
    target += gev_lpdf(to_vector(y) | mu, sigma, xi);
  }

  profile("gev_to_gauss") {
    for (i in 1:D) {
      for (j in 1:n_obs) {
        Z[i, j] = inv_Phi(gev_cdf(y[i, j] | mu, sigma, xi));
      }
    }
  }

  profile("folded") { 
    target += matern_folded_copula_lpdf(Z | dim1, dim2, rho[1], rho[2], nu) / 4;
  }

  // Priors
  target += exponential_lpdf(xi | 1);
  target += std_normal_lpdf(rho);
}
