library(stdmatern)
library(cmdstanr)
library(sparseMVN)
library(tidyverse)


# Simulate data
Q <- make_AR_prec_matrix(dim = 60, rho = -0.8)
Z <- rmvn.sparse(n = 1, mu = rep(0, nrow(Q)), CH = Cholesky(Q)) |> as.numeric()
U <- pnorm(Z)
y <- qgev(U, loc = 10, scale = 5, shape = 0.1)

ml_fit <- evd::fgev(y)
inits <- list(
  mu = ml_fit$estimate[1],
  sigma = ml_fit$estimate[2],
  xi = pmax(ml_fit$estimate[3], 1e-3),
  rho = 0.1
)

stan_data <- list(
  n_obs = length(y),
  y = y
)

model <- cmdstan_model(
  here::here("1d-ar1", "Stan", "gev_ar1.stan"),
  include_paths = here::here("stanfunctions")
)

fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  init = list(inits, inits, inits, inits)
)

fit$summary()
