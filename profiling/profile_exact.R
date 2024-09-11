library(cmdstanr)
library(stdmatern)
library(evd)
library(tidyverse)
library(loo)
exact <- cmdstan_model(here::here("profiling", "Stan", "exact.stan"))

dim <- c(12, 12)
rho <- c(0.7, 0.5)
nu <- 2
n_obs <- 1
mu <- 6
sigma <- 2
xi <- 0.1

Z <- rmatern_copula(n_obs, dim, rho, nu)
U <- pnorm(Z)
Y <- qgev(U, loc = mu, scale = sigma, shape = xi)



gev_fit <- evd::fgev(Y)$estimate
inits <- list(
  mu = gev_fit[1],
  sigma = gev_fit[2],
  xi = pmax(gev_fit[3], 1e-3),
  rho = c(0.1, 0.1)
)

results_nu0 <- exact$sample(
  data = list(
    dim1 = dim[1],
    dim2 = dim[2],
    nu = 0,
    n_obs = n_obs,
    y = Y
  ),
  chains = 4,
  parallel_chains = 4,
  init = list(inits, inits, inits, inits)
)

results_nu0$summary(c("rho", "mu", "sigma", "xi"))

results_nu1 <- exact$sample(
  data = list(
    dim1 = dim[1],
    dim2 = dim[2],
    nu = 1,
    n_obs = n_obs,
    y = Y
  ),
  chains = 4,
  parallel_chains = 4,
  init = list(inits, inits, inits, inits)
)

results_nu2 <- exact$sample(
  data = list(
    dim1 = dim[1],
    dim2 = dim[2],
    nu = 2,
    n_obs = n_obs,
    y = Y
  ),
  chains = 4,
  parallel_chains = 4,
  init = list(inits, inits, inits, inits)
)

results_nu0$summary(c("rho", "mu", "sigma", "xi"))
results_nu1$summary(c("rho", "mu", "sigma", "xi"))
results_nu2$summary(c("rho", "mu", "sigma", "xi"))

loo_nu0 <- loo(results_nu0$draws("log_lik"))
loo_nu1 <- loo(results_nu1$draws("log_lik"))
loo_nu2 <- loo(results_nu2$draws("log_lik"))

loo_compare(
  list(
    "nu0" = loo_nu0,
    "nu1" = loo_nu1,
    "nu2" = loo_nu2
  )
)
