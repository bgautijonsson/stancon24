library(cmdstanr)
library(stdmatern)
library(tidyverse)
library(evd)

box::use(
  casestudies / gev_margins / functions[
    fit_folded,
    prep_data,
    fit_circulant,
    fit_exact,
    fit_all
  ]
)

n_obs <- 10
dim <- c(10, 10)
rho <- c(0.8, 0.5)
nu <- 1
pars <- list(mu = 6, sigma = 2, xi = 0.1)

data <- prep_data(n_obs, dim, rho, nu, pars)


folded <- fit_folded(data$stan_data, data$inits, recomp = TRUE)
circulant <- fit_circulant(data$stan_data, data$inits, recomp = TRUE)
exact <- fit_exact(data$stan_data, data$inits, recomp = TRUE)

bind_rows(folded, circulant, exact)
