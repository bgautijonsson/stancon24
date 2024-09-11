library(cmdstanr)
library(stdmatern)
library(tidyverse)
library(evd)

box::use(
  benchmark / functions[
    fit_folded,
    prep_data,
    fit_circulant,
    fit_exact,
    fit_all
  ]
)

n_obs <- 1
dim <- c(10, 10)
rho <- c(0.7, 0.9)
nu <- 2
pars <- list(
  mu = 6,
  sigma = 2,
  xi = 0.1
)

data <- prep_data(n_obs, dim, rho, nu, pars)

exact <- fit_exact(data$stan_data, data$inits, recomp = TRUE)

exact

folded <- fit_folded(data$stan_data, data$ints, recomp = TRUE)

folded
