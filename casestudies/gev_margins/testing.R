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


d <- crossing(
  n_obs = 1,
  dim1 = c(10, 15, 20, 25, 30),
  dim2 = c(10, 15, 20, 25, 30),
  rho1 = c(0.25, 0.5, 0.75),
  rho2 = c(0.25, 0.5, 0.75),
  nu = c(0, 1, 2),
  mu = 6,
  sigma = 2,
  xi = 0.1
) |>
  arrange(runif(n())) |>
  mutate(
    iter = row_number()
  ) |>
  group_by(iter) |>
  group_split() |>
  map(
    ~ fit_all(
      n_obs = .$n_obs,
      dim = c(.$dim1, .$dim2),
      rho = c(.$rho1, .$rho2),
      nu = .$nu,
      pars = list(mu = .$mu, sigma = .$sigma, xi = .$xi)
    )
  )
