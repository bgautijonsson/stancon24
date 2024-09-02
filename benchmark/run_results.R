library(cmdstanr)
library(stdmatern)
library(tidyverse)
library(evd)

box::use(
  benchmarks / gev_margins / functions[
    fit_all
  ]
)


crossing(
  dim1 = 40,
  dim2 = 40,
  n_obs = 1,
  rho1 = 0.5,
  rho2 = 0.5,
  nu = c(0, 1, 2),
  mu = 6,
  sigma = 2,
  xi = 0.1
) |>
  mutate(
    grid_size = dim1 * dim2
  ) |>
  arrange(grid_size, n_obs, nu) |>
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
