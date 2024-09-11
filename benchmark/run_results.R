library(cmdstanr)
library(stdmatern)
library(tidyverse)
library(evd)

box::use(
  benchmark / functions[
    fit_all
  ]
)


crossing(
  dim = c(10),
  n_obs = c(10, 20),
  rho = seq(0.01, 0.99, length.out = 10),
  nu = c(0),
  mu = 6,
  sigma = 2,
  xi = 0.1
) |>
  mutate(
    grid_size = dim^2
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
      dim = c(.$dim, .$dim),
      rho = c(.$rho, .$rho),
      nu = .$nu,
      pars = list(mu = .$mu, sigma = .$sigma, xi = .$xi)
    )
  )
