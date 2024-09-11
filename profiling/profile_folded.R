library(cmdstanr)
library(stdmatern)
library(evd)
library(tidyverse)
folded <- cmdstan_model(here::here("profiling", "Stan", "folded.stan"))

dim <- c(20, 40)
rho <- c(0.5, 0.9)
nu <- 2
n_obs <- 5
mu <- 6
sigma <- 2
xi <- 0.1

Z <- rmatern_copula_folded_full(n_obs, dim[1], dim[2], rho[1], rho[2], nu)
U <- pnorm(Z)
Y <- qgev(U, loc = mu, scale = sigma, shape = xi)


stan_data <- list(
  dim1 = dim[1],
  dim2 = dim[2],
  nu = nu,
  n_obs = n_obs,
  y = Y
)

gev_fit <- evd::fgev(Y)$estimate
inits <- list(
  mu = gev_fit[1],
  sigma = gev_fit[2],
  xi = pmax(gev_fit[3], 1e-2),
  rho = c(0.5, 0.5)
)

results <- folded$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  init = list(inits, inits, inits, inits)
)

results$summary()
results$time()

results$profiles() |>
  purrr::map(\(x) select(x, name, total_time, forward_time, reverse_time)) |>
  purrr::list_rbind() |>
  dplyr::as_tibble() |>
  group_by(name) |>
  summarise_at(
    vars(total_time, forward_time, reverse_time),
    mean
  )
