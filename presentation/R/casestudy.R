library(stdmatern)
library(cmdstanr)
library(tidyverse)
library(evd)
library(patchwork)
library(here)

to_tibble <- function(y, value_label = "value") {
  y |>
    as.data.frame() |>
    as_tibble() |>
    mutate(
      lat = rep(seq_len(dim[1]), each = dim[2]),
      lon = rep(seq_len(dim[2]), times = dim[1])
    ) |>
    pivot_longer(
      c(-lat, -lon),
      names_to = "obs",
      values_to = value_label,
      names_prefix = "V",
      names_transform = as.integer
    )
}

theme_set(bggjphd::theme_bggj())
dim <- c(30, 20)
rho <- c(0.7, 0.9)
nu <- 2
n_obs <- 10
mu <- 6
sigma <- 2
xi <- 0.1

Z <- rmatern_copula(n_obs, dim, rho, nu)
U <- pnorm(Z)
Y <- qgev(U, loc = mu, scale = sigma, shape = xi)

model <- cmdstanr::cmdstan_model(
  here::here("casestudy", "Stan", "folded.stan"),
  include_paths = here::here("stanfunctions"),
  quiet = TRUE,
  force_recompile = TRUE
)

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
  xi = pmax(gev_fit[3], 1e-3),
  rho = c(0.5, 0.5)
)


results <- exact$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  init = list(inits, inits, inits, inits)
)
