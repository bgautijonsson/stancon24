library(cmdstanr)
library(stdmatern)
library(tidyverse)
library(evd)

box::use(
  casestudies / gev_margins / functions[run_iteration]
)




res <- purrr::map(
  seq_len(30),
  safely(run_iteration)
)

res

d <- res |>
  purrr::map(1) |>
  purrr::list_rbind()
