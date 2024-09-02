library(cmdstanr)
library(stdmatern)
library(tidyverse)
library(evd)

box::use(
  benchmarks / gev_margins / functions[
    fit_folded,
    prep_data,
    fit_circulant,
    fit_exact,
    fit_all
  ]
)
