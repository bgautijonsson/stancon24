library(stdmatern)
library(tidyverse)
library(purrr)

dim1 <- 20
dim2 <- 20
rho1 <- 0.8
rho2 <- 0.3
nu <- 2
n <- 10

Z <- rmatern_copula(n, c(dim1, dim2), c(rho1, rho2), nu)

dmatern_copula_eigen(Z, dim1, dim2, rho1, rho2, nu)
dmatern_copula_folded(Z, dim1, dim2, rho1, rho2, nu)



loglik_fun <- function(pars) {
  -sum(dmatern_copula_folded(Z, dim1, dim2, plogis(pars[1]), plogis(pars[2]), nu))
}

fit_fun <- function(iter) {

  rho <- optim(
    rnorm(n = 2),
    loglik_fun
  )$par |> 
    plogis()

  tibble(
    iter = iter,
    par = c("rho1", "rho2"),
    value = rho
  )
}


d <- map(1:50, fit_fun)

d |> 
  list_rbind() |> 
  ggplot(aes(value, fill = par)) +
  geom_histogram(position = "identity") +
  facet_wrap("par", scales = "free")
