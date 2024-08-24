#' @export
fit_exact <- function(stan_data, inits, force_recompile = FALSE) {
  exact <- cmdstanr::cmdstan_model(
    here::here("casestudies", "gev_margins", "Stan", "exact.stan"),
    include_paths = here::here("stanfunctions"),
    force_recompile = force_recompile,
    quiet = TRUE
  )


  exact_results <- exact$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    init = list(inits, inits, inits, inits)
  )

  exact_results$summary(c("mu", "xi", "sigma", "rho", "log_lik")) |>
    dplyr::select(variable, median, rhat, ess_bulk, ess_tail) |>
    dplyr::mutate_if(is.numeric, round, 2) |>
    dplyr::mutate(
      time = exact_results$time()$total,
      model = "exact"
    )
}

#' @export
fit_folded <- function(stan_data, inits, force_recompile = FALSE) {
  folded <- cmdstanr::cmdstan_model(
    here::here("casestudies", "gev_margins", "Stan", "folded.stan"),
    include_paths = here::here("stanfunctions"),
    force_recompile = force_recompile,
    quiet = TRUE
  )

  folded_results <- folded$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    init = list(inits, inits, inits, inits)
  )

  folded_results$summary(c("mu", "xi", "sigma", "rho", "log_lik")) |>
    dplyr::select(variable, median, rhat, ess_bulk, ess_tail) |>
    dplyr::mutate_if(is.numeric, round, 2) |>
    dplyr::mutate(
      time = folded_results$time()$total,
      model = "folded"
    )
}

#' @export
fit_circulant <- function(stan_data, inits, force_recompile = FALSE) {
  circulant <- cmdstanr::cmdstan_model(
    here::here("casestudies", "gev_margins", "Stan", "circulant.stan"),
    include_paths = here::here("stanfunctions"),
    force_recompile = force_recompile,
    quiet = TRUE
  )

  circulant_results <- circulant$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    init = list(inits, inits, inits, inits)
  )

  circulant_results$summary(c("mu", "xi", "sigma", "rho", "log_lik")) |>
    dplyr::select(variable, median, rhat, ess_bulk, ess_tail) |>
    dplyr::mutate_if(is.numeric, round, 2) |>
    dplyr::mutate(
      time = circulant_results$time()$total,
      model = "circulant"
    )
}

#' @export
fit_all <- function(n_obs, dim, rho, nu, pars, force_recompile = FALSE) {
  z <- stdmatern::rmatern_copula(n_obs, dim, rho, nu)
  z_test <- stdmatern::rmatern_copula(n_obs, dim, rho, nu)
  u <- stats::pnorm(z)
  u_test <- stats::pnorm(z_test)
  y <- evd::qgev(
    u,
    loc = pars$mu,
    scale = pars$sigma,
    shape = pars$xi
  )

  y_test <- evd::qgev(
    u_test,
    loc = pars$mu,
    scale = pars$sigma,
    shape = pars$xi
  )

  gev_fit <- evd::fgev(y)$estimate
  inits <- list(
    mu = gev_fit[1],
    sigma = gev_fit[2],
    xi = pmax(gev_fit[3], 0),
    rho = c(0.5, 0.5)
  )

  stan_data <- list(
    dim1 = dim[1],
    dim2 = dim[2],
    n_obs = n_obs,
    y = y,
    y_test = y_test,
    nu = nu
  )

  exact_results <- fit_exact(stan_data, inits, force_recompile)
  circulant_results <- fit_circulant(stan_data, inits, force_recompile)
  folded_results <- fit_folded(stan_data, inits, force_recompile)

  dplyr::bind_rows(
    exact = exact_results,
    circulant = circulant_results,
    folded = folded_results
  )
}

#' @export
run_iteration <- function(iter, force_recompile = FALSE) {
  n_obs <- sample(10:40, size = 1)
  dim <- sample(5:25, size = 2, replace = TRUE)
  rho <- stats::runif(n = 2, min = 0.05, max = 0.95)
  nu <- sample(c(0, 1, 2), size = 1)

  pars <- list(
    mu = 6,
    sigma = 2,
    xi = 0.1
  )



  res <- fit_all(
    n_obs,
    dim,
    rho,
    nu,
    pars,
    force_recompile = force_recompile
  ) |>
    dplyr::left_join(
      dplyr::tibble(
        variable = c("mu", "sigma", "xi", "rho[1]", "rho[2]"),
        value = c(pars$mu, pars$sigma, pars$xi, rho[1], rho[2])
      ),
      by = dplyr::join_by(variable)
    ) |>
    dplyr::mutate(
      n_obs = n_obs,
      dim1 = dim[1],
      dim2 = dim[2],
      nu = nu
    )

  save_results(res)
  res
}

#' @export
save_results <- function(res) {
  cur_results <- here::here(
    "casestudies",
    "gev_margins",
    "results"
  ) |>
    list.files() |>
    readr::parse_number()

  max_iter <- max(cur_results, default = 0)

  out_path <- here::here(
    "casestudies",
    "gev_margins",
    "results",
    glue::glue("iter={max_iter + 1}")
  )

  dir.create(out_path)

  arrow::write_parquet(
    res,
    here::here(out_path, "part-0.parquet")
  )
}
