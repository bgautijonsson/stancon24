#' Fit i.i.d. model
#'
#' This function fits an i.i.d. model using Stan,
#' assuming independence between observations
#' conditional on the GEV parameters
#'
#' @param stan_data A list containing the data for Stan
#' @param inits A list of initial values for the parameters
#' @param recomp Logical, whether to recompile the Stan model
#'
#' @return A tibble with parameter estimates and diagnostics
#' @export
fit_iid <- function(stan_data, inits, recomp = FALSE) {
  iid <- cmdstanr::cmdstan_model(
    here::here("benchmark", "Stan", "iid.stan"),
    include_paths = here::here("stanfunctions"),
    force_recompile = recomp,
    quiet = TRUE
  )

  iid_results <- iid$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    init = list(inits, inits, inits, inits)
  )

  iid_results$summary(c("mu", "xi", "sigma", "log_lik")) |>
    dplyr::select(variable, median, rhat, ess_bulk, ess_tail) |>
    dplyr::mutate_if(is.numeric, round, 2) |>
    dplyr::mutate(
      time = iid_results$time()$total,
      model = "iid"
    )
}




#' Fit exact Matérn model
#'
#' This function fits an exact Matérn model using Stan
#'
#' @param stan_data A list containing the data for Stan
#' @param inits A list of initial values for the parameters
#' @param recomp Logical, whether to recompile the Stan model
#'
#' @return A tibble with parameter estimates and diagnostics
#' @export
fit_exact <- function(stan_data, inits, recomp = FALSE) {
  exact <- cmdstanr::cmdstan_model(
    here::here("benchmark", "Stan", "exact.stan"),
    include_paths = here::here("stanfunctions"),
    force_recompile = recomp,
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

#' Fit folded Matérn model
#'
#' This function fits a folded approximation of a Matérn model using Stan
#'
#' @param stan_data A list containing the data for Stan
#' @param inits A list of initial values for the parameters
#' @param recomp Logical, whether to recompile the Stan model
#'
#' @return A tibble with parameter estimates and diagnostics
#' @export
fit_folded <- function(stan_data, inits, recomp = FALSE) {
  folded <- cmdstanr::cmdstan_model(
    here::here("benchmark", "Stan", "folded.stan"),
    include_paths = here::here("stanfunctions"),
    force_recompile = recomp,
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

#' Fit circulant Matérn model
#'
#' This function fits a circulant approximation of a Matérn model using Stan
#'
#' @param stan_data A list containing the data for Stan
#' @param inits A list of initial values for the parameters
#' @param recomp Logical, whether to recompile the Stan model
#'
#' @return A tibble with parameter estimates and diagnostics
#' @export
fit_circulant <- function(stan_data, inits, recomp = FALSE) {
  circulant <- cmdstanr::cmdstan_model(
    here::here("benchmark", "Stan", "circulant.stan"),
    include_paths = here::here("stanfunctions"),
    force_recompile = recomp,
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

#' Prepare data for Matérn models
#'
#' This function prepares data for fitting Matérn models with GEV margins
#'
#' @param n_obs Number of observations
#' @param dim A vector of length 2 specifying the dimensions of the spatial grid
#' @param rho A vector of length 2 specifying the spatial correlation parameters
#' @param nu The smoothness parameter of the Matérn covariance
#' @param pars A list containing GEV parameters (mu, sigma, xi)
#'
#' @return A list containing stan_data and initial values for the model
#' @export
prep_data <- function(n_obs, dim, rho, nu, pars) {
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
    xi = pmax(gev_fit[3], 1e-2),
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

  list(stan_data = stan_data, inits = inits)
}

#' Fit all Matérn models
#'
#' This function fits exact, circulant, and folded Matérn models
#'
#' @param n_obs Number of observations
#' @param dim A vector of length 2 specifying the dimensions of the spatial grid
#' @param rho A vector of length 2 specifying the spatial correlation parameters
#' @param nu The smoothness parameter of the Matérn covariance
#' @param pars A list containing GEV parameters (mu, sigma, xi)
#' @param recomp Logical, whether to recompile the Stan models
#'
#' @return A tibble with parameter estimates and diagnostics
#' @export
fit_all <- function(n_obs, dim, rho, nu, pars, recomp = FALSE) {
  data <- prep_data(n_obs, dim, rho, nu, pars)
  stan_data <- data$stan_data
  inits <- data$inits

  exact_results <- fit_exact(stan_data, inits, recomp)
  # circulant_results <- fit_circulant(stan_data, inits, recomp)
  # folded_results <- fit_folded(stan_data, inits, recomp)
  iid_results <- fit_iid(stan_data, inits, recomp)

  res <- dplyr::bind_rows(
    exact = exact_results,
    # circulant = circulant_results,
    # folded = folded_results,
    iid = iid_results
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
}

#' Save simulation results
#'
#' This function saves the results of a simulation iteration to a parquet file.
#'
#' @param res A tibble containing the results of a simulation iteration
#'
#' @return None (called for side effects)
#' @export
save_results <- function(res) {
  cur_results <- here::here(
    "benchmark",
    "results"
  ) |>
    list.files() |>
    readr::parse_number()

  max_iter <- max(cur_results, default = 0)

  out_path <- here::here(
    "benchmark",
    "results",
    glue::glue("iter={max_iter + 1}")
  )

  dir.create(out_path)

  arrow::write_parquet(
    res,
    here::here(out_path, "part-0.parquet")
  )
}
