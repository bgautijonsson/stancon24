library(tidyverse)
library(here)
library(arrow)
library(gt)
library(gtExtras)
library(scales)
theme_set(bggjphd::theme_bggj())

d <- here("casestudies", "gev_margins", "results") |>
  open_dataset() |>
  collect()

d |>
  select(iter, model, time, n_obs, dim1, dim2, nu) |>
  distinct() |>
  arrange(iter) |>
  mutate(
    grid_size = dim1 * dim2,
    time_obs = time / n_obs
  ) |>
  ggplot(aes(grid_size, time_obs, color = model)) +
  geom_line() +
  geom_point()

d |>
  select(iter, model, time, n_obs, dim1, dim2, nu) |>
  distinct() |>
  arrange(iter) |>
  mutate(
    grid_size = dim1 * dim2
  ) |> 
  pivot_wider(names_from = model, values_from = time) |> 
  pivot_longer(c(circulant, folded)) |> 
  mutate(
    diff = value / exact
  ) |> 
  ggplot(aes(grid_size, diff, color = name)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_smooth(method = "lm") +
  geom_point() +
  scale_x_continuous(
    trans = "log10",
    breaks = breaks_log()
  ) +
  scale_y_continuous(
    trans = "log10",
    labels = \(x) scales::percent(x - 1),
    breaks = breaks_log()
  ) +
  facet_wrap("nu", labeller = label_both)



d |>
  filter(!is.na(value)) |>
  mutate(err = median - value) |>
  summarise(
    mean = median(abs(err)),
    sd = sd(err),
    .by = c(model, variable)
  ) |>
  pivot_wider(names_from = variable, values_from = c(mean, sd)) |>
  select(
    model,
    contains("mu"),
    contains("sigma"),
    contains("xi"),
    contains("rho[1]"),
    contains("rho[2]")
  ) |>
  gt() |>
  cols_label(
    model = "",
    mean_mu = "Mean",
    mean_sigma = "Mean",
    mean_xi = "Mean",
    `mean_rho[1]` = "Mean",
    `mean_rho[2]` = "Mean",
    sd_mu = "SD",
    sd_sigma = "SD",
    sd_xi = "SD",
    `sd_rho[1]` = "SD",
    `sd_rho[2]` = "SD",
  ) |>
  tab_spanner(
    label = "Mu",
    columns = 2:3
  ) |>
  tab_spanner(
    label = "Sigma",
    columns = 4:5
  ) |>
  tab_spanner(
    label = "Xi",
    columns = 6:7
  ) |>
  tab_spanner(
    label = "Rho 1",
    columns = 8:9
  ) |>
  tab_spanner(
    label = "Rho 2",
    columns = 10:11
  ) |>
  fmt_number(decimals = 3) |>
  gt_color_rows(-1, palette = "Greys")




d |>
  filter(variable == "log_lik") |>
  mutate(
    mean_ll = median / (n_obs * dim1 * dim2)
  ) |>
  select(iter, model, mean_ll) |>
  summarise(
    mean = mean(mean_ll),
    .by = model
  )

d |>
  filter(variable %in% c("log_lik", "rho[1]", "rho[2]")) |>
  select(iter, model, variable, median, n_obs, dim1, dim2) |>
  pivot_wider(names_from = variable, values_from = median) |> 
  janitor::clean_names() |> 
  distinct() |> 
  mutate(
    mean_ll = log_lik / (n_obs * dim1 * dim2)
  ) |> 
  ggplot(aes(x = rho_1, y = mean_ll, color = model)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  labs(x = "Mean Rho[1]", y = "Mean Log-Likelihood", title = "Mean Log-Likelihood vs Rho[1]")
