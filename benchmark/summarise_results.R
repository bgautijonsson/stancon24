library(tidyverse)
library(here)
library(arrow)
library(gt)
library(gtExtras)
library(scales)
library(ggh4x)
library(geomtextpath)
theme_set(bggjphd::theme_bggj())

d <- here("benchmark", "results") |>
  open_dataset() |>
  collect()

d |>
  mutate(
    rho1 = value[variable == "rho[1]"] |> max(),
    rho2 = value[variable == "rho[2]"] |> max(),
    rho_mean = (rho1 + rho2) / 2,
    .by = c(iter)
  ) |>
  filter(
    variable == "log_lik"
  ) |>
  mutate(
    mean_ll = median / (n_obs * dim1 * dim2),
  ) |>
  select(
    iter, model, mean_ll, rho1, rho2, rho_mean, nu, dim1
  ) |>
  arrange(iter) |>
  ggplot(aes(rho_mean, mean_ll, col = model, fill = model)) +
  geom_smooth(col = NA, alpha = 0.1) +
  geom_textsmooth(aes(label = model)) +
  scale_x_continuous(
    limits = c(0, 1),
    guide = guide_axis_truncated()
  ) +
  scale_y_continuous(
    guide = guide_axis_truncated()
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(cols = vars(nu)) +
  labs(
    x = expression(rho),
    y = "elppd"
  ) +
  theme(legend.position = "none")

d |>
  mutate(
    rho1 = value[variable == "rho[1]"] |> max(),
    rho2 = value[variable == "rho[2]"] |> max(),
    rho_mean = (rho1 + rho2) / 2,
    .by = c(iter)
  ) |>
  filter(
    variable == "log_lik"
  ) |>
  mutate(
    mean_ll = median / (n_obs * dim1 * dim2),
  ) |>
  select(
    iter, model, mean_ll, rho1, rho2, rho_mean, nu, dim1
  ) |>
  arrange(iter) |>
  pivot_wider(names_from = model, values_from = mean_ll) |>
  mutate(
    diff = exact - iid,
    nu = factor(nu)
  ) |>
  ggplot(aes(rho_mean, diff, col = nu, fill = nu, group = nu)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_textsmooth(
    aes(label = glue::glue("&nu; = {nu}")),
    rich = TRUE,
    size = 5
  ) +
  geom_smooth(col = NA, alpha = 0.2) +
  scale_x_continuous(
    limits = c(0, 1),
    guide = guide_axis_truncated()
  ) +
  scale_y_continuous(
    guide = guide_axis_truncated()
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x = expression(rho),
    y = expression(paste(elppd - elppd[i.i.d.]))
  ) +
  theme(legend.position = "none")

ggsave(
  filename = "diff_elppd.png",
  width = 8, height = 0.5 * 8, scale = 1.3
)

d |>
  mutate(
    rho1 = value[variable == "rho[1]"] |> max(),
    rho2 = value[variable == "rho[2]"] |> max(),
    .by = c(iter)
  ) |>
  filter(
    variable == "log_lik"
  ) |>
  mutate(
    mean_ll = median / (n_obs * dim1 * dim2),
  ) |>
  select(
    iter, model, mean_ll, rho1, rho2, nu
  ) |>
  arrange(iter) |>
  pivot_wider(names_from = model, values_from = mean_ll) |>
  mutate(
    diff = exact - iid
  ) |>
  drop_na(diff) |>
  ggplot(aes(rho1, rho2, fill = diff)) +
  geom_raster(
    interpolate = TRUE
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(0)
  ) +
  scale_fill_distiller(palette = "Spectral", direction = -1)

d |>
  select(iter, model, time, n_obs, dim1, dim2, nu) |>
  distinct() |>
  arrange(iter) |>
  mutate(
    grid_size = dim1 * dim2,
    time_obs = time / n_obs,
    model = str_to_title(model) |>
      fct_relevel("Exact")
  ) |>
  ggplot(aes(grid_size, time_obs, color = model, fill = model, group = model)) +
  geom_smooth(method = "lm") +
  geom_point() +
  scale_x_continuous(
    trans = "log10",
    labels = label_number()
  ) +
  scale_y_continuous(
    trans = "log10",
    labels = label_timespan(),
    breaks = c(10, 60, 15 * 60, 3 * 60 * 60, 24 * 60 * 60)
  ) +
  scale_colour_manual(
    values = c(
      "#525252",
      "#e41a1c",
      "#377eb8"
    )
  ) +
  scale_fill_manual(
    values = c(
      "#525252",
      "#e41a1c",
      "#377eb8"
    )
  ) +
  theme(
    legend.position = "top"
  ) +
  labs(
    x = "Grid size",
    y = NULL,
    title = "How long does it take Stan to sample from the model?",
    fill = NULL,
    col = NULL
  )

d |>
  select(iter, model, time, n_obs, dim1, dim2, nu) |>
  distinct() |>
  arrange(iter) |>
  mutate(
    grid_size = dim1 * dim2
  ) |>
  pivot_wider(names_from = model, values_from = time) |>
  mutate(
    diff = exact / iid
  ) |>
  ggplot(aes(grid_size, diff)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_smooth(method = "lm") +
  scale_x_continuous(
    trans = "log10",
    breaks = breaks_log(6),
    guide = guide_axis_logticks()
  ) +
  scale_y_continuous(
    trans = "log10",
    labels = \(x) scales::number(x - 1, suffix = "x"),
    breaks = breaks_log(8),
    guide = guide_axis_logticks()
  ) +
  scale_colour_brewer(
    palette = "Set1"
  ) +
  scale_fill_brewer(
    palette = "Set1"
  ) +
  # facet_wrap("nu", labeller = label_both) +
  theme(
    legend.position = "top",
    plot.margin = margin(t = 5, r = 15, b = 5, l = 15)
  ) +
  labs(
    x = "Grid size",
    y = NULL,
    subtitle = "Sampling time compared to the exact method (% difference)",
    title = "Comparing the speed of approximations to exact method",
    col = NULL,
    fill = NULL
  )

ggsave(
  here("benchmarks", "gev_margins", "figures", "speed_comparison.png"),
  scale = 1.4, width = 8, height = 0.621 * 8
)


d |>
  filter(!is.na(value)) |>
  mutate(err = median - value) |>
  summarise(
    mean = median(abs(median - value)),
    sd = mean((median - value)^2),
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
    mean_mu = "MAD",
    mean_sigma = "MAD",
    mean_xi = "MAD",
    `mean_rho[1]` = "MAD",
    `mean_rho[2]` = "MAD",
    sd_mu = "MSE",
    sd_sigma = "MSE",
    sd_xi = "MSE",
    `sd_rho[1]` = "MSE",
    `sd_rho[2]` = "MSE",
  ) |>
  tab_spanner(
    label = md("$\\mu$"),
    columns = 2:3
  ) |>
  tab_spanner(
    label = md("$\\sigma$"),
    columns = 4:5
  ) |>
  tab_spanner(
    label = md("$\\xi$"),
    columns = 6:7
  ) |>
  tab_spanner(
    label = md("$\\rho_1$"),
    columns = 8:9
  ) |>
  tab_spanner(
    label = md("$\\rho_2$"),
    columns = 10:11
  ) |>
  fmt_number(decimals = 3)


d |>
  select(iter, model, variable, median, value) |>
  mutate(
    err = value / median,
    model = fct_relevel(
      model,
      "folded",
      "circulant",
      "exact"
    ) |>
      fct_recode(
        "Folded" = "folded",
        "Circulant" = "circulant",
        "Exact" = "exact"
      ),
    variable = fct_relevel(
      variable,
      "rho[2]",
      "rho[1]",
      "xi",
      "sigma",
      "mu"
    )
  ) |>
  drop_na() |>
  reframe(
    p = seq(0.025, 0.475, 0.025),
    size = 1 - 2 * p,
    lower = quantile(err, p),
    upper = quantile(err, 1 - p),
    .by = c(model, variable)
  ) |>
  ggplot() +
  geom_vline(xintercept = 1, lty = 2) +
  geom_segment(
    aes(
      x = lower, xend = upper,
      y = model,
      col = size, alpha = -size, linewidth = -size,
      group = model
    )
  ) +
  scale_x_continuous(
    trans = "log10",
    labels = \(x) percent(x - 1, style_positive = "plus"),
    breaks = breaks_log(8),
    guide = guide_axis_truncated()
  ) +
  scale_y_discrete(
    guide = guide_axis_truncated()
  ) +
  scale_colour_distiller(
    guide = "none",
    palette = "Blues"
  ) +
  scale_alpha_continuous(
    guide = "none",
    range = c(0.1, 0.6)
  ) +
  scale_linewidth_continuous(
    guide = "none",
    range = c(0.1, 5)
  ) +
  facet_grid(
    rows = vars(variable),
    labeller = label_parsed
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title = "Do the approximations capture the true parameters well enough?"
  )


ggsave(
  here("benchmarks", "gev_margins", "figures", "estimate_comparison_log.png"),
  scale = 1.4, width = 8, height = 0.4 * 8
)

d |>
  filter(
    rhat < 1.1
  ) |>
  select(iter, model, variable, median, value) |>
  mutate(
    err = value - median,
    model = fct_relevel(
      model,
      "folded",
      "circulant",
      "exact"
    ) |>
      fct_recode(
        "Folded" = "folded",
        "Circulant" = "circulant",
        "Exact" = "exact"
      ),
    variable = fct_relevel(
      variable,
      "rho[2]",
      "rho[1]",
      "xi",
      "sigma",
      "mu"
    )
  ) |>
  drop_na() |>
  reframe(
    p = seq(0.025, 0.475, 0.025),
    size = 1 - 2 * p,
    lower = quantile(err, p),
    upper = quantile(err, 1 - p),
    .by = c(model, variable)
  ) |>
  ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_segment(
    aes(
      x = lower, xend = upper,
      y = model,
      col = size, alpha = -size, linewidth = -size,
      group = model
    )
  ) +
  scale_x_continuous(
    labels = label_number(style_positive = "plus"),
    breaks = breaks_pretty(),
    guide = guide_axis_truncated()
  ) +
  scale_y_discrete(
    guide = guide_axis_truncated()
  ) +
  scale_colour_distiller(
    guide = "none",
    limits = c(0, 1),
  ) +
  scale_alpha_continuous(
    guide = "none",
    range = c(0.3, 0.7)
  ) +
  scale_linewidth_continuous(
    guide = "none",
    range = c(0.1, 5)
  ) +
  facet_grid(
    rows = vars(variable),
    labeller = label_parsed
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title = "Do the approximations capture the true parameters well enough?"
  )

ggsave(
  here("benchmarks", "gev_margins", "figures", "estimate_comparison_linear.png"),
  scale = 1.4, width = 8, height = 0.4 * 8
)
