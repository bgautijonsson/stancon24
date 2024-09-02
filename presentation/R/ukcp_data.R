library(bggjphd)
library(tidyverse)
library(tmap)
d <- stations |>
  stations_to_sf() |>
  points_to_grid()

uk <- get_uk_spatial()


tm_shape(uk) +
  tm_polygons("grey") +
  tm_shape(d) +
  tm_polygons("station", alpha = 0.1, border.col = NULL) +
  tm_legend(show = FALSE)

ggplot(uk) +
  geom_sf() +
  geom_sf(
    data = d,
    col = NA,
    alpha = 0.2,
    fill = "#3182bd"
  ) +
  theme_bggj()

ggsave(
  filename = here::here("presentation", "images", "ukcp_data.png"),
  scale = 1.1, width = 8, height = 1.2 * 8
)
