
# Required packages
require(ggspatial)
require(patchwork)
require(terra)
require(tidyverse)

# Sample metadata
data <- data.table::fread("generated_data/sample_metadata.csv")

# ── READ MAP VECTORS ──────────────────────────────────────────────────────────

# Australia vector
aus_map <- rnaturalearth::ne_countries(
  scale = 'medium', type = 'map_units', returnclass = 'sf'
) %>%
  filter(name == 'Australia') %>%
  ggplot() +
  geom_sf(color = "grey80", fill = NA) +
  coord_sf(
    xlim = c(113, 153),
    ylim = c(-43.5, -10),
    expand = FALSE
  ) +
  scale_x_continuous(
    breaks = c(120, 130, 140, 150),
    labels = c("120°E", "130°E", "140°E", "150°E")
  ) +
  scale_y_continuous(
    breaks = seq(-40, -10, 10)
  ) +
  xlab(NULL) +
  ylab(NULL) +
  theme_minimal()

# Biome vector
temperate_forest_vect <- vect(
  "../georeferenced_covariates/Ecoregions2017/Ecoregions2017.shp"
  ) %>%
  tidyterra::filter(BIOME_NAME == "Temperate Broadleaf & Mixed Forests") %>%
  # Crop to Australia extent
  terra::crop(ext(110, 155, -45, -10))

# ── CREATE MAP ────────────────────────────────────────────────────────────────

# Plot study area
insert_map <- aus_map +
  tidyterra::geom_spatvector(
    data = temperate_forest_vect,
    fill = "grey70",
    color = "grey70"
  ) +
  theme_void()

# Plot points on temperate forest map
main_map <- ggplot() +
  tidyterra::geom_spatvector(
    data = temperate_forest_vect,
    fill = "grey70",
    color = "grey70"
  ) +
  geom_point(
    data = data,
    aes(x = longitude, y = latitude),
    color = "#b2182b",
    size = 1,
    shape = 15,
    alpha = 0.6
  ) +
  scale_x_continuous(
    breaks = c(140, 145, 150),
    labels = c("140°E","145°E", "150°E")
  ) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  # Add north arrow
  annotation_north_arrow(
    location = "bl",  # bottom left
    which_north = "true",
    style = north_arrow_fancy_orienteering(
      line_width = 0.75,
      line_col = "grey80",
      fill = c("white", "grey80"),
      text_col = "grey80",
      text_size = 8,
    ),
    height = unit(1.2, "cm"),
    width = unit(1.2, "cm"),
    pad_x = unit(0.05, "cm"),
    pad_y = unit(0.15, "cm")
  ) +
  # Add scale bar
  annotation_scale(
    location = "br",
    width_hint = 0.3,
    style = "ticks",
    line_col = "grey80",
    text_col = "grey80",
    text_cex = 0.7,
    pad_x = unit(0.15, "cm"),
    pad_y = unit(0.15, "cm")
  ) +
  labs(
    x = NULL,
    y = NULL
  )

# Combine with inset map in top left corner
final_map <- main_map +
  inset_element(
    insert_map,
    left = 0.01,   # Position from left edge (0-1 scale)
    bottom = 0.625,  # Position from bottom edge (0-1 scale)
    right = 0.65,   # Position from left edge (0-1 scale)
    top = 0.99     # Position from bottom edge (0-1 scale)
  )

# Check final map
print(final_map)

# ── SAVE MAP ──────────────────────────────────────────────────────────────────

ggsave(
  filename = "output/figure_1.png",
  plot = final_map,
  width = 8,
  height = 13,
  unit = "cm",
  bg = "white",
  dpi = 600
)
ggsave(
  filename = "output/figure_1.tiff",
  plot = final_map,
  width = 8,
  height = 13,
  unit = "cm",
  bg = "white",
  dpi = 600
)

