library(patchwork)
library(ggtext)
library(paletteer)
library(tidyverse)

source("code/functions_taxon_names.R")

# Load ECM TITAN results
load("generated_data/figure_5.Rdata")

# Set nitrogen axis limits from observed data range
nitrogen_limits <- data.table::fread("generated_data/sample_covariates.txt") %>%
  mutate(mineral_nitrogen = nitrate + ammonium) %>%
  pull(mineral_nitrogen) %>%
  range()

# ─────────────────────────────────────────────────────────────────────────────
# Shared theme
# ─────────────────────────────────────────────────────────────────────────────
theme_custom <- theme_minimal() +
  theme(
    panel.border      = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
    panel.grid        = element_blank(),
    axis.ticks        = element_line(colour = "grey80", linewidth = 0.25),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.title        = element_markdown(size = 10),
    axis.text         = element_markdown(size = 8),
    axis.text.y.right = element_markdown(size = 8, hjust = 0),
    plot.margin       = margin(5, 1, 1, 1),
    legend.title      = element_markdown(size = 10),
    legend.text       = element_markdown(size = 8),
    plot.tag          = element_markdown(size = 14, margin = margin(0, 0, 0, 0))
  )

# ─────────────────────────────────────────────────────────────────────────────
# Lineage colour palette (consistent across both panels)
# ─────────────────────────────────────────────────────────────────────────────
sorted_lineages <- sort(unique(c(titan_taxa_n_spp$lineage, titan_taxa_n_lin$lineage)))
lineage_colors  <- setNames(
  paletteer::paletteer_c("grDevices::Spectral", n = length(sorted_lineages)),
  sorted_lineages
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared size scale (comparable point sizes across panels and figures 4 & 5)
# ─────────────────────────────────────────────────────────────────────────────
# Load AM data into a separate environment to avoid variable name collisions
amf_env <- new.env()
load("generated_data/figure_4.Rdata", envir = amf_env)
z_min       <- min(c(titan_taxa_n_spp$z.median, titan_taxa_n_lin$z.median,
                      amf_env$titan_taxa_n_spp$z.median, amf_env$titan_taxa_n_fam$z.median), na.rm = TRUE)
z_max       <- max(abs(c(titan_taxa_n_spp$z.median, titan_taxa_n_lin$z.median,
                          amf_env$titan_taxa_n_spp$z.median, amf_env$titan_taxa_n_fam$z.median)), na.rm = TRUE)
rm(amf_env)
size_limits <- c(
  floor(z_min), 
  ceiling(max(abs(
    c(titan_taxa_n_spp$z.median, titan_taxa_n_lin$z.median)
  ), na.rm = TRUE))
)
size_breaks <- seq(floor(z_min), ceiling(z_max), by = 3)

# ─────────────────────────────────────────────────────────────────────────────
# -- Panel (a): Lineage-level data --------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
lin_data <- titan_taxa_n_lin %>%
  mutate(
    threshold_bt       = exp(threshold_log),
    threshold_lower_bt = exp(threshold_lower),
    threshold_upper_bt = exp(threshold_upper),
    lineage_color      = unname(lineage_colors[lineage]),
    # Axis labels with colored squares
    label_left  = paste0(lineage, " <span style='color:", lineage_color, "; font-size:9pt;'>&#9608;</span>"),
    label_right = paste0("<span style='color:", lineage_color, "; font-size:9pt;'>&#9608;</span> ", lineage)
  )

# Assign y-positions: Nitrophobic on integers, Nitrophilic on half-integers (interleaved)
lin_neg <- lin_data %>%
  filter(response == "Nitrophobic") %>%
  arrange(desc(threshold_bt)) %>%
  mutate(y_pos = row_number())

lin_pos <- lin_data %>%
  filter(response == "Nitrophilic") %>%
  arrange(threshold_bt) %>%
  mutate(y_pos = row_number() + 0.5)

lin_data_pos <- bind_rows(lin_neg, lin_pos) %>%
  mutate(response = factor(response, levels = c("Nitrophobic", "Nitrophilic")))

# Add dummy data for legend (z+ point, placed outside plot limits)
lin_data_pos <- lin_data_pos %>%
  add_row(
    threshold_bt       = 1,
    threshold_lower_bt = 1,
    threshold_upper_bt = 1,
    y_pos              = -999,  # Outside y-axis limits
    response           = factor("Nitrophilic", levels = c("Nitrophobic", "Nitrophilic")),
    z.median             = 10,
    lineage_color      = unname(lineage_colors[1])
  )

# Axis breaks and labels
lin_left_breaks  <- lin_neg %>% arrange(y_pos) %>% pull(y_pos)
lin_left_labels  <- lin_neg %>% arrange(y_pos) %>% pull(label_left)
lin_right_breaks <- lin_pos %>% arrange(y_pos) %>% pull(y_pos)
lin_right_labels <- lin_pos %>% arrange(y_pos) %>% pull(label_right)

# Community threshold ranges
lin_ranges <- lin_data_pos %>%
  filter(y_pos > 0) %>%  # Exclude dummy data
  group_by(response) %>%
  summarise(min_y = min(y_pos), max_y = max(y_pos), .groups = "drop")

comm_lin <- community_thresholds_n_lin %>%
  filter(grepl("fsumz", group)) %>%
  mutate(response = if_else(grepl("Nitrophobic", response), "Nitrophobic", "Nitrophilic")) %>%
  left_join(lin_ranges, by = "response")

# ─────────────────────────────────────────────────────────────────────────────
# -- Panel (b): Species-level data --------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
plot_data <- titan_taxa_n_spp %>%
  mutate(
    threshold_bt       = exp(threshold_log),
    threshold_lower_bt = exp(threshold_lower),
    threshold_upper_bt = exp(threshold_upper),
    species_display = format_taxon_name(species, italics = TRUE),
    lineage_color = unname(lineage_colors[lineage]),
    # Axis labels with colored squares
    species_label_left  = paste0(species_display, " <span style='color:", lineage_color, "; font-size:9pt;'>&#9608;</span>"),
    species_label_right = paste0("<span style='color:", lineage_color, "; font-size:9pt;'>&#9608;</span> ", species_display)
  )

# Assign y-positions (same interleaved logic as panel a)
nitrophobic_pos <- plot_data %>%
  filter(response == "Nitrophobic") %>%
  arrange(desc(threshold_bt)) %>%
  mutate(y_pos = row_number())

nitrophilic_pos <- plot_data %>%
  filter(response == "Nitrophilic") %>%
  arrange(threshold_bt) %>%
  mutate(y_pos = row_number() + 0.5)

plot_data_pos <- bind_rows(nitrophobic_pos, nitrophilic_pos)

# Axis breaks and labels
spp_left_breaks  <- nitrophobic_pos %>% arrange(y_pos) %>% pull(y_pos)
spp_left_labels  <- nitrophobic_pos %>% arrange(y_pos) %>% pull(species_label_left)
spp_right_breaks <- nitrophilic_pos %>% arrange(y_pos) %>% pull(y_pos)
spp_right_labels <- nitrophilic_pos %>% arrange(y_pos) %>% pull(species_label_right)

# Community threshold ranges
spp_ranges <- plot_data_pos %>%
  group_by(response) %>%
  summarise(min_y = min(y_pos), max_y = max(y_pos), .groups = "drop")

comm_spp <- community_thresholds_n_spp %>%
  filter(grepl("fsumz", group)) %>%
  mutate(response = if_else(grepl("Nitrophobic", response), "Nitrophobic", "Nitrophilic")) %>%
  left_join(spp_ranges, by = "response")

# Dummy data for lineage legend (invisible points to register all lineages)
lineage_legend_data <- data.frame(lineage = names(lineage_colors), x = 1, y = 1)

# ─────────────────────────────────────────────────────────────────────────────
# Panel (a): Lineage-level plot
# ─────────────────────────────────────────────────────────────────────────────
panel_a <- lin_data_pos %>%
  ggplot(aes(x = threshold_bt, y = y_pos, colour = response,
             xmin = threshold_lower_bt, xmax = threshold_upper_bt)) +
  # Community threshold bands
  geom_rect(
    data = filter(comm_lin, response == "Nitrophilic"),
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = min_y - 0.5, ymax = max_y + 0.5),
    fill = "#d73027", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_rect(
    data = filter(comm_lin, response == "Nitrophobic"),
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = min_y - 0.5, ymax = max_y + 0.5),
    fill = "#2166ac", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_segment(
    data = comm_lin,
    aes(x = threshold_original, xend = threshold_original,
        y = min_y - 0.5, yend = max_y + 0.5, colour = response),
    linetype = "dotted", linewidth = 0.8, alpha = 0.7, inherit.aes = FALSE
  ) +
  # Taxon thresholds
  geom_errorbarh(height = 0.3, linewidth = 0.6) +
  geom_point(aes(size = abs(z.median))) +
  # Scales
  # drop = FALSE retains both factor levels in the legend even when z+ is absent
  scale_colour_manual(
    values = c("Nitrophobic" = "#2166ac", "Nitrophilic" = "#d73027"),
    labels = c("Nitrophobic" = "*z*<sup> −</sup>", "Nitrophilic" = "*z*<sup> +</sup>"),
    drop   = FALSE,
    name   = "Response direction"
  ) +
  scale_size_continuous(name = "Response magnitude", breaks = size_breaks, limits = size_limits) +
  # No limits here — limits are set via coord_cartesian so the dummy z+ row
  # at y_pos = -999 is clipped visually but still registers with the colour scale
  scale_y_continuous(
    breaks   = lin_left_breaks,
    labels   = lin_left_labels,
    sec.axis = sec_axis(~ ., breaks = lin_right_breaks, labels = lin_right_labels)
  ) +
  coord_cartesian(ylim = c(0.5, max(filter(lin_data_pos, y_pos > 0)$y_pos) + 0.5)) +
  scale_x_log10(limits = nitrogen_limits) +
  labs(x = NULL, y = NULL, tag = "(**a**)") +
  guides(
    colour = guide_legend(override.aes = list(size = 4, alpha = 1), order = 1),
    size   = guide_legend(order = 2)
  ) +
  theme_minimal(base_size = 10) +
  theme_custom +
  theme(
    legend.position      = c(1.5, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = alpha("white", 0.8), colour = NA),
    legend.key.size      = unit(0.4, "cm"),
    axis.text.x          = element_blank(),
    plot.tag.location    = "panel",
    plot.tag.position    = c(0.93, 0.88)
  )

# ─────────────────────────────────────────────────────────────────────────────
# Panel (b): Species-level plot
# ─────────────────────────────────────────────────────────────────────────────
panel_b <- plot_data_pos %>%
  ggplot(aes(x = threshold_bt, y = y_pos, colour = response,
             xmin = threshold_lower_bt, xmax = threshold_upper_bt)) +
  # Dummy points for lineage legend
  geom_point(data = lineage_legend_data, aes(x = x, y = y, fill = lineage),
             size = 0, alpha = 0, inherit.aes = FALSE) +
  # Community threshold bands
  geom_rect(
    data = filter(comm_spp, response == "Nitrophilic"),
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = min_y - 0.5, ymax = max_y + 0.5),
    fill = "#d73027", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_rect(
    data = filter(comm_spp, response == "Nitrophobic"),
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = min_y - 0.5, ymax = max_y + 0.5),
    fill = "#2166ac", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_segment(
    data = comm_spp,
    aes(x = threshold_original, xend = threshold_original,
        y = min_y - 0.5, yend = max_y + 0.5, colour = response),
    linetype = "dotted", linewidth = 0.8, alpha = 0.7, inherit.aes = FALSE
  ) +
  # Taxon thresholds
  geom_errorbarh(height = 0.3, linewidth = 0.6) +
  geom_point(aes(size = abs(z.median))) +
  # Scales
  scale_colour_manual(
    values = c("Nitrophobic" = "#2166ac", "Nitrophilic" = "#d73027"),
    name   = "Direction"
  ) +
  scale_fill_manual(values = lineage_colors, name = "Lineage") +
  scale_size_continuous(name = "Magnitude", breaks = size_breaks, limits = size_limits) +
  scale_y_continuous(
    breaks   = spp_left_breaks,
    labels   = spp_left_labels,
    expand   = expansion(add = 0.5),
    sec.axis = sec_axis(~ ., breaks = spp_right_breaks, labels = spp_right_labels)
  ) +
  scale_x_log10(limits = nitrogen_limits) +
  labs(
    x = "Log-scaled mineral nitrogen threshold (mg kg<sup>-1</sup>)",
    y = NULL, 
    tag = "(**b**)"
  ) +
  guides(
    colour = "none",
    size   = "none",
    fill   = guide_legend(ncol = 4, override.aes = list(size = 4.5, alpha = 1, shape = 22, stroke = 0.15))
  ) +
  theme_minimal(base_size = 10) +
  theme_custom +
  theme(
    legend.position      = "bottom",
    legend.title         = element_blank(),
    legend.key.size      = unit(10, "pt"),
    legend.key.spacing.y = unit(0.1, "pt"),
    plot.tag.location    = "panel",
    plot.tag.position    = c(0.93, 0.96)
  )

# ─────────────────────────────────────────────────────────────────────────────
# -- Combine and save ---------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
figure <- wrap_plots(
  panel_a, panel_b, 
  ncol = 1,
  heights = c(nrow(lin_data) + 2, nrow(plot_data))
)

print(figure)

ggsave(
  "output/figure_5.png",
  figure,
  width  = 173,
  height = 210,
  units  = "mm",
  dpi    = 600
)
ggsave(
  "output/figure_5.tiff",
  figure,
  width  = 173,
  height = 210,
  units  = "mm",
  dpi    = 600
)

