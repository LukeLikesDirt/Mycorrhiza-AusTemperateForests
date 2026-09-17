library(patchwork)
library(ggtext)
library(paletteer)
library(tidyverse)

source("code/functions_taxon_names.R")

# Load AM TITAN results
load("generated_data/figure_4.Rdata")

# Set nitrogen axis limits from observed data range
nitrogen_limits <- data.table::fread("data/sample_covariates.txt") %>%
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
    axis.text         = element_markdown(size = 7.5),
    axis.text.y.right = element_markdown(size = 7.5, hjust = 0),
    plot.margin       = margin(5, 1, 1, 1),
    legend.title      = element_markdown(size = 10),
    legend.text       = element_markdown(size = 8),
    plot.tag          = element_markdown(size = 14, margin = margin(0, 0, 0, 0))
  )

# ─────────────────────────────────────────────────────────────────────────────
# Family colour palette and display labels (consistent across both panels)
# Glomeromycota families get a blue-green (GnBu) palette, Densosporales families
# get a yellow-red (YlOrRd) palette -- matching the AM colour scheme in figure_3.R
#
# Legend order: Glomeromycota families (alphabetical), then Densosporales families
# (alphabetical). Colours run dark -> light within each block, in that same order.
# ─────────────────────────────────────────────────────────────────────────────
sorted_families_all <- sort(unique(c(titan_taxa_n_spp$family, titan_taxa_n_fam$family)))

family_lineage <- data.table::fread("data/amf/classification.txt") %>%
  mutate(lineage = if_else(am_group == "Endogonomycetes", "Densosporales", "Glomeromycota")) %>%
  distinct(family, lineage) %>%
  deframe()

gam_families <- sorted_families_all[family_lineage[sorted_families_all] == "Glomeromycota"]
dam_families <- sorted_families_all[family_lineage[sorted_families_all] == "Densosporales"]
sorted_families <- c(gam_families, dam_families)

family_colors <- c(
  setNames(rev(RColorBrewer::brewer.pal(max(3, length(gam_families)), "GnBu"))[seq_along(gam_families)], gam_families),
  setNames(paletteer::paletteer_c("grDevices::YlOrRd", length(dam_families)), dam_families)
)[sorted_families]

family_labels <- setNames(format_taxon_name(sorted_families, italics = FALSE), sorted_families)

# ─────────────────────────────────────────────────────────────────────────────
# Shared size scale (comparable point sizes across panels and figures 4 & 5)
# ─────────────────────────────────────────────────────────────────────────────
# Load ECM data into a separate environment to avoid variable name collisions
emf_env <- new.env()
load("generated_data/figure_5.Rdata", envir = emf_env)
z_min       <- min(c(titan_taxa_n_spp$z.median, titan_taxa_n_fam$z.median,
                      emf_env$titan_taxa_n_spp$z.median, emf_env$titan_taxa_n_lin$z.median), na.rm = TRUE)
z_max       <- max(abs(c(titan_taxa_n_spp$z.median, titan_taxa_n_fam$z.median,
                          emf_env$titan_taxa_n_spp$z.median, emf_env$titan_taxa_n_lin$z.median)), na.rm = TRUE)
rm(emf_env)
size_limits <- c(
  floor(z_min), 
  ceiling(max(abs(
    c(titan_taxa_n_spp$z.median, titan_taxa_n_fam$z.median)
  ), na.rm = TRUE))
)
size_breaks <- seq(5, ceiling(z_max), by = 5)

# ─────────────────────────────────────────────────────────────────────────────
# Panel (a): Family-level data
# ─────────────────────────────────────────────────────────────────────────────
fam_data <- titan_taxa_n_fam %>%
  mutate(
    threshold_bt       = exp(threshold_log),
    threshold_lower_bt = exp(threshold_lower),
    threshold_upper_bt = exp(threshold_upper),
    family_display = unname(family_labels[family]),
    family_color   = unname(family_colors[family]),
    # Axis labels with coloured squares
    label_left  = paste0(family_display, " <span style='color:", family_color, "; font-size:9pt;'>&#9608;</span>"),
    label_right = paste0("<span style='color:", family_color, "; font-size:9pt;'>&#9608;</span> ", family_display)
  )

# Assign y-positions: Nitrophobic on integers, Nitrophilic on half-integers (interleaved)
fam_neg <- fam_data %>%
  filter(response == "Nitrophobic") %>%
  arrange(desc(threshold_bt)) %>%
  mutate(y_pos = row_number())

fam_pos <- fam_data %>%
  filter(response == "Nitrophilic") %>%
  arrange(threshold_bt) %>%
  mutate(y_pos = row_number() + 0.5)

fam_data_pos <- bind_rows(fam_neg, fam_pos) %>%
  mutate(response = factor(response, levels = c("Nitrophobic", "Nitrophilic")))

# Add dummy data for complete legend (placed outside y-axis limits)
# Ensures both z- and z+ appear in legend regardless of which response types are present
fam_data_pos <- fam_data_pos %>%
  add_row(
    threshold_bt       = 1,
    threshold_lower_bt = 1,
    threshold_upper_bt = 1,
    y_pos              = -999,
    response           = factor("Nitrophobic", levels = c("Nitrophobic", "Nitrophilic")),
    z.median             = 10,
    family_color       = unname(family_colors[1])
  ) %>%
  add_row(
    threshold_bt       = 1,
    threshold_lower_bt = 1,
    threshold_upper_bt = 1,
    y_pos              = -999,
    response           = factor("Nitrophilic", levels = c("Nitrophobic", "Nitrophilic")),
    z.median             = 10,
    family_color       = unname(family_colors[1])
  )

# Axis breaks and labels
fam_left_breaks  <- fam_neg %>% arrange(y_pos) %>% pull(y_pos)
fam_left_labels  <- fam_neg %>% arrange(y_pos) %>% pull(label_left)
fam_right_breaks <- fam_pos %>% arrange(y_pos) %>% pull(y_pos)
fam_right_labels <- fam_pos %>% arrange(y_pos) %>% pull(label_right)

# Community threshold y-ranges (dummy rows excluded)
fam_ranges <- fam_data_pos %>%
  filter(y_pos > 0) %>%
  group_by(response) %>%
  summarise(min_y = min(y_pos), max_y = max(y_pos), .groups = "drop") %>%
  # Add missing response type with default range if needed
  complete(
    response = factor(c("Nitrophobic", "Nitrophilic"), 
                      levels = c("Nitrophobic", "Nitrophilic")),
    fill = list(min_y = 0.5, max_y = 0.5)
  )

comm_fam <- community_thresholds_n_fam %>%
  filter(grepl("fsumz", group)) %>%
  mutate(response = if_else(grepl("Nitrophobic", response), "Nitrophobic", "Nitrophilic")) %>%
  left_join(fam_ranges, by = "response")

# ─────────────────────────────────────────────────────────────────────────────
# Panel (b): Species-level data
# ─────────────────────────────────────────────────────────────────────────────
plot_data <- titan_taxa_n_spp %>%
  mutate(
    threshold_bt       = exp(threshold_log),
    threshold_lower_bt = exp(threshold_lower),
    threshold_upper_bt = exp(threshold_upper),
    species_display = format_taxon_name(species, italics = TRUE),
    family_display  = unname(family_labels[family]),
    family_color    = unname(family_colors[family]),
    # Axis labels with coloured squares
    species_label_left  = paste0(species_display, " <span style='color:", family_color, "; font-size:9pt;'>&#9608;</span>"),
    species_label_right = paste0("<span style='color:", family_color, "; font-size:9pt;'>&#9608;</span> ", species_display)
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

# Community threshold y-ranges
spp_ranges <- plot_data_pos %>%
  group_by(response) %>%
  summarise(min_y = min(y_pos), max_y = max(y_pos), .groups = "drop")

comm_spp <- community_thresholds_n_spp %>%
  filter(grepl("fsumz", group)) %>%
  mutate(response = if_else(grepl("Nitrophobic", response), "Nitrophobic", "Nitrophilic")) %>%
  left_join(spp_ranges, by = "response")

# Dummy data for family legend (invisible points to register all families)
family_legend_data <- data.frame(
  family         = factor(sorted_families, levels = sorted_families),
  family_display = unname(family_labels[sorted_families]),
  x = 1,
  y = 1
)

# ─────────────────────────────────────────────────────────────────────────────
# Panel (a): Family-level plot
# ─────────────────────────────────────────────────────────────────────────────
panel_a <- fam_data_pos %>%
  ggplot(aes(x = threshold_bt, y = y_pos, colour = response,
             xmin = threshold_lower_bt, xmax = threshold_upper_bt)) +
  # Community threshold bands
  geom_rect(
    data = filter(comm_fam, response == "Nitrophilic"),
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = min_y - 0.5, ymax = max_y + 0.5),
    fill = "#d73027", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_rect(
    data = filter(comm_fam, response == "Nitrophobic"),
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = min_y - 0.5, ymax = max_y + 0.5),
    fill = "#2166ac", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_segment(
    data = comm_fam,
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
    labels = c("Nitrophobic" = "*z*<sup> −</sup>", "Nitrophilic" = "*z*<sup> +</sup>"),
    name   = "Response<br>direction"
  ) +
  scale_size_continuous(name = "Response<br>magnitude", breaks = size_breaks, limits = size_limits) +
  scale_y_continuous(
    breaks   = fam_left_breaks,
    labels   = fam_left_labels,
    sec.axis = sec_axis(~ ., breaks = fam_right_breaks, labels = fam_right_labels)
  ) +
  # Clip from y = 1 to remove white space at bottom
  coord_cartesian(ylim = c(0.75, max(fam_data_pos$y_pos) + 0.25)) +
  scale_x_log10(limits = nitrogen_limits) +
  labs(x = NULL, y = NULL, tag = "(**a**)") +
  guides(
    colour = guide_legend(override.aes = list(size = 4, alpha = 1), order = 1),
    size   = guide_legend(order = 2)
  ) +
  theme_minimal(base_size = 10) +
  theme_custom +
  theme(
    legend.position      = c(-0.60, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_blank(),
    legend.key.size      = unit(0.4, "cm"),
    axis.text.x          = element_blank(),
    plot.tag.location    = "panel",
    plot.tag.position    = c(0.065, 0.875)
  )

# ─────────────────────────────────────────────────────────────────────────────
# Panel (b): Species-level plot
# ─────────────────────────────────────────────────────────────────────────────
panel_b <- plot_data_pos %>%
  mutate(family = factor(family, levels = sorted_families)) %>%
  ggplot(aes(x = threshold_bt, y = y_pos, colour = response,
             xmin = threshold_lower_bt, xmax = threshold_upper_bt)) +
  # Dummy points for family legend
  geom_point(data = family_legend_data, aes(x = x, y = y, fill = family),
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
  scale_fill_manual(
    values = family_colors,
    labels = family_labels,
    name   = "Family"
  ) +
  scale_size_continuous(name = "Magnitude", breaks = size_breaks, limits = size_limits) +
  scale_y_continuous(
    breaks   = spp_left_breaks,
    labels   = spp_left_labels,
    expand   = expansion(add = 0.5),
    sec.axis = sec_axis(~ ., breaks = spp_right_breaks, labels = spp_right_labels)
  ) +
  scale_x_log10(limits = nitrogen_limits) +
  coord_cartesian(ylim = c(0.75, max(plot_data_pos$y_pos) + 0.25)) +
  labs(
    x   = "Log-scaled mineral nitrogen threshold (mg kg<sup>-1</sup>)",
    y   = NULL,
    tag = "(**b**)"
  ) +
  guides(
    colour = "none",
    size   = "none",
    fill   = guide_legend(ncol = 4, byrow = FALSE, override.aes = list(size = 4.5, alpha = 1, shape = 22, stroke = 0.15))
  ) +
  theme_minimal(base_size = 10) +
  theme_custom +
  theme(
    legend.position      = "bottom",
    legend.title         = element_blank(),
    legend.key.size      = unit(10, "pt"),
    legend.key.spacing.y = unit(0.1, "pt"),
    plot.tag.location    = "panel",
    plot.tag.position    = c(0.07, 0.95)
  )

# ─────────────────────────────────────────────────────────────────────────────
# Combine and save
# ─────────────────────────────────────────────────────────────────────────────
figure <- wrap_plots(
  panel_a, panel_b,
  ncol    = 1,
  heights = c(nrow(fam_data) + 2, nrow(plot_data))
)

print(figure)

ggsave(
  "output/figure_4.png",
  figure,
  width  = 173,
  height = 180,
  units  = "mm",
  dpi    = 600
)

ggsave(
  "output/figure_4.tiff",
  figure,
  width  = 173,
  height = 180,
  units  = "mm",
  dpi    = 600
)

