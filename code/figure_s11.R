library(patchwork)
library(ggtext)
library(paletteer)
library(tidyverse)

source("code/functions_taxon_names.R")

# Load AM TITAN results
load("generated_data/figure_s11.Rdata")

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
# Genus colour palette and display labels (consistent across both panels)
# Glomeromycota genera get a blue-green (GnBu) palette, Densosporales genera
# get a yellow-red (YlOrRd) palette -- matching the AM colour scheme in figure_3.R
#
# Legend order: Glomeromycota genera (alphabetical), then Densosporales genera
# (alphabetical). Colours run dark -> light within each block, in that same order.
#
# GnBu is used as a continuous ramp (not RColorBrewer::brewer.pal, capped at 9
# colours) because genus counts routinely exceed the family-level counts.
# ─────────────────────────────────────────────────────────────────────────────
sorted_genera_all <- sort(unique(c(titan_taxa_n_spp$genus, titan_taxa_n_gen$genus)))

genus_lineage <- data.table::fread("data/amf/classification.txt") %>%
  mutate(lineage = if_else(am_group == "Endogonomycetes", "Densosporales", "Glomeromycota")) %>%
  distinct(genus, lineage) %>%
  deframe()

gam_genera <- sorted_genera_all[genus_lineage[sorted_genera_all] == "Glomeromycota"]
dam_genera <- sorted_genera_all[genus_lineage[sorted_genera_all] == "Densosporales"]
sorted_genera <- c(gam_genera, dam_genera)

genus_colors <- c(
  setNames(paletteer::paletteer_c("grDevices::GnBu", length(gam_genera)), gam_genera),
  setNames(paletteer::paletteer_c("grDevices::YlOrRd", length(dam_genera)), dam_genera)
)[sorted_genera]

genus_labels <- setNames(format_taxon_name(sorted_genera, italics = FALSE), sorted_genera)

# ─────────────────────────────────────────────────────────────────────────────
# Shared size scale (comparable point sizes across panels and figures 4 & 5)
# ─────────────────────────────────────────────────────────────────────────────
# Load ECM data into a separate environment to avoid variable name collisions
emf_env <- new.env()
load("generated_data/figure_5.Rdata", envir = emf_env)
z_min       <- min(c(titan_taxa_n_spp$z.median, titan_taxa_n_gen$z.median,
                      emf_env$titan_taxa_n_spp$z.median, emf_env$titan_taxa_n_lin$z.median), na.rm = TRUE)
z_max       <- max(abs(c(titan_taxa_n_spp$z.median, titan_taxa_n_gen$z.median,
                          emf_env$titan_taxa_n_spp$z.median, emf_env$titan_taxa_n_lin$z.median)), na.rm = TRUE)
rm(emf_env)
size_limits <- c(
  floor(z_min),
  ceiling(max(abs(
    c(titan_taxa_n_spp$z.median, titan_taxa_n_gen$z.median)
  ), na.rm = TRUE))
)
size_breaks <- seq(5, ceiling(z_max), by = 5)

# ─────────────────────────────────────────────────────────────────────────────
# Panel (a): Genus-level data
# ─────────────────────────────────────────────────────────────────────────────
gen_data <- titan_taxa_n_gen %>%
  mutate(
    threshold_bt       = exp(threshold_log),
    threshold_lower_bt = exp(threshold_lower),
    threshold_upper_bt = exp(threshold_upper),
    genus_display = unname(genus_labels[genus]),
    genus_color   = unname(genus_colors[genus]),
    # Axis labels with coloured squares
    label_left  = paste0(genus_display, " <span style='color:", genus_color, "; font-size:9pt;'>&#9608;</span>"),
    label_right = paste0("<span style='color:", genus_color, "; font-size:9pt;'>&#9608;</span> ", genus_display)
  )

# Assign y-positions: Nitrophobic on integers, Nitrophilic on half-integers (interleaved)
gen_neg <- gen_data %>%
  filter(response == "Nitrophobic") %>%
  arrange(desc(threshold_bt)) %>%
  mutate(y_pos = row_number())

gen_pos <- gen_data %>%
  filter(response == "Nitrophilic") %>%
  arrange(threshold_bt) %>%
  mutate(y_pos = row_number() + 0.5)

gen_data_pos <- bind_rows(gen_neg, gen_pos) %>%
  mutate(response = factor(response, levels = c("Nitrophobic", "Nitrophilic")))

# Add dummy data for complete legend (placed outside y-axis limits)
# Ensures both z- and z+ appear in legend regardless of which response types are present
gen_data_pos <- gen_data_pos %>%
  add_row(
    threshold_bt       = 1,
    threshold_lower_bt = 1,
    threshold_upper_bt = 1,
    y_pos              = -999,
    response           = factor("Nitrophobic", levels = c("Nitrophobic", "Nitrophilic")),
    z.median             = 10,
    genus_color        = unname(genus_colors[1])
  ) %>%
  add_row(
    threshold_bt       = 1,
    threshold_lower_bt = 1,
    threshold_upper_bt = 1,
    y_pos              = -999,
    response           = factor("Nitrophilic", levels = c("Nitrophobic", "Nitrophilic")),
    z.median             = 10,
    genus_color        = unname(genus_colors[1])
  )

# Axis breaks and labels
gen_left_breaks  <- gen_neg %>% arrange(y_pos) %>% pull(y_pos)
gen_left_labels  <- gen_neg %>% arrange(y_pos) %>% pull(label_left)
gen_right_breaks <- gen_pos %>% arrange(y_pos) %>% pull(y_pos)
gen_right_labels <- gen_pos %>% arrange(y_pos) %>% pull(label_right)

# Community threshold y-ranges (dummy rows excluded)
gen_ranges <- gen_data_pos %>%
  filter(y_pos > 0) %>%
  group_by(response) %>%
  summarise(min_y = min(y_pos), max_y = max(y_pos), .groups = "drop") %>%
  # Add missing response type with default range if needed
  complete(
    response = factor(c("Nitrophobic", "Nitrophilic"),
                      levels = c("Nitrophobic", "Nitrophilic")),
    fill = list(min_y = 0.5, max_y = 0.5)
  )

comm_gen <- community_thresholds_n_gen %>%
  filter(grepl("fsumz", group)) %>%
  mutate(response = if_else(grepl("Nitrophobic", response), "Nitrophobic", "Nitrophilic")) %>%
  left_join(gen_ranges, by = "response")

# ─────────────────────────────────────────────────────────────────────────────
# Panel (b): Species-level data
# ─────────────────────────────────────────────────────────────────────────────
plot_data <- titan_taxa_n_spp %>%
  mutate(
    threshold_bt       = exp(threshold_log),
    threshold_lower_bt = exp(threshold_lower),
    threshold_upper_bt = exp(threshold_upper),
    species_display = format_taxon_name(species, italics = TRUE),
    genus_display   = unname(genus_labels[genus]),
    genus_color     = unname(genus_colors[genus]),
    # Axis labels with coloured squares
    species_label_left  = paste0(species_display, " <span style='color:", genus_color, "; font-size:9pt;'>&#9608;</span>"),
    species_label_right = paste0("<span style='color:", genus_color, "; font-size:9pt;'>&#9608;</span> ", species_display)
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

# Dummy data for genus legend (invisible points to register all genera)
genus_legend_data <- data.frame(
  genus         = factor(sorted_genera, levels = sorted_genera),
  genus_display = unname(genus_labels[sorted_genera]),
  x = 1,
  y = 1
)

# ─────────────────────────────────────────────────────────────────────────────
# Panel (a): Genus-level plot
# ─────────────────────────────────────────────────────────────────────────────
panel_a <- gen_data_pos %>%
  ggplot(aes(x = threshold_bt, y = y_pos, colour = response,
             xmin = threshold_lower_bt, xmax = threshold_upper_bt)) +
  # Community threshold bands
  geom_rect(
    data = filter(comm_gen, response == "Nitrophilic"),
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = min_y - 0.5, ymax = max_y + 0.5),
    fill = "#d73027", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_rect(
    data = filter(comm_gen, response == "Nitrophobic"),
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = min_y - 0.5, ymax = max_y + 0.5),
    fill = "#2166ac", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_segment(
    data = comm_gen,
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
    breaks   = gen_left_breaks,
    labels   = gen_left_labels,
    sec.axis = sec_axis(~ ., breaks = gen_right_breaks, labels = gen_right_labels)
  ) +
  # Clip from y = 1 to remove white space at bottom
  coord_cartesian(ylim = c(0.75, max(gen_data_pos$y_pos) + 0.25)) +
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
  mutate(genus = factor(genus, levels = sorted_genera)) %>%
  ggplot(aes(x = threshold_bt, y = y_pos, colour = response,
             xmin = threshold_lower_bt, xmax = threshold_upper_bt)) +
  # Dummy points for genus legend
  geom_point(data = genus_legend_data, aes(x = x, y = y, fill = genus),
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
    values = genus_colors,
    labels = genus_labels,
    name   = "Genus"
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
  heights = c(nrow(gen_data), nrow(plot_data))
)

print(figure)

ggsave(
  "output/figure_s11.png",
  figure,
  width  = 173,
  height = 200,
  units  = "mm",
  dpi    = 600
)
