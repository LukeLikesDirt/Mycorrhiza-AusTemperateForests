
require(data.table)
require(ggtext)
require(patchwork)
require(tidyverse)

# ─────────────────────────────────────────────────────────────────────────────
# -- EMF ----------------------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
# Load taxonomy
taxa_emf <- fread("generated_data/emf/classification.txt") %>%
  left_join(
    fread("generated_data/emf/otu_table_srs_emf_genus.txt") %>%
      select(genus, lineage),
    by = "genus"
  )
# Load OTUs
otu_long_emf <- fread("generated_data/emf/otu_table_srs.txt") %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id",
    values_to = "abundance") %>%
  inner_join(
    taxa_emf %>%
      select(otu_id, lineage),
    by = "otu_id"
  )

# Build summary table
lineage_summary_emf <- list(
  abundance = otu_long_emf %>% summarise(abundance = sum(abundance), .by = lineage),
  richness = otu_long_emf %>% filter(abundance > 0) %>% summarise(abs_richness = n_distinct(otu_id), .by = lineage)
) %>%
  reduce(inner_join, by = "lineage") %>%
  mutate(
    rel_abundance = abundance / sum(abundance),
    rel_richness = abs_richness / sum(abs_richness)
  ) %>%
  select(lineage, n_otus = abs_richness, Abundance = rel_abundance, Richness = rel_richness) %>%
  arrange(desc(Richness)) %>%
  print(n = Inf) %>%
  # Mutate groups with relative abundance & richness < 2% to "Other EMF lineages"
  mutate(
    lineage = case_when(
      Abundance < 0.02 & Richness < 0.02 ~ "Other EMF lineages",
      TRUE ~ as.character(lineage)
    )) %>%
  group_by(lineage) %>%
  summarise(
    Abundance = sum(Abundance),
    Richness = sum(Richness),
    n_otus = sum(n_otus)
  ) %>%
  ungroup() %>%
  # Order by richness (ascending), with "Other EMF lineages" at the end. The "[n]"
  # richness suffix is appended last so "Other EMF lineages" shows the summed OTU
  # count of the lineages it absorbed, not an individual lineage's count.
  arrange(desc(lineage == "Other EMF lineages"), Richness) %>%
  mutate(
    lineage = paste0(lineage, " [", n_otus, "]"),
    lineage = factor(lineage, levels = lineage)
  ) %>%
  # Remove underscores from lineage names (fct_relabel preserves factor order)
  mutate(lineage = fct_relabel(lineage, str_replace_all, "_", " ")) %>%
  pivot_longer(
    cols = c(Richness, Abundance),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("Richness", "Abundance")))

# Create custom color palette with mid grey for "Other EMF lineages" -- this grey is
# the average of the light/dark "Other" greys used for the G-AMF/M-AMF split above
lineage_levels <- levels(lineage_summary_emf$lineage)
n_colors <- length(lineage_levels)

# TolRainbow11 palette from Bokeh (converted from hex codes)
# tol_rainbow11 <- rev(c("#882E72", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265",
#                        "#CAE0AB", "#F7F056", "#F4A736", "#E8601C", "#DC050C",
#                        "#72190E"))
spectral_colors <- paletteer::paletteer_c("grDevices::Spectral", 11)

# Extend if needed
if (n_colors > 11) {
  spectral_colors <- colorRampPalette(spectral_colors)(n_colors)
} else {
  spectral_colors <- spectral_colors[1:n_colors]
}

# Replace "Other EMF lineages" color with mid grey
color_palette <- spectral_colors
other_position <- which(startsWith(lineage_levels, "Other EMF lineages"))
if (length(other_position) > 0) {
  color_palette[other_position] <- "grey75"
}
names(color_palette) <- lineage_levels

# Create the plot
plot_emf_lineage <- ggplot(
  lineage_summary_emf,
  aes(x = metric, y = value, fill = lineage)
) +
  geom_col(position = "fill", width = 0.8) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(
    labels = scales::percent_format(),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y          = element_text(angle = 90, hjust = 0.5, size = 10),
    axis.ticks.y         = element_blank(),
    axis.line.y          = element_blank(),
    axis.title.x         = element_text(size = 10),
    axis.text.x          = element_text(size = 8.5),
    axis.ticks.x         = element_line(colour = "grey70"),
    axis.ticks.length    = unit(-0.1, "cm"),
    axis.line.x          = element_line(colour = "grey70"),
    legend.text          = element_text(size = 8.5),
    legend.key.size      = unit(12, "pt"),
    legend.key.spacing.y = unit(2.5, "pt"),
    legend.position      = "bottom",
    plot.tag             = element_markdown(size = 14),
    plot.tag.location    = "plot",
    plot.tag.position    = c(0.01, 1.01),
    plot.margin          = margin(t = 20, r = 10, b = 1, l = 1),
    aspect.ratio         = 1/3
  ) +
  guides(fill = guide_legend(ncol = 4, byrow = FALSE, reverse = TRUE)) +
  labs(
    x = NULL,
    y = "Relative proportion",
    fill = NULL,
    tag = "(**a**)"
  )

# Display the plot
print(plot_emf_lineage)

# ─────────────────────────────────────────────────────────────────────────────
# -- AMF ------------------------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

# Load taxonomy
taxa_amf <- fread("generated_data/amf/classification.txt")

# Load OTUs
otu_long_amf <- fread("generated_data/amf/otu_table_srs.txt") %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  inner_join(
    taxa_amf %>% select(otu_id, family, am_group),
    by = "otu_id"
  )

# Build summary table (proportions are relative to the whole AMF assemblage, i.e.
# Glomeromycota + Densosporales combined, so the two halves of the bar sum to 100%)
family_summary_amf_raw <- list(
  abundance = otu_long_amf %>% summarise(abundance = sum(abundance), .by = c(family, am_group)),
  richness = otu_long_amf %>% filter(abundance > 0) %>% summarise(abs_richness = n_distinct(otu_id), .by = c(family, am_group))
) %>%
  reduce(inner_join, by = c("family", "am_group")) %>%
  mutate(
    rel_abundance = abundance / sum(abundance),
    rel_richness = abs_richness / sum(abs_richness)
  ) %>%
  select(family, am_group, n_otus = abs_richness, Abundance = rel_abundance, Richness = rel_richness) %>%
  arrange(desc(Richness)) %>%
  print(n = Inf)

# Mutate families with relative abundance & richness < 2% to a group-specific "Other"
# bucket, then order each group ascending by richness with "Other" first in the row
# order -- for a single group this is exactly the EMF logic below; putting the two
# groups' "Other" row first/last respectively is what lands G-AMF on the left half of
# the bar, M-AMF on the right half, with each "Other" bridging out of its own half.
# The "[n]" richness suffix is appended last so the "Other" bucket shows the summed
# OTU count of the families it absorbed, not an individual family's count.
collapse_am_group <- function(data, other_label) {
  data %>%
    mutate(
      family = case_when(
        Abundance < 0.02 & Richness < 0.02 ~ other_label,
        TRUE ~ as.character(family)
      )) %>%
    group_by(family) %>%
    summarise(Abundance = sum(Abundance), Richness = sum(Richness), n_otus = sum(n_otus)) %>%
    ungroup() %>%
    arrange(desc(family == other_label), Richness) %>%
    mutate(family = paste0(family, " [", n_otus, "]"))
}

g_amf_summary <- family_summary_amf_raw %>%
  filter(am_group == "Glomeromycota") %>%
  collapse_am_group("Other G-AMF families")

m_amf_summary <- family_summary_amf_raw %>%
  filter(am_group == "Endogonomycetes") %>%
  collapse_am_group("Other M-AMF families")

# Sequential fill palettes: blue-green for G-AMF, yellow-red for M-AMF. The two "Other"
# greys (light for G-AMF, dark for M-AMF) are chosen so their average equals the mid
# grey used for the EMF plot's "Other EMF lineages" (#D3D3D3), tying the two panels
# together.
n_g_amf <- sum(!startsWith(g_amf_summary$family, "Other G-AMF families"))
n_m_amf <- sum(!startsWith(m_amf_summary$family, "Other M-AMF families"))

g_amf_summary$fill <- "grey60"
g_amf_summary$fill[!startsWith(g_amf_summary$family, "Other G-AMF families")] <-
  RColorBrewer::brewer.pal(n_g_amf, "GnBu")

m_amf_summary$fill <- "grey90"
m_amf_summary$fill[!startsWith(m_amf_summary$family, "Other M-AMF families")] <-
  rev(paletteer::paletteer_c("grDevices::YlOrRd", n_m_amf))

# Combine with M-AMF rows first so, after the row-order reversal that geom_col's
# stacking applies, G-AMF ends up as the first half of the bar and M-AMF the second
family_summary_amf <- bind_rows(m_amf_summary, g_amf_summary) %>%
  mutate(family = factor(family, levels = family)) %>%
  pivot_longer(
    cols = c(Richness, Abundance),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("Richness", "Abundance"))) %>%
  # Remove underscores from family names (fct_relabel preserves factor order)
  mutate(family = fct_relabel(family, str_replace_all, "_", " "))

# Named colour vector matching the final (relabelled) factor levels
amf_fill_colors <- family_summary_amf %>%
  distinct(family, fill) %>%
  tibble::deframe()

# Create the plot
plot_amf_family <- ggplot(
  family_summary_amf,
  aes(x = metric, y = value, fill = family)
) +
  geom_col(position = "fill", width = 0.8) +
  scale_fill_manual(values = amf_fill_colors) +
  scale_y_continuous(
    labels = scales::percent_format(),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y          = element_text(angle = 90, hjust = 0.5, size = 10),
    axis.ticks.y         = element_blank(),
    axis.line.y          = element_blank(),
    axis.title.x         = element_text(size = 10),
    axis.text.x          = element_text(size = 8.5),
    axis.ticks.x         = element_line(colour = "grey70"),
    axis.ticks.length    = unit(-0.1, "cm"),
    axis.line.x          = element_line(colour = "grey70"),
    legend.text          = element_text(size = 8.5),
    legend.key.size      = unit(12, "pt"),
    legend.key.spacing.y = unit(2.5, "pt"),
    legend.position      = "bottom",
    plot.tag             = element_markdown(size = 14),
    plot.tag.location    = "plot",
    plot.tag.position    = c(0.01, 1.01),
    plot.margin          = margin(t = 10, r = 10, b = 1, l = 1),
    aspect.ratio         = 1/3
  ) +
  guides(fill = guide_legend(ncol = 4, byrow = FALSE, reverse = TRUE)) +
  labs(
    x = NULL,
    y = "Relative proportion",
    fill = NULL,
    tag = "(**b**)"
  )

# Display the plot
print(plot_amf_family)

# ─────────────────────────────────────────────────────────────────────────────
# -- JOIN & SAVE --------------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

plot_combined <- wrap_plots(
  plot_emf_lineage,
  plot_amf_family,
  ncol = 1
)

ggsave(
  filename = "output/figure_3.png",
  plot = plot_combined,
  width = 173,
  height = 190,
  units = "mm",
  dpi = 600
)
ggsave(
  filename = "output/figure_3.tiff",
  plot = plot_combined,
  width = 173,
  height = 190,
  units = "mm",
  dpi = 600
)

