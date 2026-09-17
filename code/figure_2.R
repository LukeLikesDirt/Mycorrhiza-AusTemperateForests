# ─────────────────────────────────────────────────────────────────────────────
# Figure 2a — mineral nitrogen responses of the three mycorrhizal guilds
#
#   G-AMF  Glomeromycota
#   M-AMF  Densosporales
#   EMF    Ectomycorrhizal fungi
#
# 3 x 3 grid: guild in columns, diversity metric in rows (Hill–Shannon
# taxonomic diversity, Hill–Shannon phylogenetic diversity, SRS abundance).
#
# Input is built by code/01e_prepare_generated_data.R.
# ─────────────────────────────────────────────────────────────────────────────

library(patchwork)
library(ggtext)
library(tidyverse)

load("generated_data/figure_2_panel.Rdata")

# ─────────────────────────────────────────────────────────────────────────────
# Style
# ─────────────────────────────────────────────────────────────────────────────

# Blues carried over from Figure S3 (code/01a_alpha_diversity_amf.R): the same
# fit/ribbon pair is used in every panel, so guild identity is read from the
# column title rather than from colour.
fit_colour    <- "#3366FF"
ribbon_colour <- "#3181FF"

theme_custom <- theme_minimal() +
  theme(
    panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    axis.ticks = element_line(colour = "grey80", linewidth = 0.25),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.title = element_markdown(size = 10),
    axis.text = element_markdown(size = 8),
    plot.title = element_markdown(size = 12, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(t = 2, r = 2, b = 0, l = 2),
    aspect.ratio = 1
  )

# Common x range across all nine panels
x_limits_nitro <- c(-2.12, 2.35)

# Axis titles
axis_title_labels <- c(
  "Hill–Shannon"    = "Taxanomic diversity",
  "Hill–Shannon PD" = "Phylogenetic diversity",
  "Abundance"       = "Abundance"
)

# Set guild levels
guild_levels <- c("EMF", "G-AMF", "M-AMF")

# ─────────────────────────────────────────────────────────────────────────────
# Panels
# ─────────────────────────────────────────────────────────────────────────────

# One panel of the 3 x 3 grid. Column titles appear on the top row only, y-axis
# titles on the first column only, and x-axis text on the bottom row only, so
# the grid reads as a single block rather than nine separate plots.
make_panel <- function(guild_name, response_name, show_title, show_y_title,
                       show_x_text, show_x_title) {

  pts   <- panel_points %>% filter(guild == guild_name, response == response_name)
  prd   <- panel_pred   %>% filter(guild == guild_name, response == response_name)
  annot <- panel_annot  %>% filter(guild == guild_name, response == response_name)

  p <- ggplot(pts, aes(x = x, y = y)) +
    geom_point(alpha = 0.2, size = 2, colour = "#4A4A4A", shape = 16) +
    geom_ribbon(
      data = prd,
      aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
      alpha = 0.2, fill = ribbon_colour, inherit.aes = FALSE
    ) +
    geom_line(
      data = prd,
      aes(x = x, y = predicted),
      colour = fit_colour, linewidth = 1.2, inherit.aes = FALSE
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    coord_cartesian(ylim = c(0, NA)) +
    scale_x_continuous(limits = x_limits_nitro, breaks = scales::pretty_breaks(n = 5)) +
    annotate(
      "richtext", x = Inf, y = Inf, label = annot$label,
      hjust = 1, vjust = 1, size = 2.4, fill = NA, label.color = NA
    ) +
    labs(
      x     = if (show_x_title) "Mineral nitrogen" else NULL,
      y     = if (show_y_title) axis_title_labels[[as.character(response_name)]] else NULL,
      title = if (show_title) paste0("**", guild_name, "**") else NULL
    ) +
    theme_custom

  if (!show_x_text) p <- p + theme(axis.text.x = element_blank())

  p
}

# The x-axis title hangs off the bottom-centre panel (M-AMF / Abundance), which
# centres it under the grid without needing a separate spacer plot.
middle_column <- ceiling(length(guild_levels) / 2)

panels <- list()
for (i in seq_along(response_levels)) {
  for (j in seq_along(guild_levels)) {
    panels[[length(panels) + 1]] <- make_panel(
      guild_name    = guild_levels[j],
      response_name = response_levels[i],
      show_title    = i == 1,
      show_y_title  = j == 1,
      show_x_text   = i == length(response_levels),
      show_x_title  = i == length(response_levels) && j == middle_column
    )
  }
}

figure_2 <- wrap_plots(panels, ncol = length(guild_levels))

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────

ggsave(
  "output/figure_2.png",
  figure_2,
  width = 18, height = 18, units = "cm",
  dpi = 300
)

ggsave(
  "output/figure_2.tiff",
  figure_2,
  width = 18, height = 18, units = "cm",
  dpi = 600
)

cat("Saved output/figure_2.png and output/figure_2.tiff\n")
