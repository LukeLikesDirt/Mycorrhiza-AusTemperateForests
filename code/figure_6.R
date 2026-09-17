library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggtext)
library(tidyverse)

source("code/functions_taxon_names.R")   # format_taxon_name()

# ─────────────────────────────────────────────────────────────────────────────
# Figure 6 — Phylogenetic clustering of nitrogen-response indicators (plotting)
#
# Plotting only; all computation is in code/03_phylogenetic_analysis.R, which
# saves generated_data/figure_6.Rdata. Run that first.
#
# One row per guild showing the pre-specified metric for each marker (B1 in
# 03_phylogenetic_analysis.R): SES-MPD for G-AMF and M-AMF (SSU V4 resolves deep
# structure), SES-MNTD for EMF (ITS resolves terminal structure), computed on
# the phylogenetic basis only -- this is exactly `primary_raw`, saved by
# 03_phylogenetic_analysis.R, so the guild x metric selection lives in one place.
# Each row has two points: 2 directions (red = z+, blue = z-). Point = SES; wide
# band = jackknife range (taxon influence); thin bar = approximate bootstrap 95%
# CI from TITAN2's site-resampling replicates (membership uncertainty).
# Convention: negative = clustered, positive = overdispersed; significance is
# the ses rank p-value / empirical-ROPE classification (table 6), so no
# +/-1.96 reference lines are drawn.
#
# CAVEAT (state in caption/Methods): the bootstrap reconstructs each replicate's
# indicators by direction + IndVal p<=0.05, which is more permissive than the
# reported purity/reliability filter, so replicate sets run larger than the
# reported sets and the CI characterises a related, more-permissive estimand
# rather than a sampling CI for the reported SES. See 03_phylogenetic_analysis.R.
#
#   1. Figure 6  — SES per guild (MPD for G-AMF/M-AMF, MNTD for EMF), legend at right
#   2. Figure S  — phylogenetic trees with change-point bar plots
# ─────────────────────────────────────────────────────────────────────────────

load("generated_data/figure_6.Rdata")   # results, primary_raw, freedom, guild_levels, GUILD_LABELS, dat_*, tree_*

# Response-direction colours, consistent with figures 4 and 5
ind_cols <- c("Nitrophobic (z-)" = "#2166ac",
              "Nitrophilic (z+)" = "#d73027",
              "Not retained"     = "grey70")
dir_cols <- c("Nitrophilic (z+)" = "#d73027", "Nitrophobic (z-)" = "#2166ac")

# z-/z+ legend labels formatted as in figure_4.R
dir_labels <- c("Nitrophilic (z+)" = "*z*<sup> +</sup>", "Nitrophobic (z-)" = "*z*<sup> −</sup>")

# Plot constants
tag_size   <- 14
title_size <- 10
text_size  <- 9

# Plot theme
theme_custom <- theme_minimal() +
  theme(
    panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    axis.ticks = element_line(colour = "grey80", linewidth = 0.25),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.title = element_markdown(size = title_size),
    axis.text.y = element_markdown(size = title_size),
    axis.text.x.bottom = element_markdown(size = text_size),
    axis.text.x.top = element_text(size = 6),
    plot.title = element_markdown(size = title_size, hjust = 0.5),
    plot.tag = element_markdown(size = tag_size),
    legend.position = "none",
    plot.margin = margin(t = 2, r = 2, b = 0, l = 2),
    aspect.ratio = 1
  )


# ─────────────────────────────────────────────────────────────────────────────
# 1. Figure 6 — SES per guild (MPD for G-AMF/M-AMF, MNTD for EMF)
# ─────────────────────────────────────────────────────────────────────────────
# One row per guild, using `primary_raw` saved by 03_phylogenetic_analysis.R:
# SES-MPD for G-AMF and M-AMF, SES-MNTD for EMF, phylogenetic basis only (B1's
# PRIMARY_BASIS / PRIMARY_METRIC). Two points per row: 2 directions (red = z+,
# blue = z-), dodged vertically. Convention: negative = clustered, positive =
# overdispersed. No +/-1.96 reference lines -- significance is the ses rank
# p-value / empirical-ROPE classification, reported in table 6.
# Two nested bars per point:
#   wide translucent = jackknife range (taxon influence)
#   thin solid       = bootstrap 95% CI (ses_boot_lo/ses_boot_hi), drawn only
#                      OUTSIDE the jackknife band so the two never overlap
# Row order, top to bottom: EMF, G-AMF, M-AMF. Must be set before plot_df
# below, since gy is computed from this ordering -- reassigning guild_levels
# afterwards would leave the axis labels (drawn from this) out of sync with
# the point/bar y-positions (baked into plot_df from whatever ordering was in
# effect when gy was computed).
guild_levels <- c(GUILD_LABELS[["emf"]], GUILD_LABELS[["gamf"]], GUILD_LABELS[["mamf"]])

plot_df <- primary_raw %>%
  mutate(
    direction = factor(set, levels = c("Nitrophilic (z+)", "Nitrophobic (z-)")),
    gy = as.numeric(factor(guild, levels = rev(guild_levels))),
    y  = gy + ifelse(direction == "Nitrophilic (z+)", 0.18, -0.18),
    jack_lo = pmin(jack_min, ses, na.rm = TRUE),
    jack_hi = pmax(jack_max, ses, na.rm = TRUE)
  )

# Bootstrap-CI whiskers, drawn only where the CI extends beyond the jackknife
# band so the thin bar and the wide band never overlap. Left segment runs from
# ses_boot_lo to the band's lower edge, right from the band's upper edge to
# ses_boot_hi; either is dropped (x >= xend) when the CI sits inside the band.
whisk <- bind_rows(
  plot_df %>% transmute(y, direction, x = ses_boot_lo, xend = ses),
  plot_df %>% transmute(y, direction, x = ses, xend = ses_boot_hi)
) %>% filter(!is.na(x), !is.na(xend), x < xend)

x_lim <- c(floor(min(c(plot_df$ses_boot_lo, plot_df$jack_lo), na.rm = TRUE)),
           ceiling(max(c(plot_df$ses_boot_hi, plot_df$jack_hi), na.rm = TRUE)))

figure_6 <- ggplot(plot_df) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = 1.96, linetype = "dotted", colour = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = -1.96, linetype = "dotted", colour = "grey30", linewidth = 0.6) +
  geom_segment(data = whisk,
               aes(x = x, xend = xend, y = y, yend = y, colour = direction),
               linewidth = 0.7, na.rm = TRUE) +
  geom_errorbarh(aes(y = y, xmin = jack_lo, xmax = jack_hi, colour = direction),
                 height = 0, linewidth = 4, alpha = 0.25, na.rm = TRUE) +
  geom_point(aes(x = ses, y = y, colour = direction),
             shape = 16, size = 4, stroke = 1.1) +
  scale_colour_manual(values = dir_cols, labels = dir_labels, name = "Response<br>direction",
                      breaks = c("Nitrophilic (z+)", "Nitrophobic (z-)")) +
  scale_y_continuous(breaks = seq_along(guild_levels), labels = rev(guild_levels),
                     limits = c(0.5, length(guild_levels) + 0.5)) +
  scale_x_continuous(
    limits = x_lim,
    breaks = c(-1.96, 0, 1.96),
    # Single break at x = 0: the label is drawn at that tick's own position,
    # not centred across the full (asymmetric) axis like a title would be
    sec.axis = dup_axis(name = NULL, breaks = 0,
                        labels = "        Clustered ← Random → Overdispersed")
  ) +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 4))) +
  labs(x = "Standardised effect size (SES)", y = NULL) +
  theme_custom +
  theme(
    legend.position    = "right",
    legend.title       = element_markdown(size = title_size),
    legend.text        = element_markdown(size = text_size),
    axis.text.x.top    = element_text(size = text_size, colour = "grey30", hjust = 0.5, vjust = 1),
    axis.text.x.bottom = element_markdown(size = text_size),
    axis.ticks         = element_blank()
  )

ggsave("output/figure_6.png",  figure_6, width = 120, height = 100, units = "mm", dpi = 600, bg = "white")
ggsave("output/figure_6.tiff", figure_6, width = 120, height = 100, units = "mm", dpi = 600, bg = "white")
cat("Saved output/figure_6.png and output/figure_6.tiff\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Figure S — phylogenetic trees with change-point bar plots
# ─────────────────────────────────────────────────────────────────────────────
# One figure per guild: rectangular phylogram with tips coloured by response
# direction, and an aligned change-point bar (mg N/kg, ggtreeExtra::geom_fruit)
# for each taxon coloured by family (AMF) or FungalTraits lineage (EMF). Two
# legends -- response direction and family/lineage -- sit inside the panel at
# top-left. G-AMF and M-AMF are shown separately, matching the effect-size
# analysis.
dir.create("output/phylo_dispersion", showWarnings = FALSE, recursive = TRUE)

# Discrete, high-contrast palette for families/lineages (up to 36 categories)
fam_palette <- function(lvls) setNames(grDevices::palette.colors(length(lvls), "Polychrome 36"), lvls)

# Legend display names: AMF family / pseudo-family formatting; strip the leading
# "/" from FungalTraits lineage codes for EMF.
fam_display <- function(x, guild) {
  if (guild == GUILD_LABELS[["emf"]]) sub("^/", "", x) else format_taxon_name(x, italics = FALSE)
}

make_change_point_figure <- function(tr, dat, guild, tip_size = 1.6) {
  plot_data <- dat %>%
    filter(otu_id %in% tr$tip.label) %>%
    mutate(
      indicator_label   = factor(indicator, levels = c("neg", "pos", "none"),
                                 labels = names(ind_cols)),
      family_or_lineage = factor(family_or_lineage, levels = sort(unique(family_or_lineage)))
    )
  fam_cols  <- fam_palette(levels(plot_data$family_or_lineage))
  fam_labs  <- fam_display(levels(plot_data$family_or_lineage), guild)
  rank_name <- if (guild == GUILD_LABELS[["emf"]]) "Lineage" else "Family"

  # Tree tips coloured by response direction; change-point bars (geom_fruit)
  # coloured by family/lineage. Colour (direction) and fill (family) are separate
  # scales, so both legends appear -- placed inside the panel, top-left.
  # NB: attaching data with %<+% breaks geom_fruit's mapping, so tips are drawn
  # from the tree's own coordinates instead.
  p <- ggtree(tr, layout = "rectangular", ladderize = TRUE, linewidth = 0.3, colour = "grey55")
  tip_xy <- p$data %>%
    filter(isTip) %>%
    left_join(plot_data %>% select(otu_id, indicator_label), by = c("label" = "otu_id"))

  p +
    geom_point(data = tip_xy, aes(x = x, y = y, colour = indicator_label),
               inherit.aes = FALSE, size = tip_size) +
    scale_colour_manual(values = ind_cols, name = "Response direction",
                        guide = guide_legend(order = 1, override.aes = list(size = 2.4))) +
    geom_fruit(
      data    = plot_data,
      geom    = geom_bar,
      mapping = aes(y = otu_id, x = threshold_mgkg, fill = family_or_lineage),
      orientation = "y", stat = "identity", pwidth = 0.55, offset = 0.06, colour = NA,
      axis.params = list(axis = "x", text.size = 2, nbreak = 4)
    ) +
    scale_fill_manual(values = fam_cols, labels = fam_labs, name = rank_name,
                      guide = guide_legend(order = 2, ncol = 1)) +
    theme(
      legend.position        = "inside",
      legend.position.inside  = c(0.01, 0.99),
      legend.justification    = c(0, 1),
      legend.title            = element_text(size = 8, face = "bold"),
      legend.text             = element_text(size = 6.5),
      legend.key.size         = unit(0.30, "cm"),
      legend.spacing.y        = unit(1, "pt"),
      legend.background       = element_rect(fill = alpha("white", 0.55), colour = NA)
    )
}

fig_specs <- list(
  list(tr = tree_emf,  dat = dat_emf,  guild = GUILD_LABELS[["emf"]],  file = "figure_s7.png", h = 42, tip = 1.2),
  list(tr = tree_gamf, dat = dat_gamf, guild = GUILD_LABELS[["gamf"]], file = "figure_s8.png", h = 26, tip = 1.6),
  list(tr = tree_mamf, dat = dat_mamf, guild = GUILD_LABELS[["mamf"]], file = "figure_s9.png", h = 16, tip = 2.2)
)

for (s in fig_specs) {
  fig <- make_change_point_figure(s$tr, s$dat, s$guild, tip_size = s$tip)
  ggsave(file.path("output/", s$file), fig,
         width = 24, height = s$h, units = "cm", dpi = 300, limitsize = FALSE)
  cat("Saved output/", s$file, "\n", sep = "")
}

cat("\n=== FIGURE 6 COMPLETE ===\n")

