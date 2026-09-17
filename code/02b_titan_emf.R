require(data.table)
require(openxlsx)
require(TITAN2)
require(tidyverse)

# ─────────────────────────────────────────────────────────────────────────────
# TITAN2 — Threshold Indicator Taxa Analysis: EMF fungi, nitrogen gradient
#
# Identifies the mineral nitrogen value at which each OTU shows an abrupt
# change in abundance, classifying taxa as:
#   z+  Nitrophilic  — increase above threshold
#   z-  Nitrophobic  — decrease above threshold
#
# Two taxonomic levels are analysed:
#   Species level   — individual OTUs (faceted by lineage in plots)
#   Lineage level   — genera aggregated to EMF lineage; prevalence filtered
#                     after aggregation; lineage is the y-axis in plots
#
# Lineage assignments come from the SRS-normalised genus-level OTU table
# (data/emf/otu_table_srs_emf_genus.txt), which contains a `lineage` column.
#
# Gradient: log(nitrate + ammonium); thresholds back-transformed with exp()
# for reporting in mg/kg.
#
# Outputs saved to generated_data/figure_5.Rdata.
# ─────────────────────────────────────────────────────────────────────────────

# Guild naming — single lookup, do not scatter string literals (naming migration)
GUILD_LABELS <- c(gamf = "G-AMF", mamf = "M-AMF", emf = "EMF", amf = "AMF")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Sample metadata
# ─────────────────────────────────────────────────────────────────────────────
sample_metadata <- inner_join(
  fread("data/emf/diversity_emf.txt"),
  fread("data/sample_covariates.txt"),
  by = "sample_id"
)

# Log-transformed nitrogen gradient used as the TITAN environmental variable
titan_env_n <- sample_metadata %>%
  arrange(sample_id) %>%
  mutate(log_nitrogen = log(nitrate + ammonium)) %>%
  pull(log_nitrogen)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Taxonomy
# lineage is sourced from the genus-level OTU table (the species-level
# classification file does not include lineage directly)
# ─────────────────────────────────────────────────────────────────────────────
taxa_emf <- fread("data/emf/classification.txt") %>%
  left_join(
    fread("data/emf/otu_table_srs_emf_genus.txt") %>% select(genus, lineage),
    by = "genus"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 3. Species-level OTU table
#
# Prevalence filter: keep OTUs present in >= 5 samples, then build a
# samples × species matrix ordered to match titan_env_n.
# ─────────────────────────────────────────────────────────────────────────────
otu_long_emf <- fread("data/emf/otu_table_srs.txt") %>%
  pivot_longer(cols = -otu_id, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0)

n_sites <- n_distinct(otu_long_emf$sample_id)

otu_prev_spp <- otu_long_emf %>%
  group_by(otu_id) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop") %>%
  mutate(prevalence = n_samples / n_sites * 100) %>%
  filter(n_samples >= 5) %>%
  left_join(taxa_emf %>% select(otu_id, lineage, species), by = "otu_id") %>%
  arrange(desc(prevalence))

# Part C.1: two OTUs sharing a species string would silently produce
# list-columns in the pivot below, which TITAN would then misbehave on.
stopifnot(!any(duplicated(na.omit(otu_prev_spp$species))))

otu_table_spp <- otu_long_emf %>%
  filter(otu_id %in% otu_prev_spp$otu_id) %>%
  left_join(taxa_emf %>% select(otu_id, species), by = "otu_id") %>%
  select(-otu_id) %>%
  pivot_wider(names_from = species, values_from = abundance, values_fill = 0) %>%
  right_join(sample_metadata %>% select(sample_id), by = "sample_id") %>%
  mutate(across(-sample_id, ~ replace_na(., 0))) %>%
  arrange(sample_id) %>%
  column_to_rownames("sample_id")

stopifnot(length(titan_env_n) == nrow(otu_table_spp))
cat("Species TITAN — input OK:", length(titan_env_n), "sites,",
    ncol(otu_table_spp), "OTUs\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4. TITAN — species level, nitrogen gradient
# ─────────────────────────────────────────────────────────────────────────────
set.seed(1986)
titan_nitrogen_spp <- titan(
  env     = titan_env_n,
  txa     = otu_table_spp,
  numPerm = 1000,   # permutations for taxon-level significance
  boot    = TRUE,
  nBoot   = 1000,   # bootstrap replicates for CI on thresholds
  ivTot   = TRUE,   # include community-level sum(z) change points
  pur.cut = 0.95,   # purity: direction consistent in 95% of bootstraps
  rel.cut = 0.90,   # reliability: taxon significant in 90% of bootstraps
  ncpus   = 1       # single-core: TITAN2's parallel bootstrap is not seed-reproducible
)

# Persist species-level bootstrap replicates for the phylogenetic SES bootstrap
# (code/03_phylogenetic_analysis.R). metricArray is taxa x 4 x nBoot with no
# dimnames; slice 1 = maxgrp (1 = z-, 2 = z+), slice 4 = obsiv.prob (reliability
# thresholds this at < 0.05). Taxa order = rownames(sppmax).
if (!dir.exists("generated_data")) dir.create("generated_data")
saveRDS(list(metricArray = titan_nitrogen_spp$metricArray[, c(1, 4), , drop = FALSE],
             taxa        = rownames(titan_nitrogen_spp$sppmax),
             arguments   = titan_nitrogen_spp$arguments,
             slice_index = c(maxgrp = 1L, obsiv.prob = 4L)),
        "generated_data/titan_boot_emf.rds")

# Retain taxa meeting purity and reliability thresholds (Baker & King 2010)
titan_taxa_n_spp <- titan_nitrogen_spp$sppmax %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  filter(purity >= 0.95, reliability >= 0.9) %>%
  mutate(
    response           = case_when(
      filter == 1 ~ "Nitrophobic",   # z-: decline above threshold
      filter == 2 ~ "Nitrophilic",   # z+: increase above threshold
      TRUE        ~ "Neutral"
    ),
    threshold_original = exp(`50%`)  # back-transform log threshold to mg/kg
  ) %>%
  arrange(response, z.median) %>%
  select(
    species, response,
    threshold_log    = `50%`,   # median threshold (log scale)
    threshold_lower  = `5%`,    # 5th percentile bootstrap CI
    threshold_upper  = `95%`,   # 95th percentile bootstrap CI
    threshold_original,          # median threshold (mg/kg)
    zscore, z.median, reliability, purity
  ) %>%
  left_join(otu_prev_spp %>% select(species, lineage), by = "species")

cat("\nSignificant nitrogen TITAN taxa — species level (EMF):\n")
print(titan_taxa_n_spp)

# Community-level change points: peak in sum(z) marks where the most taxa
# change synchronously. fsumz (filtered) rows are preferred for robustness.
community_thresholds_n_spp <- titan_nitrogen_spp$sumz.cp %>%
  as.data.frame() %>%
  rownames_to_column("group") %>%
  mutate(
    response = case_when(
      group == "sumz-"  ~ "Nitrophobic community (sum z-)",
      group == "sumz+"  ~ "Nitrophilic community (sum z+)",
      group == "fsumz-" ~ "Nitrophobic community (filtered fsum z-)",
      group == "fsumz+" ~ "Nitrophilic community (filtered fsum z+)",
      TRUE ~ group
    ),
    threshold_log_median     = `0.50`,
    threshold_log_lower      = `0.05`,
    threshold_log_upper      = `0.95`,
    threshold_original       = exp(`0.50`),   # back-transform to mg/kg
    threshold_original_lower = exp(`0.05`),
    threshold_original_upper = exp(`0.95`)
  ) %>%
  select(
    group, response,
    threshold_log_median, threshold_log_lower, threshold_log_upper,
    threshold_original, threshold_original_lower, threshold_original_upper, cp
  )

cat("\nCommunity-level nitrogen thresholds — species level (EMF):\n")
print(community_thresholds_n_spp, digits = 4)

# ─── Build species-level taxa plot ───────────────────────────────────────────
# Points = OTUs; x = threshold (mg/kg, log axis); bars = 95% CI; size = |z|
plot_titan_taxa_n_spp <- titan_taxa_n_spp %>%
  filter(response != "Neutral") %>%
  mutate(
    threshold_bt       = exp(threshold_log),
    threshold_lower_bt = exp(threshold_lower),
    threshold_upper_bt = exp(threshold_upper)
  ) %>%
  arrange(lineage, threshold_bt) %>%
  mutate(species = factor(species, levels = unique(species))) %>%
  ggplot(aes(
    x = threshold_bt, y = species, colour = response,
    xmin = threshold_lower_bt, xmax = threshold_upper_bt
  )) +
  geom_errorbarh(height = 0.3, linewidth = 0.6) +
  geom_point(aes(size = abs(zscore))) +
  scale_colour_manual(
    values = c("Nitrophilic" = "#d73027", "Nitrophobic" = "#2166ac"),
    name   = NULL
  ) +
  scale_size_continuous(name = "z score", range = c(1, 5)) +
  scale_x_log10() +
  facet_grid(lineage ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Mineral nitrogen threshold (mg kg⁻¹)", y = NULL,
    title = paste0("TITAN thresholds — ", GUILD_LABELS[["emf"]], " species, nitrogen response")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border       = element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
    strip.text.y       = element_text(angle = 0, hjust = 0, face = "bold", size = 9),
    axis.text.y        = element_text(face = "italic", size = 8),
    legend.position    = "bottom"
  )

# Filtered community thresholds (fsumz) for ribbon + dashed-line overlay
comm_lines_n_spp <- community_thresholds_n_spp %>%
  filter(grepl("^fsumz", group)) %>%
  mutate(
    response_short = if_else(grepl("Nitrophobic", response),
                             "Nitrophobic community", "Nitrophilic community")
  )

plot_titan_taxa_with_comm_n_spp <- plot_titan_taxa_n_spp +
  geom_rect(
    data = comm_lines_n_spp,
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = -Inf, ymax = Inf, fill = response_short),
    alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_vline(
    data = comm_lines_n_spp,
    aes(xintercept = threshold_original, colour = response_short),
    linetype = "dotted", linewidth = 0.8, alpha = 0.7
  ) +
  scale_colour_manual(
    values = c(
      "Nitrophilic" = "#d73027", "Nitrophobic" = "#2166ac",
      "Nitrophilic community" = "#d73027", "Nitrophobic community" = "#2166ac"
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(
      "Nitrophilic community" = "#d73027", "Nitrophobic community" = "#2166ac"
    ),
    name = NULL
  ) +
  labs(
    title   = paste0("TITAN thresholds — ", GUILD_LABELS[["emf"]], " species and community nitrogen response"),
    caption = "Ribbons = 95% bootstrap CI of filtered community thresholds (fsum z)\nDashed lines = median community change point"
  ) +
  theme(legend.position = "bottom", legend.box = "vertical")

print(plot_titan_taxa_with_comm_n_spp)

# ─────────────────────────────────────────────────────────────────────────────
# 5. Lineage-level OTU table
#
# Genera are summed to lineage within each sample BEFORE prevalence filtering,
# using the SRS-normalised genus table (which carries the lineage column).
# ─────────────────────────────────────────────────────────────────────────────
lineage_long <- fread("data/emf/otu_table_srs_emf_genus.txt") %>%
  select(-exploration_type) %>%
  pivot_longer(
    cols      = -c(genus, lineage),
    names_to  = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  group_by(lineage, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# Prevalence filter applied to aggregated lineage abundances
lineage_prev <- lineage_long %>%
  group_by(lineage) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop") %>%
  mutate(prevalence = n_samples / n_sites * 100) %>%
  filter(n_samples >= 5)

otu_table_lin <- lineage_long %>%
  filter(lineage %in% lineage_prev$lineage) %>%
  pivot_wider(names_from = lineage, values_from = abundance, values_fill = 0) %>%
  right_join(sample_metadata %>% select(sample_id), by = "sample_id") %>%
  mutate(across(-sample_id, ~ replace_na(., 0))) %>%
  arrange(sample_id) %>%
  column_to_rownames("sample_id")

stopifnot(length(titan_env_n) == nrow(otu_table_lin))
cat("Lineage TITAN — input OK:", length(titan_env_n), "sites,",
    ncol(otu_table_lin), "lineages\n")

# ─────────────────────────────────────────────────────────────────────────────
# 6. TITAN — lineage level, nitrogen gradient
# ─────────────────────────────────────────────────────────────────────────────
set.seed(1986)
titan_nitrogen_lin <- titan(
  env     = titan_env_n,
  txa     = otu_table_lin,
  numPerm = 1000,
  boot    = TRUE,
  nBoot   = 1000,
  ivTot   = TRUE,
  pur.cut = 0.95,
  rel.cut = 0.90,
  ncpus   = 1       # single-core for reproducibility (matches the species run)
)

# lineage column comes from rownames of sppmax (which are lineage names)
titan_taxa_n_lin <- titan_nitrogen_lin$sppmax %>%
  as.data.frame() %>%
  rownames_to_column("lineage") %>%
  filter(purity >= 0.95, reliability >= 0.9) %>%
  mutate(
    response           = case_when(
      filter == 1 ~ "Nitrophobic",
      filter == 2 ~ "Nitrophilic",
      TRUE        ~ "Neutral"
    ),
    threshold_original = exp(`50%`)
  ) %>%
  arrange(response, z.median) %>%
  select(
    lineage, response,
    threshold_log    = `50%`,
    threshold_lower  = `5%`,
    threshold_upper  = `95%`,
    threshold_original,
    zscore, z.median, reliability, purity
  )

cat("\nSignificant nitrogen TITAN taxa — lineage level (EMF):\n")
print(titan_taxa_n_lin)

community_thresholds_n_lin <- titan_nitrogen_lin$sumz.cp %>%
  as.data.frame() %>%
  rownames_to_column("group") %>%
  mutate(
    response = case_when(
      group == "sumz-"  ~ "Nitrophobic community (sum z-)",
      group == "sumz+"  ~ "Nitrophilic community (sum z+)",
      group == "fsumz-" ~ "Nitrophobic community (filtered fsum z-)",
      group == "fsumz+" ~ "Nitrophilic community (filtered fsum z+)",
      TRUE ~ group
    ),
    threshold_log_median     = `0.50`,
    threshold_log_lower      = `0.05`,
    threshold_log_upper      = `0.95`,
    threshold_original       = exp(`0.50`),
    threshold_original_lower = exp(`0.05`),
    threshold_original_upper = exp(`0.95`)
  ) %>%
  select(
    group, response,
    threshold_log_median, threshold_log_lower, threshold_log_upper,
    threshold_original, threshold_original_lower, threshold_original_upper, cp
  )

cat("\nCommunity-level nitrogen thresholds — lineage level (EMF):\n")
print(community_thresholds_n_lin, digits = 4)

# ─── Build lineage-level taxa plot ───────────────────────────────────────────
# Each point is a lineage; no faceting since lineage is already the y-axis
plot_titan_taxa_n_lin <- titan_taxa_n_lin %>%
  filter(response != "Neutral") %>%
  mutate(
    threshold_bt       = exp(threshold_log),
    threshold_lower_bt = exp(threshold_lower),
    threshold_upper_bt = exp(threshold_upper)
  ) %>%
  arrange(threshold_bt) %>%
  mutate(lineage = factor(lineage, levels = unique(lineage))) %>%
  ggplot(aes(
    x = threshold_bt, y = lineage, colour = response,
    xmin = threshold_lower_bt, xmax = threshold_upper_bt
  )) +
  geom_errorbarh(height = 0.3, linewidth = 0.6) +
  geom_point(aes(size = abs(zscore))) +
  scale_colour_manual(
    values = c("Nitrophilic" = "#d73027", "Nitrophobic" = "#2166ac"),
    name   = NULL
  ) +
  scale_size_continuous(name = "z score", range = c(1, 5)) +
  scale_x_log10() +
  labs(
    x = "Mineral nitrogen threshold (mg kg⁻¹)", y = NULL,
    title = paste0("TITAN thresholds — ", GUILD_LABELS[["emf"]], " lineages, nitrogen response")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border       = element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
    axis.text.y        = element_text(size = 9),
    legend.position    = "bottom"
  )

comm_lines_n_lin <- community_thresholds_n_lin %>%
  filter(grepl("^fsumz", group)) %>%
  mutate(
    response_short = if_else(grepl("Nitrophobic", response),
                             "Nitrophobic community", "Nitrophilic community")
  )

plot_titan_taxa_with_comm_n_lin <- plot_titan_taxa_n_lin +
  geom_rect(
    data = comm_lines_n_lin,
    aes(xmin = threshold_original_lower, xmax = threshold_original_upper,
        ymin = -Inf, ymax = Inf, fill = response_short),
    alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_vline(
    data = comm_lines_n_lin,
    aes(xintercept = threshold_original, colour = response_short),
    linetype = "dotted", linewidth = 0.8, alpha = 0.7
  ) +
  scale_colour_manual(
    values = c(
      "Nitrophilic" = "#d73027", "Nitrophobic" = "#2166ac",
      "Nitrophilic community" = "#d73027", "Nitrophobic community" = "#2166ac"
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(
      "Nitrophilic community" = "#d73027", "Nitrophobic community" = "#2166ac"
    ),
    name = NULL
  ) +
  labs(
    title   = paste0("TITAN thresholds — ", GUILD_LABELS[["emf"]], " lineages and community nitrogen response"),
    caption = "Ribbons = 95% bootstrap CI of filtered community thresholds (fsum z)\nDashed lines = median community change point"
  ) +
  theme(legend.position = "bottom", legend.box = "vertical")

print(plot_titan_taxa_with_comm_n_lin)

# ─────────────────────────────────────────────────────────────────────────────
# Save figure data — only the objects required by code/figure_5.R:
#   titan_taxa_n_spp / _lin  : filtered indicator taxa for threshold plots
#   community_thresholds_n_*  : community-level change points for plot overlays
# ─────────────────────────────────────────────────────────────────────────────
save(
  titan_taxa_n_spp, titan_taxa_n_lin,
  community_thresholds_n_spp, community_thresholds_n_lin,
  file = "generated_data/figure_5.Rdata"
)

cat("\nSaved generated_data/figure_5.Rdata\n")

# ─────────────────────────────────────────────────────────────────────────────
# Export full (unfiltered) TITAN results to Excel
# Each sheet contains the complete sppmax table before purity/reliability
# filtering, allowing inspection of all taxa regardless of significance.
# ─────────────────────────────────────────────────────────────────────────────

wb <- createWorkbook()

addWorksheet(wb, "species")
writeData(wb, "species",
          titan_nitrogen_spp$sppmax %>%
            as.data.frame() %>%
            rownames_to_column("species"))

addWorksheet(wb, "lineage")
writeData(wb, "lineage",
          titan_nitrogen_lin$sppmax %>%
            as.data.frame() %>%
            rownames_to_column("lineage"))

saveWorkbook(wb, "output/titan_emf.xlsx", overwrite = TRUE)

cat("Saved output/titan_emf.xlsx\n")
