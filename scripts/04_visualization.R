#!/usr/bin/env Rscript
# Script 04: Visualization
# Purpose: Generate publication-quality figures
# Author: Luke Florence
# Date: 2025

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)  # For combining plots

# Set random seed for reproducibility
set.seed(42)

# Define paths
output_tables <- "outputs/tables"
output_figures <- "outputs/figures"

# Create output directory if it doesn't exist
dir.create(output_figures, showWarnings = FALSE, recursive = TRUE)

# Set ggplot theme
theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.minor = blank(),
                  legend.position = "right"))

cat("===== Visualization =====\n\n")

# ============================================================================
# 1. Load Data
# ============================================================================
cat("Step 1: Loading data for visualization...\n")

# Load alpha diversity
if (file.exists(file.path(output_tables, "alpha_diversity.csv"))) {
  alpha_div <- read_csv(file.path(output_tables, "alpha_diversity.csv"),
                       show_col_types = FALSE)
  cat("  - Loaded alpha diversity data\n")
} else {
  cat("  - WARNING: alpha_diversity.csv not found\n")
  alpha_div <- NULL
}

# Load NMDS scores
if (file.exists(file.path(output_tables, "nmds_scores.csv"))) {
  nmds_data <- read_csv(file.path(output_tables, "nmds_scores.csv"),
                       show_col_types = FALSE)
  cat("  - Loaded NMDS scores\n")
} else {
  cat("  - WARNING: nmds_scores.csv not found\n")
  nmds_data <- NULL
}

# Load model results
if (file.exists(file.path(output_tables, "alpha_diversity_nutrient_models.csv"))) {
  model_results <- read_csv(file.path(output_tables, "alpha_diversity_nutrient_models.csv"),
                           show_col_types = FALSE)
  cat("  - Loaded model results\n")
} else {
  cat("  - WARNING: model results not found\n")
  model_results <- NULL
}

# Load PERMANOVA results
if (file.exists(file.path(output_tables, "permanova_results.csv"))) {
  permanova <- read_csv(file.path(output_tables, "permanova_results.csv"),
                       show_col_types = FALSE)
  cat("  - Loaded PERMANOVA results\n")
} else {
  cat("  - WARNING: PERMANOVA results not found\n")
  permanova <- NULL
}

# Load phylum abundance
if (file.exists(file.path(output_tables, "phylum_relative_abundance.csv"))) {
  phylum_abund <- read_csv(file.path(output_tables, "phylum_relative_abundance.csv"),
                          show_col_types = FALSE)
  cat("  - Loaded phylum abundance data\n")
} else {
  cat("  - WARNING: phylum abundance not found\n")
  phylum_abund <- NULL
}

cat("\n")

# ============================================================================
# 2. Alpha Diversity Distributions
# ============================================================================
cat("Step 2: Creating alpha diversity distribution plots...\n")

if (!is.null(alpha_div)) {
  # Reshape data for faceting
  alpha_long <- alpha_div %>%
    select(sample_id, richness, shannon, simpson) %>%
    pivot_longer(-sample_id, names_to = "metric", values_to = "value")
  
  # Create violin/box plots
  p1 <- ggplot(alpha_long, aes(x = metric, y = value, fill = metric)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    facet_wrap(~metric, scales = "free") +
    labs(title = "Distribution of Alpha Diversity Metrics",
         x = "Diversity Metric", y = "Value") +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "grey90"))
  
  ggsave(file.path(output_figures, "alpha_diversity_distributions.png"),
         plot = p1, width = 10, height = 4, dpi = 300)
  cat("  - Saved alpha_diversity_distributions.png\n")
}

cat("\n")

# ============================================================================
# 3. Nutrient-Diversity Relationships
# ============================================================================
cat("Step 3: Creating nutrient-diversity relationship plots...\n")

if (!is.null(alpha_div)) {
  # Check for available nutrient variables
  nutrient_vars <- c("N_total", "P_total", "NP_ratio")
  available_nutrients <- intersect(nutrient_vars, names(alpha_div))
  
  if (length(available_nutrients) > 0) {
    # Create scatter plots for each nutrient-diversity combination
    plots <- list()
    
    for (nutrient in available_nutrients) {
      # Richness
      p_rich <- ggplot(alpha_div, aes_string(x = nutrient, y = "richness")) +
        geom_point(alpha = 0.6, size = 2) +
        geom_smooth(method = "lm", se = TRUE, color = "blue") +
        labs(title = paste("Richness vs", nutrient),
             x = nutrient, y = "OTU Richness") +
        theme_bw()
      
      # Shannon
      p_shan <- ggplot(alpha_div, aes_string(x = nutrient, y = "shannon")) +
        geom_point(alpha = 0.6, size = 2) +
        geom_smooth(method = "lm", se = TRUE, color = "blue") +
        labs(title = paste("Shannon vs", nutrient),
             x = nutrient, y = "Shannon Diversity") +
        theme_bw()
      
      # Combine plots
      combined <- p_rich + p_shan + plot_layout(ncol = 2)
      
      ggsave(file.path(output_figures, paste0("diversity_vs_", nutrient, ".png")),
             plot = combined, width = 10, height = 4, dpi = 300)
      cat("  - Saved diversity_vs_", nutrient, ".png\n", sep = "")
    }
  } else {
    cat("  - No nutrient variables available for plotting\n")
  }
}

cat("\n")

# ============================================================================
# 4. NMDS Ordination Plot
# ============================================================================
cat("Step 4: Creating NMDS ordination plot...\n")

if (!is.null(nmds_data)) {
  # Check for available environmental variables
  env_vars <- c("N_total", "P_total", "pH")
  available_env <- intersect(env_vars, names(nmds_data))
  
  if (length(available_env) > 0) {
    # Create NMDS plot colored by first environmental variable
    color_var <- available_env[1]
    
    p4 <- ggplot(nmds_data, aes_string(x = "NMDS1", y = "NMDS2", color = color_var)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = median(nmds_data[[color_var]], na.rm = TRUE)) +
      labs(title = "NMDS Ordination of Fungal Communities",
           subtitle = paste("Colored by", color_var),
           x = "NMDS1", y = "NMDS2") +
      theme_bw() +
      theme(legend.position = "right")
    
    ggsave(file.path(output_figures, "nmds_ordination.png"),
           plot = p4, width = 8, height = 6, dpi = 300)
    cat("  - Saved nmds_ordination.png\n")
    
    # Create multiple NMDS plots for different variables
    if (length(available_env) > 1) {
      nmds_plots <- list()
      for (var in available_env) {
        nmds_plots[[var]] <- ggplot(nmds_data, aes_string(x = "NMDS1", y = "NMDS2", color = var)) +
          geom_point(size = 2.5, alpha = 0.7) +
          scale_color_gradient2(low = "blue", mid = "white", high = "red",
                               midpoint = median(nmds_data[[var]], na.rm = TRUE)) +
          labs(title = var, x = "NMDS1", y = "NMDS2") +
          theme_bw() +
          theme(legend.position = "bottom")
      }
      
      combined_nmds <- wrap_plots(nmds_plots, ncol = min(3, length(available_env)))
      
      ggsave(file.path(output_figures, "nmds_multiple_variables.png"),
             plot = combined_nmds, width = 12, height = 4 * ceiling(length(available_env)/3), 
             dpi = 300)
      cat("  - Saved nmds_multiple_variables.png\n")
    }
  } else {
    # Basic NMDS without environmental coloring
    p4 <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
      geom_point(size = 3, alpha = 0.7) +
      labs(title = "NMDS Ordination of Fungal Communities",
           x = "NMDS1", y = "NMDS2") +
      theme_bw()
    
    ggsave(file.path(output_figures, "nmds_ordination.png"),
           plot = p4, width = 8, height = 6, dpi = 300)
    cat("  - Saved nmds_ordination.png\n")
  }
}

cat("\n")

# ============================================================================
# 5. Taxonomic Composition
# ============================================================================
cat("Step 5: Creating taxonomic composition plots...\n")

if (!is.null(phylum_abund)) {
  # Get top phyla
  top_phyla <- phylum_abund %>%
    group_by(Phylum) %>%
    summarise(total_abund = sum(rel_abundance), .groups = "drop") %>%
    arrange(desc(total_abund)) %>%
    head(10) %>%
    pull(Phylum)
  
  # Filter for top phyla
  phylum_plot_data <- phylum_abund %>%
    mutate(Phylum = ifelse(Phylum %in% top_phyla, Phylum, "Other"))
  
  # Create stacked bar plot
  p5 <- ggplot(phylum_plot_data, aes(x = sample_id, y = rel_abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Paired") +
    labs(title = "Fungal Community Composition at Phylum Level",
         x = "Sample", y = "Relative Abundance") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "right")
  
  ggsave(file.path(output_figures, "phylum_composition.png"),
         plot = p5, width = 12, height = 6, dpi = 300)
  cat("  - Saved phylum_composition.png\n")
}

cat("\n")

# ============================================================================
# 6. Model Results Summary Plot
# ============================================================================
cat("Step 6: Creating model results summary plot...\n")

if (!is.null(model_results)) {
  # Create coefficient plot
  p6 <- ggplot(model_results, 
               aes(x = predictor, y = estimate, color = significance)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = estimate - std_error, 
                     ymax = estimate + std_error),
                 width = 0.2, position = position_dodge(width = 0.5)) +
    facet_wrap(~response, scales = "free_y") +
    scale_color_manual(values = c("***" = "red", "**" = "orange", 
                                   "*" = "yellow", "ns" = "grey60")) +
    labs(title = "Effects of Nutrients on Alpha Diversity",
         subtitle = "Error bars show standard error",
         x = "Predictor Variable", 
         y = "Model Coefficient",
         color = "Significance") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey90"))
  
  ggsave(file.path(output_figures, "model_coefficients.png"),
         plot = p6, width = 10, height = 6, dpi = 300)
  cat("  - Saved model_coefficients.png\n")
}

cat("\n")

# ============================================================================
# 7. PERMANOVA Results Visualization
# ============================================================================
cat("Step 7: Creating PERMANOVA results visualization...\n")

if (!is.null(permanova)) {
  # Create bar plot of R2 values
  p7 <- ggplot(permanova, aes(x = reorder(predictor, R2), y = R2, 
                              fill = significance)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("***" = "red", "**" = "orange", 
                                  "*" = "yellow", "ns" = "grey60")) +
    coord_flip() +
    labs(title = "PERMANOVA Results: Variance Explained by Nutrients",
         subtitle = "Effect on community composition (Bray-Curtis distance)",
         x = "Predictor Variable", 
         y = "R² (Variance Explained)",
         fill = "Significance") +
    theme_bw()
  
  ggsave(file.path(output_figures, "permanova_results.png"),
         plot = p7, width = 8, height = 6, dpi = 300)
  cat("  - Saved permanova_results.png\n")
}

cat("\n")

# ============================================================================
# 8. Summary Figure (Multi-panel)
# ============================================================================
cat("Step 8: Creating multi-panel summary figure...\n")

# Combine key figures into a summary panel
# This requires that individual plots were created above
if (!is.null(alpha_div) && !is.null(nmds_data)) {
  # Recreate simplified versions for the panel
  
  # Panel A: Richness distribution
  pa <- ggplot(alpha_div, aes(x = "", y = richness)) +
    geom_violin(fill = "lightblue", alpha = 0.6) +
    geom_boxplot(width = 0.2) +
    labs(title = "A) OTU Richness", x = "", y = "Richness") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))
  
  # Panel B: Shannon distribution
  pb <- ggplot(alpha_div, aes(x = "", y = shannon)) +
    geom_violin(fill = "lightgreen", alpha = 0.6) +
    geom_boxplot(width = 0.2) +
    labs(title = "B) Shannon Diversity", x = "", y = "Shannon") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))
  
  # Panel C: NMDS
  pc <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(title = "C) NMDS Ordination", x = "NMDS1", y = "NMDS2") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))
  
  # Combine panels
  summary_fig <- (pa + pb) / pc + plot_layout(heights = c(1, 1.2))
  
  ggsave(file.path(output_figures, "summary_figure.png"),
         plot = summary_fig, width = 10, height = 8, dpi = 300)
  cat("  - Saved summary_figure.png\n")
}

cat("\n===== Visualization Complete =====\n")
cat("Figures saved to:", output_figures, "\n")

# Save session info for reproducibility
sink(file.path(output_figures, "04_session_info.txt"))
sessionInfo()
sink()
