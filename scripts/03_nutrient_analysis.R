#!/usr/bin/env Rscript
# Script 03: Nutrient-Diversity Relationship Analysis
# Purpose: Analyze relationships between soil N/P and fungal diversity
# Author: Luke Florence
# Date: 2025

# Load required packages
library(vegan)
library(dplyr)
library(tidyr)
library(readr)
library(MASS)  # For glm.nb (negative binomial GLM)

# Set random seed for reproducibility
set.seed(42)

# Define paths
processed_data_path <- "data/processed"
output_path <- "outputs/tables"

# Create output directory if it doesn't exist
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

cat("===== Nutrient-Diversity Relationship Analysis =====\n\n")

# ============================================================================
# 1. Load Data
# ============================================================================
cat("Step 1: Loading data...\n")

# Load alpha diversity results
if (file.exists(file.path(output_path, "alpha_diversity.csv"))) {
  alpha_div <- read_csv(file.path(output_path, "alpha_diversity.csv"),
                       show_col_types = FALSE)
  cat("  - Loaded alpha diversity data:", nrow(alpha_div), "samples\n")
} else {
  stop("ERROR: alpha_diversity.csv not found. Run 02_diversity_analysis.R first.")
}

# Load NMDS scores
if (file.exists(file.path(output_path, "nmds_scores.csv"))) {
  nmds_data <- read_csv(file.path(output_path, "nmds_scores.csv"),
                       show_col_types = FALSE)
  cat("  - Loaded NMDS scores:", nrow(nmds_data), "samples\n")
} else {
  stop("ERROR: nmds_scores.csv not found. Run 02_diversity_analysis.R first.")
}

# Load OTU table for multivariate analyses
if (file.exists(file.path(processed_data_path, "fungal_counts_filtered.csv"))) {
  otu_table <- read_csv(file.path(processed_data_path, "fungal_counts_filtered.csv"),
                       show_col_types = FALSE)
  
  # Transpose for vegan
  otu_id_col <- names(otu_table)[1]
  sample_cols <- names(otu_table)[-1]
  otu_matrix <- t(otu_table[, sample_cols])
  colnames(otu_matrix) <- otu_table[[otu_id_col]]
  
  cat("  - Loaded OTU table:", ncol(otu_matrix), "OTUs\n")
} else {
  stop("ERROR: fungal_counts_filtered.csv not found.")
}

cat("\n")

# ============================================================================
# 2. Check for Nutrient Variables
# ============================================================================
cat("Step 2: Checking for nutrient variables...\n")

# Identify available nutrient variables
nutrient_vars <- c("N_total", "P_total", "N_available", "P_available", 
                   "NP_ratio", "CN_ratio", "pH")
available_vars <- intersect(nutrient_vars, names(alpha_div))

if (length(available_vars) == 0) {
  cat("  - WARNING: No nutrient variables found in data\n")
  cat("  - Please ensure soil chemistry data was merged in step 01\n")
  cat("  - Available columns:", paste(names(alpha_div), collapse = ", "), "\n")
  stop("Cannot proceed without nutrient data")
}

cat("  - Available nutrient variables:", paste(available_vars, collapse = ", "), "\n")

# Check for collinearity
if (length(available_vars) > 1) {
  cor_matrix <- cor(alpha_div[, available_vars], use = "pairwise.complete.obs")
  cat("  - Correlation among predictors:\n")
  print(round(cor_matrix, 2))
  
  # Flag high correlations (>0.8)
  high_cor <- which(abs(cor_matrix) > 0.8 & cor_matrix != 1, arr.ind = TRUE)
  if (nrow(high_cor) > 0) {
    cat("  - WARNING: High correlations detected (>0.8)\n")
    cat("  - Consider using principal components or selecting fewer predictors\n")
  }
}

cat("\n")

# ============================================================================
# 3. Alpha Diversity ~ Nutrients Models
# ============================================================================
cat("Step 3: Testing nutrient effects on alpha diversity...\n")

# Function to test nutrient effects
test_nutrient_effects <- function(response_var, data, predictors) {
  results <- list()
  
  for (pred in predictors) {
    # Skip if too many NAs
    complete_cases <- complete.cases(data[[response_var]], data[[pred]])
    if (sum(complete_cases) < 10) {
      cat("    - Skipping", pred, "(insufficient data)\n")
      next
    }
    
    # Fit GLM
    formula_str <- paste(response_var, "~", pred)
    
    if (response_var == "richness") {
      # Richness: count data, use Poisson or negative binomial
      model <- tryCatch({
        glm.nb(as.formula(formula_str), data = data[complete_cases, ])
      }, error = function(e) {
        glm(as.formula(formula_str), data = data[complete_cases, ], family = poisson)
      })
    } else {
      # Shannon/Simpson: continuous, use Gaussian
      model <- glm(as.formula(formula_str), data = data[complete_cases, ], 
                  family = gaussian)
    }
    
    # Extract results
    model_summary <- summary(model)
    coef_table <- coef(model_summary)
    
    results[[pred]] <- data.frame(
      response = response_var,
      predictor = pred,
      estimate = coef_table[2, 1],
      std_error = coef_table[2, 2],
      p_value = coef_table[2, 4],
      n = sum(complete_cases),
      R2 = 1 - (model$deviance / model$null.deviance)
    )
  }
  
  return(do.call(rbind, results))
}

# Test effects on richness
cat("  - Testing effects on richness...\n")
richness_results <- test_nutrient_effects("richness", alpha_div, available_vars)

# Test effects on Shannon diversity
cat("  - Testing effects on Shannon diversity...\n")
shannon_results <- test_nutrient_effects("shannon", alpha_div, available_vars)

# Test effects on Simpson diversity
cat("  - Testing effects on Simpson diversity...\n")
simpson_results <- test_nutrient_effects("simpson", alpha_div, available_vars)

# Combine results
alpha_models <- rbind(richness_results, shannon_results, simpson_results)

# Adjust p-values for multiple testing (FDR correction)
alpha_models$p_adjusted <- p.adjust(alpha_models$p_value, method = "fdr")

# Add significance stars
alpha_models$significance <- cut(alpha_models$p_adjusted, 
                                 breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                 labels = c("***", "**", "*", "ns"))

# Save results
write_csv(alpha_models, 
          file.path(output_path, "alpha_diversity_nutrient_models.csv"))
cat("  - Saved alpha diversity model results\n")

# Print summary
cat("\n  - Significant relationships (p < 0.05 after FDR correction):\n")
sig_results <- alpha_models %>% filter(p_adjusted < 0.05)
if (nrow(sig_results) > 0) {
  print(sig_results[, c("response", "predictor", "estimate", "p_adjusted", "significance")])
} else {
  cat("    - No significant relationships detected\n")
}

cat("\n")

# ============================================================================
# 4. Beta Diversity ~ Nutrients (PERMANOVA)
# ============================================================================
cat("Step 4: Testing nutrient effects on community composition (PERMANOVA)...\n")

# Calculate distance matrix
dist_matrix <- vegdist(otu_matrix, method = "bray")

# Prepare environmental data
env_data <- nmds_data %>%
  select(sample_id, all_of(available_vars)) %>%
  filter(complete.cases(.))

# Match samples
matched_samples <- rownames(otu_matrix) %in% env_data$sample_id
dist_matrix_matched <- as.dist(as.matrix(dist_matrix)[matched_samples, matched_samples])
env_data_matched <- env_data %>% arrange(sample_id)

if (nrow(env_data_matched) < 10) {
  cat("  - WARNING: Insufficient samples with complete environmental data\n")
  cat("  - Skipping PERMANOVA\n")
} else {
  # Run PERMANOVA for each nutrient variable
  permanova_results <- list()
  
  for (pred in available_vars) {
    if (sum(!is.na(env_data_matched[[pred]])) >= 10) {
      formula_str <- paste("dist_matrix_matched ~", pred)
      perm_test <- adonis2(as.formula(formula_str), 
                          data = env_data_matched, 
                          permutations = 999,
                          method = "bray")
      
      permanova_results[[pred]] <- data.frame(
        predictor = pred,
        F_statistic = perm_test$F[1],
        R2 = perm_test$R2[1],
        p_value = perm_test$`Pr(>F)`[1]
      )
    }
  }
  
  # Combine PERMANOVA results
  if (length(permanova_results) > 0) {
    permanova_table <- do.call(rbind, permanova_results)
    permanova_table$p_adjusted <- p.adjust(permanova_table$p_value, method = "fdr")
    permanova_table$significance <- cut(permanova_table$p_adjusted,
                                        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                        labels = c("***", "**", "*", "ns"))
    
    write_csv(permanova_table, 
              file.path(output_path, "permanova_results.csv"))
    cat("  - Saved PERMANOVA results\n")
    
    # Print summary
    cat("\n  - PERMANOVA results:\n")
    print(permanova_table)
  }
}

cat("\n")

# ============================================================================
# 5. Multivariate Analysis (RDA/CCA)
# ============================================================================
cat("Step 5: Multivariate ordination with environmental constraints...\n")

if (nrow(env_data_matched) >= 10 && length(available_vars) > 0) {
  # Prepare data
  otu_hellinger <- decostand(otu_matrix[matched_samples, ], method = "hellinger")
  
  # Build formula with all available nutrients
  rda_formula <- paste("otu_hellinger ~", paste(available_vars, collapse = " + "))
  
  # Run RDA
  rda_model <- rda(as.formula(rda_formula), data = env_data_matched)
  
  cat("  - RDA complete\n")
  cat("    Variance explained by constraints:", 
      round(RsquareAdj(rda_model)$adj.r.squared * 100, 1), "%\n")
  
  # Test significance
  rda_anova <- anova.cca(rda_model, permutations = 999)
  cat("    Overall significance: p =", rda_anova$`Pr(>F)`[1], "\n")
  
  # Test individual axes
  rda_axes <- anova.cca(rda_model, by = "axis", permutations = 999)
  cat("    RDA1 significance: p =", rda_axes$`Pr(>F)`[1], "\n")
  
  # Extract results
  rda_summary <- summary(rda_model)
  
  # Save RDA results
  rda_results <- list(
    variance_explained = RsquareAdj(rda_model)$adj.r.squared,
    overall_pvalue = rda_anova$`Pr(>F)`[1],
    axis1_pvalue = rda_axes$`Pr(>F)`[1]
  )
  
  saveRDS(rda_results, file.path(output_path, "rda_results.rds"))
  cat("  - Saved RDA results\n")
  
} else {
  cat("  - Insufficient data for RDA analysis\n")
}

cat("\n")

# ============================================================================
# 6. Summary Report
# ============================================================================
cat("Step 6: Generating summary report...\n")

# Create report
sink(file.path(output_path, "nutrient_analysis_summary.txt"))
cat("===== Nutrient-Diversity Analysis Summary =====\n")
cat("Generated:", as.character(Sys.time()), "\n\n")

cat("Available nutrient variables:\n")
cat(paste(" ", available_vars, collapse = "\n"), "\n\n")

if (exists("alpha_models")) {
  cat("Alpha diversity models - Significant results:\n")
  print(sig_results[, c("response", "predictor", "estimate", "p_adjusted")])
  cat("\n")
}

if (exists("permanova_table")) {
  cat("PERMANOVA results:\n")
  print(permanova_table)
  cat("\n")
}

if (exists("rda_results")) {
  cat("RDA results:\n")
  cat("  Adjusted R-squared:", round(rda_results$variance_explained, 3), "\n")
  cat("  Overall p-value:", round(rda_results$overall_pvalue, 3), "\n")
  cat("\n")
}

cat("===== End of Summary =====\n")
sink()

cat("  - Saved summary report\n")

cat("\n===== Nutrient Analysis Complete =====\n")
cat("Results saved to:", output_path, "\n")

# Save session info for reproducibility
sink(file.path(output_path, "03_session_info.txt"))
sessionInfo()
sink()
