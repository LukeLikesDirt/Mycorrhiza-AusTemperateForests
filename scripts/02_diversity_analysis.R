#!/usr/bin/env Rscript
# Script 02: Diversity Analysis
# Purpose: Calculate alpha and beta diversity metrics
# Author: Luke Florence
# Date: 2025

# Load required packages
library(vegan)
library(dplyr)
library(tidyr)
library(readr)

# Set random seed for reproducibility
set.seed(42)

# Define paths
processed_data_path <- "data/processed"
output_path <- "outputs/tables"

# Create output directory if it doesn't exist
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

cat("===== Diversity Analysis =====\n\n")

# ============================================================================
# 1. Load Processed Data
# ============================================================================
cat("Step 1: Loading processed data...\n")

# Load filtered OTU table
if (file.exists(file.path(processed_data_path, "fungal_counts_filtered.csv"))) {
  otu_table <- read_csv(file.path(processed_data_path, "fungal_counts_filtered.csv"),
                        show_col_types = FALSE)
  cat("  - Loaded OTU table:", nrow(otu_table), "OTUs\n")
} else {
  stop("ERROR: fungal_counts_filtered.csv not found. Run 01_data_processing.R first.")
}

# Load complete metadata
if (file.exists(file.path(processed_data_path, "metadata_complete.csv"))) {
  metadata <- read_csv(file.path(processed_data_path, "metadata_complete.csv"),
                      show_col_types = FALSE)
  cat("  - Loaded metadata:", nrow(metadata), "samples\n")
} else {
  stop("ERROR: metadata_complete.csv not found. Run 01_data_processing.R first.")
}

# Load taxonomy
if (file.exists(file.path(processed_data_path, "fungal_taxonomy_clean.csv"))) {
  taxonomy <- read_csv(file.path(processed_data_path, "fungal_taxonomy_clean.csv"),
                      show_col_types = FALSE)
  cat("  - Loaded taxonomy:", nrow(taxonomy), "OTUs\n")
} else {
  cat("  - WARNING: fungal_taxonomy_clean.csv not found\n")
  taxonomy <- NULL
}

cat("\n")

# ============================================================================
# 2. Prepare Data for Vegan
# ============================================================================
cat("Step 2: Preparing data for diversity calculations...\n")

# Convert OTU table to matrix format (samples as rows, OTUs as columns)
# Assuming first column is OTU_ID
otu_id_col <- names(otu_table)[1]
sample_cols <- names(otu_table)[-1]

# Transpose OTU table for vegan (samples as rows)
otu_matrix <- t(otu_table[, sample_cols])
colnames(otu_matrix) <- otu_table[[otu_id_col]]

cat("  - OTU matrix dimensions:", nrow(otu_matrix), "samples x", 
    ncol(otu_matrix), "OTUs\n")

# Check if rarefaction is needed
read_depth <- rowSums(otu_matrix)
cat("  - Read depth range:", min(read_depth), "-", max(read_depth), "\n")

# Rarefy to minimum depth if variation is high
if (max(read_depth) / min(read_depth) > 2) {
  min_depth <- min(read_depth)
  cat("  - Rarefying to", min_depth, "reads per sample...\n")
  otu_matrix_rare <- rrarefy(otu_matrix, sample = min_depth)
  cat("  - Rarefaction complete\n")
} else {
  cat("  - Sequencing depth is relatively even, using unrarefied data\n")
  otu_matrix_rare <- otu_matrix
}

cat("\n")

# ============================================================================
# 3. Calculate Alpha Diversity
# ============================================================================
cat("Step 3: Calculating alpha diversity metrics...\n")

# Calculate diversity indices
alpha_diversity <- data.frame(
  sample_id = rownames(otu_matrix_rare),
  richness = specnumber(otu_matrix_rare),  # Number of species
  shannon = diversity(otu_matrix_rare, index = "shannon"),  # Shannon index
  simpson = diversity(otu_matrix_rare, index = "simpson"),  # Simpson index
  evenness = diversity(otu_matrix_rare, index = "shannon") / log(specnumber(otu_matrix_rare))  # Pielou's evenness
)

cat("  - Alpha diversity summary:\n")
cat("    Richness: mean =", round(mean(alpha_diversity$richness), 1), 
    ", SD =", round(sd(alpha_diversity$richness), 1), "\n")
cat("    Shannon: mean =", round(mean(alpha_diversity$shannon), 2), 
    ", SD =", round(sd(alpha_diversity$shannon), 2), "\n")
cat("    Simpson: mean =", round(mean(alpha_diversity$simpson), 2), 
    ", SD =", round(sd(alpha_diversity$simpson), 2), "\n")

# Merge with metadata
alpha_diversity_full <- alpha_diversity %>%
  left_join(metadata, by = "sample_id")

# Save alpha diversity results
write_csv(alpha_diversity_full, 
          file.path(output_path, "alpha_diversity.csv"))
cat("  - Saved alpha diversity results\n")

cat("\n")

# ============================================================================
# 4. Calculate Beta Diversity
# ============================================================================
cat("Step 4: Calculating beta diversity...\n")

# Calculate dissimilarity matrices
bray_dist <- vegdist(otu_matrix_rare, method = "bray")  # Bray-Curtis
jaccard_dist <- vegdist(otu_matrix_rare, method = "jaccard", binary = TRUE)  # Jaccard

cat("  - Calculated Bray-Curtis dissimilarity\n")
cat("  - Calculated Jaccard distance\n")

# Perform NMDS ordination
set.seed(42)
nmds_bray <- metaMDS(otu_matrix_rare, distance = "bray", k = 2, 
                     trymax = 100, trace = FALSE)
cat("  - NMDS ordination complete (stress =", round(nmds_bray$stress, 3), ")\n")

# Extract NMDS scores
nmds_scores <- scores(nmds_bray, display = "sites") %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(.))

# Merge with metadata
nmds_results <- nmds_scores %>%
  left_join(metadata, by = "sample_id")

# Save NMDS results
write_csv(nmds_results, 
          file.path(output_path, "nmds_scores.csv"))
cat("  - Saved NMDS ordination scores\n")

# Save dissimilarity matrices
write.csv(as.matrix(bray_dist), 
          file.path(output_path, "bray_curtis_dissimilarity.csv"))
cat("  - Saved Bray-Curtis dissimilarity matrix\n")

cat("\n")

# ============================================================================
# 5. Taxonomic Summary
# ============================================================================
cat("Step 5: Generating taxonomic summaries...\n")

if (!is.null(taxonomy)) {
  # Calculate relative abundance for each sample
  otu_rel_abund <- sweep(otu_matrix_rare, 1, rowSums(otu_matrix_rare), "/")
  
  # Aggregate by phylum (or another taxonomic level)
  if ("Phylum" %in% names(taxonomy)) {
    # Match OTU IDs
    taxonomy_matched <- taxonomy %>%
      filter(.data[[otu_id_col]] %in% colnames(otu_rel_abund))
    
    # Sum abundances by phylum for each sample
    phylum_abund <- data.frame(otu_rel_abund) %>%
      mutate(otu_id = colnames(otu_rel_abund)) %>%
      pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance") %>%
      left_join(taxonomy_matched %>% select(all_of(c(otu_id_col, "Phylum"))), 
                by = setNames(otu_id_col, "otu_id")) %>%
      group_by(sample_id, Phylum) %>%
      summarise(rel_abundance = sum(abundance), .groups = "drop")
    
    write_csv(phylum_abund, 
              file.path(output_path, "phylum_relative_abundance.csv"))
    cat("  - Saved phylum-level relative abundance\n")
  }
  
  # Overall taxonomic composition
  top_phyla <- taxonomy %>%
    count(Phylum, sort = TRUE) %>%
    head(10)
  
  write_csv(top_phyla, 
            file.path(output_path, "top_phyla.csv"))
  cat("  - Saved top 10 phyla\n")
}

cat("\n")

# ============================================================================
# 6. Summary Statistics
# ============================================================================
cat("Step 6: Generating summary statistics...\n")

# Create comprehensive summary
summary_stats <- list(
  "Alpha Diversity Summary" = summary(alpha_diversity[, -1]),
  "NMDS Stress" = nmds_bray$stress,
  "Number of Samples" = nrow(otu_matrix_rare),
  "Number of OTUs" = ncol(otu_matrix_rare),
  "Total Read Depth" = sum(otu_matrix_rare)
)

# Save summary
sink(file.path(output_path, "diversity_summary.txt"))
cat("===== Diversity Analysis Summary =====\n")
cat("Generated:", as.character(Sys.time()), "\n\n")
print(summary_stats)
sink()

cat("  - Saved diversity summary\n")

cat("\n===== Diversity Analysis Complete =====\n")
cat("Results saved to:", output_path, "\n")

# Save session info for reproducibility
sink(file.path(output_path, "02_session_info.txt"))
sessionInfo()
sink()
