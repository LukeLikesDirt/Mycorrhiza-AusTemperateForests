#!/usr/bin/env Rscript
# =============================================================================
# Extract Unique Sequences by Taxonomic Rank
# Description: Create multiple FASTA files with unique sequences:
#              1. Glomeromycota phylum only
#              2. All taxa with identified species
#              3. All taxa with identified genus
#              4. All taxa with identified family
# =============================================================================

# Load required packages
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(Biostrings)
})

# Input and output file paths
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 13) {
  stop("Usage: Rscript extract_unique_sequences.R <input_fasta> <classification_file> <glom_out> <glom_prev_species> <glom_prev_genus> <glom_prev_family> <glom_prev_order> <glom_prev_class> <species_out> <genus_out> <family_out> <all_out> <prevalence_cutoff>")
}

fasta_in <- args[1]
classification_in <- args[2]
glom_out <- args[3]
glom_prev_species <- args[4]
glom_prev_genus <- args[5]
glom_prev_family <- args[6]
glom_prev_order <- args[7]
glom_prev_class <- args[8]
species_out <- args[9]
genus_out <- args[10]
family_out <- args[11]
all_out <- args[12]
prevalence_cutoff <- as.numeric(args[13])

# Validate cutoff
if (is.na(prevalence_cutoff) || prevalence_cutoff <= 0 || prevalence_cutoff >= 1) {
  stop("prevalence_cutoff must be a number between 0 and 1 (e.g., 0.5 for 50%, 0.66 for 66%)")
}

cat(sprintf("Using prevalence cutoff: %.0f%%\n\n", prevalence_cutoff * 100))

# =============================================================================
# HELPER FUNCTION: Calculate taxonomic completeness score
# =============================================================================
calc_tax_score <- function(kingdom, phylum, class, order, family, genus, species) {
  score <- 0
  if (!is.na(kingdom) & kingdom != "" & kingdom != "unidentified") score <- score + 1
  if (!is.na(phylum) & phylum != "" & phylum != "unidentified") score <- score + 1
  if (!is.na(class) & class != "" & class != "unidentified") score <- score + 1
  if (!is.na(order) & order != "" & order != "unidentified") score <- score + 1
  if (!is.na(family) & family != "" & family != "unidentified") score <- score + 1
  if (!is.na(genus) & genus != "" & genus != "unidentified") score <- score + 1
  if (!is.na(species) & species != "" & species != "unidentified") score <- score + 1
  return(score)
}

# =============================================================================
# HELPER FUNCTION: Process sequences to keep best taxonomic information
# =============================================================================
process_unique_seqs <- function(seq_ids, classification_df, fasta_seqs) {
  # Subset classification data
  subset_df <- classification_df %>%
    filter(id %in% seq_ids)
  
  # Extract sequences - only keep IDs that exist in FASTA
  subset_seqs <- fasta_seqs[names(fasta_seqs) %in% seq_ids]
  
  # Filter classification to only include IDs present in FASTA
  subset_df <- subset_df %>%
    filter(id %in% names(subset_seqs))
  
  # Add sequences to dataframe
  subset_df <- subset_df %>%
    mutate(sequence = as.character(subset_seqs[match(id, names(subset_seqs))]))
  
  # Calculate taxonomic completeness score
  subset_df <- subset_df %>%
    rowwise() %>%
    mutate(tax_score = calc_tax_score(kingdom, phylum, class, order, family, genus, species)) %>%
    ungroup()
  
  # Create a taxonomy string (excluding id and sequence)
  subset_df <- subset_df %>%
    mutate(taxonomy = paste(kingdom, phylum, class, order, family, genus, species, sep = "|"))

  # Calculate sequence length
  subset_df <- subset_df %>%
    mutate(seq_length = nchar(sequence))

  # Group by sequence, keeping the one with highest taxonomic score
  # If tied, keep first occurrence
  filtered_df <- subset_df %>%
    group_by(sequence) %>%
    arrange(desc(tax_score)) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  # Create DNAStringSet with filtered IDs
  filtered_ids <- filtered_df$id
  filtered_seqs <- subset_seqs[names(subset_seqs) %in% filtered_ids]
  
  # Return results
  list(
    sequences = filtered_seqs,
    n_original = length(seq_ids),
    n_unique = length(filtered_seqs),
    n_removed = length(seq_ids) - length(filtered_seqs),
    n_missing = length(seq_ids) - nrow(subset_df)
  )
}

# =============================================================================
# HELPER FUNCTION: Prevalence Filtering function with configurable cutoff
# =============================================================================
# Filters sequences so that no single taxon at the specified rank accounts 
# for more than the specified proportion of sequences.
# For example, at genus level with cutoff=0.5, no single genus can account 
# for >50% of sequences.

prevalence_filter <- function(df, rank, cutoff = 0.5) {
  # Ensure we only work with sequences that have valid classification at this rank
  # Exclude unidentified, unclassified, and incertae sedis taxa
  df <- df %>%
    filter(
      !is.na(!!sym(rank)) & 
        !!sym(rank) != "" &
        !tolower(!!sym(rank)) %in% c("unidentified", "unclassified") &
        !grepl("incertae sedis", !!sym(rank), ignore.case = TRUE) &
        !grepl("sp\\.", !!sym(rank))
    )
  
  # Count sequences per taxon at this rank
  rank_counts <- df %>%
    group_by(!!sym(rank)) %>%
    summarise(n = n(), .groups = 'drop') %>%
    arrange(desc(n))  # Sort from largest to smallest
  
  # If only one taxon, return all sequences
  if (nrow(rank_counts) == 1) {
    cat(sprintf("  Only one %s found, including all %d sequences\n\n", 
                rank, nrow(df)))
    return(df)
  }
  
  # Get the largest taxon
  largest_taxon <- rank_counts[[1, rank]]
  largest_count <- rank_counts$n[1]
  
  # Sum all sequences from smaller taxa (non-dominant)
  smaller_taxa_count <- sum(rank_counts$n[-1])
  
  # Calculate max allowed from largest taxon
  # If cutoff = 0.66, then:
  # largest / (largest + smaller) = 0.66
  # largest = 0.66 * (largest + smaller)
  # largest = 0.66 * largest + 0.66 * smaller
  # largest - 0.66 * largest = 0.66 * smaller
  # 0.34 * largest = 0.66 * smaller
  # largest = (0.66 / 0.34) * smaller
  # largest = (cutoff / (1 - cutoff)) * smaller
  
  max_from_largest <- floor((cutoff / (1 - cutoff)) * smaller_taxa_count)
  
  selected_ids <- c()
  
  # First, take ALL sequences from smaller taxa
  for (i in 2:nrow(rank_counts)) {
    current_taxon <- rank_counts[[i, rank]]
    taxon_ids <- df %>%
      filter(!!sym(rank) == current_taxon) %>%
      pull(id)
    selected_ids <- c(selected_ids, taxon_ids)
  }
  
  # Now handle the largest taxon
  largest_ids <- df %>%
    filter(!!sym(rank) == largest_taxon) %>%
    pull(id)
  
  if (largest_count <= max_from_largest) {
    # Largest taxon is already <= max allowed, take all
    selected_ids <- c(selected_ids, largest_ids)
    final_pct <- (largest_count / (largest_count + smaller_taxa_count)) * 100
    cat(sprintf("  Largest %s '%s': included all %d sequences (%.1f%% of total)\n",
                rank, largest_taxon, largest_count, final_pct))
  } else {
    # Randomly sample from largest to match max allowed
    sampled_ids <- sample(largest_ids, size = max_from_largest)
    selected_ids <- c(selected_ids, sampled_ids)
    final_pct <- (max_from_largest / (max_from_largest + smaller_taxa_count)) * 100
    cat(sprintf("  Largest %s '%s': sampled %d of %d sequences (%.1f%% of final dataset)\n",
                rank, largest_taxon, max_from_largest, largest_count, final_pct))
  }
  
  # Filter the dataframe to selected IDs
  filtered_df <- df %>%
    filter(id %in% selected_ids)
  
  cat(sprintf("  Total: filtered from %d to %d sequences at %s rank\n", 
              nrow(df), nrow(filtered_df), rank))
  cat(sprintf("  Number of %s represented: %d\n",
              rank, length(unique(filtered_df[[rank]]))))
  
  # Verify the actual percentage
  actual_pct <- (max(table(filtered_df[[rank]])) / nrow(filtered_df)) * 100
  cat(sprintf("  Largest %s now accounts for %.1f%% of filtered sequences (target: %.0f%%)\n\n",
              rank, actual_pct, cutoff * 100))
  
  return(filtered_df)
}

# =============================================================================
# READ INPUT FILES
# =============================================================================
cat("Reading classification file...\n")
classification_df <- read_tsv(classification_in, show_col_types = FALSE)

cat("Reading FASTA file...\n")
fasta_seqs <- readDNAStringSet(fasta_in)

cat("Total sequences in FASTA:", length(fasta_seqs), "\n\n")

# =============================================================================
# 1. GLOMEROMYCOTA UNIQUE SEQUENCES
# =============================================================================
cat("=== PROCESSING GLOMEROMYCOTA ===\n")
glom_df <- classification_df %>%
  filter(phylum == "Glomeromycota")

cat("Found", nrow(glom_df), "Glomeromycota sequences\n")

if (nrow(glom_df) > 0) {
  glom_results <- process_unique_seqs(glom_df$id, classification_df, fasta_seqs)
  
  writeXStringSet(glom_results$sequences, glom_out)
  
  cat("Original sequences:", glom_results$n_original, "\n")
  if (glom_results$n_missing > 0) {
    cat("Missing from FASTA:", glom_results$n_missing, "\n")
  }
  cat("Unique sequences (after deduplication):", glom_results$n_unique, "\n")
  cat("Removed (duplicates + lower taxonomic info):", glom_results$n_removed, "\n")
  cat("Output:", glom_out, "\n\n")
} else {
  cat("WARNING: No Glomeromycota sequences found\n\n")
}

# =============================================================================
# 2. APPLY PREVALENCE FILTERING TO GLOMEROMYCOTA AT EACH RANK
# =============================================================================

cat("=== PREVALENCE FILTERING GLOMEROMYCOTA ===\n\n")

# Read the unique sequences (output from previous step)
glom_unique <- readDNAStringSet(glom_out)
glom_unique_ids <- names(glom_unique)

# Get classification info for these unique sequences
glom_unique_df <- classification_df %>%
  filter(id %in% glom_unique_ids)

# Define the taxonomic ranks to filter
rank_outputs <- list(
  list(rank = "class", output = glom_prev_class),
  list(rank = "order", output = glom_prev_order),
  list(rank = "family", output = glom_prev_family),
  list(rank = "genus", output = glom_prev_genus),
  list(rank = "species", output = glom_prev_species)
)

# Apply prevalence filtering at each rank
for (level in rank_outputs) {
  cat(sprintf("--- Filtering at %s rank ---\n", level$rank))
  
  # Filter to sequences that have valid, identified classifications at this rank
  # Exclude unidentified, unclassified, and incertae sedis taxa
  rank_df <- glom_unique_df %>%
    filter(
      !is.na(!!sym(level$rank)) & 
        !!sym(level$rank) != "" &
        !tolower(!!sym(level$rank)) %in% c("unidentified", "unclassified") &
        !grepl("incertae sedis", !!sym(level$rank), ignore.case = TRUE) &
        !grepl("sp\\.", !!sym(level$rank))
    )
  
  cat(sprintf("  Starting with %d sequences with valid %s classification\n", 
              nrow(rank_df), level$rank))
  
  if (nrow(rank_df) > 0) {
    # Apply prevalence filter
    filtered_df <- prevalence_filter(rank_df, level$rank, cutoff = prevalence_cutoff)
    
    # Extract sequences
    filtered_seqs <- glom_unique[names(glom_unique) %in% filtered_df$id]
    
    # Write output
    writeXStringSet(filtered_seqs, level$output)
    cat(sprintf("  Written %d sequences to: %s\n\n", 
                length(filtered_seqs), level$output))
  } else {
    cat(sprintf("  WARNING: No sequences with valid %s classification\n\n", 
                level$rank))
  }
}

cat("=== PREVALENCE FILTERING COMPLETE ===\n")

# =============================================================================
# 3. IDENTIFIED SPECIES (ALL TAXA)
# =============================================================================
cat("=== PROCESSING IDENTIFIED SPECIES ===\n")
species_df <- classification_df %>%
  filter(species != "unidentified" & species != "" & !is.na(species))

cat("Found", nrow(species_df), "sequences with identified species\n")

if (nrow(species_df) > 0) {
  species_results <- process_unique_seqs(species_df$id, classification_df, fasta_seqs)
  
  writeXStringSet(species_results$sequences, species_out)
  
  cat("Original sequences:", species_results$n_original, "\n")
  if (species_results$n_missing > 0) {
    cat("Missing from FASTA:", species_results$n_missing, "\n")
  }
  cat("Unique sequences (after deduplication):", species_results$n_unique, "\n")
  cat("Removed (duplicates + lower taxonomic info):", species_results$n_removed, "\n")
  cat("Output:", species_out, "\n\n")
} else {
  cat("WARNING: No sequences with identified species found\n\n")
}

# =============================================================================
# 4. IDENTIFIED GENUS (ALL TAXA)
# =============================================================================
cat("=== PROCESSING IDENTIFIED GENUS ===\n")
genus_df <- classification_df %>%
  filter(genus != "unidentified" & genus != "" & !is.na(genus))

cat("Found", nrow(genus_df), "sequences with identified genus\n")

if (nrow(genus_df) > 0) {
  genus_results <- process_unique_seqs(genus_df$id, classification_df, fasta_seqs)
  
  writeXStringSet(genus_results$sequences, genus_out)
  
  cat("Original sequences:", genus_results$n_original, "\n")
  if (genus_results$n_missing > 0) {
    cat("Missing from FASTA:", genus_results$n_missing, "\n")
  }
  cat("Unique sequences (after deduplication):", genus_results$n_unique, "\n")
  cat("Removed (duplicates + lower taxonomic info):", genus_results$n_removed, "\n")
  cat("Output:", genus_out, "\n\n")
} else {
  cat("WARNING: No sequences with identified genus found\n\n")
}

# =============================================================================
# 5. IDENTIFIED FAMILY (ALL TAXA)
# =============================================================================
cat("=== PROCESSING IDENTIFIED FAMILY ===\n")
family_df <- classification_df %>%
  filter(family != "unidentified" & family != "" & !is.na(family))

cat("Found", nrow(family_df), "sequences with identified family\n")

if (nrow(family_df) > 0) {
  family_results <- process_unique_seqs(family_df$id, classification_df, fasta_seqs)
  
  writeXStringSet(family_results$sequences, family_out)
  
  cat("Original sequences:", family_results$n_original, "\n")
  cat("Unique sequences (after deduplication):", family_results$n_unique, "\n")
  cat("Removed (duplicates + lower taxonomic info):", family_results$n_removed, "\n")
  cat("Output:", family_out, "\n\n")
} else {
  cat("WARNING: No sequences with identified family found\n\n")
}

# =============================================================================
# 6. UNIQUE SEQUENCES FOR ALL TAXA
# =============================================================================
cat("=== PROCESSING ALL UNIQUE SEQUENCES ===\n")
all_df <- classification_df

cat("Found", nrow(all_df), "total sequences\n")

if (nrow(all_df) > 0) {
  all_results <- process_unique_seqs(all_df$id, classification_df, fasta_seqs)
  
  writeXStringSet(all_results$sequences, all_out)
  
  cat("Original sequences:", all_results$n_original, "\n")
  if (all_results$n_missing > 0) {
    cat("Missing from FASTA:", all_results$n_missing, "\n")
  }
  cat("Unique sequences (after deduplication):", all_results$n_unique, "\n")
  cat("Removed (duplicates + lower taxonomic info):", all_results$n_removed, "\n")
  cat("Output:", all_out, "\n\n")
} else {
  cat("WARNING: No sequences found\n\n")
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== SUMMARY OF PREVELENCE FILTERED DATASETS ===\n\n")
for (level in rank_outputs) {
  if (file.exists(level$output)) {
    seqs <- readDNAStringSet(level$output)
    seq_ids <- names(seqs)
    
    # Get classification for these sequences
    seq_df <- classification_df %>%
      filter(id %in% seq_ids)
    
    # Count taxa at this rank
    rank_dist <- table(seq_df[[level$rank]])
    max_taxon <- names(which.max(rank_dist))
    max_count <- max(rank_dist)
    max_pct <- (max_count / length(seqs)) * 100
    
    cat(sprintf("%s level: %d sequences\n", 
                toupper(level$rank), length(seqs)))
    cat(sprintf("  Largest %s: %s (%d sequences, %.1f%%)\n", 
                level$rank, max_taxon, max_count, max_pct))
    cat(sprintf("  Number of %s: %d\n\n", 
                level$rank, length(unique(seq_df[[level$rank]]))))
  }
}
cat("=== FINAL SUMMARY ===\n")
cat("Input FASTA:", fasta_in, "\n")
cat("Classification file:", classification_in, "\n")
cat("Total sequences processed:", length(fasta_seqs), "\n")
cat("\nDeduplication strategy:\n")
cat("  - For identical sequences: kept highest taxonomic completeness\n")
cat("\nOutput files created:\n")
cat("  1. Glomeromycota unique:", glom_out, "\n")
cat("  2. Glomeromycota class prevalence filtered:", glom_prev_class, "\n")
cat("  3. Glomeromycota order prevalence filtered:", glom_prev_order, "\n")
cat("  4. Glomeromycota family prevalence filtered:", glom_prev_family, "\n")
cat("  5. Glomeromycota genus prevalence filtered:", glom_prev_genus, "\n")
cat("  6. Glomeromycota species prevalence filtered:", glom_prev_species, "\n")
cat("  7. Eukaryote all unique:", all_out, "\n")
cat("  8. Eukaryote family unique:", family_out, "\n")
cat("  9. Eukaryote genus unique:", genus_out, "\n")
cat(" 10. Eukaryote species unique:", species_out, "\n")

cat("\nProcessing complete!\n")