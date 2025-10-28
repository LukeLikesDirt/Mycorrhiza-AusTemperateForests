
# Run from the main project directory 

# Required packages and functions
require(Biostrings)
require(SRS)
require(data.table)
require(tidyverse)
source("./code/helper_functions/abundance_filters.R")

# Read in OTU table and sequences
otu_table <- fread("./data/ASVs.txt") %>%
  rename(otu_id = OTU_ID)
otu_seqs <- readDNAStringSet("./data/ASVs.fasta")

# (1) Remove puutative index-switching artifacts -------------------------------

# Remove rare occurrences of abundant otus using a 0.1% threshold (1 in 1000)
# Computes the total read count per OTU, then removes any occurrence of that OTU
# from a sample if it is less than 0.1% of the total reads for that OTU.
otu_table_filtered <- filter_library(otu_table, threshold = 0.1)

# Compute the number of reads before and after filtering
cat("Total reads before filtering:", sum(otu_table[, -1]), "\n")
cat("Total reads after filtering:", sum(otu_table_filtered[, -1]), "\n")

# Compute the proportion of reads removed
cat("Percentage of reads removed by filtering:", ((sum(otu_table[, -1]) - sum(otu_table_filtered[, -1])) / sum(otu_table[, -1])) * 100, "\n")

# (1) Assess minimum read depth ------------------------------------------------

# Generate a tibble
otu_tibble <- otu_table_filtered %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id", 
    values_to = "abundance"
    ) %>%
  filter(abundance > 0)

# Compute the read depth per sample_id
sample_id_depth <- otu_tibble %>%
  # Count total reads per sample_id
  group_by(sample_id) %>%
  summarise(n_seqs = sum(abundance)) %>%
  arrange(n_seqs) %>%
  print(n = 20)

# Visualise sample_id depth and range
sample_id_depth %>%
  ggplot(aes(x = 1:nrow(.), y = n_seqs)) +
  scale_y_log10() +
  geom_line() +
  geom_point()
sample_id_depth %>%
  ggplot(aes(x = 1, y = n_seqs)) +
  scale_y_log10() +
  geom_jitter()

# Set the minimum read depth threshold
min_depth <- 3750
max_depth <- max(sample_id_depth$n_seqs)

# Create a vector of low abundance sample_ids
low_abundance_sample_ids <- sample_id_depth %>%
  filter(n_seqs < min_depth) %>%
  pull(sample_id)

# OTU richness per sample_id
otu_richness = otu_tibble %>%
  group_by(sample_id) %>%
  summarise(n_otus = n_distinct(otu_id)) %>%
  arrange(n_otus)

# Evaluate the relationship between read depth and OTU richness
sample_id_depth %>%
  left_join(otu_richness, by = "sample_id") %>%
  # Remove the low abundance sample_ids
  filter(n_seqs >= min_depth) %>%
  ggplot(aes(x = n_seqs, y = n_otus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(
    # Add the r2 label
    aes(label = after_stat(rr.label))
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Read depth", y = "OTU richness")

# (3) SRS normalisation --------------------------------------------------------

# Remove the low abundance sample_ids from the OTU table
otu_table_sample_filtered <- otu_table_filtered %>%
  # Remove the low abundance sample_ids
  select(-all_of(low_abundance_sample_ids)) %>%
  # Remove any lost OTUs
  filter(rowSums(select(., -otu_id)) > 0)

# SRS expects columns are sample_ids and rows are OTUs
otus_srs <- SRS(
  # OTU ID's as rownames
  data = otu_table_sample_filtered %>% column_to_rownames(var = 'otu_id'),
  Cmin = min_depth,  # <-- change to the minimum depth of the prevalence filtered OTU table for beta-diversity analyses
  set_seed = TRUE,  # <-- for reproducibility
  seed = 1986
) %>%
  # Convert the SRS output to a data frame
  as_tibble() %>%
  # Bind the OTU_ID's and SRS normalised counts
  bind_cols(
    tibble(otu_id = otu_table_sample_filtered[["otu_id"]]),
    .
  ) %>%
  # Remove any OTUs that are now zero abundance
  filter(rowSums(select(., -otu_id)) > 0)

# Compare OTU richness before and after SRS normalisation
cat("Number of OTUs before SRS normalisation:", nrow(otu_table_sample_filtered), "\n")
cat("Number of OTUs after SRS normalisation:", nrow(otus_srs), "\n")
cat("Percentage of OTUs removed by SRS normalisation:", ((nrow(otu_table_sample_filtered) - nrow(otus_srs)) / nrow(otu_table_sample_filtered)) * 100, "\n")

# Re-compute richness after SRS normalisation
otu_tibble_srs <- otus_srs %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id", 
    values_to = "abundance"
  ) %>%
  filter(abundance > 0)

# Regraph the relationship between read depth and OTU richness after SRS
otu_tibble_srs %>%
  group_by(sample_id) %>%
  summarise(n_otus = n_distinct(otu_id)) %>%
  left_join(sample_id_depth, by = "sample_id") %>%
  ggplot(aes(x = n_seqs, y = n_otus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(
    # Add the r2 label
    aes(label = after_stat(rr.label))
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Read depth", y = "OTU richness")

# (4) Filter ASV sequences -----------------------------------------------------

# Remove size annotations from sequence names
otu_seqs_names_clean <- str_remove(names(otu_seqs), ";size=.*$")
names(otu_seqs) <- otu_seqs_names_clean

# Filter sequences to match SRS-normalized OTU table
otu_seqs_filtered <- otu_seqs[names(otu_seqs) %in% otus_srs$otu_id]

# Recalculate size annotations from SRS-normalized data
otu_seqs_sizes <- otu_tibble_srs %>%
  group_by(otu_id) %>%
  summarise(size = sum(abundance)) %>%
  mutate(otu_id_size = paste0(otu_id, ";size=", size)) %>%
  select(otu_id, otu_id_size)

# Apply new size annotations to filtered sequences
otu_seqs_filtered_size <- otu_seqs_filtered
names(otu_seqs_filtered_size) <- otu_seqs_sizes$otu_id_size[
  match(names(otu_seqs_filtered_size), otu_seqs_sizes$otu_id)
]

# Check that sequences and OTU IDs match between otu_seqs_filtered and otu_seqs_filtered_size
# Remove size info from both
otu_ids_original <- str_remove(names(otu_seqs_filtered), ";size=.*$")
otu_ids_renamed <- str_remove(names(otu_seqs_filtered_size), ";size=.*$")

# Check sequences and OTU IDs match
sequences_match <- all(as.character(otu_seqs_filtered) == as.character(otu_seqs_filtered_size))
otu_ids_match <- all(otu_ids_original == otu_ids_renamed)

# Report results
if (!sequences_match || !otu_ids_match) {
  error_msg <- c()
  
  if (!sequences_match) {
    n_seq_mismatches <- sum(as.character(otu_seqs_filtered) != as.character(otu_seqs_filtered_size))
    error_msg <- c(error_msg, sprintf("- %d sequences do not match", n_seq_mismatches))
  }
  
  if (!otu_ids_match) {
    n_id_mismatches <- sum(otu_ids_original != otu_ids_renamed)
    error_msg <- c(error_msg, sprintf("- %d OTU IDs do not match", n_id_mismatches))
  }
  
  stop(paste(c("Error:", error_msg), collapse = "\n"))
} else {
  message("✓ All sequences and OTU IDs match between otu_seqs_filtered and otu_seqs_filtered_size")
}

# Save the filtered ASV table and sequences ------------------------------------

# Write the asv table
fwrite(otus_srs, "./data/ASVs_filtered.txt", sep = "\t")

# Write the asv fasta
writeXStringSet(otu_seqs_filtered_size, "./data/ASVs_filtered.fasta")

