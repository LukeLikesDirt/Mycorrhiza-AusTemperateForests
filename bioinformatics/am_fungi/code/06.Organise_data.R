# Arbuscular Mycorrhizal Fungi Pipeline to organise bioinformatic outputs for
# downstram stream statistical analyses

# Run from the am_fungi project root directory

# 1. SETUP ---------------------------------------------------------------------

## 1.1 Load required libraries -------------------------------------------------
require(data.table)
require(vegan)
require(SRS)
require(ape)
require(picante)
require(phangorn)
require(DECIPHER)
require(Biostrings)
require(parallel)
require(foreach)
require(doParallel)
require(tidyverse)

## 1.2 Helper functions --------------------------------------------------------

# Compute Shannon diversity
shannon_diversity <- function(x) {
  # Remove zeros before calculating proportions
  x <- x[x > 0]
  
  # If no observations, return 0
  if (length(x) == 0) return(0)
  
  p <- x / sum(x)
  -sum(p * log(p))
}

# Index switching filter
filter_library <- function(library_specific_otu_table, threshold = 0.1) {
  
  # Remove rare occurrences of abundant otus using a threshold
  # This filters based on each OTU's distribution across samples
  
  filtered_otu_table <- library_specific_otu_table %>%
    # Convert to long format for grouped operations
    pivot_longer(-otu_id, names_to = "sample",
                 values_to = "otu_sample_abundance") %>%
    
    # Group by OTU to calculate library-wide statistics
    group_by(otu_id) %>%
    mutate(
      # Calculate total abundance of this OTU across all samples
      otu_library_abundance = sum(otu_sample_abundance),
      # Calculate relative abundance of this OTU in each sample
      # (as a percentage of the OTU's total library abundance)
      rel_abundance = (otu_sample_abundance / otu_library_abundance) * 100
    ) %>%
    
    # Filter: set to zero if relative abundance within OTU <= threshold (default 0.1%)
    mutate(otu_sample_abundance = replace(
      otu_sample_abundance, rel_abundance <= threshold, 0)) %>%
    ungroup() %>%
    
    # Remove temporary calculation columns
    select(-c(otu_library_abundance, rel_abundance)) %>%
    
    # Convert back to wide format
    pivot_wider(
      names_from = "sample",
      values_from = "otu_sample_abundance",
      values_fill = 0
    )
  
  return(filtered_otu_table)
  
}

# Parallel rarefaction
rrarefy_parallel <- function(x, sample, n_cores = NULL) {
  x <- as.matrix(x)
  
  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  
  # Split samples into chunks for parallel processing
  n_samples <- nrow(x)
  sample_indices <- split(1:n_samples, cut(1:n_samples, n_cores, labels = FALSE))
  
  # Parallel rarefaction
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("x", "sample"), envir = environment())
  clusterEvalQ(cl, library(vegan))
  
  rarefied_chunks <- parLapply(cl, sample_indices, function(indices) {
    rrarefy(x[indices, , drop = FALSE], sample[indices])
  })
  
  stopCluster(cl)
  
  # Combine results
  do.call(rbind, rarefied_chunks)
}

# 2. DATA IMPORT ---------------------------------------------------------------

## 2.1 Import glomeromycota data -----------------------------------------------
# Read in the classification table for glomeromycota
otu_classification <- fread("./code/tmp_clusters/otu_glom_classification.txt")

# Read sequences for glomeromycota OTUs
otu_seqs <- readDNAStringSet("./code/tmp_clusters/otu_glom_sequences.fasta")

# Read in OTU table for glomeromycota
otu_table_glom <- fread("./code/tmp_clusters/otu_glom_table.txt")

## 2.2 Import non-glomeromycota data -------------------------------------------
# Read OTU table for all non-glomeromycota otus
non_glom_ids <- fread("./code/tmp_clusters/family_non_glomeromycota_clusters.txt") %>%
  mutate(otu_id = gsub(";size=.*", "", otu_id)) %>%
  pull(otu_id)

unid_ids <- fread("./code/tmp_clusters/family_unidentified_clusters.txt") %>%
  mutate(otu_id = gsub(";size=.*", "", otu_id)) %>%
  pull(otu_id)

otu_table_non_glom <- fread("./data/ASVs.txt") %>%
  rename(otu_id = OTU_ID) %>%
  filter(otu_id %in% non_glom_ids | otu_id %in% unid_ids)

## 2.3 Combine OTU tables ------------------------------------------------------
# Generate otu tibble for glomeromycota
otu_tibble_glom <- otu_table_glom %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id", 
    values_to = "abundance"
  )

# Generate otu tibble for non-glomeromycota
otu_tibble_non_glom <- otu_table_non_glom %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id", 
    values_to = "abundance"
  )

# Join the glomeromycota and non-glomeromycota OTU tables
otu_table_all <- bind_rows(otu_tibble_glom, otu_tibble_non_glom) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = abundance,
    values_fill = 0
  )

# 3. QUALITY CONTROL -----------------------------------------------------------

## 3.1 Remove putative index-switching artefacts -------------------------------
otu_table_library_filtered <- filter_library(otu_table_all)

message("---- Library-wide filter (index-switching artefacts) ----")
message("Total reads before filtering: ", sum(otu_table_all[, -1]))
message("Total reads after filtering:  ", sum(otu_table_library_filtered[, -1]))
message("Percentage of reads removed:  ",
        round(((sum(otu_table_all[, -1]) - sum(otu_table_library_filtered[, -1])) /
                 sum(otu_table_all[, -1])) * 100, 3), "%")

## 3.2 Assess minimum read depth -----------------------------------------------

# Generate a tibble
otu_tibble <- otu_table_library_filtered  %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id", 
    values_to = "abundance"
  ) %>%
  filter(abundance > 0)

# Compute the read depth per sample
sample_id_depth <- otu_tibble %>%
  # Count total reads per sample_id
  group_by(sample_id) %>%
  summarise(
    n_seqs = sum(abundance),
    n_otus = n_distinct(otu_id)
  ) %>%
  arrange(n_seqs) %>%
  print(n = Inf)

# Set minimum read depth threshold
min_depth <- 3493

# Identify low abundance samples to exclude
low_abundance_sample_ids <- sample_id_depth %>%
  filter(n_seqs < min_depth) %>%
  pull(sample_id)

## 3.3 Visualise sample depth and richness relationships -----------------------

# Visualise sample depth distribution
sample_id_depth %>%
  ggplot(aes(x = 1:nrow(.), y = n_seqs)) +
  geom_hline(yintercept = min_depth, linetype = "dashed", color = "red") +
  scale_y_log10() +
  geom_line() +
  geom_point() +
  theme(aspect.ratio = 1)

sample_id_depth %>%
  ggplot(aes(x = 1, y = n_seqs)) +
  geom_hline(yintercept = min_depth, linetype = "dashed", color = "red") +
  scale_y_log10() +
  geom_jitter() +
  theme(aspect.ratio = 1)

# Evaluate the relationship between read depth and OTU richness
sample_id_depth %>%
  filter(n_seqs > min_depth) %>%
  ggplot(aes(x = n_seqs, y = n_otus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label))) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Read depth", y = "OTU richness") +
  theme_classic() +
  theme(aspect.ratio = 1)

# 4. PREPARE RAW AM DATA -------------------------------------------------------

## 4.1 Identify AM samples and OTUs --------------------------------------------

# AM sample identifiers (exclude low abundance samples)
am_sample_ids <- otu_table_glom %>%
  select(-otu_id, -all_of(low_abundance_sample_ids)) %>%
  names()

# AM OTU identifiers
am_otu_ids <- otu_table_library_filtered %>%
  select(-all_of(low_abundance_sample_ids)) %>%
  filter(otu_id %in% otu_table_glom$otu_id) %>%
  # Filter any lost OTUs
  filter(rowSums(select(., -otu_id)) > 0) %>%
  pull(otu_id)

## 4.2 Create and save raw AM OTU table ----------------------------------------

# Create output directory if it doesn't exist
if (!dir.exists("./output")) {
  dir.create("./output")
}

# Save the raw AM OTU table
otu_table_glom_raw <- otu_table_library_filtered %>%
  filter(otu_id %in% am_otu_ids) %>%
  select(otu_id, all_of(am_sample_ids))

otu_table_glom_raw %>%
  fwrite("./output/otu_table.txt", sep = "\t")

# Save the AM classification table
otu_classification %>%
  filter(otu_id %in% am_otu_ids) %>%
  fwrite("./output/otu_classification.txt", sep = "\t")

# Filter and save the AM sequences
otu_seqs_am <- otu_seqs[names(otu_seqs) %in% am_otu_ids]
writeXStringSet(otu_seqs_am, "./output/otu_sequences.fasta")

# 5. PHYLOGENETIC TREE ---------------------------------------------------------

## 5.1 Generate phylogenetic tree ----------------------------------------------

# Convert to DNAStringSet for alignment
dna_seqs <- readDNAStringSet("output/otu_sequences_am.fasta")

# Multiple sequence alignment
aligned_seqs <- AlignSeqs(dna_seqs, verbose = TRUE, processors = 10)

# Convert to phyDat format
phyDat_seqs <- as.phyDat(as.matrix(aligned_seqs), type = "DNA")

# Create distance matrix
dm <- dist.ml(phyDat_seqs)

# Build neighbor-joining tree
nj_tree <- NJ(dm)

# Root the tree (using midpoint rooting)
rooted_tree <- midpoint(nj_tree)

# Verify tree structure
print(paste("Tree is rooted:", ape::is.rooted(rooted_tree)))
print(paste("Number of tips:", length(rooted_tree$tip.label)))

# Write rooted tree to file
ape::write.tree(rooted_tree, file = "output/tree.newick")

# 6. AITCHISON DISTANCE (PREVALENCE FILTERED) ----------------------------------

## 6.1 Setup gemelli environment -----------------------------------------------

# Create a directory for gemelli output
system("mkdir -p data/distances")

# Convert to biom
system('conda run -n gemelli_env biom convert -i output/otu_table.txt -o data/distances/otus.biom --table-type="OTU table" --to-json')

## 6.2 Run gemelli: Identity-based (prevalence filtered) -----------------------

# Run gemelli for taxa that occur in at least 2 samples
# With jst over 100 samples, min frequency of 1 requires presence in at least 2 
# samples
system("conda run -n gemelli_env gemelli rpca --in-biom data/distances/otus.biom --output-dir data/distances/ --min-feature-frequency 1")

# Rename output files
system("mv data/distances/ordination.txt data/distances/ordination_identity.txt")
system("mv data/distances/distance-matrix.tsv output/dist_identity.txt")

## 6.3 Extract ordination results: Identity-based ------------------------------

# Read the ordination file
ordination_lines <- readLines("data/distances/ordination_identity.txt")

# Find the line indices for different sections
proportion_line <- which(grepl("Proportion explained", ordination_lines))[1]
site_line <- which(grepl("^Site", ordination_lines))[1]

# Extract the proportion explained values
proportions <- as.numeric(strsplit(ordination_lines[proportion_line + 1], "\t")[[1]])

# Convert to percentages and round
prop1 <- round(proportions[1] * 100, 0)
prop2 <- round(proportions[2] * 100, 0)

# Read the site coordinates
site_start <- site_line + 1
site_end <- length(ordination_lines)

# Find where site data ends
for (i in site_start:length(ordination_lines)) {
  if (ordination_lines[i] == "" || grepl("^Biplot", ordination_lines[i]) || 
      grepl("^Site constraints", ordination_lines[i])) {
    site_end <- i - 1
    break
  }
}
site_data <- ordination_lines[site_start:site_end]

# Parse site data into a data frame
site_coords <- data.frame(
  sample_id = character(),
  PC1 = numeric(),
  PC2 = numeric(),
  stringsAsFactors = FALSE
)

for (line in site_data) {
  if (line != "" && !grepl("^\t", line)) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      site_coords <- rbind(site_coords, data.frame(
        sample_id = parts[1],
        PC1 = as.numeric(parts[2]),
        PC2 = as.numeric(parts[3]),
        stringsAsFactors = FALSE
      )) %>%
        filter(!sample_id %in% c("Biplot", "Site constraints"))
    }
  }
}

# Create ordination results with dynamic column names
ordination_results <- site_coords %>%
  rename(
    !!paste0("pco1_identity_", prop1) := PC1,
    !!paste0("pco2_identity_", prop2) := PC2
  ) %>%
  select(sample_id, everything())

## 6.4 Run gemelli: Phylogeny-based (prevalence filtered) ----------------------

# Run phylogenetic-rpca with prevalence filtering
system("conda run -n gemelli_env gemelli phylogenetic-rpca --in-biom data/distances/otus.biom --in-phylogeny output/tree_am.newick --output-dir data/distances/ --min-feature-frequency 1")

# Rename output files
system("mv data/distances/ordination.txt data/distances/ordination_phylogeny.txt")
system("mv data/distances/distance-matrix.tsv output/dist_phylogeny.txt")

## 6.5 Extract ordination results: Phylogeny-based -----------------------------

# Read the ordination file
ordination_lines <- readLines("data/distances/ordination_phylogeny.txt")

# Find the line indices for different sections
proportion_line <- which(grepl("Proportion explained", ordination_lines))[1]
site_line <- which(grepl("^Site", ordination_lines))[1]

# Extract the proportion explained values
proportions <- as.numeric(strsplit(ordination_lines[proportion_line + 1], "\t")[[1]])

# Convert to percentages and round
prop1 <- round(proportions[1] * 100, 0)
prop2 <- round(proportions[2] * 100, 0)

# Read the site coordinates
site_start <- site_line + 1
site_end <- length(ordination_lines)

# Find where site data ends
for (i in site_start:length(ordination_lines)) {
  if (ordination_lines[i] == "" || grepl("^Biplot", ordination_lines[i]) || 
      grepl("^Site constraints", ordination_lines[i])) {
    site_end <- i - 1
    break
  }
}
site_data <- ordination_lines[site_start:site_end]

# Parse site data into a data frame
site_coords <- data.frame(
  sample_id = character(),
  PC1 = numeric(),
  PC2 = numeric(),
  stringsAsFactors = FALSE
)

for (line in site_data) {
  if (line != "" && !grepl("^\t", line)) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      site_coords <- rbind(site_coords, data.frame(
        sample_id = parts[1],
        PC1 = as.numeric(parts[2]),
        PC2 = as.numeric(parts[3]),
        stringsAsFactors = FALSE
      )) %>%
        filter(!sample_id %in% c("Biplot", "Site constraints"))
    }
  }
}

# Add phylogeny-based ordination to results
ordination_results <- ordination_results %>%
  left_join(
    site_coords %>%
      rename(
        !!paste0("pco1_phylogeny_", prop1) := PC1,
        !!paste0("pco2_phylogeny_", prop2) := PC2
      ) %>%
      select(sample_id, everything()),
    by = "sample_id"
  )

## 6.7 Save ordination results -------------------------------------------------

fwrite(ordination_results, "data/distances/ordination_combined.txt", sep = "\t")

# 7. ALPHA DIVERSITY -----------------------------------------------------------

## 7.1 Rarefied AM richness ----------------------------------------------------

# Set parameters
set.seed(1986)
n_sims <- 1000

# Filter otu table to only samples with sufficient depth
otu_table_abundant <- otu_table_library_filtered %>%
  select(-all_of(low_abundance_sample_ids)) %>%
  # Remove any lost OTUs
  filter(rowSums(select(., -otu_id)) > 0)

# Prepare matrix
otu_matrix <- otu_table_abundant %>%
  column_to_rownames(var = "otu_id") %>%
  t()

# Identify AM OTU column indices
am_otu_cols <- which(colnames(otu_matrix) %in% am_otu_ids)
n_samples <- nrow(otu_matrix)

# Pre-allocate results matrix
rarefaction_matrix <- matrix(NA_integer_, 
                             nrow = n_samples, 
                             ncol = n_sims,
                             dimnames = list(rownames(otu_matrix), NULL))

# Setup parallel processing
n_cores <- detectCores() - 1
message(paste0("Using ", n_cores, " cores for parallel processing"))
message(paste0("Performing ", n_sims, " rarefactions..."))

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary objects
clusterExport(cl, c("otu_matrix", "am_otu_cols", "min_depth"))
clusterEvalQ(cl, library(vegan))

# Parallel rarefaction
rarefaction_matrix <- foreach(
  sim = 1:n_sims,
  .combine = 'cbind',
  .packages = 'vegan'
) %dopar% {
  # Rarefy entire matrix at once
  rare_mat <- rrarefy(otu_matrix, sample = min_depth)
  
  # Extract AM OTUs and count richness per sample
  am_rare <- rare_mat[, am_otu_cols, drop = FALSE]
  rowSums(am_rare > 0)  # Richness = number of non-zero OTUs
}

stopCluster(cl)

# Calculate summary statistics
rare_am_richness <- tibble(
  sample_id = rownames(rarefaction_matrix),
  mean_richness = rowMeans(rarefaction_matrix),
  sd_richness = apply(rarefaction_matrix, 1, sd),
  ci_lower = apply(rarefaction_matrix, 1, quantile, probs = 0.025),
  ci_upper = apply(rarefaction_matrix, 1, quantile, probs = 0.975),
  median_richness = apply(rarefaction_matrix, 1, median)
) %>%
  arrange(mean_richness)

# Evaluate rarefied AM richness to sequence depth relationship
rare_am_richness %>%
  left_join(sample_id_depth, by = "sample_id") %>%
  ggplot(aes(x = n_seqs, y = mean_richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label))) +
  scale_x_log10() +
  labs(x = "Read depth", y = "Mean rarefied AM OTU richness") +
  theme_classic() +
  theme(aspect.ratio = 1)

## 7.2 SRS normalisation -------------------------------------------------------

# SRS expects columns are sample_ids and rows are OTUs
otus_srs <- SRS(
  data = otu_table_abundant %>% column_to_rownames(var = 'otu_id'),
  Cmin = min_depth,
  set_seed = TRUE,
  seed = 1986
) %>%
  as_tibble() %>%
  bind_cols(
    tibble(otu_id = otu_table_abundant[["otu_id"]]),
    .
  ) %>%
  # Remove any OTUs that are now zero abundance
  filter(rowSums(select(., -otu_id)) > 0)

# Compare OTU richness before and after SRS normalisation
cat("Number of OTUs before SRS normalisation:", nrow(otu_table_abundant), "\n")
cat("Number of OTUs after SRS normalisation:", nrow(otus_srs), "\n")
cat("Percentage of OTUs removed by SRS normalisation:", 
    ((nrow(otu_table_abundant) - nrow(otus_srs)) / nrow(otu_table_abundant)) * 100, "\n")

## 7.3 Re-compute richness after SRS normalisation -----------------------------

# Generate tibble
otu_tibble_srs <- otus_srs %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id", 
    values_to = "abundance"
  ) %>%
  filter(abundance > 0)

# Total SRS richness
srs_richness <- otu_tibble_srs %>%
  group_by(sample_id) %>%
  summarise(n_otus = n_distinct(otu_id)) %>%
  ungroup() %>%
  arrange(n_otus) %>%
  print(n = Inf)

# AM SRS richness
srs_am_richness <- otu_tibble_srs %>%
  filter(otu_id %in% am_otu_ids) %>%
  group_by(sample_id) %>%
  summarise(n_otus = n_distinct(otu_id)) %>%
  ungroup() %>%
  arrange(n_otus) %>%
  print(n = Inf)

# Visualise relationship between read depth and SRS richness
srs_richness %>%
  left_join(
    sample_id_depth %>%
      select(sample_id, n_seqs),
    by = "sample_id") %>%
  ggplot(aes(x = n_seqs, y = n_otus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label))) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Read depth", y = "OTU richness") +
  theme_classic() +
  theme(aspect.ratio = 1)

# Visualise AM richness after SRS normalisation
otu_tibble_srs %>%
  filter(otu_id %in% am_otu_ids) %>%
  group_by(sample_id) %>%
  summarise(n_otus = n_distinct(otu_id)) %>%
  left_join(
    sample_id_depth %>%
      select(sample_id, n_seqs),
    by = "sample_id") %>%
  ggplot(aes(x = n_seqs, y = n_otus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label))) +
  scale_x_log10() +
  labs(x = "Read depth", y = "AM OTU richness") +
  theme_classic() +
  theme(aspect.ratio = 1)

# Compare rarefied and SRS normalised AM OTU richness
otu_tibble_srs %>%
  filter(otu_id %in% am_otu_ids) %>%
  group_by(sample_id) %>%
  summarise(n_otus = n_distinct(otu_id)) %>%
  ungroup() %>%
  full_join(
    rare_am_richness %>%
      select(sample_id, median_richness),
    by = "sample_id"
  ) %>%
  ggplot(aes(x = median_richness, y = n_otus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label))) +
  labs(x = "Median rarefied AM OTU richness", y = "SRS normalised AM OTU richness") +
  theme_classic() +
  theme(aspect.ratio = 1)

# 8. GENERATE FINAL OTU TABLES -------------------------------------------------

# Identify samples in glom vs filtered
samples_glom <- otu_table_glom %>%
  select(-otu_id, -all_of(low_abundance_sample_ids)) %>%
  names()

## 8.1 Non-normalised AM OTU table ---------------------------------------------

# Generate AM fungi non-normalised table
otu_raw_AM <- otu_table_glom %>%
  # Remove low abundance sample_ids
  select(-all_of(low_abundance_sample_ids)) %>%
  # Filter where row sums = 0
  filter(rowSums(select(., -otu_id)) > 0)

## 8.2 SRS-normalised AM OTU table ---------------------------------------------

# Generate AM fungi SRS normalised table
otus_srs_AM <- otus_srs %>%
  # Remove low abundance sample_ids
  select(otu_id, all_of(samples_glom)) %>%
  # Filter to AM OTUs only
  filter(otu_id %in% am_otu_ids) %>%
  # Filter where row sums = 0
  filter(rowSums(select(., -otu_id)) > 0)

## 8.3 SRS-normalised family-level OTU table -----------------------------------

# Load non-glom classification
non_glom_classification <- fread("./code/tmp_clusters/family_non_glomeromycota_clusters.txt") %>%
  mutate(otu_id = gsub(";size=.*", "", otu_id))

# Get glom IDs
glom_ids <- am_otu_ids

# Subset the OTU table to identified families only
otu_raw_families <- otu_table_library_filtered %>%
  filter(otu_id %in% non_glom_ids | otu_id %in% glom_ids)

otu_srs_families <- otus_srs %>%
  filter(otu_id %in% non_glom_ids | otu_id %in% glom_ids)

# Add phylum, family and genus levels to the classification table
family_classification_all <- bind_rows(
  otu_classification %>%
    mutate(phylum = "Glomeromycota") %>%
    select(otu_id, phylum, family, genus),
  non_glom_classification %>%
    select(otu_id, phylum, family, genus)
)

# Add the classification to the OTU tables
otu_raw_families <- otu_raw_families %>%
  left_join(family_classification_all, by = "otu_id") %>%
  select(otu_id, phylum, family, genus, everything())

otu_srs_families <- otu_srs_families %>%
  left_join(family_classification_all, by = "otu_id") %>%
  select(otu_id, phylum, family, genus, everything())

## 8.5 Shannon diversity and evenness ------------------------------------------

# Calculate diversity metrics for glomeromycota OTU table
glomeromycota_alpha_diversity <- otus_srs_AM %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  group_by(sample_id) %>%
  summarise(
    abundance_am_srs = sum(abundance),
    richness = n_distinct(otu_id[abundance > 0]),
    shannon_div = shannon_diversity(abundance)
  ) %>%
  # Calculate evenness (Pielou's J')
  mutate(
    evenness = case_when(
      richness == 0 ~ 0,           # No OTUs = no evenness
      richness == 1 ~ 0,           # Only 1 OTU = no evenness (complete dominance)
      shannon_div == 0 ~ 0,        # Shannon = 0 means complete dominance
      TRUE ~ shannon_div / log(richness)
    )
  ) %>%
  select(sample_id, shannon_div, evenness)

## 8.6 Phylogenetic alpha diversity --------------------------------------------

# Convert to samples x OTUs with sample as rownames
otu_matrix <- otus_srs_AM %>%
  column_to_rownames("otu_id") %>%
  t()

# Calculate Faith's PD
pd_result <- pd(
  samp = otu_matrix,
  tree = rooted_tree,
  include.root = TRUE
)

# 9. GENERATE FINAL METADATA FILE ----------------------------------------------

## 9.1 Compile final metadata --------------------------------------------------

diversity_estimate <- otu_table_library_filtered %>%
  filter(otu_id %in% am_otu_ids) %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  group_by(sample_id) %>%
  summarise(
    richness_observed = n_distinct(otu_id)
    ) %>%
  ungroup() %>%
  # Add absolute sample abundance
  inner_join(
    sample_id_depth %>%
      select(sample_id, total_sample_abundance = n_seqs),
    by = "sample_id"
  ) %>%
  inner_join(
    rare_am_richness %>%
      select(sample_id, richness_rare = median_richness),
    by = "sample_id"
  ) %>% 
  inner_join(
    srs_am_richness %>%
      rename(richness_srs = n_otus),
    by = "sample_id"
  ) %>%
  inner_join(
    glomeromycota_alpha_diversity,
    by = "sample_id"
  ) %>%
  inner_join(
    pd_result %>%
      rownames_to_column("sample_id") %>%
      select(sample_id, PD),
    by = "sample_id"
  ) %>%
  left_join(
    ordination_results,
    by = "sample_id"
  ) %>%
  select(
    sample_id,
    total_sample_abundance,
    richness_observed,
    richness_rare,
    richness_srs,
    evenness,
    shannon_div,
    faith_phylo_div = PD,
    pco1_identity_58,
    pco2_identity_p28,
    pco1_phylogeny_47,
    pco2_phylogeny_36
  ) %>%
  arrange(richness_rare) %>%
  print(n = 20)

## 9.4 Save final metadata -----------------------------------------------------

fwrite(diversity_estimate, "./output/diversity_estimate.txt", sep = "\t")

# 10. SAVE FILTERED OTU TABLES ----------------------------------------------------

# Save SRS-normalised AM OTU table
fwrite(otus_srs_AM, "./output/otu_table_srs.txt", sep = "\t")

# Save raw family-level OTU table
fwrite(otu_raw_families, "./output/otu_table_families.txt", sep = "\t")

# Save SRS-normalised family-level OTU table
fwrite(otu_srs_families, "./output/otu_table_families_srs.txt", sep = "\t")

message("---- Pipeline complete ----")
message("All output files saved to ./output/")
message("Distance matrices saved to ./data/distances/")

# Move the output contents tha am statistics directory:
# Create an "am" directory in the main data folder
dir.create("../../data/am", showWarnings = FALSE)

# Move all output files to the new "am" directory
file.copy(from = list.files("./output/", full.names = TRUE),
          to = "../../data/am/",
          recursive = TRUE)
