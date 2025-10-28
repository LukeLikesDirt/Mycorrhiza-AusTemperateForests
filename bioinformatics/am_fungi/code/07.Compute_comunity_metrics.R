
# Required packages and functions
require(phangorn)
require(DECIPHER)
require(Biostrings)
require(data.table)
require(tidyverse)

# Generate phylogenetic tree ---------------------------------------------------

# Convert to DNAStringSet for alignment
dna_seqs <- readDNAStringSet("output/otu_sequences.fasta")

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

# Compute Aitchinson distance --------------------------------------------------

# Create a directory for gemelli output
system("mkdir -p data/distances")

# Convert to biom
system('conda run -n gemelli_env biom convert -i output/otu_table.txt -o data/distances/otus.biom --table-type="OTU table" --to-json')

#### Run gemeilli based on identity ####

# Run gemelli for all taxa
system("conda run -n gemelli_env gemelli rpca --in-biom data/distances/otus.biom --output-dir data/distances/")

# Rename "ordination.txt" and "distance-matrix.txt" to avoid overwriting
system("mv data/distances/ordination.txt data/distances/ordination_identity_all.txt")
system("mv data/distances/distance-matrix.tsv output/dist_identity_all.txt")

# Run gemelli for taxa that occur in at least 2 samples: Gemelli works on 
# percentage presence, so with 122 samples min frequency of 1 will require
# presence in at least 2 samples
system("conda run -n gemelli_env gemelli rpca --in-biom data/otus.biom --output-dir data/distances/ --min-feature-frequency 1")

# Rename "ordination.txt" and "distance-matrix.txt" to avoid overwriting
system("mv data/distances/ordination.txt data/distances/ordination_identity_prev.txt")
system("mv data/distances/distance-matrix.tsv output/dist_identity_prev.txt")

# List the output files
system("ls data/distances")

#### * All: ordination results * ####

# Read the ordination file
ordination_lines <- readLines("data/distances/ordination_identity_all.txt")

# Find the line indices for different sections
proportion_line <- which(grepl("Proportion explained", ordination_lines))[1]  # Get first match
site_line <- which(grepl("^Site", ordination_lines))[1]  # Get first match, use ^ for start of line

# Extract the proportion explained values
proportions <- as.numeric(strsplit(ordination_lines[proportion_line + 1], "\t")[[1]])

# Convert to percentages and round
prop1 <- round(proportions[1] * 100, 0)
prop2 <- round(proportions[2] * 100, 0)

# Read the site coordinates (starting from the line after "Site")
site_start <- site_line + 1
site_end <- length(ordination_lines)

# Find where site data ends (look for empty lines or "Biplot"/"Site constraints")
for (i in site_start:length(ordination_lines)) {
  if (ordination_lines[i] == "" || grepl("^Biplot", ordination_lines[i]) || grepl("^Site constraints", ordination_lines[i])) {
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
  if (line != "" && !grepl("^\t", line)) {  # Skip empty lines and lines starting with tab
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


# Create the final data frame with your specified column names
ordination_results <- site_coords %>%
  rename(
    !!paste0("pco1_identity_all_", prop1) := PC1,
    !!paste0("pco2_identity_all_", prop2) := PC2
  ) %>%
  select(sample_id, everything())

#### * Prevelance filtered: ordination results * ####

# Read the ordination file
ordination_lines <- readLines("data/distances/ordination_identity_prev.txt")

# Find the line indices for different sections
proportion_line <- which(grepl("Proportion explained", ordination_lines))[1]  # Get first match
site_line <- which(grepl("^Site", ordination_lines))[1]  # Get first match, use ^ for start of line

# Extract the proportion explained values
proportions <- as.numeric(strsplit(ordination_lines[proportion_line + 1], "\t")[[1]])

# Convert to percentages and round
prop1 <- round(proportions[1] * 100, 0)
prop2 <- round(proportions[2] * 100, 0)

# Read the site coordinates (starting from the line after "Site")
site_start <- site_line + 1
site_end <- length(ordination_lines)

# Find where site data ends (look for empty lines or "Biplot"/"Site constraints")
for (i in site_start:length(ordination_lines)) {
  if (ordination_lines[i] == "" || grepl("^Biplot", ordination_lines[i]) || grepl("^Site constraints", ordination_lines[i])) {
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
  if (line != "" && !grepl("^\t", line)) {  # Skip empty lines and lines starting with tab
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

# Add to the final data frame with your specified column names
ordination_results <- ordination_results %>%
  left_join(
    site_coords %>%
      rename(
        !!paste0("pco1_identity_prev_", prop1) := PC1,
        !!paste0("pco2_identity_prev_", prop2) := PC2
      ) %>%
      select(sample_id, everything()),
    by = "sample_id"
  )

# Regress the pco1 against each other
ordination_results %>%
  ggplot(aes(x = pco2_identity_all_19,
             y = pco2_identity_prev_22)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = paste0("PCoA1 Identity All (", round(proportions[1] * 100, 0), "%)"),
    y = paste0("PCoA1 Identity Prevalence Filtered (", round(proportions[1] * 100, 0), "%)")
  ) +
  theme_minimal()

#### Run gemeilli based on phylogeny ####

# Run gemelli for fungi first and then rename
system("conda run conda run -n gemelli_env gemelli phylogenetic-rpca --in-biom data/otus.biom --in-phylogeny output/tree.newick --output-dir data/distances/")

# List the output files
system("ls data/distances")

# Rename "ordination.txt" and "distance-matrix.txt" to avoid overwriting
system("mv data/distances/ordination.txt data/distances/ordination_phylogeny_all.txt")
system("mv data/distances/distance-matrix.tsv output/dist_phylogeny_all.txt")

# List the output files
system("ls data/distances")

# Run gemelli for prevalence filtered fungi
system("conda run -n gemelli_env gemelli phylogenetic-rpca --in-biom data/otus.biom --in-phylogeny output/tree.newick --output-dir data/distances/ --min-feature-frequency 1")

# Rename "ordination.txt" and "distance-matrix.txt" to avoid overwriting
system("mv data/distances/ordination.txt data/distances/ordination_phylogeny_prev.txt")
system("mv data/distances/distance-matrix.tsv output/dist_phylogeny_prev.txt")

#### * All: ordination results * ####

# Read the ordination file
ordination_lines <- readLines("data/distances/ordination_phylogeny_all.txt")

# Find the line indices for different sections
proportion_line <- which(grepl("Proportion explained", ordination_lines))[1]  # Get first match
site_line <- which(grepl("^Site", ordination_lines))[1]  # Get first match, use ^ for start of line

# Extract the proportion explained values
proportions <- as.numeric(strsplit(ordination_lines[proportion_line + 1], "\t")[[1]])

# Convert to percentages and round
prop1 <- round(proportions[1] * 100, 0)
prop2 <- round(proportions[2] * 100, 0)

# Read the site coordinates (starting from the line after "Site")
site_start <- site_line + 1
site_end <- length(ordination_lines)

# Find where site data ends (look for empty lines or "Biplot"/"Site constraints")
for (i in site_start:length(ordination_lines)) {
  if (ordination_lines[i] == "" || grepl("^Biplot", ordination_lines[i]) || grepl("^Site constraints", ordination_lines[i])) {
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
  if (line != "" && !grepl("^\t", line)) {  # Skip empty lines and lines starting with tab
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

# Add to the final data frame with your specified column names
ordination_results <- ordination_results %>%
  left_join(
    site_coords %>%
      rename(
        !!paste0("pco1_phylogeny_all_", prop1) := PC1,
        !!paste0("pco2_phylogeny_all_", prop2) := PC2
      ) %>%
      select(sample_id, everything()),
    by = "sample_id"
  )

#### * Prevelance filtered: ordination results * ####

# Read the ordination file
ordination_lines <- readLines("data/distances/ordination_phylogeny_prev.txt")

# Find the line indices for different sections
proportion_line <- which(grepl("Proportion explained", ordination_lines))[1]  # Get first match
site_line <- which(grepl("^Site", ordination_lines))[1]  # Get first match, use ^ for start of line

# Extract the proportion explained values
proportions <- as.numeric(strsplit(ordination_lines[proportion_line + 1], "\t")[[1]])

# Convert to percentages and round
prop1 <- round(proportions[1] * 100, 0)
prop2 <- round(proportions[2] * 100, 0)

# Read the site coordinates (starting from the line after "Site")
site_start <- site_line + 1
site_end <- length(ordination_lines)

# Find where site data ends (look for empty lines or "Biplot"/"Site constraints")
for (i in site_start:length(ordination_lines)) {
  if (ordination_lines[i] == "" || grepl("^Biplot", ordination_lines[i]) || grepl("^Site constraints", ordination_lines[i])) {
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
  if (line != "" && !grepl("^\t", line)) {  # Skip empty lines and lines starting with tab
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

# Add to the final data frame with your specified column names
ordination_results <- ordination_results %>%
  left_join(
    site_coords %>%
      rename(
        !!paste0("pco1_phylogeny_prev_", prop1) := PC1,
        !!paste0("pco2_phylogeny_prev_", prop2) := PC2
      ) %>%
      select(sample_id, everything()),
    by = "sample_id"
  )

# Save the final ordination results --------------------------------------------

fwrite(ordination_results, "data/distances/ordination_combined.txt", sep = "\t")

# The end ----------------------------------------------------------------------
