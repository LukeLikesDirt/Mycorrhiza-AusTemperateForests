
# Required packages and functions
require(ape)
require(picante)
require(data.table)
require(tidyverse)
source("./code/helper_functions/abundance_filters.R")

# Helper function for computing Shannon diversity
shannon_diversity <- function(x) {
  p <- x / sum(x)
  -sum(p * log(p))
}

# Required files
asv_table <- fread("./data/ASVs.txt")         # ASV table abundance for all ASVs
otu_table <- fread("./output/otu_table.txt")    # Final abundance filtered and clustered OTU table
phylo_tree <- read.tree("./output/tree.newick")    # Phylogenetic tree of clustered OTUs
community_composition <- fread("./data/distances/ordination_combined.txt")  # Community composition file

# Compute abundance filtered ASV table 
asv_table_filtered <- filter_library(
  fread("./data/ASVs.txt") %>% rename(otu_id = OTU_ID),
  0.1
)

# Total read abundance per sample ----------------------------------------------
total_reads_per_sample <- asv_table_filtered %>%
  # Remove samples not in the OTU sample matrix
  select(otu_id, names(otu_table)[-1]) %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  group_by(sample_id) %>%
  summarise(total_sample_abundance = sum(abundance), .groups = "drop")

# Alpha diversity metrics ------------------------------------------------------

# Glomeromycota richness, evenness and Shannon diversity per sample
glomeromycota_reads_per_sample <- otu_table %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  group_by(sample_id) %>%
  summarise(
    richness = n_distinct(otu_id[abundance > 0]),
    shannon_div = shannon_diversity(abundance)
  ) %>%
  # Calculate evenness (Pielou's J')
  mutate(
    evenness = if_else(
      richness > 1,
      shannon_div / log(richness),
      1  # If only 1 OTU, evenness = 1 (perfectly even by definition)
    )
  )

# Phylogenetic diversity metrics -----------------------------------------------

# Convert to samples x OTUs with sample as rowsnames
otu_matrix <- otu_table %>%
  column_to_rownames("otu_id") %>%
  t()

# Calculate Faith's PD
pd_result <- pd(
  samp = otu_matrix,
  tree = phylo_tree,
  include.root = TRUE
)

# Combine all diversity metrics into metadata ---------------------------------

# Combine metadata
metadata <- left_join(
  total_reads_per_sample, 
  glomeromycota_reads_per_sample,
  by = "sample_id"
  ) %>%
  left_join(
    pd_result %>%
      rownames_to_column("sample_id") %>%
      select(sample_id, PD),
    by = "sample_id"
  ) %>%
  rename(faith_pylo_div = PD) %>%
  inner_join(
    community_composition,
    by = "sample_id"
  )

# Visualise abundance against diversity measures -------------------------------

# Ranges of metrics
range(metadata$total_sample_abundance)
range(metadata$richness)
range(metadata$evenness)
range(metadata$shannon_div)
range(metadata$faith_pylo_div)

# Richness vs Tttal read abundance
ggplot(metadata, aes(x = total_sample_abundance, y = richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, ceiling(max(metadata$richness)))) +
  labs(
    x = "Total read abundance",
    y = "Richness"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1)

# Evenness vs total read abundance
ggplot(metadata, aes(x = total_sample_abundance, y = evenness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Total read abundance",
    y = "Evenness"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1)

# Shannon diversity vs total read abundance
ggplot(metadata, aes(x = total_sample_abundance, y = shannon_div)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, ceiling(max(metadata$shannon_div)))) +
  labs(
    x = "Total read abundance",
    y = "Shannon diversity"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1)


# Faith's phylogenetic diversity vs total read abundance
ggplot(metadata, aes(x = total_sample_abundance, y = faith_pylo_div)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, ceiling(max(metadata$faith_pylo_div)))) +
  labs(
    x = "Total read abundance",
    y = "Faith's phylogenetic diversity"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1)

# Save data to the statistics data folder --------------------------------------

# Create an "am" directory in the main data folder
dir.create("../../data/am", showWarnings = FALSE)

# Save metadata file
metadata %>%
  fwrite("../../data/am/diversity_metrics.txt", sep = "\t")

# Move contents from the distance folder to the am folder:
# All files with "dist_"
file.copy(
  from = list.files(
    path = "./output",
    pattern = "dist_",
    full.names = TRUE
  ),
  to = "../../data/am/",
  overwrite = TRUE
)
