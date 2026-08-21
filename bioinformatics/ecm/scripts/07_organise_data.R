# ECM pipeline: organise bioinformatic outputs for downstream statistical analyses
# Run from the ecm project root directory

# 1. SETUP -----------------------------------------------------------------------

require(data.table)
require(SRS)
require(hillR)
require(ape)
require(phangorn)
require(DECIPHER)
require(Biostrings)
require(tidyverse)

# Index-switching filter: zero out reads below threshold% of each OTU's library total
filter_library <- function(otu_table, threshold = 0.1) {
  otu_table %>%
    pivot_longer(-otu_id, names_to = "sample", values_to = "reads") %>%
    group_by(otu_id) %>%
    mutate(reads = replace(reads, (reads / sum(reads)) * 100 <= threshold, 0)) %>%
    ungroup() %>%
    pivot_wider(names_from = "sample", values_from = "reads", values_fill = 0)
}

# Low abundance filter: zero out reads below threshold% of sample total OR below min_count (default 1 read)
filter_samples <- function(otu_tab, threshold = 0.1, min_count = 1) {
  otu_tab %>%
    pivot_longer(-otu_id, names_to = 'sample') %>%
    filter(value > 0) %>%
    group_by(sample) %>%
    mutate(rel_abund = value / sum(value) * 100) %>%
    ungroup() %>%
    mutate(value = if_else(rel_abund <= threshold | value <= min_count, 0, value)) %>%
    select(-rel_abund) %>%
    filter(value > 0) %>%
    pivot_wider(names_from = "sample", values_from = "value", values_fill = 0)
}

# Parse a gemelli ordination file; returns site coords with named PCoA columns
parse_ordination <- function(path, pco1_prefix, pco2_prefix) {
  lines      <- readLines(path)
  prop_line  <- which(grepl("Proportion explained", lines))[1]
  props      <- as.numeric(strsplit(lines[prop_line + 1], "\t")[[1]])
  prop1      <- round(props[1] * 100, 0)
  prop2      <- round(props[2] * 100, 0)
  site_start <- which(grepl("^Site", lines))[1] + 1
  site_end   <- length(lines)
  for (i in site_start:length(lines)) {
    if (lines[i] == "" || grepl("^Biplot|^Site constraints", lines[i])) {
      site_end <- i - 1; break
    }
  }
  lines[site_start:site_end] %>%
    keep(~ nchar(.) > 0 && !grepl("^\t", .)) %>%
    map_dfr(~ {
      p <- strsplit(., "\t")[[1]]
      if (length(p) >= 3) tibble(sample_id = p[1], PC1 = as.numeric(p[2]), PC2 = as.numeric(p[3]))
    }) %>%
    filter(!sample_id %in% c("Biplot", "Site constraints")) %>%
    rename(!!paste0(pco1_prefix, "_", prop1) := PC1,
           !!paste0(pco2_prefix, "_", prop2) := PC2)
}

# 2. DATA IMPORT -----------------------------------------------------------------

otu_classification <- fread("./tmp_clusters/otu_ecm_classification.txt")
otu_seqs           <- readDNAStringSet("./tmp_clusters/otu_ecm_sequences.fasta")
otu_table_ecm      <- fread("./tmp_clusters/otu_ecm_table.txt")

ecm_lineage <- fread("./data/ref_seqs/ecm_genera.txt") %>%
  select(genus, lineage, exploration_type)

non_ecm_ids <- fread("./tmp_clusters/non_ecm_clusters.txt") %>%
  mutate(otu_id = gsub(";size=.*", "", otu_id)) %>%
  pull(otu_id)

otu_table_all <- bind_rows(
  otu_table_ecm %>%
    pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance"),
  fread("./data/asv_table.txt") %>%
    rename(otu_id = OTU_ID) %>%
    filter(otu_id %in% non_ecm_ids) %>%
    pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance")
) %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0)

# 3. QUALITY CONTROL -------------------------------------------------------------

otu_table_library_filtered <- filter_library(otu_table_all, threshold = 0.1) %>%
  filter_samples(., threshold = 0.05, min_count = 1)

message("---- Library-wide and sample-specific filter ----")
message("Total reads before filtering: ", sum(otu_table_all[, -1]))
message("Total reads after filtering:  ", sum(otu_table_library_filtered[, -1]))
message("Percentage of reads removed:  ",
        round(((sum(otu_table_all[, -1]) - sum(otu_table_library_filtered[, -1])) /
                 sum(otu_table_all[, -1])) * 100, 3), "%")
message("Total OTUs before filtering: ", nrow(otu_table_all))
message("Total OTUs after filtering:  ", nrow(otu_table_library_filtered))
message("Percentage of OTUs removed:  ",
        round(((nrow(otu_table_all) - nrow(otu_table_library_filtered)) /
                 nrow(otu_table_all)) * 100, 3), "%")

# Check taxanomic composition of the "moc1" and "moc2" samples
moc_ids <- otu_table_library_filtered %>%
  select(otu_id, moc1, moc2) %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id", 
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  group_by(otu_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  arrange(desc(abundance)) %>%
  # Add genus information
  left_join(
    otu_classification %>%
      select(otu_id, family, genus, species),
    by = "otu_id"
  ) %>%
  print(n = Inf) %>%
  pull(otu_id)

# Compute the proprtion of reads in each sample that are from moc OTUs
otu_table_library_filtered %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id", 
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  # Compute total sample abundance for each sample
  group_by(sample_id) %>%
  mutate(
    sample_abundance = sum(abundance),
    rel_abundance = (abundance / sample_abundance) * 100
  ) %>%
  ungroup() %>%
  filter(otu_id %in% moc_ids) %>%
  arrange(otu_id, desc(abundance)) %>%
  print(n = Inf)

# Calculate sequencing depth
sample_id_depth <- otu_table_library_filtered %>%
  pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  group_by(sample_id) %>%
  summarise(n_seqs = sum(abundance), .groups = "drop") %>%
  arrange(n_seqs) %>%
  print()

min_depth <- 3600

low_abundance_sample_ids <- sample_id_depth %>%
  filter(n_seqs < min_depth | grepl("moc|neg", sample_id)) %>%
  pull(sample_id)

# 4. PREPARE ECM DATA ------------------------------------------------------------

dir.create("./output", showWarnings = FALSE)

# ECM sample IDs
ecm_sample_ids <- otu_table_library_filtered %>%
  select(-otu_id, -any_of(low_abundance_sample_ids)) %>%
  names()

# ECM OTU IDs and abundances
ecm_otu_ids <- otu_table_library_filtered %>%
  select(-any_of(low_abundance_sample_ids)) %>%
  filter(otu_id %in% otu_table_ecm$otu_id) %>%
  pivot_longer(
    cols = -otu_id, 
    names_to = "sample_id", 
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  group_by(otu_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# SRS normalisation
otus_srs <- SRS(
  data     = otu_table_library_filtered %>%
                select(-any_of(low_abundance_sample_ids)) %>%
                column_to_rownames("otu_id"),
  Cmin     = min_depth,
  set_seed = TRUE,
  seed     = 1986
) %>%
  as_tibble() %>%
  bind_cols(tibble(otu_id = otu_table_library_filtered$otu_id), .) %>%
  filter(rowSums(select(., -otu_id)) > 0)

srs_sample_ids <- setdiff(names(otus_srs), "otu_id")

# Raw ECM OTU table filtered to SRS-surviving samples (ordination input)
otu_table_ecm_raw <- otu_table_library_filtered %>%
  filter(otu_id %in% ecm_otu_ids$otu_id) %>%
  select(otu_id, any_of(srs_sample_ids)) %>%
  filter(rowSums(select(., -otu_id)) > 0)

fwrite(otu_table_ecm_raw, "./output/otu_table.txt", sep = "\t")

# SRS ECM OTU table
otu_table_ecm_srs <- otus_srs %>%
  filter(otu_id %in% ecm_otu_ids$otu_id) %>%
  select(otu_id, any_of(srs_sample_ids)) %>%
  filter(rowSums(select(., -otu_id)) > 0)

fwrite(otu_table_ecm_srs, "./output/otu_table_srs.txt", sep = "\t")

# ECM classification
otu_classification %>%
  select(-abundance) %>%
  inner_join(
    ecm_otu_ids, by = "otu_id"
  ) %>%
  arrange(desc(abundance)) %>%
  fwrite("./output/classification.txt", sep = "\t")

# ECM sequences
otu_seqs_ecm <- otu_seqs[names(otu_seqs) %in% ecm_otu_ids$otu_id]
writeXStringSet(otu_seqs_ecm, "./output/sequences.fasta")

# 5. PHYLOGENETIC TREE -----------------------------------------------------------

aligned_seqs <- readDNAStringSet("output/sequences.fasta") %>%
  AlignSeqs(verbose = TRUE, processors = 10)

rooted_tree <- as.phyDat(as.matrix(aligned_seqs), type = "DNA") %>%
  dist.ml() %>%
  NJ() %>%
  midpoint()

ape::write.tree(rooted_tree, file = "output/tree.newick")

# 6. ORDINATION (Aitchison RPCA via gemelli) -------------------------------------

# Note: --min-feature-frequency =  Minimum percentage of samples a feature must appear with a value greater than zero.

system("mkdir -p data/distances")
system('conda run -n gemelli_env biom convert -i output/otu_table.txt -o data/distances/otus.biom --table-type="OTU table" --to-json')

system("conda run -n gemelli_env gemelli rpca --in-biom data/distances/otus.biom --output-dir data/distances/ --min-feature-frequency 5")
system("mv data/distances/ordination.txt data/distances/ordination_identity.txt")
system("mv data/distances/distance-matrix.tsv output/dist_identity.txt")

system("conda run -n gemelli_env gemelli phylogenetic-rpca --in-biom data/distances/otus.biom --in-phylogeny output/tree.newick --output-dir data/distances/ --min-feature-frequency 5")
system("mv data/distances/ordination.txt data/distances/ordination_phylogeny.txt")
system("mv data/distances/distance-matrix.tsv output/dist_phylogeny.txt")

ordination_results <-
  parse_ordination("data/distances/ordination_identity.txt",  "pco1_tax", "pco2_tax") %>%
  left_join(
    parse_ordination("data/distances/ordination_phylogeny.txt", "pco1_phy", "pco2_phy"),
    by = "sample_id"
  )

fwrite(ordination_results, "data/distances/ordination_combined.txt", sep = "\t")

# 7. OTU TABLES FOR DIFFERENTIAL ABUNDANCE ANALYSIS --------------------------------------

# ECM SRS-normalised table with lineage and exploration type information
otu_table_ecm_srs %>%
  inner_join(
    otu_classification %>% select(otu_id, genus), 
    by = "otu_id"
  ) %>%
  pivot_longer(
    cols = -c(otu_id, genus), 
    names_to = "sample_id", 
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  group_by(genus, sample_id) %>%
  summarise(
    abundance = sum(abundance),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = abundance,
    values_fill = 0
  ) %>%
  inner_join(
    ecm_lineage,
    by = "genus"
  ) %>%
  select(genus, lineage, exploration_type, everything()) %>%
  fwrite("./output/otu_table_srs_ecm_genus.txt", sep = "\t")

# 8. DIVERSITY METRICS -----------------------------------------------------------

# Sample abundance and richness (SRS-normalised ECM)
ecm_alpha_diversity <- otu_table_ecm_srs %>%
  pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance_srs = sum(abundance),
    richness_srs  = n_distinct(otu_id),
    .groups = "drop"
  )

# Hill taxonomic and phylogenetic diversity + evenness (SRS table)
comm_matrix_srs_ECM <- otu_table_ecm_srs %>% column_to_rownames("otu_id") %>% t()
srs_ecm_tree        <- keep.tip(rooted_tree, tip = colnames(comm_matrix_srs_ECM))

hill_diversity <- list(
  hill_taxa(comm_matrix_srs_ECM, q = 0)                       %>% enframe("sample_id", "hill_tax_div_q0"),
  hill_taxa(comm_matrix_srs_ECM, q = 1)                       %>% enframe("sample_id", "hill_tax_div_q1"),
  hill_taxa(comm_matrix_srs_ECM, q = 2)                       %>% enframe("sample_id", "hill_tax_div_q2"),
  hill_phylo(comm_matrix_srs_ECM, tree = srs_ecm_tree, q = 0) %>% enframe("sample_id", "hill_phy_div_q0"),
  hill_phylo(comm_matrix_srs_ECM, tree = srs_ecm_tree, q = 1) %>% enframe("sample_id", "hill_phy_div_q1"),
  hill_phylo(comm_matrix_srs_ECM, tree = srs_ecm_tree, q = 2) %>% enframe("sample_id", "hill_phy_div_q2")
) %>%
  reduce(left_join, by = "sample_id") %>%
  mutate(
    hill_tax_even_q1 = hill_tax_div_q1 / hill_tax_div_q0,
    hill_tax_even_q2 = hill_tax_div_q2 / hill_tax_div_q0,
    hill_phy_even_q1 = hill_phy_div_q1 / hill_phy_div_q0,
    hill_phy_even_q2 = hill_phy_div_q2 / hill_phy_div_q0
  )

diversity_estimate <- ecm_alpha_diversity %>%
  inner_join(sample_id_depth %>% select(sample_id, total_sample_abundance = n_seqs),
             by = "sample_id") %>%
  left_join(hill_diversity,     by = "sample_id") %>%
  left_join(ordination_results, by = "sample_id") %>%
  select(
    sample_id, total_sample_abundance, abundance_srs, richness_srs,
    hill_tax_div_q0, hill_tax_div_q1, hill_tax_div_q2,
    hill_tax_even_q1, hill_tax_even_q2,
    hill_phy_div_q0, hill_phy_div_q1, hill_phy_div_q2,
    hill_phy_even_q1, hill_phy_even_q2,
    starts_with("pco1_tax_"), starts_with("pco2_tax_"),
    starts_with("pco1_phy_"), starts_with("pco2_phy_")
  )

fwrite(diversity_estimate, "./output/diversity_ecm.txt", sep = "\t")

# 9. LINEAGE-SPECIFIC DIVERSITY --------------------------------------------------

# Compute the same metrics as diversity_ecm but restricted to OTUs of one lineage.
# All ECM samples are retained; samples absent from the lineage receive 0s.
lineage_diversity <- function(lin_otu_ids, all_sample_ids, otu_srs, tree, sample_depth) {

  # Abundance and richness for samples with >0 reads in this lineage
  abund_rich <- otu_srs %>%
    filter(otu_id %in% lin_otu_ids) %>%
    pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance") %>%
    filter(abundance > 0) %>%
    group_by(sample_id) %>%
    summarise(abundance_srs = sum(abundance), richness_srs = n_distinct(otu_id), .groups = "drop")

  present_samples <- abund_rich$sample_id

  # Hill diversity — computed only on present samples, then zeros are joined back
  if (length(lin_otu_ids) > 0 && length(present_samples) > 0) {
    comm <- otu_srs %>%
      filter(otu_id %in% lin_otu_ids) %>%
      select(otu_id, all_of(present_samples)) %>%
      column_to_rownames("otu_id") %>%
      t()

    tax_hills <- list(
      hill_taxa(comm, q = 0) %>% enframe("sample_id", "hill_tax_div_q0"),
      hill_taxa(comm, q = 1) %>% enframe("sample_id", "hill_tax_div_q1"),
      hill_taxa(comm, q = 2) %>% enframe("sample_id", "hill_tax_div_q2")
    ) %>% reduce(left_join, by = "sample_id")

    # Phylogenetic Hill requires >=2 tree tips
    tree_tips <- intersect(lin_otu_ids, colnames(comm))
    if (length(tree_tips) >= 2) {
      lin_tree  <- keep.tip(tree, tip = tree_tips)
      phy_hills <- list(
        hill_phylo(comm, tree = lin_tree, q = 0) %>% enframe("sample_id", "hill_phy_div_q0"),
        hill_phylo(comm, tree = lin_tree, q = 1) %>% enframe("sample_id", "hill_phy_div_q1"),
        hill_phylo(comm, tree = lin_tree, q = 2) %>% enframe("sample_id", "hill_phy_div_q2")
      ) %>% reduce(left_join, by = "sample_id")
      hill_res <- left_join(tax_hills, phy_hills, by = "sample_id")
    } else {
      hill_res <- mutate(tax_hills,
        hill_phy_div_q0 = NA_real_, hill_phy_div_q1 = NA_real_, hill_phy_div_q2 = NA_real_)
    }
  } else {
    hill_res <- tibble(
      sample_id       = character(),
      hill_tax_div_q0 = numeric(), hill_tax_div_q1 = numeric(), hill_tax_div_q2 = numeric(),
      hill_phy_div_q0 = numeric(), hill_phy_div_q1 = numeric(), hill_phy_div_q2 = numeric()
    )
  }

  tibble(sample_id = all_sample_ids) %>%
    left_join(sample_depth %>% select(sample_id, total_sample_abundance = n_seqs), by = "sample_id") %>%
    left_join(abund_rich, by = "sample_id") %>%
    left_join(hill_res,   by = "sample_id") %>%
    mutate(
      abundance_srs = replace_na(abundance_srs, 0),
      richness_srs  = replace_na(richness_srs,  0),
      across(starts_with("hill_"), ~ replace_na(., 0)),
      hill_tax_even_q1 = ifelse(hill_tax_div_q0 > 0, hill_tax_div_q1 / hill_tax_div_q0, 0),
      hill_tax_even_q2 = ifelse(hill_tax_div_q0 > 0, hill_tax_div_q2 / hill_tax_div_q0, 0),
      hill_phy_even_q1 = ifelse(hill_phy_div_q0 > 0, hill_phy_div_q1 / hill_phy_div_q0, 0),
      hill_phy_even_q2 = ifelse(hill_phy_div_q0 > 0, hill_phy_div_q2 / hill_phy_div_q0, 0)
    ) %>%
    select(
      sample_id, total_sample_abundance, abundance_srs, richness_srs,
      hill_tax_div_q0, hill_tax_div_q1, hill_tax_div_q2,
      hill_tax_even_q1, hill_tax_even_q2,
      hill_phy_div_q0, hill_phy_div_q1, hill_phy_div_q2,
      hill_phy_even_q1, hill_phy_even_q2
    )
}

# OTU-to-lineage lookup (ECM OTUs only)
otu_lineage_map <- otu_classification %>%
  select(otu_id, genus) %>%
  inner_join(ecm_lineage %>% select(genus, lineage), by = "genus") %>%
  filter(otu_id %in% ecm_otu_ids$otu_id)

# --- Smoke test ---------------------------------------------------------------
message("otu_lineage_map rows: ", nrow(otu_lineage_map))
message("Unique lineages in otu_lineage_map:")
print(sort(unique(otu_lineage_map$lineage)))

message("Genera in otu_classification not matched in ecm_lineage:")
print(setdiff(
  otu_classification %>% filter(otu_id %in% ecm_otu_ids$otu_id) %>% pull(genus) %>% unique(),
  ecm_lineage$genus
))

target_lineages <- c("/russula-lactarius", "/cortinarius", "/tomentella-thelephora",
                     "/laccaria", "/inocybe", "/sebacina")

message("OTUs per target lineage:")
for (lin in target_lineages) {
  n <- otu_lineage_map %>% filter(lineage == lin) %>% nrow()
  message("  ", lin, ": ", n, " OTUs")
}
# --- End smoke test -----------------------------------------------------------

for (lin in target_lineages) {
  lin_otu_ids <- otu_lineage_map %>% filter(lineage == lin) %>% pull(otu_id) %>% unique()
  lineage_diversity(lin_otu_ids, diversity_estimate$sample_id,
                    otu_table_ecm_srs, rooted_tree, sample_id_depth) %>%
    fwrite(paste0("./output/diversity_", sub("^/", "", lin), ".txt"), sep = "\t")
}

# 10. COPY OUTPUTS TO DATA FOLDER -------------------------------------------------

dir.create("../../data/emf", showWarnings = FALSE)
file.copy(from = list.files("./output/", full.names = TRUE),
          to   = "../../data/emf/",
          recursive = TRUE)

# Rename the guild-specific outputs in data/emf/ to the EMF terminology used
# downstream (alpha-diversity scripts, TITAN) -- this renames the copies in
# data/emf/ only, not the pipeline-local files under ./output/ (or, for
# otu_classification_non_ecm.txt / otu_table_ecm_genus.txt, the upstream
# utils/cluster_otus.R outputs they're copied from), so re-running this script
# is idempotent with respect to that local naming. Files with no guild-specific
# token (classification.txt, otu_table.txt, otu_table_srs.txt, dist_identity.txt,
# dist_phylogeny.txt, sequences.fasta, tree.newick, otu_table_all_genus.txt,
# and the lineage-specific diversity_*.txt files) are left as-is.
file.rename("../../data/emf/diversity_ecm.txt", "../../data/emf/diversity_emf.txt")

file.rename("../../data/emf/otu_table_srs_ecm_genus.txt", "../../data/emf/otu_table_srs_emf_genus.txt")
file.rename("../../data/emf/otu_table_ecm_genus.txt",      "../../data/emf/otu_table_emf_genus.txt")

# From utils/cluster_otus.R, swept into ./output/ before this script runs
file.rename("../../data/emf/otu_classification_non_ecm.txt", "../../data/emf/otu_classification_non_emf.txt")

message("---- Pipeline complete ----")
