# AM pipeline: organise bioinformatic outputs for downstream statistical analyses
# Run from the am project root directory

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

# Low abundance filter: zero out reads below threshold% of sample total OR below min_count
filter_samples <- function(otu_tab, threshold = 0.01, min_count = 1) {
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

# AM fungi evaluated here span two lineages: phylum Glomeromycota, and order
# Densosporales (class Endogonomycetes, Mucoromycota) — both are AM fungi. The
# otu_am_* files are the AM pool (see --am_orders in 02_classify_asvs.sh);
# otu_classification carries an `am_group` column (Glomeromycota /
# Endogonomycetes) that propagates to output/classification.txt.
otu_classification <- fread("./tmp_clusters/otu_am_classification.txt")
otu_seqs           <- readDNAStringSet("./tmp_clusters/otu_am_sequences.fasta")
otu_table_am       <- fread("./tmp_clusters/otu_am_table.txt")

non_am_ids <- fread("./tmp_clusters/non_am_clusters.txt") %>%
  mutate(otu_id = gsub(";size=.*", "", otu_id)) %>%
  pull(otu_id)

otu_table_all <- bind_rows(
  otu_table_am %>%
    pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance"),
  fread("./data/asv_table.txt") %>%
    rename(otu_id = OTU_ID) %>%
    filter(otu_id %in% non_am_ids) %>%
    pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance")
) %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0)

# 3. QUALITY CONTROL -------------------------------------------------------------

otu_table_library_filtered <- filter_library(otu_table_all, threshold = 0.1) %>%
  filter_samples(., threshold = 0.01, min_count = 1)

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

sample_id_depth <- otu_table_library_filtered %>%
  pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  group_by(sample_id) %>%
  summarise(n_seqs = sum(abundance), .groups = "drop") %>%
  arrange(n_seqs) %>%
  print()

min_depth <- 2500

low_abundance_sample_ids <- sample_id_depth %>%
  filter(n_seqs < min_depth | grepl("moc|neg", sample_id)) %>%
  pull(sample_id)

# 4. PREPARE AM DATA -------------------------------------------------------------

dir.create("./output", showWarnings = FALSE)

am_sample_ids <- otu_table_am %>%
  select(-otu_id, -any_of(low_abundance_sample_ids)) %>%
  names()

am_otu_ids <- otu_table_library_filtered %>%
  select(-any_of(low_abundance_sample_ids)) %>%
  filter(otu_id %in% otu_table_am$otu_id,
         rowSums(select(., -otu_id)) > 0) %>%
  pivot_longer(cols = -otu_id, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  group_by(otu_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

otu_classification %>%
  select(-abundance) %>%
  inner_join(am_otu_ids, by = "otu_id") %>%
  arrange(desc(abundance)) %>%
  fwrite("./output/classification.txt", sep = "\t")

otu_seqs_am <- otu_seqs[names(otu_seqs) %in% am_otu_ids$otu_id]
writeXStringSet(otu_seqs_am, "./output/sequences.fasta")

# SRS normalisation on all OTUs (run before ordination to identify surviving samples)
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

# Raw AM OTU table filtered to SRS-surviving samples (ordination input)
otu_table_am_raw <- otu_table_library_filtered %>%
  filter(otu_id %in% am_otu_ids$otu_id) %>%
  select(otu_id, any_of(srs_sample_ids)) %>%
  filter(rowSums(select(., -otu_id)) > 0)

fwrite(otu_table_am_raw, "./output/otu_table.txt", sep = "\t")

# SRS AM OTU table
otu_table_am_srs <- otus_srs %>%
  filter(otu_id %in% am_otu_ids$otu_id) %>%
  select(otu_id, any_of(srs_sample_ids)) %>%
  filter(rowSums(select(., -otu_id)) > 0)

fwrite(otu_table_am_srs, "./output/otu_table_srs.txt", sep = "\t")

# 5. PHYLOGENETIC TREE -----------------------------------------------------------

aligned_seqs <- readDNAStringSet("output/sequences.fasta") %>%
  AlignSeqs(verbose = TRUE, processors = 10)

rooted_tree <- as.phyDat(as.matrix(aligned_seqs), type = "DNA") %>%
  dist.ml() %>%
  NJ() %>%
  midpoint()

ape::write.tree(rooted_tree, file = "output/tree.newick")

# 6. ORDINATION (Aitchison RPCA via gemelli) -------------------------------------

system("mkdir -p data/distances")
system('conda run -n gemelli_env biom convert -i output/otu_table.txt -o data/distances/otus.biom --table-type="OTU table" --to-json')

system("conda run -n gemelli_env gemelli rpca --in-biom data/distances/otus.biom --output-dir data/distances/ --min-feature-frequency 1")
system("mv data/distances/ordination.txt data/distances/ordination_identity.txt")
system("mv data/distances/distance-matrix.tsv output/dist_identity.txt")

system("conda run -n gemelli_env gemelli phylogenetic-rpca --in-biom data/distances/otus.biom --in-phylogeny output/tree.newick --output-dir data/distances/ --min-feature-frequency 1")
system("mv data/distances/ordination.txt data/distances/ordination_phylogeny.txt")
system("mv data/distances/distance-matrix.tsv output/dist_phylogeny.txt")

ordination_results <-
  parse_ordination("data/distances/ordination_identity.txt",  "pco1_tax", "pco2_tax") %>%
  left_join(
    parse_ordination("data/distances/ordination_phylogeny.txt", "pco1_phy", "pco2_phy"),
    by = "sample_id"
  )

fwrite(ordination_results, "data/distances/ordination_combined.txt", sep = "\t")

# 7. DIVERSITY METRICS -----------------------------------------------------------

# Sample abundance and richness (SRS-normalised AM)
am_alpha_diversity <- otu_table_am_srs %>%
  pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  group_by(sample_id) %>%
  summarise(
    abundance_srs = sum(abundance),
    richness_srs  = n_distinct(otu_id),
    .groups = "drop"
  )

# Hill taxonomic and phylogenetic diversity + evenness (SRS table)
comm_matrix_srs_AM <- otu_table_am_srs %>% column_to_rownames("otu_id") %>% t()
srs_am_tree        <- keep.tip(rooted_tree, tip = colnames(comm_matrix_srs_AM))

hill_diversity <- list(
  hill_taxa(comm_matrix_srs_AM, q = 0)                      %>% enframe("sample_id", "hill_tax_div_q0"),
  hill_taxa(comm_matrix_srs_AM, q = 1)                      %>% enframe("sample_id", "hill_tax_div_q1"),
  hill_taxa(comm_matrix_srs_AM, q = 2)                      %>% enframe("sample_id", "hill_tax_div_q2"),
  hill_phylo(comm_matrix_srs_AM, tree = srs_am_tree, q = 0) %>% enframe("sample_id", "hill_phy_div_q0"),
  hill_phylo(comm_matrix_srs_AM, tree = srs_am_tree, q = 1) %>% enframe("sample_id", "hill_phy_div_q1"),
  hill_phylo(comm_matrix_srs_AM, tree = srs_am_tree, q = 2) %>% enframe("sample_id", "hill_phy_div_q2")
) %>%
  reduce(left_join, by = "sample_id") %>%
  mutate(
    hill_tax_even_q1 = hill_tax_div_q1 / hill_tax_div_q0,
    hill_tax_even_q2 = hill_tax_div_q2 / hill_tax_div_q0,
    hill_phy_even_q1 = hill_phy_div_q1 / hill_phy_div_q0,
    hill_phy_even_q2 = hill_phy_div_q2 / hill_phy_div_q0
  )

# 8. FINAL OUTPUTS ---------------------------------------------------------------

diversity_estimate <- am_alpha_diversity %>%
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

fwrite(diversity_estimate, "./output/diversity_am.txt", sep = "\t")

# 9. GROUP-SPECIFIC (GLOMEROMYCOTA / ENDOGONOMYCETES) DIVERSITY METRICS --------
#
# Same methodology as diversity_am.txt (SRS-normalised alpha diversity, Hill
# numbers, Aitchison RPCA ordination), computed separately per am_group so
# within-lineage diversity/ordination isn't diluted by the other lineage's
# presence/absence pattern. Reuses the SRS normalisation, raw AM table and
# tree already computed above (sections 4-6) rather than re-running QC/SRS.

compute_group_diversity <- function(group_label, file_suffix) {

  group_otu_ids <- otu_classification %>%
    filter(am_group == group_label) %>%
    pull(otu_id)

  group_otu_table_raw <- otu_table_am_raw %>% filter(otu_id %in% group_otu_ids)
  group_otu_table_srs <- otu_table_am_srs %>% filter(otu_id %in% group_otu_ids)

  raw_table_path <- paste0("output/otu_table_", file_suffix, ".txt")
  biom_path       <- paste0("data/distances/otus_", file_suffix, ".biom")
  fwrite(group_otu_table_raw, raw_table_path, sep = "\t")

  system(paste0(
    "conda run -n gemelli_env biom convert -i ", raw_table_path,
    " -o ", biom_path, ' --table-type="OTU table" --to-json'
  ))

  system(paste0(
    "conda run -n gemelli_env gemelli rpca --in-biom ", biom_path,
    " --output-dir data/distances/ --min-feature-frequency 1"
  ))
  system(paste0("mv data/distances/ordination.txt data/distances/ordination_identity_", file_suffix, ".txt"))
  system(paste0("mv data/distances/distance-matrix.tsv output/dist_identity_", file_suffix, ".txt"))

  comm_matrix_group <- group_otu_table_srs %>% column_to_rownames("otu_id") %>% t()
  # hill_phylo() needs a tree pruned to the SRS community's OTU set; gemelli's
  # phylogeny input must instead be a superset of the (raw-table-derived) biom's
  # OTUs, so it reuses the full combined tree (output/tree.newick) exactly like
  # the combined-AM ordination above does — never a group-pruned tree.
  group_tree <- keep.tip(rooted_tree, tip = colnames(comm_matrix_group))

  system(paste0(
    "conda run -n gemelli_env gemelli phylogenetic-rpca --in-biom ", biom_path,
    " --in-phylogeny output/tree.newick --output-dir data/distances/ --min-feature-frequency 1"
  ))
  system(paste0("mv data/distances/ordination.txt data/distances/ordination_phylogeny_", file_suffix, ".txt"))
  system(paste0("mv data/distances/distance-matrix.tsv output/dist_phylogeny_", file_suffix, ".txt"))

  group_ordination <-
    parse_ordination(paste0("data/distances/ordination_identity_", file_suffix, ".txt"), "pco1_tax", "pco2_tax") %>%
    left_join(
      parse_ordination(paste0("data/distances/ordination_phylogeny_", file_suffix, ".txt"), "pco1_phy", "pco2_phy"),
      by = "sample_id"
    )

  group_alpha_diversity <- group_otu_table_srs %>%
    pivot_longer(-otu_id, names_to = "sample_id", values_to = "abundance") %>%
    filter(abundance > 0) %>%
    group_by(sample_id) %>%
    summarise(abundance_srs = sum(abundance), richness_srs = n_distinct(otu_id), .groups = "drop")

  group_hill_diversity <- list(
    hill_taxa(comm_matrix_group, q = 0)                        %>% enframe("sample_id", "hill_tax_div_q0"),
    hill_taxa(comm_matrix_group, q = 1)                        %>% enframe("sample_id", "hill_tax_div_q1"),
    hill_taxa(comm_matrix_group, q = 2)                        %>% enframe("sample_id", "hill_tax_div_q2"),
    hill_phylo(comm_matrix_group, tree = group_tree, q = 0)    %>% enframe("sample_id", "hill_phy_div_q0"),
    hill_phylo(comm_matrix_group, tree = group_tree, q = 1)    %>% enframe("sample_id", "hill_phy_div_q1"),
    hill_phylo(comm_matrix_group, tree = group_tree, q = 2)    %>% enframe("sample_id", "hill_phy_div_q2")
  ) %>%
    reduce(left_join, by = "sample_id") %>%
    mutate(
      hill_tax_even_q1 = hill_tax_div_q1 / hill_tax_div_q0,
      hill_tax_even_q2 = hill_tax_div_q2 / hill_tax_div_q0,
      hill_phy_even_q1 = hill_phy_div_q1 / hill_phy_div_q0,
      hill_phy_even_q2 = hill_phy_div_q2 / hill_phy_div_q0
    )

  group_alpha_diversity %>%
    inner_join(sample_id_depth %>% select(sample_id, total_sample_abundance = n_seqs),
               by = "sample_id") %>%
    left_join(group_hill_diversity, by = "sample_id") %>%
    left_join(group_ordination,     by = "sample_id") %>%
    select(
      sample_id, total_sample_abundance, abundance_srs, richness_srs,
      hill_tax_div_q0, hill_tax_div_q1, hill_tax_div_q2,
      hill_tax_even_q1, hill_tax_even_q2,
      hill_phy_div_q0, hill_phy_div_q1, hill_phy_div_q2,
      hill_phy_even_q1, hill_phy_even_q2,
      starts_with("pco1_tax_"), starts_with("pco2_tax_"),
      starts_with("pco1_phy_"), starts_with("pco2_phy_")
    )
}

diversity_g_am <- compute_group_diversity("Glomeromycota", "g_am")
fwrite(diversity_g_am, "./output/diversity_g_am.txt", sep = "\t")

diversity_e_am <- compute_group_diversity("Endogonomycetes", "e_am")
fwrite(diversity_e_am, "./output/diversity_e_am.txt", sep = "\t")

dir.create("../../data/amf", showWarnings = FALSE)
file.copy(from = list.files("./output/", full.names = TRUE),
          to   = "../../data/amf/",
          recursive = TRUE)

# Rename the group-specific outputs in data/amf/ to the AMF/G-AMF/M-AMF
# terminology used downstream (alpha-diversity scripts, TITAN) -- this renames
# the copies in data/amf/ only, not the pipeline-local files under ./output/
# (or, for otu_classification_non_am.txt / otu_table_am_family.txt, the
# upstream utils/cluster_otus.R outputs they're copied from), so re-running
# this script is idempotent with respect to that local naming. Files with no
# guild-specific token (classification.txt, otu_table.txt, otu_table_srs.txt,
# dist_identity.txt, dist_phylogeny.txt, sequences.fasta, tree.newick,
# otu_table_all_family.txt) are left as-is.
file.rename("../../data/amf/diversity_am.txt",   "../../data/amf/diversity_amf.txt")
file.rename("../../data/amf/diversity_g_am.txt", "../../data/amf/diversity_g_amf.txt")
file.rename("../../data/amf/diversity_e_am.txt", "../../data/amf/diversity_m_amf.txt")

file.rename("../../data/amf/otu_table_g_am.txt",      "../../data/amf/otu_table_g_amf.txt")
file.rename("../../data/amf/otu_table_e_am.txt",      "../../data/amf/otu_table_m_amf.txt")
file.rename("../../data/amf/dist_identity_g_am.txt",  "../../data/amf/dist_identity_g_amf.txt")
file.rename("../../data/amf/dist_identity_e_am.txt",  "../../data/amf/dist_identity_m_amf.txt")
file.rename("../../data/amf/dist_phylogeny_g_am.txt", "../../data/amf/dist_phylogeny_g_amf.txt")
file.rename("../../data/amf/dist_phylogeny_e_am.txt", "../../data/amf/dist_phylogeny_m_amf.txt")

# From utils/cluster_otus.R, swept into ./output/ before this script runs
file.rename("../../data/amf/otu_classification_non_am.txt", "../../data/amf/otu_classification_non_amf.txt")
file.rename("../../data/amf/otu_table_am_family.txt",        "../../data/amf/otu_table_amf_family.txt")

message("---- Pipeline complete ----")
