#!/usr/bin/env Rscript
# classify_ecm.R — Classify ASVs and split into ECM / non-ECM pools
#
# Reads vsearch --usearch_global output files (one per rank-specific reference
# subset plus one for the full reference), classifies ASVs using dynamic
# similarity cutoffs, then writes three output files:
#   1) data/asv_classification.txt      — full ASV classification table
#   2) tmp_clusters/ecm_clusters.txt    — genus %in% ecm_genera list
#   3) tmp_clusters/non_ecm_clusters.txt — all other ASVs
#
# ECM designation is based on genus membership in the provided ecm_genera list.
# Downstream OTU clustering is handled by scripts/06_cluster_otus.sh.
#
# Usage:
#   Rscript utils/classify_ecm.R \
#     --cutoffs         cutoffs_ecm_ITS.txt \
#     --classification  eukaryome_ITS.classification \
#     --vsearch_species vsearch_species.txt \
#     --vsearch_genus   vsearch_genus.txt \
#     --vsearch_all     vsearch_all.txt \
#     --input           asv_sequences.fasta \
#     --output          asv_classification.txt \
#     --ecm_genera      ecm_genera.txt

# ── Packages ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
  library(tidyverse)
})

source("./utils/lca.R")
taxonomy_ranks <- c("kingdom", "phylum", "class", "order", "family", "lineage", "genus", "species")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-c", "--cutoffs"),
              type = "character", metavar = "FILE",
              help = "Tab-delimited cutoff file (columns: Rank, Dataset, cut-off) [required]"),
  make_option("--classification",
              type = "character", metavar = "FILE",
              help = "Tab-delimited reference classification file (must have 'id' column) [required]"),
  make_option("--vsearch_species",
              type = "character", default = "", metavar = "FILE",
              help = "vsearch --usearch_global output against species-level reference [required]"),
  make_option("--vsearch_genus",
              type = "character", default = "", metavar = "FILE",
              help = "vsearch --usearch_global output against genus-level reference [required]"),
  make_option("--vsearch_all",
              type = "character", default = "", metavar = "FILE",
              help = "vsearch --usearch_global output against full reference [required]"),
  make_option(c("-i", "--input"),
              type = "character", metavar = "FILE",
              help = "Input ASV FASTA file [required]"),
  make_option(c("-o", "--output"),
              type = "character", metavar = "FILE",
              help = "Output classification file [required]"),
  make_option("--ecm_genera",
              type = "character", metavar = "FILE",
              help = "Plain-text file listing ECM genera (one per line) [required]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage = "%prog [options]",
    description = paste(
      "Classify ASVs from vsearch --usearch_global output and split into",
      "ECM / non-ECM pools using dynamic similarity cutoffs and an ECM genera list."
    )
  )
)

# ── Validation ────────────────────────────────────────────────────────────────

if (is.null(opt$cutoffs))            stop("--cutoffs is required.")
if (is.null(opt$classification))     stop("--classification is required.")
if (nchar(opt$vsearch_species) == 0) stop("--vsearch_species is required.")
if (nchar(opt$vsearch_genus)   == 0) stop("--vsearch_genus is required.")
if (nchar(opt$vsearch_all)     == 0) stop("--vsearch_all is required.")
if (is.null(opt$input))              stop("--input is required.")
if (is.null(opt$output))             stop("--output is required.")
if (is.null(opt$ecm_genera))         stop("--ecm_genera is required.")

# ── Assign parameters ─────────────────────────────────────────────────────────

cutoff_file           <- opt$cutoffs
classification_file   <- opt$classification
vsearch_species_file  <- opt$vsearch_species
vsearch_genus_file    <- opt$vsearch_genus
vsearch_all_file      <- opt$vsearch_all
asv_sequence_file     <- opt$input
classification_output <- opt$output
ecm_genera_file       <- opt$ecm_genera

# ── Helper functions ──────────────────────────────────────────────────────────

# Parse vsearch --usearch_global output (--userfields query+target+id+ql+tl+alnlen+qcov+tcov).
# Filters to qcov > 90%, selects the hit with the highest id then highest qcov per query,
# and resolves ties (identical id and qcov) with the LCA approach.
# Returns a data.frame ready for classify_taxonomy().
format_vsearch <- function(vsearch_file, ref_classification) {
  empty <- data.frame(
    otu_id = character(), reference_id = character(),
    kingdom = character(), phylum = character(), class = character(),
    order = character(), family = character(), lineage = character(),
    genus = character(), species = character(), score = numeric(),
    abundance = integer()
  )

  if (!file.exists(vsearch_file) || file.info(vsearch_file)$size == 0L) {
    return(empty)
  }

  hits <- fread(vsearch_file) %>%
    mutate(
      otu_id    = V1,
      abundance = as.integer(str_extract(V1, "(?<=;size=)\\d+"))
    ) %>%
    select(
      otu_id,
      reference_id = V2,
      pident       = V3,
      qcov         = V7,
      abundance
    ) %>%
    filter(qcov > 90) %>%
    mutate(score = pident / 100) %>%
    left_join(ref_classification, by = c("reference_id" = "id"))

  if (nrow(hits) == 0) return(empty)

  # Select best hit(s) per query: highest id, then highest qcov
  best <- hits %>%
    group_by(otu_id) %>%
    filter(pident == max(pident)) %>%
    filter(qcov   == max(qcov)) %>%
    ungroup()

  # Resolve ties (same max id AND max qcov) with LCA
  resolved <- best %>%
    group_by(otu_id) %>%
    group_modify(~ {
      if (nrow(.x) == 1) return(.x)
      lca_result <- resolve_lca_rows(.x)
      out <- .x[1, , drop = FALSE]
      for (r in taxonomy_ranks) out[[r]] <- lca_result[[r]]
      out
    }) %>%
    ungroup()

  resolved %>%
    select(otu_id, reference_id, kingdom, phylum, class, order, family, lineage, genus,
           species, score, abundance)
}

# Classify ASVs at a single rank.
# vsearch_data: pre-read data.frame from format_vsearch() — NOT a file path.
classify_taxonomy <- function(vsearch_data, cutoff_file, rank = "species",
                              excluded_otus = NULL) {

  valid_ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")
  if (!rank %in% valid_ranks) {
    stop("Invalid rank. Must be one of: ", paste(valid_ranks, collapse = ", "))
  }

  rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "lineage", "genus", "species")
  rank_index     <- which(rank_hierarchy == rank)
  higher_ranks   <- rank_hierarchy[1:(rank_index - 1)]
  lower_ranks    <- if (rank_index < length(rank_hierarchy)) {
    rank_hierarchy[(rank_index + 1):length(rank_hierarchy)]
  } else {
    character(0)
  }

  vsearch_data <- vsearch_data %>%
    filter(!otu_id %in% excluded_otus)

  all_otu_ids <- distinct(vsearch_data, otu_id)
  cutoff_rank <- rank

  # Get cutoffs for OTUs with specific taxonomic matches.
  # For genus rank, look up by genus name first (genus-specific cutoffs), then
  # fall back to higher ranks. For all other ranks, use higher ranks only.
  lookup_cols <- if (rank == "genus") c("genus", higher_ranks) else higher_ranks

  if (length(lookup_cols) > 0) {
    cutoffs_matched <- vsearch_data %>%
      mutate(rank = cutoff_rank) %>%
      select(otu_id, rank, all_of(lookup_cols)) %>%
      pivot_longer(-c(otu_id, rank), names_to = "level", values_to = "taxa") %>%
      left_join(cutoff_file, by = c("rank", "taxa")) %>%
      filter(!is.na(cutoff)) %>%
      group_by(otu_id) %>%
      slice(1) %>%
      ungroup() %>%
      select(otu_id, cutoff)
  } else {
    cutoffs_matched <- data.frame(otu_id = character(), cutoff = numeric())
  }

  # Find OTUs without specific cutoffs and assign "All" cutoff
  cutoffs_unmatched <- all_otu_ids %>%
    anti_join(cutoffs_matched, by = "otu_id") %>%
    mutate(cutoff = cutoff_file$cutoff[cutoff_file$taxa == "All" & cutoff_file$rank == cutoff_rank][1])

  all_cutoffs <- bind_rows(cutoffs_matched, cutoffs_unmatched)

  result <- vsearch_data %>%
    left_join(all_cutoffs, by = "otu_id") %>%
    filter(score >= cutoff) %>%
    mutate(rank = !!rank) %>%
    group_by(otu_id) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    ungroup()

  # Set all ranks below the classification rank to "unclassified"
  if (length(lower_ranks) > 0) {
    for (lower_rank in lower_ranks) {
      result <- result %>%
        mutate(!!lower_rank := "unclassified")
    }
  }

  return(result)
}

# Get cutoffs for any rank
get_rank_cutoff <- function(taxa_file, unique_taxa_cutoffs, cutoff_file, rank, superranks) {

  # Handle kingdom (no superranks)
  if (length(superranks) == 0) {
    default_cutoff <- cutoff_file %>%
      filter(rank == !!rank & taxa == "All") %>%
      pull(cutoff)

    return(taxa_file %>% mutate(!!paste0(rank, "_cutoff") := default_cutoff))
  }

  cutoffs_all <- taxa_file %>%
    select(rank, all_of(superranks), kingdom) %>%
    filter(rank == !!rank) %>%
    pivot_longer(-rank, values_to = "taxa") %>%
    left_join(unique_taxa_cutoffs %>% filter(rank == !!rank), by = c("rank", "taxa")) %>%
    filter(!is.na(cutoff)) %>%
    select(-rank, rank = name) %>%
    unique()

  cutoffs_list <- setNames(
    lapply(superranks, function(sr) cutoffs_all %>% filter(rank == sr) %>% select(-rank)),
    superranks
  )

  # Build case_when conditions dynamically
  conditions <- lapply(superranks, function(sr) {
    expr(!!sym(sr) %in% cutoffs_list[[!!sr]][["taxa"]] ~
           cutoffs_list[[!!sr]][["cutoff"]][match(!!sym(sr), cutoffs_list[[!!sr]][["taxa"]])])
  })

  default_cutoff <- cutoff_file %>%
    filter(rank == !!rank & taxa == "All") %>%
    pull(cutoff)

  default_cutoff <- if (length(default_cutoff) > 0) default_cutoff else NA_real_

  taxa_file %>%
    mutate(!!paste0(rank, "_cutoff") := case_when(!!!conditions, TRUE ~ default_cutoff))
}

# ── Read reference classification, cutoff file and ECM genera list ─────────────

# Read reference classification once — reused by every format_vsearch() call
ref_classification <- fread(classification_file)

taxon_cutoffs <- fread(cutoff_file) %>%
  select(rank = rank, taxa = dataset, cutoff = "cut-off")

ecm_genera <- fread(ecm_genera_file)$genus
ecm_genera <- ecm_genera[nchar(trimws(ecm_genera)) > 0]

# ── Pre-read vsearch files (each file read exactly once) ─────────────────────

message("Reading vsearch outputs...")

vsearch_data_species <- format_vsearch(vsearch_species_file, ref_classification)
vsearch_data_genus   <- format_vsearch(vsearch_genus_file,   ref_classification)
vsearch_data_all     <- format_vsearch(vsearch_all_file,     ref_classification)

# ── Classify vsearch results ──────────────────────────────────────────────────
#
# Steps run in order; each step excludes OTUs already classified in prior steps.
# Genus is attempted twice: first against the genus-specific reference subset,
# then against the full reference (vsearch_all) using genus-specific cutoffs.

classification_steps <- list(
  list(data = vsearch_data_species, rank = "species", label = "species (vsearch_species)"),
  list(data = vsearch_data_genus,   rank = "genus",   label = "genus (vsearch_genus)"),
  list(data = vsearch_data_all,     rank = "genus",   label = "genus (vsearch_all)"),
  list(data = vsearch_data_all,     rank = "family",  label = "family"),
  list(data = vsearch_data_all,     rank = "order",   label = "order"),
  list(data = vsearch_data_all,     rank = "class",   label = "class"),
  list(data = vsearch_data_all,     rank = "phylum",  label = "phylum"),
  list(data = vsearch_data_all,     rank = "kingdom", label = "kingdom")
)

classified_otus     <- character()
all_classifications <- list()

for (step in classification_steps) {
  rank <- step$rank
  message("Classifying ", step$label, "...")

  hits <- classify_taxonomy(
    vsearch_data  = step$data,
    cutoff_file   = taxon_cutoffs,
    rank          = rank,
    excluded_otus = classified_otus
  )

  if (nrow(hits) > 0) {
    hits_identified <- hits %>%
      filter(!!sym(rank) != "unclassified" & !!sym(rank) != "")

    if (nrow(hits_identified) > 0) {
      all_classifications[[length(all_classifications) + 1]] <- hits_identified
      classified_otus <- c(classified_otus, hits_identified$otu_id)
    }
  }
}

all_classifications <- bind_rows(all_classifications)

# ── Get unclassified ASVs ─────────────────────────────────────────────────────

all_unclassified_ids <- if (nrow(vsearch_data_all) > 0) {
    distinct(vsearch_data_all, otu_id)
  } else {
    data.frame(otu_id = character())
  }

all_unclassified <- all_unclassified_ids %>%
  anti_join(all_classifications, by = "otu_id") %>%
  mutate(
    reference_id = "",
    kingdom      = "unclassified",
    phylum       = "unclassified",
    class        = "unclassified",
    order        = "unclassified",
    family       = "unclassified",
    lineage      = "unclassified",
    genus        = "unclassified",
    species      = "unclassified",
    score        = as.numeric(NA),
    abundance    = as.integer(str_extract(otu_id, "(?<=;size=)\\d+")),
    cutoff       = as.numeric(NA),
    rank         = ""
  )

# ── Save the complete classification file ─────────────────────────────────────

all_classifications %>%
  bind_rows(all_unclassified) %>%
  select(otu_id, reference_id, kingdom, phylum, class, order, family, lineage, genus,
         species, rank, cutoff, score, abundance) %>%
  fwrite(classification_output, sep = "\t")

message("Saved classification: ", classification_output)

# ── Compute rank-specific cutoffs ─────────────────────────────────────────────

taxa_file <- bind_rows(all_classifications, all_unclassified)

unique_taxa_cutoffs <- taxa_file %>%
  select(kingdom, phylum, class, order, family, lineage, genus, species, rank) %>%
  pivot_longer(-rank, names_to = "level", values_to = "taxa") %>%
  filter(!taxa %in% c("unclassified", "") & !is.na(taxa)) %>%
  select(-level) %>%
  unique() %>%
  left_join(taxon_cutoffs, by = c("rank", "taxa")) %>%
  filter(!is.na(cutoff))

superranks_list <- list(
  species = c("genus", "lineage", "family", "order", "class", "phylum", "kingdom"),
  genus   = c("lineage", "family", "order", "class", "phylum", "kingdom"),
  family  = c("order", "class", "phylum", "kingdom"),
  order   = c("class", "phylum", "kingdom"),
  class   = c("phylum", "kingdom"),
  phylum  = c("kingdom"),
  kingdom = NULL
)

taxa_cutoffs <- taxa_file
for (rank in names(superranks_list)) {
  taxa_cutoffs <- get_rank_cutoff(taxa_cutoffs, unique_taxa_cutoffs, taxon_cutoffs, rank, superranks_list[[rank]])
}

# ── Write ECM and non-ECM pools ───────────────────────────────────────────────
#
# A sequence enters the ECM pool if its genus is in the ECM genera list.
# All other sequences enter the non-ECM pool.

dir.create("./tmp_clusters", showWarnings = FALSE)

ecm_clusters <- taxa_cutoffs %>%
  filter(genus %in% ecm_genera)

non_ecm_clusters <- taxa_cutoffs %>%
  filter(!genus %in% ecm_genera)

fwrite(ecm_clusters,     "./tmp_clusters/ecm_clusters.txt",     sep = "\t")
fwrite(non_ecm_clusters, "./tmp_clusters/non_ecm_clusters.txt", sep = "\t")

message("ECM pool:     ", nrow(ecm_clusters),     " sequences -> tmp_clusters/ecm_clusters.txt")
message("Non-ECM pool: ", nrow(non_ecm_clusters), " sequences -> tmp_clusters/non_ecm_clusters.txt")
