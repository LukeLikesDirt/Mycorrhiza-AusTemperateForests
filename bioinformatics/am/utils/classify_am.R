#!/usr/bin/env Rscript
# classify_am.R — Classify ASVs and split into AM / non-AM pools
#
# AM fungi evaluated here span two lineages: phylum Glomeromycota, and order
# Densosporales (class Endogonomycetes, Mucoromycota) — both are AM fungi.
# All of order Densosporales is included (not just the named families
# Densosporaceae/Planticonsortiaceae): Lutz et al. 2025 (Fungal Ecology 74,
# 101407) — co-authored by the same EUKARYOME curators (Mikryukov, Tedersoo)
# who name the Densosporales families in Tedersoo et al. 2024 (MycoKeys 107)
# — define putative E-AMF at order rank ("OTUs assigned 'Densosporales' at
# order level"), citing strong root-colonisation evidence for the order as a
# whole while flagging finer within-order resolution as future work.
#
# Reads vsearch --usearch_global output files (one per rank-specific reference
# subset plus one for the full reference), classifies ASVs using dynamic
# similarity cutoffs, then writes three output files:
#   1) data/asv_classification.txt      — full ASV classification table
#   2) tmp_clusters/am_clusters.txt     — the AM pool: phylum == "Glomeromycota"
#      (score >= min_sim_glom) PLUS any ASV whose order is in --am_orders
#      (default: Densosporales), which are AM fungi outside Glomeromycota.
#      Each AM ASV carries an `am_group` label (Glomeromycota / Endogonomycetes).
#   3) tmp_clusters/non_am_clusters.txt — all other ASVs
#
# Downstream OTU clustering is handled by scripts/03_cluster_otus.sh.
#
# Usage:
#   Rscript utils/classify_am.R \
#     --cutoffs         cutoffs_V4_am.txt \
#     --classification  eukaryome_V4.classification \
#     --vsearch_species vsearch_species.txt \
#     --vsearch_genus   vsearch_genus.txt \
#     --vsearch_family  vsearch_family.txt \
#     --vsearch_all     vsearch_all.txt \
#     --input           asv_sequences.fasta \
#     --output          asv_classification.txt \
#     --min_sim_glom    0.95 \
#     --am_orders       Densosporales

# ── Packages ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
  library(tidyverse)
})

source("./utils/lca.R")
taxonomy_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

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
  make_option("--vsearch_family",
              type = "character", default = "", metavar = "FILE",
              help = "vsearch --usearch_global output against family-level reference [required]"),
  make_option("--vsearch_all",
              type = "character", default = "", metavar = "FILE",
              help = "vsearch --usearch_global output against full reference [required]"),
  make_option(c("-i", "--input"),
              type = "character", metavar = "FILE",
              help = "Input ASV FASTA file [required]"),
  make_option(c("-o", "--output"),
              type = "character", metavar = "FILE",
              help = "Output classification file [required]"),
  make_option("--min_sim_glom",
              type = "double", default = 0.95, metavar = "FLOAT",
              help = paste("Minimum similarity (0-1) to any Glomeromycota reference",
                           "required for a phylum-Glomeromycota sequence to remain in",
                           "the AM pool (Densosporales are not subject to this filter)",
                           "[default: %default]")),
  make_option("--am_orders",
              type = "character", default = "Densosporales", metavar = "STR",
              help = paste("Comma-separated orders outside phylum Glomeromycota to",
                           "include in the AM pool as AM fungi (default: Densosporales).",
                           "Pass \"\" to include Glomeromycota only [default: %default]"))
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage = "%prog [options]",
    description = paste(
      "Classify ASVs from vsearch --usearch_global output and split into",
      "AM / non-AM pools using dynamic similarity cutoffs."
    )
  )
)

# ── Validation ────────────────────────────────────────────────────────────────

if (is.null(opt$cutoffs))            stop("--cutoffs is required.")
if (is.null(opt$classification))     stop("--classification is required.")
if (nchar(opt$vsearch_species) == 0) stop("--vsearch_species is required.")
if (nchar(opt$vsearch_genus)   == 0) stop("--vsearch_genus is required.")
if (nchar(opt$vsearch_family)  == 0) stop("--vsearch_family is required.")
if (nchar(opt$vsearch_all)     == 0) stop("--vsearch_all is required.")
if (is.null(opt$input))              stop("--input is required.")
if (is.null(opt$output))             stop("--output is required.")

# ── Assign parameters ─────────────────────────────────────────────────────────

cutoff_file           <- opt$cutoffs
classification_file   <- opt$classification
vsearch_species_file  <- opt$vsearch_species
vsearch_genus_file    <- opt$vsearch_genus
vsearch_family_file   <- opt$vsearch_family
vsearch_all_file      <- opt$vsearch_all
asv_sequence_file     <- opt$input
classification_output <- opt$output
min_sim_glom          <- opt$min_sim_glom

# Orders (outside Glomeromycota) to treat as AM fungi; empty string -> none
am_orders <- trimws(str_split(opt$am_orders, ",")[[1]])
am_orders <- am_orders[nzchar(am_orders)]

# ── Helper functions ──────────────────────────────────────────────────────────

# Parse vsearch --usearch_global output (--userfields query+target+id+ql+tl+alnlen+qcov+tcov).
# Filters to qcov > 90%, selects the hit with the highest id then highest qcov per query,
# and resolves ties (identical id and qcov) with the LCA approach.
# Returns a data.frame ready for classify_taxonomy().
format_vsearch <- function(vsearch_file, ref_classification) {
  empty <- data.frame(
    otu_id = character(), reference_id = character(),
    kingdom = character(), phylum = character(), class = character(),
    order = character(), family = character(), genus = character(),
    species = character(), score = numeric(), abundance = integer()
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
    select(otu_id, reference_id, kingdom, phylum, class, order, family, genus, species,
           score, abundance)
}

# Classify ASVs at a single rank.
# vsearch_data: pre-read data.frame from format_vsearch() — NOT a file path.
classify_taxonomy <- function(vsearch_data, cutoff_file, rank = "species",
                              excluded_otus = NULL) {

  valid_ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")
  if (!rank %in% valid_ranks) {
    stop("Invalid rank. Must be one of: ", paste(valid_ranks, collapse = ", "))
  }

  rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
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

  # Get cutoffs for OTUs with specific taxonomic matches
  if (length(higher_ranks) > 0) {
    cutoffs_matched <- vsearch_data %>%
      mutate(rank = cutoff_rank) %>%
      select(otu_id, rank, all_of(higher_ranks)) %>%
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

# ── Read reference classification and cutoff file ─────────────────────────────

# Read reference classification once — reused by every format_vsearch() call
ref_classification <- fread(classification_file)

taxon_cutoffs <- fread(cutoff_file, fill = TRUE) %>%
  select(rank = rank, taxa = dataset, cutoff = "cut-off")

# ── Pre-read vsearch files (each file read exactly once) ─────────────────────

message("Reading vsearch outputs...")

vsearch_data_species <- format_vsearch(vsearch_species_file, ref_classification)
vsearch_data_genus   <- format_vsearch(vsearch_genus_file,   ref_classification)
vsearch_data_family  <- format_vsearch(vsearch_family_file,  ref_classification)
vsearch_data_all     <- format_vsearch(vsearch_all_file,     ref_classification)

# Map each classification rank to its pre-read data
vsearch_data_by_rank <- list(
  species = vsearch_data_species,
  genus   = vsearch_data_genus,
  family  = vsearch_data_family,
  order   = vsearch_data_all,
  class   = vsearch_data_all,
  phylum  = vsearch_data_all,
  kingdom = vsearch_data_all
)

# ── Classify vsearch results ──────────────────────────────────────────────────

classified_otus     <- character()
all_classifications <- list()

ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")

for (i in seq_along(ranks)) {
  rank <- ranks[i]
  message("Classifying ", rank, "...")

  hits <- classify_taxonomy(
    vsearch_data  = vsearch_data_by_rank[[rank]],
    cutoff_file   = taxon_cutoffs,
    rank          = rank,
    excluded_otus = classified_otus
  )

  if (nrow(hits) > 0) {
    hits_identified <- hits %>%
      filter(!!sym(rank) != "unclassified" & !!sym(rank) != "")

    if (nrow(hits_identified) > 0) {
      all_classifications[[rank]] <- hits_identified
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
  select(otu_id, reference_id, kingdom, phylum, class, order, family, genus, species,
         rank, cutoff, score, abundance) %>%
  fwrite(classification_output, sep = "\t")

message("Saved classification: ", classification_output)

# ── Compute rank-specific cutoffs ─────────────────────────────────────────────

taxa_file <- bind_rows(all_classifications, all_unclassified)

unique_taxa_cutoffs <- taxa_file %>%
  select(kingdom, phylum, class, order, family, genus, species, rank) %>%
  pivot_longer(-rank, names_to = "level", values_to = "taxa") %>%
  filter(!taxa %in% c("unclassified", "") & !is.na(taxa)) %>%
  select(-level) %>%
  unique() %>%
  left_join(taxon_cutoffs, by = c("rank", "taxa")) %>%
  filter(!is.na(cutoff))

superranks_list <- list(
  species = c("genus", "family", "order", "class", "phylum", "kingdom"),
  genus   = c("family", "order", "class", "phylum", "kingdom"),
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

# ── Pre-filter: demote Glomeromycota sequences below min_sim_glom ─────────────
#
# A phylum-Glomeromycota sequence is kept in the AM pool only if the score from
# the classification step is >= min_sim_glom (Densosporales are not subject to
# this filter). Sequences below this threshold are demoted to
# "Glomeromycota(unclassified)" and enter the non-AM pool.

demote_flag <- taxa_cutoffs$phylum == "Glomeromycota" &
               (is.na(taxa_cutoffs$score) | taxa_cutoffs$score < min_sim_glom)

message("Glomeromycota pre-filter (min_sim_glom = ", min_sim_glom, "): ",
        sum(demote_flag), " sequence(s) demoted to non-AM pool.")

taxa_cutoffs <- taxa_cutoffs %>%
  mutate(
    phylum  = if_else(demote_flag, "Glomeromycota(unclassified)", phylum),
    class   = if_else(demote_flag, "unclassified", class),
    order   = if_else(demote_flag, "unclassified", order),
    family  = if_else(demote_flag, "unclassified", family),
    genus   = if_else(demote_flag, "unclassified", genus),
    species = if_else(demote_flag, "unclassified", species),
    rank    = if_else(demote_flag, "unclassified", rank),
    cutoff  = if_else(demote_flag, NA_real_, cutoff)
  )

# ── Write AM (Glomeromycota + configured orders) and non-AM pools ──────────────
#
# The AM pool contains all Glomeromycota ASVs plus any ASV classified to an
# order listed in --am_orders (default: Densosporales) — AM fungi outside
# phylum Glomeromycota. Each AM ASV is labelled in `am_group` as
# "Glomeromycota" or "Endogonomycetes" (the class Densosporales belongs to).
# The two pools are an exact partition, so no ASV is duplicated or dropped.

dir.create("./tmp_clusters", showWarnings = FALSE)

is_am <- (taxa_cutoffs$phylum == "Glomeromycota") |
         (taxa_cutoffs$order %in% am_orders)   # %in% is NA-safe
is_am[is.na(is_am)] <- FALSE

taxa_cutoffs$am_group <- case_when(
  taxa_cutoffs$order %in% am_orders      ~ "Endogonomycetes",
  taxa_cutoffs$phylum == "Glomeromycota" ~ "Glomeromycota",
  TRUE                                   ~ NA_character_
)

am_clusters     <- taxa_cutoffs[is_am, ]
non_am_clusters <- taxa_cutoffs[!is_am, ]

fwrite(am_clusters,     "./tmp_clusters/am_clusters.txt",     sep = "\t")
fwrite(non_am_clusters, "./tmp_clusters/non_am_clusters.txt", sep = "\t")

message("AM pool:     ", nrow(am_clusters), " sequences (Glomeromycota",
        if (length(am_orders) > 0) paste0(" + ", paste(am_orders, collapse = ", ")) else "",
        ") -> tmp_clusters/am_clusters.txt")
message("Non-AM pool: ", nrow(non_am_clusters), " sequences -> tmp_clusters/non_am_clusters.txt")
