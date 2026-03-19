#!/usr/bin/env Rscript
# prepare_glomeromycota_clusters.R — Classify ASVs and prepare Glomeromycota clusters
#
# Reads vsearch --usearch_global output files (one per rank-specific reference
# subset plus one for the full reference), classifies ASVs using dynamic
# similarity cutoffs, then performs reference-based clustering for unclassified
# ASVs within Glomeromycota and non-Glomeromycota groups.
#
# Usage:
#   Rscript utils/prepare_glomeromycota_clusters.R \
#     --cutoffs       cutoffs_glom_V4.txt \
#     --classification eukaryome_V4.classification \
#     --vsearch_species vsearch_species.txt \
#     --vsearch_genus   vsearch_genus.txt \
#     --vsearch_family  vsearch_family.txt \
#     --vsearch_all     vsearch_all.txt \
#     --input           asv_sequences.fasta \
#     --output          asv_classification.txt \
#     --threads         64

# ── Packages ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
  library(tidyverse)
})

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
  make_option("--threads",
              type = "integer", default = 1L, metavar = "INT",
              help = "Number of threads for vsearch [default: %default]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage = "%prog [options]",
    description = paste(
      "Classify ASVs from vsearch --usearch_global output and prepare",
      "Glomeromycota clusters using dynamic similarity cutoffs."
    )
  )
)

# ── Validation ────────────────────────────────────────────────────────────────

if (is.null(opt$cutoffs))          stop("--cutoffs is required.")
if (is.null(opt$classification))   stop("--classification is required.")
if (nchar(opt$vsearch_species) == 0) stop("--vsearch_species is required.")
if (nchar(opt$vsearch_genus)   == 0) stop("--vsearch_genus is required.")
if (nchar(opt$vsearch_family)  == 0) stop("--vsearch_family is required.")
if (nchar(opt$vsearch_all)     == 0) stop("--vsearch_all is required.")
if (is.null(opt$input))            stop("--input is required.")
if (is.null(opt$output))           stop("--output is required.")

# ── Assign parameters ─────────────────────────────────────────────────────────

cutoff_file          <- opt$cutoffs
classification_file  <- opt$classification
vsearch_species_file <- opt$vsearch_species
vsearch_genus_file   <- opt$vsearch_genus
vsearch_family_file  <- opt$vsearch_family
vsearch_all_file     <- opt$vsearch_all
asv_sequence_file    <- opt$input
classification_output <- opt$output
threads              <- opt$threads

# ── Helper functions ──────────────────────────────────────────────────────────

# Parse vsearch --usearch_global output (--userfields query+target+id+ql+tl+alnlen)
format_vsearch <- function(vsearch_file, classification_file) {
  fread(vsearch_file) %>%
    mutate(
      otu_id    = V1,
      abundance = as.integer(str_extract(V1, "(?<=;size=)\\d+"))
    ) %>%
    select(
      otu_id, reference_id = V2,
      pident = V3,
      abundance
    ) %>%
    mutate(
      score = round(pident / 100 / 0.005) * 0.005
    ) %>%
    left_join(fread(classification_file), by = c("reference_id" = "id")) %>%
    select(otu_id, reference_id, kingdom, phylum, class, order, family, genus, species,
           score, abundance)
}

# Main classification function
classify_taxonomy <- function(vsearch_file, classification_file, cutoff_file, rank = "species",
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

  vsearch_data <- format_vsearch(vsearch_file, classification_file) %>%
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

  # Set all ranks below the classification rank to "unidentified"
  if (length(lower_ranks) > 0) {
    for (lower_rank in lower_ranks) {
      result <- result %>%
        mutate(!!lower_rank := "unidentified")
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

# Clustering function (uses vsearch --usearch_global internally)
cluster_at_rank <- function(this_rank, this_subrank, this_supertaxon, this_superrank = NULL,
                            taxa_cutoffs, asv_sequences, threads) {

  # Initialise or filter taxa file
  if (!is.null(this_superrank)) {
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt"))

    # Only filter if not processing all taxa
    if (this_supertaxon != "all") {
      taxa_cutoffs <- taxa_cutoffs %>%
        filter(!!sym(this_superrank) == this_supertaxon)
    }
  }

  # Ensure abundance column exists
  if (!"abundance" %in% names(taxa_cutoffs)) {
    taxa_cutoffs <- taxa_cutoffs %>%
      mutate(abundance = as.integer(str_extract(otu_id, "(?<=;size=)\\d+")))
  }

  fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

  # Get the cutoff for this rank
  this_cutoff <- taxa_cutoffs %>%
    pull(!!sym(paste0(this_rank, "_cutoff"))) %>%
    unique() %>%
    na.omit() %>%
    first()

  if (is.na(this_cutoff) || length(this_cutoff) == 0) {
    stop("No valid cutoff found for ", this_rank, " in ", this_supertaxon)
  }

  message("Starting clustering for ", this_rank, " in ", this_supertaxon,
          " at cutoff ", this_cutoff, " !!!")

  # Define rank hierarchy for proper taxonomy handling
  rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  rank_index     <- which(rank_hierarchy == this_rank)
  lower_ranks    <- if (rank_index < length(rank_hierarchy)) {
    rank_hierarchy[(rank_index + 1):length(rank_hierarchy)]
  } else {
    character(0)
  }

  ### Reference-based clustering loop ###
  repeat {
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

    # Find identified ASVs at this rank
    identified_asvs <- taxa_cutoffs %>%
      filter(!!sym(this_rank) != "unidentified", !!sym(this_rank) != "") %>%
      pull(otu_id)

    if (length(identified_asvs) == 0) {
      message("No identified ASVs to use as references for ", this_supertaxon)
      break
    }

    # Find unidentified ASVs at this rank
    unidentified_asvs <- taxa_cutoffs %>%
      filter(!!sym(this_rank) == "unidentified" | !!sym(this_rank) == "") %>%
      pull(otu_id)

    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      message("No unidentified ASVs remaining for ", this_supertaxon)
      break
    }

    # Filter sequences
    identified_sequences   <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]

    writeXStringSet(identified_sequences,   "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")

    message("Clustering ", initial_unidentified_count, " unidentified ASVs to ",
            this_rank, " in ", this_supertaxon, " at ", this_cutoff * 100, "% similarity...")

    # Run vsearch --usearch_global
    system2("vsearch", args = c(
      "--usearch_global", "./tmp/unidentified_fasta",
      "--db",             "./tmp/identified_fasta",
      "--id",             "0.5",
      "--strand",         "both",
      "--maxaccepts",     "1",
      "--maxrejects",     "0",
      "--userout",        "./tmp/unidentified_vsearch.out",
      "--userfields",     "query+target+id+ql+tl+alnlen",
      "--threads",        as.character(threads)
    ), stdout = FALSE, stderr = FALSE)

    # Check if vsearch output is empty
    vsearch_output_size <- file.info("./tmp/unidentified_vsearch.out")$size

    if (is.na(vsearch_output_size) || vsearch_output_size == 0) {
      message("No vsearch hits found for ", this_supertaxon, " - stopping clustering")
      break
    }

    # Process vsearch results
    new_clusters <- fread("./tmp/unidentified_vsearch.out") %>%
      setnames(c("otu_id", "reference_id", "pident", "qlen", "tl", "alnlen")) %>%
      group_by(otu_id) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        score     = round(pident / 100 / 0.005) * 0.005,
        abundance = as.integer(str_extract(otu_id, "(?<=;size=)\\d+"))
      ) %>%
      filter(score >= this_cutoff)

    # Check if any sequences passed the similarity filter
    if (nrow(new_clusters) == 0) {
      message("No sequences passed similarity filter for ", this_supertaxon)
      break
    }

    # Join with reference taxonomy
    new_clusters <- new_clusters %>%
      left_join(
        taxa_cutoffs %>%
          select(reference_id = otu_id, all_of(rank_hierarchy), ends_with("_cutoff")),
        by = "reference_id"
      ) %>%
      mutate(
        rank   = this_rank,
        cutoff = this_cutoff
      )

    # Set lower ranks to unidentified (after inheriting from reference)
    if (length(lower_ranks) > 0) {
      for (lower_rank in lower_ranks) {
        new_clusters <- new_clusters %>%
          mutate(!!lower_rank := "unidentified")
      }
    }

    # Clean up vsearch-specific columns before binding
    new_clusters <- new_clusters %>%
      select(-pident, -alnlen, -qlen, -tl)

    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!otu_id %in% new_clusters$otu_id)
    )

    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(!!sym(this_rank) == "unidentified" | !!sym(this_rank) == "") %>%
      nrow()

    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

    message("Assigned ", initial_unidentified_count - remaining_unidentified_count,
            " ASVs. ", remaining_unidentified_count, " remain unidentified.")

    # Check stopping conditions
    if (remaining_unidentified_count == 0) {
      message("No more unidentified ASVs for ", this_supertaxon)
      break
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message("No new assignments for ", this_supertaxon, " - stopping")
      break
    } else {
      message("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon)
    }
  }

  # Standardize column order for output files
  standard_cols <- c("otu_id", "reference_id", "kingdom", "phylum", "class", "order",
                     "family", "genus", "species", "rank", "cutoff", "score", "abundance",
                     "kingdom_cutoff", "phylum_cutoff", "class_cutoff", "order_cutoff",
                     "family_cutoff", "genus_cutoff", "species_cutoff")

  final_data <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

  # Ensure abundance exists
  if (!"abundance" %in% names(final_data)) {
    final_data <- final_data %>%
      mutate(abundance = as.integer(str_extract(otu_id, "(?<=;size=)\\d+")))
  }

  # Select only columns that exist
  existing_cols <- intersect(standard_cols, names(final_data))

  # Write outputs with standardized columns
  final_data %>%
    filter(!!sym(this_rank) == "unidentified" | !!sym(this_rank) == "") %>%
    select(all_of(existing_cols)) %>%
    fwrite(paste0("./tmp_clusters/", this_rank, "_unidentified.txt"), sep = "\t")

  final_data %>%
    filter(!!sym(this_rank) != "unidentified" & !!sym(this_rank) != "") %>%
    select(all_of(existing_cols)) %>%
    fwrite(paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

  message("Completed clustering for ", this_rank, " in ", this_supertaxon, " !!!")
}

# ── Define vsearch files per rank ─────────────────────────────────────────────

vsearch_files <- list(
  species = vsearch_species_file,
  genus   = vsearch_genus_file,
  family  = vsearch_family_file,
  order   = vsearch_all_file,
  class   = vsearch_all_file,
  phylum  = vsearch_all_file,
  kingdom = vsearch_all_file
)

# ── Read cutoff file ──────────────────────────────────────────────────────────

taxon_cutoffs <- fread(cutoff_file) %>%
  select(rank = rank, taxa = dataset, cutoff = "cut-off")

# ── Classify vsearch results ──────────────────────────────────────────────────

classified_otus    <- character()
all_classifications <- list()

# Classify each rank
ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")

for (i in seq_along(ranks)) {
  rank <- ranks[i]
  message("Classifying ", rank, "...")

  hits <- classify_taxonomy(
    vsearch_file        = vsearch_files[[rank]],
    classification_file = classification_file,
    cutoff_file         = taxon_cutoffs,
    rank                = rank,
    excluded_otus       = classified_otus
  )

  if (nrow(hits) > 0) {
    hits_identified <- hits %>%
      filter(!!sym(rank) != "unidentified" & !!sym(rank) != "")

    if (nrow(hits_identified) > 0) {
      all_classifications[[rank]] <- hits_identified
      classified_otus <- c(classified_otus, hits_identified$otu_id)
    }
  }
}

all_classifications <- bind_rows(all_classifications)

# ── Get unclassified ASVs ─────────────────────────────────────────────────────

all_unclassified <- fread(vsearch_all_file) %>%
  distinct(otu_id = V1) %>%
  anti_join(all_classifications, by = "otu_id") %>%
  mutate(
    reference_id = "",
    kingdom      = "unidentified",
    phylum       = "unidentified",
    class        = "unidentified",
    order        = "unidentified",
    family       = "unidentified",
    genus        = "unidentified",
    species      = "unidentified",
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

# ── Cluster unclassified ──────────────────────────────────────────────────────

# Create temporary directories
dir.create("tmp",          showWarnings = FALSE)
dir.create("tmp_clusters", showWarnings = FALSE)

taxa_file     <- bind_rows(all_classifications, all_unclassified)
asv_sequences <- readDNAStringSet(asv_sequence_file)

# Match unique taxa to cut-offs
unique_taxa_cutoffs <- taxa_file %>%
  select(kingdom, phylum, class, order, family, genus, species, rank) %>%
  pivot_longer(-rank, names_to = "level", values_to = "taxa") %>%
  filter(!taxa %in% c("unidentified", "") & !is.na(taxa)) %>%
  select(-level) %>%
  unique() %>%
  left_join(taxon_cutoffs, by = c("rank", "taxa")) %>%
  filter(!is.na(cutoff))

# Define superranks for each rank
superranks_list <- list(
  species = c("genus", "family", "order", "class", "phylum", "kingdom"),
  genus   = c("family", "order", "class", "phylum", "kingdom"),
  family  = c("order", "class", "phylum", "kingdom"),
  order   = c("class", "phylum", "kingdom"),
  class   = c("phylum", "kingdom"),
  phylum  = c("kingdom"),
  kingdom = NULL
)

# Get cutoffs for all ranks
taxa_cutoffs <- taxa_file
for (rank in names(superranks_list)) {
  taxa_cutoffs <- get_rank_cutoff(taxa_cutoffs, unique_taxa_cutoffs, taxon_cutoffs, rank, superranks_list[[rank]])
}

# ── Run clustering ────────────────────────────────────────────────────────────

# (1) GLOMEROMYCOTA: Filter and cluster at family level
standard_cols <- c("otu_id", "reference_id", "kingdom", "phylum", "class", "order",
                   "family", "genus", "species", "rank", "cutoff", "score", "abundance",
                   "kingdom_cutoff", "phylum_cutoff", "class_cutoff", "order_cutoff",
                   "family_cutoff", "genus_cutoff", "species_cutoff")

glom_seqs <- taxa_cutoffs %>% filter(phylum == "Glomeromycota")

if (nrow(glom_seqs) > 0) {
  fwrite(glom_seqs, "./tmp_clusters/order_glom_clusters.txt", sep = "\t")
  cluster_at_rank("family", "genus", "all", "order_glom", taxa_cutoffs, asv_sequences, threads)

  # Save Glomeromycota results immediately before they get overwritten
  glom_parts <- list()
  if (file.exists("./tmp_clusters/family_clusters.txt")) {
    dt <- fread("./tmp_clusters/family_clusters.txt")
    if (nrow(dt) > 0) glom_parts[[length(glom_parts) + 1]] <- dt
  }
  if (file.exists("./tmp_clusters/family_unidentified.txt")) {
    dt <- fread("./tmp_clusters/family_unidentified.txt")
    if (nrow(dt) > 0) glom_parts[[length(glom_parts) + 1]] <- dt
  }
  glom_clusters <- bind_rows(glom_parts) %>% select(any_of(standard_cols))
} else {
  message("No Glomeromycota sequences found — skipping Glomeromycota clustering.")
  glom_clusters <- glom_seqs %>% select(any_of(standard_cols))
}

fwrite(glom_clusters, "./tmp_clusters/glomeromycota_clusters.txt", sep = "\t")
glom_otu_ids <- glom_clusters %>% pull(otu_id)

# (2) NON-GLOMEROMYCOTA: Filter and cluster at family level
nonglom_seqs <- taxa_cutoffs %>%
  filter(phylum != "Glomeromycota" | is.na(phylum) | phylum == "unidentified")

if (nrow(nonglom_seqs) > 0) {
  fwrite(nonglom_seqs, "./tmp_clusters/order_nonglom_clusters.txt", sep = "\t")
  cluster_at_rank("family", "genus", "all", "order_nonglom", taxa_cutoffs, asv_sequences, threads)
} else {
  message("No non-Glomeromycota sequences found — skipping non-Glomeromycota clustering.")
}

# ── Write final outputs ───────────────────────────────────────────────────────

# (2) Non-Glomeromycota clusters (only family-identified)
non_glom_clusters <- fread("./tmp_clusters/family_clusters.txt") %>%
  filter(!otu_id %in% glom_otu_ids) %>%
  select(any_of(standard_cols))

fwrite(non_glom_clusters, "./tmp_clusters/family_non_glomeromycota_clusters.txt", sep = "\t")
non_glom_otu_ids <- non_glom_clusters %>% pull(otu_id)

# (3) Unidentified clusters (everything else)
taxa_cutoffs %>%
  filter(!otu_id %in% c(glom_otu_ids, non_glom_otu_ids)) %>%
  select(any_of(standard_cols)) %>%
  fwrite("./tmp_clusters/family_unidentified_clusters.txt", sep = "\t")
