#!/usr/bin/env Rscript
# 03_cluster_otus.R — Dynamic OTU clustering for Glomeromycota ASVs
#
# Reference-based clustering uses vsearch --usearch_global (global alignment).
#
# De novo clustering uses vsearch --allpairs_global to compute all pairwise
# similarities, then finds connected components of the threshold graph via
# union-find. This produces exact single-linkage clusters: sequences A and B
# are in the same cluster iff there is a chain A → x₁ → ... → B where every
# consecutive pair has global similarity ≥ the cutoff. The connected components
# of the threshold graph have this property by definition.
#
# Usage:
#   Rscript scripts/03_cluster_otus.R \
#     --asv_sequences    data/asv_sequences.fasta \
#     --asv_table        data/asv_table.txt \
#     --asv_classification data/asv_classification.txt \
#     --ref_classification data/ref_seqs/eukaryome_V4.classification \
#     --glom_clusters    tmp_clusters/glomeromycota_clusters.txt \
#     --threads          10

# ── Packages ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
  library(tidyverse)
})

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-s", "--asv_sequences"),
              type = "character", metavar = "FILE",
              help = "ASV FASTA file [required]"),
  make_option(c("-t", "--asv_table"),
              type = "character", metavar = "FILE",
              help = "ASV sample matrix file [required]"),
  make_option(c("-a", "--asv_classification"),
              type = "character", metavar = "FILE",
              help = "ASV classification file from 02_classify_asvs.sh [required]"),
  make_option(c("-r", "--ref_classification"),
              type = "character", metavar = "FILE",
              help = "Reference sequence classification file (for taxonomy string formatting) [required]"),
  make_option(c("-g", "--glom_clusters"),
              type = "character", default = "./tmp_clusters/glomeromycota_clusters.txt",
              metavar = "FILE",
              help = "Glomeromycota clusters file produced by 02_classify_asvs.sh [default: %default]"),
  make_option("--threads",
              type = "integer", default = 10L, metavar = "INT",
              help = "Number of threads for vsearch [default: %default]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage = "%prog [options]",
    description = paste(
      "Dynamic OTU clustering for Glomeromycota ASVs.",
      "Reference-based clustering uses vsearch --usearch_global.",
      "De novo clustering uses vsearch --allpairs_global + union-find (exact single-linkage)."
    )
  )
)

# ── Validation ────────────────────────────────────────────────────────────────

if (is.null(opt$asv_sequences))     stop("--asv_sequences is required.")
if (is.null(opt$asv_table))         stop("--asv_table is required.")
if (is.null(opt$asv_classification)) stop("--asv_classification is required.")
if (is.null(opt$ref_classification)) stop("--ref_classification is required.")

if (!file.exists(opt$glom_clusters)) {
  stop("Glomeromycota clusters file not found: ", opt$glom_clusters,
       "\nPlease run 02_classify_asvs.sh first.")
}

# ── Parameters ────────────────────────────────────────────────────────────────

asv_sequence_file           <- opt$asv_sequences
asv_sample_matrix_file      <- opt$asv_table
asv_classification_file     <- opt$asv_classification
ref_seq_classification_file <- opt$ref_classification
taxa_cutoffs_file           <- opt$glom_clusters
threads                     <- opt$threads

# ── Setup ─────────────────────────────────────────────────────────────────────

dir.create("tmp",          showWarnings = FALSE)
dir.create("tmp_clusters", showWarnings = FALSE)

# ── Read input files ──────────────────────────────────────────────────────────

# Reference sequence taxonomy (for generating taxonomy strings with rank prefixes)
ref_seq_classification <- fread(ref_seq_classification_file) %>%
  select(reference_id = id, kingdom, phylum, class, order, family, genus, species) %>%
  mutate(
    kingdom = if_else(!str_starts(kingdom, "k__"), paste0("k__", kingdom), kingdom),
    phylum  = if_else(!str_starts(phylum,  "p__"), paste0("p__", phylum),  phylum),
    class   = if_else(!str_starts(class,   "c__"), paste0("c__", class),   class),
    order   = if_else(!str_starts(order,   "o__"), paste0("o__", order),   order),
    family  = if_else(!str_starts(family,  "f__"), paste0("f__", family),  family),
    genus   = if_else(!str_starts(genus,   "g__"), paste0("g__", genus),   genus),
    species = if_else(!str_starts(species, "s__"), paste0("s__", species), species),
    species = str_replace_all(species, " ", "_"),
    reference_taxonomy = paste(kingdom, phylum, class, order, family, genus, species, sep = ";")
  ) %>%
  select(reference_id, reference_taxonomy)

asv_sequences     <- readDNAStringSet(asv_sequence_file)
asv_sample_matrix <- fread(asv_sample_matrix_file)
asv_classification <- fread(asv_classification_file) %>%
  left_join(ref_seq_classification, by = "reference_id")

# ── De novo clustering via union-find (exact single-linkage) ──────────────────
#
# Runs vsearch --allpairs_global on the supplied sequences to compute all
# pairwise global similarities, then builds connected components of the
# threshold graph using union-find with path halving and union by rank.
#
# Connected components of a threshold graph ARE single-linkage clusters:
# sequences are co-clustered iff they are connected by a chain of pairwise
# similarities each ≥ cutoff. This is the definition of SLC and of graph
# connectivity simultaneously, so no approximation is involved.
#
# Returns a data.frame with columns: otu_id, pseudo_name.

denovo_cluster <- function(seq_ids, cutoff, label, rank, asv_sequences, threads) {

  n <- length(seq_ids)

  # Single sequence — trivially its own cluster
  if (n == 1L) {
    return(data.frame(otu_id = seq_ids,
                      pseudo_name = paste0(label, "_pseudo_", rank, "_0001"),
                      stringsAsFactors = FALSE))
  }

  # Write unidentified sequences
  sub_seqs <- asv_sequences[names(asv_sequences) %in% seq_ids]
  writeXStringSet(sub_seqs, "./tmp/denovo_seqs.fasta")

  # Compute all pairwise global similarities
  out_file <- "./tmp/denovo_pairs.txt"
  system2("vsearch", args = c(
    "--allpairs_global", "./tmp/denovo_seqs.fasta",
    "--acceptall",
    "--userout",   out_file,
    "--userfields", "query+target+id",
    "--threads",   as.character(threads)
  ), stdout = FALSE, stderr = FALSE)

  # ── Union-Find with path halving and union by rank ─────────────────────────
  id_to_int <- setNames(seq_len(n), seq_ids)
  parent    <- seq_len(n)
  uf_rank   <- integer(n)

  uf_find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]   # path halving
      x <- parent[x]
    }
    x
  }

  uf_union <- function(a, b) {
    ra <- uf_find(a); rb <- uf_find(b)
    if (ra == rb) return()
    if (uf_rank[ra] < uf_rank[rb]) { tmp <- ra; ra <- rb; rb <- tmp }
    parent[rb] <<- ra
    if (uf_rank[ra] == uf_rank[rb]) uf_rank[ra] <<- uf_rank[ra] + 1L
  }

  # Add edges for all pairs meeting the cutoff
  pair_size <- file.info(out_file)$size
  if (!is.na(pair_size) && pair_size > 0L) {
    pairs <- fread(out_file, col.names = c("i", "j", "pident")) %>%
      filter(i != j) %>%
      mutate(score = round(pident / 100 / 0.005) * 0.005) %>%
      filter(score >= cutoff)

    for (k in seq_len(nrow(pairs))) {
      ai <- id_to_int[pairs$i[k]]
      bj <- id_to_int[pairs$j[k]]
      if (!is.na(ai) && !is.na(bj)) uf_union(ai, bj)
    }
  }

  # Build cluster assignments — number components in order of first appearance
  component  <- vapply(seq_len(n), uf_find, integer(1L))
  comp_order <- match(component, unique(component))

  data.frame(
    otu_id      = seq_ids,
    pseudo_name = paste0(label, "_pseudo_", rank, "_", sprintf("%04d", comp_order)),
    stringsAsFactors = FALSE
  )
}

# ── Reference-based + de novo clustering ─────────────────────────────────────
#
# For each unique supertaxon:
#   1. Iteratively cluster unidentified ASVs against identified references
#      using vsearch --usearch_global until no new assignments are made.
#   2. Apply de novo clustering (vsearch --allpairs_global + union-find) to
#      any ASVs that remain unidentified after the reference step.

cluster_with_denovo <- function(this_superrank, this_rank, this_subrank,
                                taxa_cutoffs, asv_sequences, threads) {

  message(paste0("Starting clustering for ", this_rank, " !!!"))

  # Get unique supertaxa and their cutoffs
  supertaxa_cutoffs <- taxa_cutoffs %>%
    select(!!sym(this_superrank), paste0(this_rank, "_cutoff")) %>%
    filter(!get(this_superrank) %in% c("unidentified", "")) %>%
    unique()

  # Check for NA cutoffs
  if (any(is.na(taxa_cutoffs[[paste0(this_rank, "_cutoff")]]))) {
    taxa_cutoffs %>% filter(is.na(!!sym(paste0(this_rank, "_cutoff")))) %>% print()
    stop("Error: ", this_rank, "_cutoff contains NA values. Exiting...")
  }

  # ── Loop over each supertaxon ──────────────────────────────────────────────
  for (i in seq_len(nrow(supertaxa_cutoffs))) {
    taxa_cutoffs    <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
    this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
    this_cutoff     <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]

    message(paste0("Starting clustering for ", this_supertaxon, " at rank ",
                   this_rank, " using cutoff ", this_cutoff, "..."))

    ### Reference-based clustering loop ###
    repeat {
      taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

      # Find identified ASVs at this rank (used as references)
      identified_asvs <- taxa_cutoffs %>%
        filter(!!sym(this_superrank) == this_supertaxon,
               !!sym(this_rank) != "unidentified", !!sym(this_rank) != "") %>%
        pull(otu_id)

      if (length(identified_asvs) == 0) {
        message("No identified ASVs to use as references for ", this_supertaxon)
        break
      }

      # Find unidentified ASVs at this rank
      unidentified_asvs <- taxa_cutoffs %>%
        filter(!!sym(this_superrank) == this_supertaxon,
               !!sym(this_rank) == "unidentified" | !!sym(this_rank) == "") %>%
        pull(otu_id)

      initial_unidentified_count <- length(unidentified_asvs)
      if (initial_unidentified_count == 0) {
        message("No unidentified ASVs remaining for ", this_supertaxon)
        break
      }

      # Write sequences to temp files
      identified_sequences   <- asv_sequences[names(asv_sequences) %in% identified_asvs]
      unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]
      writeXStringSet(identified_sequences,   "./tmp/identified_fasta")
      writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")

      message("Clustering ", initial_unidentified_count, " unidentified ASVs to ",
              this_supertaxon, " at ", this_cutoff * 100, "% similarity...")

      # Run vsearch --usearch_global
      system2("vsearch", args = c(
        "--usearch_global", "./tmp/unidentified_fasta",
        "--db",             "./tmp/identified_fasta",
        "--id",             "0.5",
        "--strand",         "both",
        "--maxaccepts",     "1",
        "--maxrejects",     "0",
        "--userout",        "./tmp/unidentified_vsearch.out",
        "--userfields",     "query+target+id",
        "--threads",        as.character(threads)
      ), stdout = FALSE, stderr = FALSE)

      # Check if vsearch output is empty
      vsearch_size <- file.info("./tmp/unidentified_vsearch.out")$size
      if (is.na(vsearch_size) || vsearch_size == 0) {
        message("No vsearch hits found for ", this_supertaxon, " - stopping reference clustering")
        break
      }

      # Process vsearch results
      rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
      rank_index     <- which(rank_hierarchy == this_rank)

      new_clusters <- fread("./tmp/unidentified_vsearch.out",
                            col.names = c("otu_id", "reference_id", "pident")) %>%
        group_by(otu_id) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(
          score     = round(pident / 100 / 0.005) * 0.005,
          abundance = as.integer(str_extract(otu_id, "(?<=;size=)\\d+"))
        ) %>%
        filter(score >= this_cutoff)

      if (nrow(new_clusters) == 0) {
        message("No sequences passed similarity filter for ", this_supertaxon)
        break
      }

      # Inherit taxonomy from matched reference
      new_clusters <- new_clusters %>%
        left_join(
          taxa_cutoffs %>%
            select(reference_id = otu_id, all_of(rank_hierarchy), ends_with("_cutoff")),
          by = "reference_id"
        ) %>%
        mutate(rank = this_rank, cutoff = this_cutoff)

      # Set lower ranks to unidentified
      if (rank_index < length(rank_hierarchy)) {
        lower_ranks <- rank_hierarchy[(rank_index + 1):length(rank_hierarchy)]
        for (lower_rank in lower_ranks) {
          new_clusters <- new_clusters %>% mutate(!!lower_rank := "unidentified")
        }
      }

      # Remove vsearch-specific columns before binding
      new_clusters <- new_clusters %>% select(-pident)

      # Update taxa table
      taxa_cutoffs <- bind_rows(
        new_clusters,
        taxa_cutoffs %>% filter(!otu_id %in% new_clusters$otu_id)
      )

      fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

      remaining_unidentified_count <- taxa_cutoffs %>%
        filter(!!sym(this_superrank) == this_supertaxon,
               !!sym(this_rank) == "unidentified" | !!sym(this_rank) == "") %>%
        nrow()

      message("Assigned ", initial_unidentified_count - remaining_unidentified_count,
              " ASVs. ", remaining_unidentified_count, " remain unidentified.")

      if (remaining_unidentified_count == 0) {
        message("No more unidentified ASVs for ", this_supertaxon)
        break
      } else if (remaining_unidentified_count == initial_unidentified_count) {
        message("No new assignments for ", this_supertaxon, " - stopping reference clustering")
        break
      } else {
        message("Continuing reference clustering for remaining ASVs in ", this_supertaxon)
      }
    }

    ### De novo clustering (exact single-linkage via union-find) ###
    remaining_unidentified_asvs <- taxa_cutoffs %>%
      filter(!!sym(this_superrank) == this_supertaxon,
             !!sym(this_rank) == "unidentified" | !!sym(this_rank) == "") %>%
      pull(otu_id)

    if (length(remaining_unidentified_asvs) == 0) {
      message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon))
    } else {
      message(paste0("Commencing de novo clustering of ", length(remaining_unidentified_asvs),
                     " unidentified ASVs in ", this_supertaxon, "..."))

      cluster_map <- denovo_cluster(
        seq_ids       = remaining_unidentified_asvs,
        cutoff        = this_cutoff,
        label         = this_supertaxon,
        rank          = this_rank,
        asv_sequences = asv_sequences,
        threads       = threads
      )

      taxa_cutoffs <- cluster_map %>%
        left_join(taxa_cutoffs, by = "otu_id") %>%
        mutate(!!sym(this_rank) := pseudo_name) %>%
        select(-pseudo_name) %>%
        bind_rows(taxa_cutoffs %>% filter(!otu_id %in% cluster_map$otu_id))

      fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

      message("De novo clustering complete for ", this_supertaxon, ": ",
              n_distinct(cluster_map$pseudo_name), " cluster(s) formed.")
    }
  }

  message(paste0("Clustering for ", this_rank, " complete!!!"))
  unlink(c("./tmp/identified_fasta", "./tmp/unidentified_fasta",
           "./tmp/unidentified_vsearch.out",
           "./tmp/denovo_seqs.fasta", "./tmp/denovo_pairs.txt"))
  return(taxa_cutoffs)
}

# ── Helper: propagate cutoffs within a rank ───────────────────────────────────

get_cutoff_for_rank <- function(taxa_file, rank, rank_cutoff_col) {
  taxa_file %>%
    group_by(!!sym(rank)) %>%
    mutate(!!paste0(rank, "_cutoff") := first(!!sym(rank_cutoff_col))) %>%
    ungroup()
}

# ── Main clustering workflow ──────────────────────────────────────────────────

taxa_cutoffs <- fread(taxa_cutoffs_file)

# Clustering hierarchy for Glomeromycota
clustering_hierarchy <- list(
  list(superrank = "family", rank = "genus",   subrank = "species"),
  list(superrank = "genus",  rank = "species",  subrank = NULL)
)

for (i in seq_along(clustering_hierarchy)) {
  level          <- clustering_hierarchy[[i]]
  this_superrank <- level$superrank
  this_rank      <- level$rank
  this_subrank   <- level$subrank

  message(paste0("\n### Processing rank: ", this_rank, " ###\n"))

  if (!is.null(this_subrank)) {
    rank_cutoff_col <- paste0(this_rank, "_cutoff")
    taxa_cutoffs <- get_cutoff_for_rank(taxa_cutoffs, this_rank, rank_cutoff_col)
  } else {
    if (!paste0(this_rank, "_cutoff") %in% names(taxa_cutoffs)) {
      stop("Species cutoffs not found in taxa_cutoffs")
    }
  }

  fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

  taxa_cutoffs <- cluster_with_denovo(
    this_superrank = this_superrank,
    this_rank      = this_rank,
    this_subrank   = this_subrank,
    taxa_cutoffs   = taxa_cutoffs,
    asv_sequences  = asv_sequences,
    threads        = threads
  )
}

# ── Final OTU clusters ────────────────────────────────────────────────────────

pre_clustered_otus <- fread("./tmp_clusters/species_clusters.txt") %>%
  filter(family != "unidentified" & family != "") %>%
  select(
    otu_id, reference_id, kingdom, phylum, class, order, family, genus, species,
    rank, cutoff, score, asv_abundance = abundance
  )

# Rename pseudo taxa at each rank
pre_clustered_otus <- pre_clustered_otus %>%
  group_by(phylum) %>%
  mutate(
    phylum = case_when(
      str_detect(phylum, "_pseudo_") ~ paste0(
        str_extract(phylum, "^[^_]+"), "_pseudo_phylum_", sprintf("%04d", cur_group_id())),
      str_detect(phylum, "unidentified") ~ paste0(
        str_extract(kingdom, "^[^_]+"), "_pseudo_phylum_", sprintf("%04d", cur_group_id())),
      (str_detect(phylum, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(phylum, "Incerate") & str_detect(phylum, "unidentified")) ~ paste0(
          str_extract(kingdom, "^[^_]+"), "_pseudo_phylum_", sprintf("%04d", cur_group_id())),
      TRUE ~ phylum
    )
  ) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(
    class = case_when(
      str_detect(class, "_pseudo_") ~ paste0(
        str_extract(class, "^[^_]+"), "_pseudo_class_", sprintf("%04d", cur_group_id())),
      str_detect(class, "unidentified") ~ paste0(
        str_extract(phylum, "^[^_]+"), "_pseudo_class_", sprintf("%04d", cur_group_id())),
      (str_detect(class, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(class, "Incerate") & str_detect(class, "unidentified")) ~ paste0(
          str_extract(phylum, "^[^_]+"), "_pseudo_class_", sprintf("%04d", cur_group_id())),
      TRUE ~ class
    )
  ) %>%
  ungroup() %>%
  group_by(order) %>%
  mutate(
    order = case_when(
      str_detect(order, "_pseudo_") ~ paste0(
        str_extract(order, "^[^_]+"), "_pseudo_order_", sprintf("%04d", cur_group_id())),
      str_detect(order, "unidentified") ~ paste0(
        str_extract(class, "^[^_]+"), "_pseudo_order_", sprintf("%04d", cur_group_id())),
      (str_detect(order, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(order, "Incerate") & str_detect(order, "unidentified")) ~ paste0(
          str_extract(class, "^[^_]+"), "_pseudo_order_", sprintf("%04d", cur_group_id())),
      TRUE ~ order
    )
  ) %>%
  ungroup() %>%
  group_by(family) %>%
  mutate(
    family = case_when(
      str_detect(family, "_pseudo_") ~ paste0(
        str_extract(family, "^[^_]+"), "_pseudo_family_", sprintf("%04d", cur_group_id())),
      str_detect(family, "unidentified") ~ paste0(
        str_extract(order, "^[^_]+"), "_pseudo_family_", sprintf("%04d", cur_group_id())),
      (str_detect(family, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(family, "Incerate") & str_detect(family, "unidentified")) ~ paste0(
          str_extract(order, "^[^_]+"), "_pseudo_family_", sprintf("%04d", cur_group_id())),
      TRUE ~ family
    )
  ) %>%
  ungroup() %>%
  group_by(genus) %>%
  mutate(
    genus = case_when(
      str_detect(genus, "_pseudo_") ~ paste0(
        str_extract(genus, "^[^_]+"), "_pseudo_genus_", sprintf("%04d", cur_group_id())),
      str_detect(genus, "unidentified") ~ paste0(
        str_extract(family, "^[^_]+"), "_pseudo_genus_", sprintf("%04d", cur_group_id())),
      (str_detect(genus, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(genus, "Incerate") & str_detect(genus, "unidentified")) ~ paste0(
          str_extract(family, "^[^_]+"), "_pseudo_genus_", sprintf("%04d", cur_group_id())),
      TRUE ~ genus
    )
  ) %>%
  ungroup() %>%
  group_by(species) %>%
  mutate(
    species = case_when(
      str_detect(species, "_pseudo_") ~ paste0(
        str_extract(species, "^[^_]+"), "_pseudo_species_", sprintf("%04d", cur_group_id())),
      str_detect(species, "unidentified") ~ paste0(
        str_extract(genus, "^[^_]+"), "_pseudo_species_", sprintf("%04d", cur_group_id())),
      (str_detect(species, "Incerate") & str_detect(species, "_pseudo_")) |
        (str_detect(species, "Incerate") & str_detect(species, "unidentified")) ~ paste0(
          str_extract(genus, "^[^_]+"), "_pseudo_species_", sprintf("%04d", cur_group_id())),
      TRUE ~ species
    )
  ) %>%
  ungroup() %>%
  mutate(asv_id = otu_id) %>%
  arrange(species) %>%
  group_by(species) %>%
  mutate(
    otu_abundance = sum(asv_abundance),
    # Prioritise de novo cluster cores over reference-based assignments.
    # Reference IDs from the database start with uppercase letters; de novo
    # cluster representatives (ASV IDs) start with lowercase hex or digits.
    score_temp = case_when(
      str_detect(reference_id, "^[0-9a-z]") ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(desc(score_temp), desc(asv_abundance), desc(score)) %>%
  mutate(
    otu_id       = first(otu_id),
    reference_id = first(reference_id),
    score        = first(score),
    cutoff       = first(cutoff)
  ) %>%
  ungroup() %>%
  arrange(desc(otu_abundance)) %>%
  select(
    otu_id, asv_id, reference_id,
    class, order, family, genus, species,
    rank, score, cutoff, abundance = otu_abundance
  ) %>%
  print(n = 100)

# ── Format and save OTU classification ───────────────────────────────────────

otu_classification <- pre_clustered_otus %>%
  filter(rank %in% c("family", "genus", "species")) %>%
  group_by(otu_id) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(desc(abundance)) %>%
  select(-reference_id) %>%
  left_join(
    asv_classification %>% select(otu_id, reference_id, reference_taxonomy),
    by = "otu_id"
  ) %>%
  select(otu_id, reference_id, reference_taxonomy,
         class, order, family, genus, species, rank, score, abundance) %>%
  mutate(
    reference_id = if_else(
      is.na(reference_id) | reference_id == "",
      pre_clustered_otus$reference_id[match(otu_id, pre_clustered_otus$otu_id)],
      reference_id
    ),
    reference_taxonomy = if_else(
      is.na(reference_taxonomy) | reference_taxonomy == "",
      "closed_reference_based_cluster",
      reference_taxonomy
    ),
    otu_id = str_replace(otu_id, ";size=\\d+$", "")
  )

# Check for empty taxonomies
otu_classification %>%
  filter(is.na(class) | class == "" | is.na(order) | order == "" |
         is.na(family) | family == "" | is.na(genus) | genus == "" |
         is.na(species) | species == "")

# ── Renumber pseudo taxa sequentially within each rank by abundance ───────────

ranks_to_renumber <- c("class", "order", "family", "genus", "species")

for (this_rank in ranks_to_renumber) {
  pseudo_pattern <- paste0("_pseudo_", this_rank, "_")
  otu_classification <- otu_classification %>%
    mutate(
      is_pseudo    = str_detect(!!sym(this_rank), pseudo_pattern),
      pseudo_prefix = if_else(
        is_pseudo,
        str_extract(!!sym(this_rank), paste0("^.+_pseudo_", this_rank, "_")),
        NA_character_
      )
    ) %>%
    group_by(pseudo_prefix) %>%
    arrange(desc(abundance)) %>%
    mutate(
      new_number = if_else(is_pseudo, row_number(), NA_integer_),
      !!sym(this_rank) := if_else(
        is_pseudo,
        paste0(pseudo_prefix, sprintf("%04d", new_number)),
        !!sym(this_rank)
      )
    ) %>%
    ungroup() %>%
    select(-is_pseudo, -pseudo_prefix, -new_number)
}

# ── Format and save OTU sample matrix ────────────────────────────────────────

otu_sample_matrix <- asv_sample_matrix %>%
  rename(asv_id = OTU_ID) %>%
  pivot_longer(cols = -asv_id, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  inner_join(
    pre_clustered_otus %>%
      filter(rank %in% c("family", "genus", "species")) %>%
      mutate(
        asv_id = str_replace(asv_id, ";size=\\d+$", ""),
        otu_id = str_replace(otu_id, ";size=\\d+$", "")
      ) %>%
      select(asv_id, otu_id) %>%
      unique(),
    by = "asv_id"
  ) %>%
  group_by(otu_id, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  arrange(desc(abundance)) %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0)

# ── Consistency checks ────────────────────────────────────────────────────────

message(paste0("Final number of OTUs: ", nrow(otu_sample_matrix)))
message(paste0("Final number of OTU classifications: ", nrow(otu_classification)))

if (nrow(otu_sample_matrix) != nrow(otu_classification)) {
  stop("Mismatch in number of OTUs between sample matrix and classification")
}

otu_ids_matrix         <- sort(otu_sample_matrix$otu_id)
otu_ids_classification <- sort(otu_classification$otu_id)

missing_in_classification <- setdiff(otu_ids_matrix, otu_ids_classification)
if (length(missing_in_classification) > 0) {
  message("OTU IDs in sample matrix but NOT in classification:")
  print(missing_in_classification)
}

missing_in_matrix <- setdiff(otu_ids_classification, otu_ids_matrix)
if (length(missing_in_matrix) > 0) {
  message("OTU IDs in classification but NOT in sample matrix:")
  print(missing_in_matrix)
}

if (!all(otu_ids_matrix %in% otu_ids_classification)) {
  stop("Some OTU IDs in sample matrix do not match classification")
}
if (!all(otu_ids_classification %in% otu_ids_matrix)) {
  stop("Some OTU IDs in classification do not match sample matrix")
}
message("All OTU IDs match between sample matrix and classification!")

# ── OTU representative sequences ─────────────────────────────────────────────

otu_representative_names <- unique(
  pre_clustered_otus %>%
    filter(rank %in% c("family", "genus", "species")) %>%
    pull(otu_id)
)

otu_representative_sequences <- asv_sequences[names(asv_sequences) %in% otu_representative_names]
names(otu_representative_sequences) <- str_replace(
  names(otu_representative_sequences), ";size=\\d+$", ""
)

message(paste0("Final number of OTU representative sequences: ", length(otu_representative_sequences)))

if (length(otu_representative_sequences) != nrow(otu_classification)) {
  stop("Mismatch in number of OTU representative sequences and classification")
}

otu_rep_names   <- sort(names(otu_representative_sequences))
otu_class_names <- sort(otu_classification$otu_id)

if (!all(otu_rep_names %in% otu_class_names)) {
  stop("Some OTU IDs in representative sequences do not match classification")
}
if (!all(otu_class_names %in% otu_rep_names)) {
  stop("Some OTU IDs in classification do not match representative sequences")
}
message("All OTU IDs match between representative sequences and classification!")

# ── Save Glomeromycota outputs ────────────────────────────────────────────────

fwrite(otu_sample_matrix, "./tmp_clusters/otu_glom_table.txt", sep = "\t")
fwrite(otu_classification, "./tmp_clusters/otu_glom_classification.txt", sep = "\t")
writeXStringSet(otu_representative_sequences, "./tmp_clusters/otu_glom_sequences.fasta")

# ── Save final combined outputs ───────────────────────────────────────────────

non_glomeromycota_clusters <- fread("./tmp_clusters/family_non_glomeromycota_clusters.txt") %>%
  mutate(otu_id = str_replace(otu_id, ";size=\\d+$", ""))

non_glomeromycota_clusters %>%
  bind_rows(
    fread("./tmp_clusters/species_clusters.txt") %>%
      filter(family == "unidentified" | family == "") %>%
      mutate(otu_id = str_replace(otu_id, ";size=\\d+$", ""))
  ) %>%
  bind_rows(
    fread("./tmp_clusters/family_unidentified_clusters.txt") %>%
      mutate(otu_id = str_replace(otu_id, ";size=\\d+$", ""))
  ) %>%
  fwrite("../output/otu_classification_non_glom.txt", sep = "\t")

# Family-level OTU tables
non_glomeromycota_family_clusters <- asv_sample_matrix %>%
  rename(otu_id = OTU_ID) %>%
  pivot_longer(cols = -otu_id, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  inner_join(non_glomeromycota_clusters %>% select(otu_id, family), by = "otu_id") %>%
  group_by(family, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  arrange(desc(abundance))

glomeromycota_family_clusters <- otu_sample_matrix %>%
  pivot_longer(cols = -otu_id, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  inner_join(otu_classification %>% select(otu_id, family), by = "otu_id") %>%
  group_by(family, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  arrange(desc(abundance))

bind_rows(non_glomeromycota_family_clusters, glomeromycota_family_clusters) %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
  fwrite("../output/otu_table_all_family.txt", sep = "\t")

glomeromycota_family_clusters %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
  fwrite("../output/otu_table_glom_family.txt", sep = "\t")
