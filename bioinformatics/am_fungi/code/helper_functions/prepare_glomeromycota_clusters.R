# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
cutoff_file <- args[1]
classification_file <- args[2]
blast_species_file <- args[3]
blast_genus_file <- args[4]
blast_family_file <- args[5]
blast_all_file <- args[6]
asv_sequence_file <- args[7]
classification_output <- args[8]
minlen <- as.integer(args[9])
threads <- as.integer(args[10])

# Require packages
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(tidyverse))

# Helper functions -------------------------------------------------------------

# Format BLAST results
format_blast <- function(blast_file, classification_file, minlen) {
  fread(blast_file) %>%
    mutate(
      otu_id = V1,
      abundance = as.integer(str_extract(V1, "(?<=;size=)\\d+"))
    ) %>%
    select(
      otu_id, reference_id = V2, 
      pident = V3, length = V4, qlen = V5, slen = V6,
      abundance
    ) %>%
    mutate(
      subject_coverage = pmin(length / slen, 1.0),
      query_coverage = pmin(length / qlen, 1.0),
      score = pident / 100,
      score = if_else(length < minlen, (score * length) / minlen, score)
    ) %>%
    left_join(fread(classification_file), by = c("reference_id" = "id")) %>%
    select(otu_id, reference_id, kingdom, phylum, class, order, family, genus, species, 
           score, subject_coverage, query_coverage, abundance)
}

# Main classification function
classify_taxonomy <- function(blast_file, classification_file, cutoff_file, rank = "species",
                              minlen = 450, min_subject_coverage = 0.9, min_query_coverage = 0.9,
                              apply_coverage_filter = TRUE, excluded_otus = NULL) {
  
  valid_ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")
  if (!rank %in% valid_ranks) {
    stop("Invalid rank. Must be one of: ", paste(valid_ranks, collapse = ", "))
  }
  
  rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  rank_index <- which(rank_hierarchy == rank)
  higher_ranks <- rank_hierarchy[1:(rank_index - 1)]
  lower_ranks <- if (rank_index < length(rank_hierarchy)) {
    rank_hierarchy[(rank_index + 1):length(rank_hierarchy)]
  } else {
    character(0)
  }
  
  blast_data <- format_blast(blast_file, classification_file, minlen) %>%
    filter(!otu_id %in% excluded_otus)
  
  all_otu_ids <- distinct(blast_data, otu_id)
  cutoff_rank <- rank
  
  # Get cutoffs for OTUs with specific taxonomic matches
  if (length(higher_ranks) > 0) {
    cutoffs_matched <- blast_data %>%
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
  
  # Filter and classify hits
  hits <- blast_data
  if (apply_coverage_filter) {
    hits <- filter(hits, subject_coverage > min_subject_coverage & query_coverage > min_query_coverage)
  }
  
  result <- hits %>%
    left_join(all_cutoffs, by = "otu_id") %>%
    filter(score >= cutoff) %>%
    mutate(rank = !!rank) %>%
    group_by(otu_id) %>%
    arrange(desc(score), desc(query_coverage)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-c(subject_coverage, query_coverage))
  
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
  
  default_cutoff <- if(length(default_cutoff) > 0) default_cutoff else NA_real_
  
  taxa_file %>%
    mutate(!!paste0(rank, "_cutoff") := case_when(!!!conditions, TRUE ~ default_cutoff))
}

# Clustering function
cluster_at_rank <- function(this_rank, this_subrank, this_supertaxon, this_superrank = NULL, 
                           taxa_cutoffs, asv_sequences, minlen, threads) {
  
  # Initialise or filter taxa file
  if (!is.null(this_superrank)) {
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
      filter(!!sym(this_superrank) == this_supertaxon)
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
  rank_index <- which(rank_hierarchy == this_rank)
  lower_ranks <- if (rank_index < length(rank_hierarchy)) {
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
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]
    
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")
    
    message("Clustering ", initial_unidentified_count, " unidentified ASVs to ", 
            this_rank, " in ", this_supertaxon, " at ", this_cutoff * 100, "% similarity...")
    
    # Build BLAST database and run BLAST
    system2("makeblastdb", args = c("-in", "./tmp/identified_fasta", "-dbtype", "nucl"),
            stdout = FALSE, stderr = FALSE)
    system2("blastn", args = c(
      "-query", "./tmp/unidentified_fasta",
      "-db", "./tmp/identified_fasta",
      "-outfmt", "'6 qseqid sseqid pident length qlen slen'",
      "-task", "blastn-short",
      "-num_threads", as.character(threads),
      "-out", "./tmp/unidentified_blast.out"
    ), stdout = FALSE, stderr = FALSE)
    
    # Check if BLAST output is empty
    blast_output_size <- file.info("./tmp/unidentified_blast.out")$size
    
    if (is.na(blast_output_size) || blast_output_size == 0) {
      message("No BLAST hits found for ", this_supertaxon, " - stopping clustering")
      break
    }
    
    # Process BLAST results
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      setnames(c("otu_id", "reference_id", "sim", "alignment_length", "qlen", "slen")) %>%
      group_by(otu_id) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        query_coverage = if_else(qlen > 0, pmin(alignment_length / qlen, 1.0) * 100, 0),
        reference_coverage = if_else(slen > 0, pmin(alignment_length / slen, 1.0) * 100, 0),
        score = sim / 100,
        score = if_else(alignment_length < minlen, (score * alignment_length) / minlen, score),
        abundance = as.integer(str_extract(otu_id, "(?<=;size=)\\d+"))
      ) %>%
      filter(score >= this_cutoff, query_coverage >= 90 | reference_coverage >= 90)
    
    # Check if any sequences passed the filters
    if (nrow(new_clusters) == 0) {
      message("No sequences passed similarity/coverage filters for ", this_supertaxon)
      break
    }
    
    # Join with reference taxonomy FIRST (this is the key fix!)
    new_clusters <- new_clusters %>%
      left_join(
        taxa_cutoffs %>% 
          select(reference_id = otu_id, all_of(rank_hierarchy), ends_with("_cutoff")),
        by = "reference_id"
      ) %>%
      mutate(
        rank = this_rank,
        cutoff = this_cutoff
      )
    
    # NOW set lower ranks to unidentified (after inheriting from reference)
    if (length(lower_ranks) > 0) {
      for (lower_rank in lower_ranks) {
        new_clusters <- new_clusters %>%
          mutate(!!lower_rank := "unidentified")
      }
    }
    
    # Clean up columns before binding
    new_clusters <- new_clusters %>%
      select(-query_coverage, -reference_coverage, -sim, -alignment_length, -qlen, -slen)
    
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

# Define blast files for each rank ---------------------------------------------

# Use blast_all_file if not defined
blast_files <- list(
  species = if(exists("blast_species_file")) blast_species_file else blast_all_file,
  genus = if(exists("blast_genus_file")) blast_genus_file else blast_all_file,
  family = if(exists("blast_family_file")) blast_family_file else blast_all_file,
  order = if(exists("blast_order_file")) blast_order_file else blast_all_file,
  class = if(exists("blast_class_file")) blast_class_file else blast_all_file,
  phylum = if(exists("blast_phylum_file")) blast_phylum_file else blast_all_file,
  kingdom = blast_all_file
)

# Read cutoff file -------------------------------------------------------------

taxon_cutoffs <- fread(cutoff_file) %>%
  select(rank = Rank, taxa = Dataset, cutoff = "cut-off")

# Classify blast results -------------------------------------------------------

classified_otus <- character()
all_classifications <- list()

# Classify each rank (requires 90% coverage between query and subject)
ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")
coverage_filter <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)

for (i in seq_along(ranks)) {
  rank <- ranks[i]
  message("Classifying ", rank, "...")
  
  hits <- classify_taxonomy(
    blast_file = blast_files[[rank]],
    classification_file = classification_file,
    cutoff_file = taxon_cutoffs,
    rank = rank,
    minlen = minlen,
    apply_coverage_filter = coverage_filter[i],
    excluded_otus = classified_otus  # This prevents re-classifying already classified OTUs
  )
  
  # Only add OTUs that were successfully classified at this rank
  # (i.e., where the rank column is not "unidentified")
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

# Get unclassified ASVs --------------------------------------------------------

all_unclassified <- fread(blast_all_file) %>%
  distinct(otu_id = V1) %>%
  anti_join(all_classifications, by = "otu_id") %>%
  mutate(
    reference_id = "",
    kingdom = "unidentified",
    phylum = "unidentified",
    class = "unidentified",
    order = "unidentified",
    family = "unidentified",
    genus = "unidentified",
    species = "unidentified",
    score = as.numeric(NA),
    abundance = as.integer(str_extract(otu_id, "(?<=;size=)\\d+")),
    cutoff = as.numeric(NA),
    rank = ""
  )

# Save the complete classification file ----------------------------------------

all_classifications %>%
  bind_rows(all_unclassified) %>%
  select(otu_id, reference_id, kingdom, phylum, class, order, family, genus, species,
         rank, cutoff, score, abundance) %>%
  fwrite(classification_output, sep = "\t")

# Cluster unclassified ---------------------------------------------------------

# Create temporary directories
dir.create("tmp", showWarnings = FALSE)
dir.create("tmp_clusters", showWarnings = FALSE)

taxa_file <- bind_rows(all_classifications, all_unclassified)
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
  genus = c("family", "order", "class", "phylum", "kingdom"),
  family = c("order", "class", "phylum", "kingdom"),
  order = c("class", "phylum", "kingdom"),
  class = c("phylum", "kingdom"),
  phylum = c("kingdom"),
  kingdom = NULL
)

# Get cutoffs for all ranks
taxa_cutoffs <- taxa_file
for (rank in names(superranks_list)) {
  taxa_cutoffs <- get_rank_cutoff(taxa_cutoffs, unique_taxa_cutoffs, taxon_cutoffs, rank, superranks_list[[rank]])
}

# Run clustering ---------------------------------------------------------------

# (1) Eukarya clustering at kingdom level
cluster_at_rank("kingdom", "phylum", "Eukarya", NULL, taxa_cutoffs, asv_sequences, minlen, threads)

# (2) Fungi clustering at phylum level
cluster_at_rank("phylum", "class", "Fungi", "kingdom", taxa_cutoffs, asv_sequences, minlen, threads)

# (3) Glomeromycota clustering at class level
cluster_at_rank("class", "order", "Glomeromycota", "phylum", taxa_cutoffs, asv_sequences, minlen, threads)

# Write final outputs ----------------------------------------------------------

# Standard column order
standard_cols <- c("otu_id", "reference_id", "kingdom", "phylum", "class", "order", 
                   "family", "genus", "species", "rank", "cutoff", "score", "abundance",
                   "kingdom_cutoff", "phylum_cutoff", "class_cutoff", "order_cutoff", 
                   "family_cutoff", "genus_cutoff", "species_cutoff")

# Glomeromycota
fread("./tmp_clusters/class_clusters.txt") %>%
  select(any_of(standard_cols)) %>%
  fwrite("./tmp_clusters/glomeromycota_clusters.txt", sep = "\t")

# Non-glomeromycota
bind_rows(
  fread("./tmp_clusters/kingdom_unidentified.txt"),
  fread("./tmp_clusters/phylum_unidentified.txt"),
  fread("./tmp_clusters/phylum_clusters.txt"),
  fread("./tmp_clusters/phylum_unidentified.txt"),
  fread("./tmp_clusters/class_unidentified.txt")
) %>%
  select(any_of(standard_cols)) %>%
  # Remove any OTUs already classified as Glomeromycota
  anti_join(
    fread("./tmp_clusters/glomeromycota_clusters.txt") %>% select(otu_id),
    by = "otu_id"
  ) %>%
  fwrite("tmp_clusters/non_glomeromycota_clusters.txt", sep = "\t")

# Cleanup
unlink("tmp", recursive = TRUE)
file.remove(c(
  "tmp_clusters/kingdom_unidentified.txt", "tmp_clusters/phylum_unidentified.txt",
  "tmp_clusters/phylum_clusters.txt", "tmp_clusters/kingdom_clusters.txt",
  "tmp_clusters/class_clusters.txt", "tmp_clusters/class_unidentified.txt"
  ))

# Summary ----------------------------------------------------------------------
message("\nClassification complete! Check tmp_clusters/ for results.")
