
# Run from the source code directory:
setwd("./code")

# Required files:
asv_sequence_file <- "../data/ASVs.fasta"
asv_sample_matrix_file <- "../data/ASVs.txt"
asv_classification_file <- "../data/asv_classification.txt"
ref_seq_classification_file <- "./dnabarcoder/data/ref_seqs_V4.classification"
taxa_cutoffs_file <- "./tmp_clusters/glomeromycota_clusters.txt"

# Check if the glomeromycota clusters file exists
if (!file.exists(taxa_cutoffs_file)) {
  stop("glomeromycota_clusters.txt not found in ./tmp_clusters/. Please run the phylum clustering step first.")
}

# Constants
minlen <- 400
threads <- 10

# Require packages
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(tidyverse))

# Create temporary directory for blast db
dir.create("tmp", showWarnings = FALSE)

# Read reference sequence taxonomy
ref_seq_classification <- fread(ref_seq_classification_file) %>%
  select(reference_id = id, kingdom, phylum, class, order, family, genus, species) %>%
  mutate(
    # Add k__ for kingdom, p__ for phylum, etc.
    kingdom = if_else(!str_starts(kingdom, "k__"), paste0("k__", kingdom), kingdom),
    phylum = if_else(!str_starts(phylum, "p__"), paste0("p__", phylum), phylum),
    class = if_else(!str_starts(class, "c__"), paste0("c__", class), class),
    order = if_else(!str_starts(order, "o__"), paste0("o__", order), order),
    family = if_else(!str_starts(family, "f__"), paste0("f__", family), family),
    genus = if_else(!str_starts(genus, "g__"), paste0("g__", genus), genus),
    species = if_else(!str_starts(species, "s__"), paste0("s__", species), species),
    # Add and underscore in species names
    species = str_replace_all(species, " ", "_"),
    # Separate ranks by ;
    reference_taxonomy = paste(kingdom, phylum, class, order, family, genus, species, sep = ";")
  ) %>%
  select(reference_id, reference_taxonomy)

# Read ASVs
asv_sequences <- readDNAStringSet(asv_sequence_file)
asv_sample_matrix <- fread(asv_sample_matrix_file)
asv_classification <- fread(asv_classification_file) %>%
  # Add reference taxonomy
  left_join(ref_seq_classification, by = "reference_id")

# Generalised clustering function with de novo capability ----------------------

cluster_with_denovo <- function(this_superrank, this_rank, this_subrank, 
                                taxa_cutoffs, asv_sequences, minlen, threads) {
  
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
  
  # Loop over each supertaxon
  for (i in 1:nrow(supertaxa_cutoffs)) {
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
    this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
    this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]
    
    message(paste0("Starting clustering for ", this_supertaxon, " at rank ", 
                   this_rank, " using cutoff ", this_cutoff, "..."))
    
    ### Reference-based clustering loop ###
    repeat {
      taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
      
      # Find identified ASVs at this rank
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
      
      # Filter sequences
      identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
      unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]
      
      writeXStringSet(identified_sequences, "./tmp/identified_fasta")
      writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")
      
      message("Clustering ", initial_unidentified_count, " unidentified ASVs to ", 
              this_supertaxon, " at ", this_cutoff * 100, "% similarity...")
      
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
          score = round(score / 0.005) * 0.005,
          abundance = as.integer(str_extract(otu_id, "(?<=;size=)\\d+"))
        ) %>%
        filter(score >= this_cutoff, query_coverage >= 95 | reference_coverage >= 95)
      
      # Check if any sequences passed the filters
      if (nrow(new_clusters) == 0) {
        message("No sequences passed similarity/coverage filters for ", this_supertaxon)
        break
      }
      
      # Define rank hierarchy
      rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
      rank_index <- which(rank_hierarchy == this_rank)
      
      # Join with reference taxonomy to inherit ALL taxonomic ranks
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
      
      # Set lower ranks to unidentified (after inheriting from reference)
      if (rank_index < length(rank_hierarchy)) {
        lower_ranks <- rank_hierarchy[(rank_index + 1):length(rank_hierarchy)]
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
      
      fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
      
      # Re-count remaining unidentified after update
      remaining_unidentified_count <- taxa_cutoffs %>%
        filter(!!sym(this_superrank) == this_supertaxon,
               !!sym(this_rank) == "unidentified" | !!sym(this_rank) == "") %>%
        nrow()
      
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
    
    ### De novo clustering ###
    remaining_unidentified_asvs <- taxa_cutoffs %>%
      filter(!!sym(this_superrank) == this_supertaxon,
             !!sym(this_rank) == "unidentified" | !!sym(this_rank) == "") %>%
      pull(otu_id)
    
    if (length(remaining_unidentified_asvs) == 0) {
      message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon))
    } else {
      message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))
      
      # Write unidentified sequences
      unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
      writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")
      
      # Run blastclust
      identity_cutoff <- this_cutoff * 100
      system2("blastclust", args = c(
        "-i", paste0("./tmp/remaining_unidentified.fasta"),
        "-S", as.character(identity_cutoff),
        "-a", as.character(threads),
        "-p", "F",
        "-o", paste0("./tmp/de_novo_clusters.txt")
      ))
      
      # Check if blastclust output is empty
      blastclust_output_size <- file.info("./tmp/de_novo_clusters.txt")$size
      
      if (is.na(blastclust_output_size) || blastclust_output_size == 0) {
        message(paste0("No clusters formed by blastclust for ", this_supertaxon, 
                       " - sequences too divergent or too few"))
        next  # Skip to next supertaxon
      }
      
      # Process de novo clusters
      pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, 
                               sep = "\n", col.names = "cluster")
      
      taxa_cutoffs <- pseudo_clusters %>%
        mutate(ASVs = str_split(cluster, " ")) %>%
        unnest(ASVs) %>%
        group_by(cluster) %>%
        mutate(pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", 
                                    sprintf("%04d", cur_group_id()))) %>%
        ungroup() %>%
        mutate(otu_id = ASVs) %>%
        select(otu_id, pseudo_name) %>%
        full_join(taxa_cutoffs, by = "otu_id") %>%
        mutate(!!sym(this_rank) := if_else(!is.na(pseudo_name), pseudo_name, 
                                           !!sym(this_rank))) %>%
        select(-pseudo_name)
      
      fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
    }
  }
  
  message(paste0("Clustering for ", this_rank, " complete!!!"))
  
  # Cleanup temporary files
  unlink(list.files("tmp", full.names = TRUE))
  
  return(taxa_cutoffs)
}

# Helper function to get cutoff for a rank ---------------------------------

get_cutoff_for_rank <- function(taxa_file, rank, rank_cutoff_col) {
  taxa_file %>%
    group_by(!!sym(rank)) %>%
    mutate(!!paste0(rank, "_cutoff") := first(!!sym(rank_cutoff_col))) %>%
    ungroup()
}

# Main clustering workflow -------------------------------------------------

# Read the phylum clusters file (output from previous step)
taxa_cutoffs <- fread(taxa_cutoffs_file)

# Define the clustering hierarchy for Glomeromycota
clustering_hierarchy <- list(
  #list(superrank = "phylum", rank = "class", subrank = "order"),
  #list(superrank = "class", rank = "order", subrank = "family"),
  #list(superrank = "order", rank = "family", subrank = "genus"),
  list(superrank = "family", rank = "genus", subrank = "species"),
  list(superrank = "genus", rank = "species", subrank = NULL)  # species has no subrank
)

# Process each rank in the hierarchy
for (i in seq_along(clustering_hierarchy)) {
  level <- clustering_hierarchy[[i]]
  this_superrank <- level$superrank
  this_rank <- level$rank
  this_subrank <- level$subrank
  
  message(paste0("\n### Processing rank: ", this_rank, " ###\n"))
  
  # Update cutoffs using subrank cutoffs
  if (!is.null(this_subrank)) {
    rank_cutoff_col <- paste0(this_rank, "_cutoff")
    taxa_cutoffs <- get_cutoff_for_rank(taxa_cutoffs, this_rank, rank_cutoff_col)
  } else {
    # For species (no subrank), cutoffs should already exist
    if (!paste0(this_rank, "_cutoff") %in% names(taxa_cutoffs)) {
      stop("Species cutoffs not found in taxa_cutoffs")
    }
  }
  
  # Write the initial file for this rank
  fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  
  # Run clustering
  taxa_cutoffs <- cluster_with_denovo(
    this_superrank = this_superrank,
    this_rank = this_rank,
    this_subrank = this_subrank,
    taxa_cutoffs = taxa_cutoffs,
    asv_sequences = asv_sequences,
    minlen = minlen,
    threads = threads
  )
}


# Final OTU clusters -----------------------------------------------------------

pre_clustered_otus <- fread("./tmp_clusters/species_clusters.txt") %>%
  filter(family != "unidentified" & family != "") %>%
  select(
    otu_id, reference_id, kingdom, phylum, class, order, family, genus, species,
    rank, cutoff, score, asv_abundance = abundance)

# Rename pseudo taxa for each rank
pre_clustered_otus <- pre_clustered_otus %>%
  group_by(phylum) %>%
  mutate(
    phylum = case_when(
      str_detect(phylum, "_pseudo_") ~ paste0(
        str_extract(phylum, "^[^_]+"), 
        "_pseudo_phylum_", 
        sprintf("%04d", cur_group_id())
      ),
      str_detect(phylum, "unidentified") ~ paste0(
        str_extract(kingdom, "^[^_]+"), 
        "_pseudo_phylum_", 
        sprintf("%04d", cur_group_id())
      ),
      (str_detect(phylum, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(phylum, "Incerate") & str_detect(phylum, "unidentified")) ~ paste0(
          str_extract(kingdom, "^[^_]+"), 
          "_pseudo_phylum_", 
          sprintf("%04d", cur_group_id())
        ),
      TRUE ~ phylum
    )
  ) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(
    class = case_when(
      str_detect(class, "_pseudo_") ~ paste0(
        str_extract(class, "^[^_]+"), 
        "_pseudo_class_", 
        sprintf("%04d", cur_group_id())
      ),
      str_detect(class, "unidentified") ~ paste0(
        str_extract(phylum, "^[^_]+"), 
        "_pseudo_class_", 
        sprintf("%04d", cur_group_id())
      ),
      (str_detect(class, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(class, "Incerate") & str_detect(class, "unidentified")) ~ paste0(
          str_extract(phylum, "^[^_]+"), 
          "_pseudo_class_", 
          sprintf("%04d", cur_group_id())
        ),
      TRUE ~ class
    )
  ) %>%
  ungroup() %>%
  group_by(order) %>%
  mutate(
    order = case_when(
      str_detect(order, "_pseudo_") ~ paste0(
        str_extract(order, "^[^_]+"), 
        "_pseudo_order_", 
        sprintf("%04d", cur_group_id())
      ),
      str_detect(order, "unidentified") ~ paste0(
        str_extract(class, "^[^_]+"), 
        "_pseudo_order_", 
        sprintf("%04d", cur_group_id())
      ),
      (str_detect(order, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(order, "Incerate") & str_detect(order, "unidentified")) ~ paste0(
          str_extract(class, "^[^_]+"), 
          "_pseudo_order_", 
          sprintf("%04d", cur_group_id())
        ),
      TRUE ~ order
    )
  ) %>%
  ungroup() %>%
  group_by(family) %>%
  mutate(
    family = case_when(
      str_detect(family, "_pseudo_") ~ paste0(
        str_extract(family, "^[^_]+"), 
        "_pseudo_family_", 
        sprintf("%04d", cur_group_id())
      ),
      str_detect(family, "unidentified") ~ paste0(
        str_extract(order, "^[^_]+"), 
        "_pseudo_family_", 
        sprintf("%04d", cur_group_id())
      ),
      (str_detect(family, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(family, "Incerate") & str_detect(family, "unidentified")) ~ paste0(
          str_extract(order, "^[^_]+"), 
          "_pseudo_family_", 
          sprintf("%04d", cur_group_id())
        ),
      TRUE ~ family
    )
  ) %>%
  ungroup() %>%
  group_by(genus) %>%
  mutate(
    genus = case_when(
      str_detect(genus, "_pseudo_") ~ paste0(
        str_extract(genus, "^[^_]+"), 
        "_pseudo_genus_", 
        sprintf("%04d", cur_group_id())
      ),
      str_detect(genus, "unidentified") ~ paste0(
        str_extract(family, "^[^_]+"), 
        "_pseudo_genus_", 
        sprintf("%04d", cur_group_id())
      ),
      (str_detect(genus, "Incerate") & str_detect(genus, "_pseudo_")) |
        (str_detect(genus, "Incerate") & str_detect(genus, "unidentified")) ~ paste0(
          str_extract(family, "^[^_]+"), 
          "_pseudo_genus_", 
          sprintf("%04d", cur_group_id())
        ),
      TRUE ~ genus
    )
  ) %>%
  ungroup() %>%
  group_by(species) %>%
  mutate(
    species = case_when(
      str_detect(species, "_pseudo_") ~ paste0(
        str_extract(species, "^[^_]+"), 
        "_pseudo_species_", 
        sprintf("%04d", cur_group_id())
      ),
      str_detect(species, "unidentified") ~ paste0(
        str_extract(genus, "^[^_]+"), 
        "_pseudo_species_", 
        sprintf("%04d", cur_group_id())
      ),
      (str_detect(species, "Incerate") & str_detect(species, "_pseudo_")) |
        (str_detect(species, "Incerate") & str_detect(species, "unidentified")) ~ paste0(
          str_extract(genus, "^[^_]+"), 
          "_pseudo_species_", 
          sprintf("%04d", cur_group_id())
        ),
      TRUE ~ species
    )
  ) %>%
  ungroup() %>%
  # Generate unique OTU_IDs and sum the abundances
  mutate(
    asv_id = otu_id
  ) %>%
  arrange(species) %>%
  group_by(species) %>%
  mutate(
    otu_abundance = sum(asv_abundance),
    # Prioritise cluster core sequences over reference-based clusters.
    # Reference-based clusters have reference_id starting lowercase letters
    # or numbers, db reference_id start with uppercase letters.
    # Assign a score of 0 to reference_ids starting with a number or lowercase
    # letter to prevent them from being selected as representative sequences.
    score_temp = case_when(
      str_detect(reference_id, "^[0-9a-z]") ~ 0,
      TRUE ~ 1
    )
  ) %>%
  # Select representative ASV sequence based on the most abundant ASVs followed by best score to reference
  arrange(desc(score_temp), desc(asv_abundance), desc(score)) %>%
  mutate(
    otu_id = first(otu_id),
    reference_id = first(reference_id),
    score = first(score),
    cutoff = first(cutoff)
  ) %>%
  ungroup() %>%
  arrange(desc(otu_abundance)) %>%
  select(
    otu_id, asv_id, reference_id, 
    class, order, family, genus, species,
    rank, score, cutoff, abundance = otu_abundance
  ) %>%
  print(n = 100)

# Format and save OTU classification file --------------------------------------

# Format OTU classification table
otu_classification <- pre_clustered_otus %>%
  filter(rank %in% c("family", "genus", "species")) %>%
  # Select otu representative sequences
  group_by(otu_id) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(desc(abundance)) %>%
  select(-reference_id) %>%
  # Add reference taxonomy
  left_join(
    asv_classification %>%
      select(otu_id, reference_id, reference_taxonomy),
    by = "otu_id"
  ) %>%
  select(otu_id, reference_id, reference_taxonomy, class, order, family, genus, species, rank, score, abundance) %>%
  # When reference_id is blank, get the reference_id from pre_clustered_otus
  # that corresponds to this otu_id, and set reference_taxonomy to "closed_reference_based_cluster"
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
    # Remove size information from otu_id
    otu_id = str_replace(otu_id, ";size=\\d+$", "")
  )

# Check for empty taxonomies
otu_classification %>%
  filter(is.na(class) | class == "" |
           is.na(order) | order == "" |
           is.na(family) | family == "" |
           is.na(genus) | genus == "" |
           is.na(species) | species == ""
  )


# Renumber pseudo taxa sequentially within each rank based on abundance
ranks_to_renumber <- c("class", "order", "family", "genus", "species")

for (this_rank in ranks_to_renumber) {
  # Filter to only pseudo taxa at this rank
  pseudo_pattern <- paste0("_pseudo_", this_rank, "_")
  
  otu_classification <- otu_classification %>%
    mutate(
      # Check if this is a pseudo taxon
      is_pseudo = str_detect(!!sym(this_rank), pseudo_pattern),
      # Extract the prefix (e.g., "Paraglomeraceae_pseudo_genus_")
      pseudo_prefix = if_else(
        is_pseudo,
        str_extract(!!sym(this_rank), paste0("^.+_pseudo_", this_rank, "_")),
        NA_character_
      )
    ) %>%
    # Group by the pseudo prefix and renumber by abundance
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

# Format and save OTU sample matrix --------------------------------------------

# Format OTU sample matrix
otu_sample_matrix <- asv_sample_matrix %>%
  rename(asv_id = OTU_ID) %>%
  pivot_longer(
    cols = -asv_id,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  inner_join(
    pre_clustered_otus %>%
      filter(rank %in% c("family", "genus", "species")) %>%
      # Remove size information from asv_id and otu_id
      mutate(
        asv_id = str_replace(asv_id, ";size=\\d+$", ""),
        otu_id = str_replace(otu_id, ";size=\\d+$", "")
      ) %>%
      select(asv_id, otu_id) %>%
      unique(),  # Get unique asv_id to otu_id mappings
    by = "asv_id"
  ) %>%
  group_by(otu_id, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  ungroup() %>%
  arrange(desc(abundance)) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = abundance,
    values_fill = 0
  )

# Check the number of OTUs in the OTU sample matrix and classification
message(paste0("Final number of OTUs: ", nrow(otu_sample_matrix)))
message(paste0("Final number of OTU classifications: ", nrow(otu_classification)))

if (nrow(otu_sample_matrix) != nrow(otu_classification)) {
  stop("Mismatch in number of OTUs between sample matrix and classification")
}

# Find OTU IDs that don't match
otu_ids_matrix <- sort(otu_sample_matrix$otu_id)
otu_ids_classification <- sort(otu_classification$otu_id)

# OTUs in sample matrix but NOT in classification
missing_in_classification <- setdiff(otu_ids_matrix, otu_ids_classification)
if (length(missing_in_classification) > 0) {
  message("OTU IDs in sample matrix but NOT in classification:")
  print(missing_in_classification)
}

# OTUs in classification but NOT in sample matrix
missing_in_matrix <- setdiff(otu_ids_classification, otu_ids_matrix)
if (length(missing_in_matrix) > 0) {
  message("OTU IDs in classification but NOT in sample matrix:")
  print(missing_in_matrix)
}

# Check if all OTU IDs match (both directions)
if (!all(otu_ids_matrix %in% otu_ids_classification)) {
  stop("Some OTU IDs in sample matrix do not match classification")
}

if (!all(otu_ids_classification %in% otu_ids_matrix)) {
  stop("Some OTU IDs in classification do not match sample matrix")
}

message("All OTU IDs match between sample matrix and classification!")

# Format and OTU representative sequences --------------------------------------

# Get representative sequences
otu_representative_names <- unique(
  pre_clustered_otus %>%
    filter(rank %in% c("family", "genus", "species")) %>%
    pull(otu_id)
)

# Grab representative sequences
otu_representative_sequences <- asv_sequences[names(asv_sequences) %in% otu_representative_names]

# Rename sequences to match otu_id without size information
names(otu_representative_sequences) <- str_replace(names(otu_representative_sequences), ";size=\\d+$", "")

# Check number of sequences and names match
message(paste0("Final number of OTU representative sequences: ", length(otu_representative_sequences)))
if (length(otu_representative_sequences) != nrow(otu_classification)) {
  stop("Mismatch in number of OTU representative sequences and classification")
}

# Check if all names match (both directions)
otu_rep_names <- sort(names(otu_representative_sequences))
otu_class_names <- sort(otu_classification$otu_id)
if (!all(otu_rep_names %in% otu_class_names)) {
  stop("Some OTU IDs in representative sequences do not match classification")
}

if (!all(otu_class_names %in% otu_rep_names)) {
  stop("Some OTU IDs in classification do not match representative sequences")
}

message("All OTU IDs match between representative sequences and classification!")

# Save the temporary classification, matrix and sequence files -----------------

# Save OTU sample matrix
otu_sample_matrix %>%
  fwrite("./tmp_clusters/otu_glom_table.txt", sep = "\t")

# Save OTU classification file
otu_classification %>%
  fwrite("./tmp_clusters/otu_glom_classification.txt", sep = "\t")

# Save representative sequences
writeXStringSet(otu_representative_sequences, "./tmp_clusters/otu_glom_sequences.fasta")

# Save the final outputs -------------------------------------------------------

# Read clusters for non-Glomeromycota identified at family level
non_glomeromycota_clusters <- fread("./tmp_clusters/family_non_glomeromycota_clusters.txt") %>%
  # Remove size information from otu_id
  mutate(
    otu_id = str_replace(otu_id, ";size=\\d+$", "")
  )

# Save clusters non-Glomeromycota clusters
non_glomeromycota_clusters %>%
  bind_rows(
    fread("./tmp_clusters/species_clusters.txt") %>%
      filter(family == "unidentified" | family == "") %>%
      # Remove size information from otu_id
      mutate(
        otu_id = str_replace(otu_id, ";size=\\d+$", "")
        )
    ) %>%
  bind_rows(
    fread("./tmp_clusters/family_unidentified_clusters.txt") %>%
      # Remove size information from otu_id
      mutate(
        otu_id = str_replace(otu_id, ";size=\\d+$", "")
        )
    ) %>%
  fwrite("../output/otu_classification_non_glom.txt", sep = "\t")

# Generate a family-level otu table for non-Glomeromycota identified at family level
non_glomeromycota_family_clusters <- asv_sample_matrix %>%
  rename(otu_id = OTU_ID) %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  inner_join(
    non_glomeromycota_clusters %>%
      select(otu_id, family),
    by = "otu_id"
  ) %>%
  group_by(family, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  ungroup() %>%
  arrange(desc(abundance)) 

# Generate family-level otu table for Glomeromycota identified at family level
glomeromycota_family_clusters <- otu_sample_matrix %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  inner_join(
    otu_classification %>%
      select(otu_id, family),
    by = "otu_id"
  ) %>%
  group_by(family, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  ungroup() %>%
  arrange(desc(abundance))

# Join and save final family-level otu table for all taxa
bind_rows(
    non_glomeromycota_family_clusters,
    glomeromycota_family_clusters
  ) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = abundance,
    values_fill = 0
  ) %>%
  fwrite("../output/otu_table_all_family.txt", sep = "\t")

# Save final family-level otu table for Glomeromycota
glomeromycota_family_clusters %>%
  pivot_wider(
    names_from = sample_id,
    values_from = abundance,
    values_fill = 0
  ) %>%
  fwrite("../output/otu_table_glom_family.txt", sep = "\t")
