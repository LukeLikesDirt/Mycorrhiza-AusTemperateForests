#!/usr/bin/env Rscript
# =============================================================================
# FASTA Header Formatting for Classification Script
# Description: Format FASTA headers to create a classification file
#              and modify headers to match DNABARCODER requirements
# =============================================================================

# Load required packages
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
})

# Input and output file paths
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments provided
if (length(args) != 3) {
  stop("Usage: Rscript format_fasta_classification.R <input_fasta> <output_fasta> <classification_output>")
}

fasta_in <- args[1]
fasta_out <- args[2]
classification_out <- args[3]

# =============================================================================
# Read the FASTA file
fasta_lines <- readLines(fasta_in)

# Extract headers and sequences
header_indices <- which(grepl("^>", fasta_lines))

# Create dataframe with headers and sequences
fasta_df <- tibble(
  header = character(),
  sequence = character()
)

for (i in seq_along(header_indices)) {
  current_header_idx <- header_indices[i]
  
  # Determine sequence end
  if (i < length(header_indices)) {
    seq_end_idx <- header_indices[i + 1] - 1
  } else {
    seq_end_idx <- length(fasta_lines)
  }
  
  # Get header and sequence
  header <- fasta_lines[current_header_idx]
  seq_lines <- fasta_lines[(current_header_idx + 1):seq_end_idx]
  sequence <- paste(seq_lines, collapse = "")
  
  fasta_df <- bind_rows(fasta_df, tibble(header = header, sequence = sequence))
}

# Process headers to create classification data
classification_df <- fasta_df %>%
  # Extract ID (everything before the first semicolon)
  mutate(id = str_extract(header, "(?<=^>)[^;]+")) %>%
  # Extract taxonomic information
  mutate(
    kingdom = str_extract(header, "k__([^;]+)") %>% str_remove("k__"),
    phylum = str_extract(header, "p__([^;]+)") %>% str_remove("p__"),
    class = str_extract(header, "c__([^;]+)") %>% str_remove("c__"),
    order = str_extract(header, "o__([^;]+)") %>% str_remove("o__"),
    family = str_extract(header, "f__([^;]+)") %>% str_remove("f__"),
    genus = str_extract(header, "g__([^;]+)") %>% str_remove("g__"),
    species_epithet = str_extract(header, "s__([^;]+)$") %>% str_remove("s__")
  ) %>%
  # Replace "unclassified" with "unidentified" across all taxonomic ranks
  mutate(across(c(kingdom, phylum, class, order, family, genus, species_epithet), 
                ~ifelse(.x == "unclassified", "unidentified", .x))) %>%
  # Replace EKARYOME genus predictions to "unidentified"
  mutate(genus = ifelse(grepl("_gen\\d+$", genus), "unidentified", genus)) %>%
  # Create proper species names (genus + species epithet for identified species)
  mutate(
    species = case_when(
      species_epithet == "unidentified" ~ "unidentified",
      is.na(species_epithet) | species_epithet == "" ~ "unidentified",
      genus == "unidentified" ~ "unidentified",
      !is.na(genus) & !is.na(species_epithet) & 
        species_epithet != "unidentified" & genus != "unidentified" ~ 
        paste(genus, species_epithet, sep = " "),
      TRUE ~ "unidentified"
    )
  ) %>%
  # Reformat genus names that contain kingdom name in parentheses
  mutate(
    genus = ifelse(
      grepl("\\(", genus),
      paste0(
        str_extract(genus, "^[^\\(]+"),               # Genus name before "("
        "_kng_",
        str_extract(genus, "(?<=\\()[^\\)]+")         # Kingdom name inside "()"
      ),
      genus
    )
  ) %>%
  # Calculate taxonomic completeness score
  mutate(
      tax_score = 
        (!is.na(kingdom) & kingdom != "" & kingdom != "unidentified") +
        (!is.na(phylum) & phylum != "" & phylum != "unidentified") +
        (!is.na(class) & class != "" & class != "unidentified") +
        (!is.na(order) & order != "" & order != "unidentified") +
        (!is.na(family) & family != "" & family != "unidentified") +
        (!is.na(genus) & genus != "" & genus != "unidentified") +
        (!is.na(species) & species != "" & species != "unidentified")
    ) %>%
  # Select unique IDs (keeping highest taxonomic score)
  group_by(id) %>%
  arrange(desc(tax_score), id) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  # Replace any NA values with empty strings if needed
  mutate(across(everything(), ~ifelse(is.na(.x), "", .x)))

# Check for duplicated IDs and exit with a warning if found
if (any(duplicated(classification_df$id))) {
  dup_ids <- classification_df$id[duplicated(classification_df$id)]
  stop(paste("Error: Duplicate IDs found in classification data:", paste(dup_ids, collapse = ", ")))
}

# Create new FASTA file using header and sequence columns
fasta_content <- classification_df %>%
  mutate(fasta_header = paste0(">", id)) %>%
  select(fasta_header, sequence) %>%
  tidyr::pivot_longer(everything(), names_to = NULL, values_to = "line") %>%
  pull(line)

writeLines(fasta_content, fasta_out)

# Write the classification file (without header, sequence, and tax_score columns)
classification_df %>%
  select(id, kingdom, phylum, class, order, family, genus, species) %>%
  write_tsv(classification_out)

# Print summary
cat("Processing complete!\n")
cat("Input FASTA:", fasta_in, "\n")
cat("Output FASTA:", fasta_out, "\n")
cat("Classification file:", classification_out, "\n")
cat("Number of sequences processed:", nrow(classification_df), "\n")
cat("Species identified:", sum(classification_df$species != "unidentified"), "\n")
cat("Species unidentified:", sum(classification_df$species == "unidentified"), "\n")