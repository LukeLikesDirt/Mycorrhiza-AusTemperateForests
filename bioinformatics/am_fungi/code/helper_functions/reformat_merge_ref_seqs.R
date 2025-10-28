#!/usr/bin/env Rscript
# =============================================================================
# FASTA Header Reformatting and Merging Script
# Description: Reformat MAARJAM headers to match EUKARYOME taxonomy format
#              and merge both databases into a single reference file
# =============================================================================

# Load required library
suppressPackageStartupMessages(
    library(Biostrings)
)

# =============================================================================
# FUNCTIONS
# =============================================================================

# Function to parse eukaryome headers
parse_eukaryome_header <- function(header) {
  header <- sub("^>", "", header)
  parts <- strsplit(header, ";")[[1]]
  id <- parts[1]

  clean_taxon <- function(taxon_string, prefix) {
    taxon <- sub(prefix, "", taxon_string)
    taxon <- trimws(taxon)  # remove spaces
    if (taxon == "unclassified" || 
        grepl("_?(phy|cl|ord|fam|gen|sp)[0-9]+$", taxon)) {
      return("unidentified")
    }
    return(taxon)
  }

  taxonomy <- list()
  for (part in parts[-1]) {
    if (grepl("^k__", part)) taxonomy$kingdom <- clean_taxon(part, "k__")
    if (grepl("^p__", part)) taxonomy$phylum  <- clean_taxon(part, "p__")
    if (grepl("^c__", part)) taxonomy$class   <- clean_taxon(part, "c__")
    if (grepl("^o__", part)) taxonomy$order   <- clean_taxon(part, "o__")
    if (grepl("^f__", part)) taxonomy$family  <- clean_taxon(part, "f__")
    if (grepl("^g__", part)) taxonomy$genus   <- clean_taxon(part, "g__")
    if (grepl("^s__", part)) taxonomy$species <- clean_taxon(part, "s__")
  }

  # Make lower ranks 'unidentified' if they duplicate the higher rank
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  for (i in seq_along(ranks)[-1]) {
    higher <- taxonomy[[ranks[i - 1]]]
    lower  <- taxonomy[[ranks[i]]]
    if (!is.null(higher) && !is.null(lower) &&
        trimws(tolower(lower)) == trimws(tolower(higher))) {
      taxonomy[[ranks[i]]] <- "unidentified"
    }
  }

  return(list(id = id, taxonomy = taxonomy))
}

# Function to parse maarjam headers
parse_maarjam_header <- function(header) {
  # Remove leading ">"
  header <- sub("^>", "", header)
  
  # Extract ID (pattern: gb|ID|)
  id_match <- regmatches(header, regexpr("(?<=gb\\|)[^|]+(?=\\|)", header, perl = TRUE))
  id <- if (length(id_match) > 0) id_match else NA
  
  # Split the remaining part by spaces
  parts <- strsplit(header, "\\s+")[[1]]
  
  # Find family, genus, and species
  # Family is typically after the ID section
  family_idx <- which(grepl("aceae$", parts))
  family <- if (length(family_idx) > 0) parts[family_idx[1]] else NA
  
  # Update taxonomy: Claroideoglomeraceae -> Entrophosporaceae
  if (!is.na(family) && family == "Claroideoglomeraceae") {
    family <- "Entrophosporaceae"
  }
  
  # Genus is typically after family
  genus_idx <- family_idx + 1
  genus <- if (length(genus_idx) > 0 && genus_idx <= length(parts)) parts[genus_idx] else NA
  
  # Species is after genus
  species_idx <- genus_idx + 1
  species <- NA
  
  if (length(species_idx) > 0 && species_idx <= length(parts)) {
    potential_species <- parts[species_idx]
    
    # Check if it's "sp." or starts with uppercase/numbers (unclassified)
    if (potential_species == "sp." || grepl("^[A-Z0-9]", potential_species)) {
      species <- "unclassified"
    } else if (grepl("^[a-z]", potential_species)) {
      # Classified species (starts with lowercase)
      species <- potential_species
    } else {
      # Default to unclassified
      species <- "unclassified"
    }
  } else {
    species <- "unclassified"
  }
  
  return(list(id = id, family = family, genus = genus, species = species))
}

# Function to produce cleaned eukaryome fasta lines
clean_eukaryome_fasta <- function(eukaryome_file) {
  lines <- readLines(eukaryome_file)
  out <- character(length(lines))
  for (i in seq_along(lines)) {
    line <- lines[i]
    if (grepl("^>", line)) {
      parsed <- parse_eukaryome_header(line)
      tax <- parsed$taxonomy
      # build cleaned header (always include id and available ranks)
      newh <- paste0(">", parsed$id)
      if (!is.null(tax$kingdom)) newh <- paste0(newh, ";k__", tax$kingdom)
      if (!is.null(tax$phylum))  newh <- paste0(newh, ";p__", tax$phylum)
      if (!is.null(tax$class))   newh <- paste0(newh, ";c__", tax$class)
      if (!is.null(tax$order))   newh <- paste0(newh, ";o__", tax$order)
      if (!is.null(tax$family))  newh <- paste0(newh, ";f__", tax$family)
      if (!is.null(tax$genus))   newh <- paste0(newh, ";g__", tax$genus)
      if (!is.null(tax$species)) newh <- paste0(newh, ";s__", tax$species)
      out[i] <- newh
    } else {
      out[i] <- line
    }
  }
  return(out)
}

# Cleaned headers and avoid 'unidentified' keys
create_taxonomy_lookup <- function(eukaryome_file) {
  cleaned <- clean_eukaryome_fasta(eukaryome_file)
  headers <- cleaned[grepl("^>", cleaned)]
  lookup <- list()
  for (header in headers) {
    parsed <- parse_eukaryome_header(header)
    family <- parsed$taxonomy$family
    # skip unidentified family keys and avoid overwriting existing family entry
    if (!is.null(family) && !is.na(family) && family != "unidentified") {
      if (is.null(lookup[[family]])) {
        lookup[[family]] <- parsed$taxonomy
      }
    }
  }
  return(lookup)
}


# Function to merge and reformat fasta files
merge_and_reformat_fasta <- function(maarjam_file, eukaryome_file, output_file) {
  # Create taxonomy lookup from eukaryome file
  cat("Creating taxonomy lookup from eukaryome file...\n")
  tax_lookup <- create_taxonomy_lookup(eukaryome_file)
  
  # Read maarjam fasta file
  cat("Reading and reformatting maarjam fasta file...\n")
  maarjam_fasta <- readLines(maarjam_file)
  
  # Process maarjam lines
  reformatted_maarjam <- character(length(maarjam_fasta))
  matched <- 0
  unmatched <- 0
  
  for (i in seq_along(maarjam_fasta)) {
    line <- maarjam_fasta[i]
    
    if (grepl("^>", line)) {
      # Parse the header
      parsed <- parse_maarjam_header(line)
      
      # Try to match family in lookup
      if (!is.na(parsed$family) && parsed$family %in% names(tax_lookup)) {
        tax <- tax_lookup[[parsed$family]]
        
        # Build new header in eukaryome format
        new_header <- paste0(">", parsed$id)
        
        if (!is.null(tax$kingdom)) new_header <- paste0(new_header, ";k__", tax$kingdom)
        if (!is.null(tax$phylum)) new_header <- paste0(new_header, ";p__", tax$phylum)
        if (!is.null(tax$class)) new_header <- paste0(new_header, ";c__", tax$class)
        if (!is.null(tax$order)) new_header <- paste0(new_header, ";o__", tax$order)
        if (!is.null(tax$family)) new_header <- paste0(new_header, ";f__", tax$family)
        
        # Add genus from maarjam (or from lookup if not available)
        genus <- if (!is.na(parsed$genus)) parsed$genus else tax$genus
        if (!is.null(genus) && !is.na(genus)) new_header <- paste0(new_header, ";g__", genus)
        
        # Add species from maarjam
        species <- if (!is.na(parsed$species)) parsed$species else "unclassified"
        new_header <- paste0(new_header, ";s__", species)
        
        reformatted_maarjam[i] <- new_header
        matched <- matched + 1
      } else {
        # No match found, keep original or create minimal header
        cat("Warning: No taxonomy match for family:", parsed$family, "in line", i, "\n")
        reformatted_maarjam[i] <- paste0(">", parsed$id, ";f__", parsed$family, 
                              ";g__", parsed$genus, ";s__", parsed$species)
        unmatched <- unmatched + 1
      }
    } else {
      # Keep sequence lines as-is
      reformatted_maarjam[i] <- line
    }
  }
  
  # Read and clean eukaryome fasta file
  cat("Reading and cleaning eukaryome fasta file...\n")
  eukaryome_fasta_clean <- clean_eukaryome_fasta(eukaryome_file)

  # Merge: reformatted maarjam + cleaned eukaryome
  cat("Merging files...\n")
  merged_fasta <- c(reformatted_maarjam, eukaryome_fasta_clean)
  
  # Write output
  cat("Writing merged fasta file...\n")
  writeLines(merged_fasta, output_file)
  
  cat("\n=== SUMMARY ===\n")
  cat("Maarjam headers processed:", matched + unmatched, "\n")
  cat("  Successfully matched:", matched, "\n")
  cat("  Unmatched (no taxonomy):", unmatched, "\n")
  cat("Eukaryome sequences added:", sum(grepl("^>", eukaryome_fasta)), "\n")
  cat("Total sequences in merged file:", sum(grepl("^>", merged_fasta)), "\n")
  cat("Output written to:", output_file, "\n")
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage: Rscript reformat_merge_fasta.R <maarjam_file> <eukaryome_file> <output_file>\n")
  cat("\nExample:\n")
  cat("  Rscript reformat_merge_fasta.R \\\n")
  cat("    ../data/ref_seqs/maarjam_SSU_NS31_2021.fasta \\\n")
  cat("    ../data/ref_seqs/eukaryome_SSU_NS31_AML2_v1.9.4.fasta \\\n")
  cat("    ../data/ref_seqs/ref_seqs_V4.fasta\n")
  quit(status = 1)
}

maarjam_file <- args[1]
eukaryome_file <- args[2]
output_file <- args[3]

# Check if input files exist
if (!file.exists(maarjam_file)) {
  cat("ERROR: Maarjam file not found:", maarjam_file, "\n")
  quit(status = 1)
}

if (!file.exists(eukaryome_file)) {
  cat("ERROR: Eukaryome file not found:", eukaryome_file, "\n")
  quit(status = 1)
}

# Run the merging and reformatting
cat("=== STARTING FASTA REFORMATTING AND MERGING ===\n")
cat("Maarjam file:", maarjam_file, "\n")
cat("Eukaryome file:", eukaryome_file, "\n")
cat("Output file:", output_file, "\n\n")

merge_and_reformat_fasta(maarjam_file, eukaryome_file, output_file)

cat("\n=== REFORMATTING COMPLETED SUCCESSFULLY ===\n")