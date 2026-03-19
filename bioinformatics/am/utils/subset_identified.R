#!/usr/bin/env Rscript
# subset_identified.R — Subset reference FASTA to sequences with identified taxonomy at a given rank
#
# Usage:
#   Rscript utils/subset_identified.R \
#     --input          eukaryome_V4.fasta \
#     --classification eukaryome_V4.classification \
#     --rank           species \
#     --output         eukaryome_V4_species.fasta

# ── Packages ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
})

# ── Taxonomy filtering ────────────────────────────────────────────────────────
#
# is_identified(x)
#
# Returns a logical vector: TRUE where x represents a valid, resolved
# taxonomic name; FALSE for missing, empty, or ambiguous entries.
#
# Patterns treated as unidentified (case-insensitive):
#   - NA or empty string
#   - "unidentified" or "unclassified" (exact, full word)
#   - starts with "uncultured"
#   - contains "incertae sedis" (with space, underscore, dot, or dash separator)
#   - ends with " sp.", "_sp.", ".sp.", or "-sp." (species placeholder)

is_identified <- function(x) {
  !is.na(x) &
  nzchar(x) &
  !grepl("^unidentified$|^unclassified$|^uncultured", x, ignore.case = TRUE) &
  !grepl("incertae[ _.-]sedis",                       x, ignore.case = TRUE) &
  !grepl("[ _.-]sp\\.",                               x, ignore.case = TRUE) &
  !grepl("Unispike1|Unispike2|Unispike3",             x, ignore.case = TRUE) &
  !grepl("Archaea|Bacteria",                          x, ignore.case = TRUE) &
  !grepl("mitochondrion|nucleomorph|plastid",         x, ignore.case = TRUE)
}

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", metavar = "FILE",
              help = "Input FASTA file [required]"),
  make_option(c("-c", "--classification"),
              type = "character", metavar = "FILE",
              help = "Tab-delimited classification file (must have 'id' column) [required]"),
  make_option(c("-r", "--rank"),
              type = "character", metavar = "STR",
              help = "Taxonomic rank to filter by (e.g. species, genus, family) [required]"),
  make_option(c("-o", "--output"),
              type = "character", metavar = "FILE",
              help = "Output FASTA file [required]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── Validation ────────────────────────────────────────────────────────────────

if (is.null(opt$input))          stop("--input is required.")
if (!file.exists(opt$input))     stop("FASTA file not found: ", opt$input)
if (is.null(opt$classification)) stop("--classification is required.")
if (is.null(opt$rank))           stop("--rank is required.")
if (is.null(opt$output))         stop("--output is required.")

# ── Subset sequences ──────────────────────────────────────────────────────────

seqs <- readDNAStringSet(opt$input)
cls  <- fread(opt$classification)

if (!opt$rank %in% colnames(cls)) {
  stop("Column '", opt$rank, "' not found in classification file.")
}

identified_ids  <- cls[is_identified(cls[[opt$rank]]), id]
seqs_filtered   <- seqs[names(seqs) %in% identified_ids]

writeXStringSet(seqs_filtered, opt$output)

message(
  "Wrote ", length(seqs_filtered), " / ", length(seqs),
  " sequences identified to rank '", opt$rank, "' → ", opt$output
)
