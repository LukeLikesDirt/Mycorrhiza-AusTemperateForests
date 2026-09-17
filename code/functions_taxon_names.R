# Taxon name formatting ########################################################
#
# Classification names in data/*/classification.txt are underscore-joined paths
# built from three kinds of segments:
#   - formal Latin names       e.g. Glomeraceae, Densosporales, Acaulospora
#   - Tedersoo et al. (2024)   alphanumeric, phylogenetically-placed but informal
#     alphanumeric codes       ranks, e.g. Glomeromycota_cl01, Densosporales_fam02
#   - pseudo-taxa              taxa defined by clustering in this study, always
#                              containing "pseudo_<rank>_<number>"
#
# format_taxon_name() renders a raw classification string for display:
#   - underscores between segments become spaces
#   - pseudo-taxon suffixes become "pseudo-<rank><NN>" (e.g. pseudo-sp01)
#   - genus/species formal Latin names are italicised (markdown *asterisks*,
#     for use with ggtext::element_markdown()); family/order/class/phylum
#     names and Tedersoo alphanumeric codes are not
#
# Examples:
#   Diversispora epigaea                       -> *Diversispora epigaea*
#   Acaulospora_pseudo_sp_0002                 -> *Acaulospora* pseudo-sp02
#   Glomeraceae_pseudo_sp_0011                 -> Glomeraceae pseudo-sp11
#   Glomeromycota_cl01_pseudo_fam_0001         -> Glomeromycota cl01 pseudo-fam01
#   Densosporales_fam02_gen01_pseudo_sp_0004   -> Densosporales fam02 gen01 pseudo-sp04
# ─────────────────────────────────────────────────────────────────────────────
format_taxon_name <- function(x, italics = TRUE) {
  vapply(x, function(name) {
    has_pseudo <- grepl("_pseudo_[a-z]+_[0-9]+$", name)

    if (has_pseudo) {
      pseudo_rank   <- sub(".*_pseudo_([a-z]+)_[0-9]+$", "\\1", name)
      pseudo_number <- as.integer(sub(".*_pseudo_[a-z]+_([0-9]+)$", "\\1", name))
      prefix        <- sub("_pseudo_[a-z]+_[0-9]+$", "", name)
    } else {
      prefix <- name
    }

    # Formal binomial species (e.g. "Diversispora epigaea"): literal space, no underscore
    if (!has_pseudo && grepl(" ", prefix, fixed = TRUE) && !grepl("_", prefix, fixed = TRUE)) {
      return(if (italics) paste0("*", prefix, "*") else prefix)
    }

    # Otherwise the prefix is an underscore-joined path of already-atomic segments
    tokens <- strsplit(prefix, "_", fixed = TRUE)[[1]]
    leaf   <- tokens[length(tokens)]
    lead   <- if (length(tokens) > 1) paste(tokens[-length(tokens)], collapse = " ") else ""

    # Italicise only a genuine formal Latin genus name: alphabetic-only, and not a
    # family/order/class/phylum-suffixed name or a Tedersoo alphanumeric code
    is_genus <- grepl("^[A-Za-z]+$", leaf) &&
      !grepl("(aceae|ales|mycetes|mycota)$", leaf) &&
      leaf != "Glomeromycota"

    leaf_display   <- if (italics && is_genus) paste0("*", leaf, "*") else leaf
    prefix_display <- if (nzchar(lead)) paste(lead, leaf_display) else leaf_display

    if (has_pseudo) {
      paste0(prefix_display, " pseudo-", pseudo_rank, sprintf("%02d", pseudo_number))
    } else {
      prefix_display
    }
  }, character(1), USE.NAMES = FALSE)
}
