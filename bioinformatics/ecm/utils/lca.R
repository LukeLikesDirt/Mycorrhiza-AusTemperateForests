# ── LCA helper ────────────────────────────────────────────────────────────────

resolve_lca_rows <- function(rows) {
  result <- character(length(taxonomy_ranks))
  names(result) <- taxonomy_ranks
  stop_propagating <- FALSE
  for (rank in taxonomy_ranks) {
    if (stop_propagating) {
      result[rank] <- "unclassified"
      next
    }
    vals <- unique(rows[[rank]])
    vals <- vals[!vals %in% c("unclassified", "unidentified", NA, "")]
    if (length(vals) == 1) {
      result[rank] <- vals
    } else {
      result[rank] <- "unclassified"
      stop_propagating <- TRUE
    }
  }
  result
}