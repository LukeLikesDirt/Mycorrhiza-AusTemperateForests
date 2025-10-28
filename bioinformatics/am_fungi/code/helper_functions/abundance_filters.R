# Filters to remove low abundance OTUs from samples -----------------------------

# Removes otus sample-wise with an abundance <= threshold

filter_samples <- function(otu_tab, threshold) {
  # Filter out zero values and calculate relative abundance
  otu_tab %>%
    pivot_longer(-otu_id, names_to = 'sample') %>%
    filter(value > 0) %>%
    group_by(sample) %>%
    mutate(rel_abund = value / sum(value) * 100) %>%
    ungroup() %>%
    
    # Set values below the threshold to zero
    mutate(value = replace(value, rel_abund <= threshold, 0)) %>%
    
    # Remove rows where all values are zero 
    select(-rel_abund) %>%
    filter(value > 0) %>%
    
    # Reshape the table
    pivot_wider(names_from = "sample", values_from = "value", values_fill = 0)
  
}

# Filter rare occurrences of abundant otus -------------------------------------

# Proportional version
filter_library <- function(library_specific_otu_table, threshold) {
  
  # Remove rare occurrences of abundant otus using a 0.5% threshold.
  filtered_otu_table <- library_specific_otu_table %>%
    pivot_longer(-otu_id, names_to = "sample",
                 values_to = "otu_sample_abundance") %>%
    group_by(otu_id) %>%
    mutate(
      otu_library_abundance = sum(otu_sample_abundance),
      rel_abundance = (otu_sample_abundance / otu_library_abundance) * 100
    ) %>%
    mutate(otu_sample_abundance = replace(
      otu_sample_abundance, rel_abundance <= threshold, 0)) %>%
    ungroup() %>%
    select(-c(otu_library_abundance, rel_abundance)) %>%
    pivot_wider(
      names_from = "sample",
      values_from = "otu_sample_abundance",
      values_fill = 0
      )
  
  return(filtered_otu_table)
  
}
