# Filters to remove low abundance OTUs from samples -----------------------------
# Removes otus sample-wise with an abundance <= threshold
filter_samples <- function(otu_tab, threshold = 0.1, min_count = 1) {
  
  # Filter OTUs from samples based on within-sample relative abundance
  # Calculates relative abundance using only the OTUs present in otu_tab
  
  otu_tab %>%
    # Convert to long format for calculations
    pivot_longer(-otu_id, names_to = 'sample') %>%
    
    # Remove zeros before calculating relative abundance
    filter(value > 0) %>%
    
    # Group by sample to calculate sample-level statistics
    group_by(sample) %>%
    # Calculate relative abundance as percentage of sample total
    # (based only on OTUs in this table)
    mutate(rel_abund = value / sum(value) * 100) %>%
    ungroup() %>%
    
    # Filter: set to zero if EITHER condition is true:
    # 1. Relative abundance within sample <= threshold
    # 2. Absolute count <= min_count (default 1 read)
    mutate(value = if_else(rel_abund <= threshold | value <= min_count, 0, value)) %>%
    
    # Remove temporary calculation column
    select(-rel_abund) %>%
    
    # Remove filtered entries (now zeros)
    filter(value > 0) %>%
    
    # Convert back to wide format
    pivot_wider(names_from = "sample", values_from = "value", values_fill = 0)
  
}


filter_samples <- function(otu_tab, threshold = 0.1, 
                           rare_otu_cutoff = 100,
                           rare_min_count = 1,
                           abundant_otu_cutoff = 1000) {
  
  otu_tab %>%
    pivot_longer(-otu_id, names_to = 'sample') %>%
    filter(value > 0) %>%
    
    # Calculate total abundance per OTU across all samples
    group_by(otu_id) %>%
    mutate(otu_total = sum(value)) %>%
    ungroup() %>%
    
    group_by(sample) %>%
    mutate(rel_abund = value / sum(value) * 100) %>%
    ungroup() %>%
    
    # Filtering strategy based on OTU total abundance:
    mutate(value = case_when(
      # Rare OTUs (< 500 reads total): keep all occurrences
      otu_total < rare_otu_cutoff ~ value,
      
      # Intermediate OTUs (500-999 reads): remove singletons only
      otu_total >= rare_otu_cutoff & otu_total < abundant_otu_cutoff ~ 
        if_else(value <= rare_min_count, 0, value),
      
      # Abundant OTUs (≥1000 reads): apply proportional filter
      otu_total >= abundant_otu_cutoff ~ 
        if_else(rel_abund <= threshold, 0, value),
      
      TRUE ~ value
    )) %>%
    
    select(otu_id, sample, value) %>%
    filter(value > 0) %>%
    pivot_wider(names_from = "sample", values_from = "value", values_fill = 0)
}

# Filter rare occurrences of abundant otus -------------------------------------
# Proportional version
filter_library <- function(library_specific_otu_table, threshold = 0.1) {
  
  # Remove rare occurrences of abundant otus using a threshold
  # This filters based on each OTU's distribution across samples
  
  filtered_otu_table <- library_specific_otu_table %>%
    # Convert to long format for grouped operations
    pivot_longer(-otu_id, names_to = "sample",
                 values_to = "otu_sample_abundance") %>%
    
    # Group by OTU to calculate library-wide statistics
    group_by(otu_id) %>%
    mutate(
      # Calculate total abundance of this OTU across all samples
      otu_library_abundance = sum(otu_sample_abundance),
      # Calculate relative abundance of this OTU in each sample
      # (as a percentage of the OTU's total library abundance)
      rel_abundance = (otu_sample_abundance / otu_library_abundance) * 100
    ) %>%
    
    # Filter: set to zero if relative abundance within OTU <= threshold (default 0.1%)
    mutate(otu_sample_abundance = replace(
      otu_sample_abundance, rel_abundance <= threshold, 0)) %>%
    ungroup() %>%
    
    # Remove temporary calculation columns
    select(-c(otu_library_abundance, rel_abundance)) %>%
    
    # Convert back to wide format
    pivot_wider(
      names_from = "sample",
      values_from = "otu_sample_abundance",
      values_fill = 0
    )
  
  return(filtered_otu_table)
  
}

# Filter samples based on a subset otu table ------------------------------
filter_subset_samples <- function(otu_tab, sample_abundance, threshold = 0.01, min_count = 1) {
  
  # Filter OTUs from samples based on sample-level relative abundance
  # Use this when working with a subset of OTUs (e.g., specific taxonomic group)
  # Requires total sample abundances from the full dataset for accurate percentages
  
  otu_tab %>%
    # Convert to long format for calculations
    pivot_longer(-otu_id, names_to = 'sample') %>%
    
    # Remove zeros before joining (efficiency)
    filter(value > 0) %>%
    
    # Join with total sample abundances from the full dataset
    left_join(sample_abundance, by = c("sample" = "sample_id")) %>%
    
    # Calculate relative abundance as percentage of total sample abundance
    # This accounts for the fact that we're working with a taxonomic subset
    mutate(rel_abund = (value / total_abundance) * 100) %>%
    
    # Filter: set to zero if EITHER condition is true:
    # 1. Relative abundance within sample <= threshold (default 0.01%)
    # 2. Absolute count <= min_count (default 1 read)
    mutate(value = if_else(rel_abund <= threshold | value <= min_count, 0, value)) %>%
    
    # Remove temporary calculation columns
    select(-rel_abund, -total_abundance) %>%
    
    # Remove filtered entries (now zeros)
    filter(value > 0) %>%
    
    # Convert back to wide format
    pivot_wider(names_from = "sample", values_from = "value", values_fill = 0)
}