# Raw Data Directory

This directory contains unprocessed data files for the mycorrhizal diversity study.

## Expected Files

### Soil Chemistry Data
- `soil_chemistry.csv` - Soil nutrient concentrations and properties
  - Required columns: site_id, N_total, P_total, N_available, P_available, pH, C_total, moisture
  
### Fungal Sequence Data
- `fungal_otus.csv` or `fungal_asv.csv` - OTU/ASV abundance table
  - Rows: OTUs/ASVs
  - Columns: Samples
  
- `fungal_taxonomy.csv` - Taxonomic assignments for each OTU/ASV
  - Columns: OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species
  
- `sample_metadata.csv` - Sample information and environmental variables
  - Required columns: sample_id, site_id, date_collected, coordinates, forest_type, elevation

### Optional Files
- `fungal_sequences.fasta` - Representative sequences for each OTU/ASV
- `phylogenetic_tree.tre` - Phylogenetic tree (if using phylogenetic metrics)

## Data Format Guidelines

- Use CSV format with UTF-8 encoding
- First row should contain column headers
- Use consistent missing data codes (NA or blank)
- Date format: YYYY-MM-DD
- Coordinate format: Decimal degrees
- Avoid special characters in column names (use underscores instead of spaces)

## Data Privacy and Ethics

- Ensure data collection complied with ethics approvals
- Obtain necessary permissions for sampling on protected lands
- Respect intellectual property rights of collaborators
- Consider data sharing agreements before publication
