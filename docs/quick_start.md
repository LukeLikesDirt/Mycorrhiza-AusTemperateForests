# Quick Start Guide

This guide will help you get started with the Mycorrhizal Diversity Analysis project.

## Prerequisites

- R version 4.0.0 or higher
- RStudio (recommended but not required)

## Step 1: Install Required Packages

Run the setup script to install all required R packages:

```R
source("setup.R")
```

This will:
- Check your R version
- Install required packages (vegan, dplyr, tidyr, readr, ggplot2, patchwork, MASS)
- Optionally install additional packages (phyloseq, lme4, indicspecies)
- Verify the project directory structure

## Step 2: Prepare Your Data

Place your data files in the `data/raw/` directory:

1. **soil_chemistry.csv** - Soil nutrient data
   - Required columns: site_id, N_total, P_total, pH
   - Optional: N_available, P_available, C_total, moisture
   - See `data/raw/soil_chemistry_EXAMPLE.csv` for format

2. **fungal_otus.csv** - OTU/ASV abundance table
   - First column: OTU_ID
   - Subsequent columns: Sample abundances (one column per sample)

3. **fungal_taxonomy.csv** - Taxonomic assignments
   - Required columns: OTU_ID, Phylum
   - Optional: Kingdom, Class, Order, Family, Genus, Species

4. **sample_metadata.csv** - Sample information
   - Required columns: sample_id, site_id
   - Optional: date_collected, latitude, longitude, forest_type, elevation
   - See `data/raw/sample_metadata_EXAMPLE.csv` for format

Note: Example files are provided for soil chemistry and sample metadata. Use these as templates for formatting your fungal sequence and taxonomy data.

## Step 3: Run the Analysis Pipeline

Execute the analysis scripts in order:

```R
# Step 1: Data processing and cleaning
source("scripts/01_data_processing.R")

# Step 2: Calculate diversity metrics
source("scripts/02_diversity_analysis.R")

# Step 3: Test nutrient-diversity relationships
source("scripts/03_nutrient_analysis.R")

# Step 4: Generate visualizations
source("scripts/04_visualization.R")
```

Alternatively, run all scripts at once:

```R
source("scripts/01_data_processing.R")
source("scripts/02_diversity_analysis.R")
source("scripts/03_nutrient_analysis.R")
source("scripts/04_visualization.R")
```

## Step 4: Review Results

Outputs are saved to:
- `data/processed/` - Cleaned and processed data
- `outputs/tables/` - Statistical results and summaries
- `outputs/figures/` - Publication-ready figures

Key output files:
- `alpha_diversity.csv` - Alpha diversity metrics
- `alpha_diversity_nutrient_models.csv` - Model results
- `permanova_results.csv` - Beta diversity test results
- `nmds_ordination.png` - Community composition plot
- `diversity_vs_*.png` - Nutrient-diversity relationships

## Step 5: Interpret Results

1. Check `data_quality_report.txt` for data processing summary
2. Review `diversity_summary.txt` for overview of diversity metrics
3. Examine `nutrient_analysis_summary.txt` for statistical test results
4. View figures in `outputs/figures/` for visual interpretation

## Troubleshooting

### Common Issues

**Issue**: Script fails with "file not found" error
- **Solution**: Ensure data files are in `data/raw/` with correct names

**Issue**: Package installation fails
- **Solution**: Install packages manually: `install.packages("package_name")`

**Issue**: Insufficient sequencing depth warning
- **Solution**: Check `data/processed/sequencing_depth_summary.csv` and consider removing low-coverage samples

**Issue**: No significant relationships detected
- **Solution**: Check if you have sufficient sample size and environmental variation

### Getting Help

1. Check `docs/methods.md` for detailed methodology
2. Review error messages carefully
3. Check session info files for version compatibility issues
4. Open an issue on GitHub for bugs or questions

## Advanced Usage

### Customizing Analyses

You can modify the scripts to:
- Change diversity metrics
- Add additional environmental variables
- Adjust filtering thresholds
- Modify visualization parameters
- Include phylogenetic analyses (requires phylogenetic tree)

### Working with Large Datasets

For datasets with >10,000 OTUs or >100 samples:
- Consider increasing filtering stringency in script 01
- Use more permutations in PERMANOVA (999 is standard)
- Adjust figure dimensions and point sizes in script 04

### Using Different Distance Metrics

Edit `02_diversity_analysis.R` to use alternative distance metrics:
- Jaccard (presence/absence)
- UniFrac (requires phylogenetic tree)
- Weighted UniFrac (abundance-weighted, requires tree)

## Next Steps

After completing the basic analysis:
1. Examine indicator species (requires indicspecies package)
2. Test for guild-specific patterns (if guild annotations available)
3. Perform variance partitioning to separate N vs P effects
4. Consider temporal or spatial structure (if applicable)
5. Integrate with other datasets (e.g., plant diversity, climate)

## Citation

If you use this workflow, please cite:
- R Core Team (2024) for R
- Oksanen et al. (2024) for vegan package
- Any other packages used (listed in session info)
- This repository: LukeLikesDirt/Mycorrhiza-AusTemperateForests

## License

This project is licensed under the MIT License - see LICENSE file for details.
