# Mycorrhizal Fungal Diversity in Australian Temperate Forests

## Overview

This project assesses the impacts of soil nitrogen (N) and phosphorus (P) on mycorrhizal fungal diversity in Australian temperate forests. Mycorrhizal fungi form symbiotic relationships with plant roots and play critical roles in nutrient cycling, ecosystem functioning, and plant community dynamics.

## Research Questions

1. How do soil nitrogen and phosphorus concentrations affect mycorrhizal fungal diversity?
2. What are the relationships between soil nutrient ratios (N:P) and fungal community composition?
3. Are there specific fungal taxa that respond differently to nitrogen vs. phosphorus availability?

## Project Structure

```
.
├── data/
│   ├── raw/          # Raw data files (soil chemistry, fungal sequences, environmental variables)
│   └── processed/    # Cleaned and processed data ready for analysis
├── scripts/
│   ├── 01_data_processing.R    # Data cleaning and preparation
│   ├── 02_diversity_analysis.R # Alpha and beta diversity analyses
│   ├── 03_nutrient_analysis.R  # N and P relationship analyses
│   └── 04_visualization.R      # Figure generation
├── outputs/
│   ├── figures/      # Generated plots and figures
│   └── tables/       # Summary statistics and model outputs
├── docs/
│   └── methods.md    # Detailed methodology and analysis notes
├── README.md         # This file
└── LICENSE           # MIT License

```

## Methods

### Study System
- **Location**: Australian temperate forests
- **Mycorrhizal types**: Arbuscular mycorrhizal (AM) and ectomycorrhizal (ECM) fungi
- **Soil nutrients**: Total and available nitrogen (N) and phosphorus (P)

### Data Collection
- Soil samples collected across forest sites with varying nutrient availability
- DNA extraction and high-throughput sequencing (ITS region for fungi)
- Soil chemical analysis for N, P, C, pH, and other nutrients

### Statistical Analyses
- Alpha diversity metrics (richness, Shannon, Simpson indices)
- Beta diversity and community composition (NMDS, PERMANOVA)
- Generalized linear models (GLMs) to test nutrient effects
- Multivariate analyses (RDA, CCA) to explore nutrient-diversity relationships

## Dependencies

This project uses R (version ≥ 4.0.0) with the following packages:
- `vegan` - Community ecology analyses
- `phyloseq` - Microbiome data handling
- `ggplot2` - Data visualization
- `dplyr` & `tidyr` - Data manipulation
- `lme4` - Mixed effects models
- `indicspecies` - Indicator species analysis

## Usage

1. Place raw data files in `data/raw/`
2. Run scripts in numerical order:
   ```R
   source("scripts/01_data_processing.R")
   source("scripts/02_diversity_analysis.R")
   source("scripts/03_nutrient_analysis.R")
   source("scripts/04_visualization.R")
   ```
3. View outputs in `outputs/figures/` and `outputs/tables/`

## Data Availability

Raw sequencing data and soil chemistry data will be made available upon publication.

## Citation

If you use this code or data, please cite:
[Publication details to be added upon manuscript acceptance]

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or collaborations, please open an issue or contact the repository maintainer.

## Acknowledgments

This research acknowledges the traditional custodians of the lands where this research was conducted.