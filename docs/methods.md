# Detailed Methods and Analysis Notes

## Study Design

### Site Selection
- Forest sites selected across a gradient of soil nitrogen and phosphorus availability
- Minimum of 3 replicate sites per nutrient level
- Sites matched for other environmental variables (elevation, aspect, forest type)

### Sampling Protocol
- Soil cores collected from 0-10 cm depth (organic + mineral layer)
- Multiple cores composited per site
- Samples divided for:
  - DNA extraction (stored at -80°C)
  - Soil chemistry analysis (air-dried or fresh)

## Laboratory Methods

### Soil Chemistry
- **Total N and C**: Dry combustion (CN analyzer)
- **Total P**: Acid digestion followed by colorimetric analysis
- **Available P**: Olsen P or Bray P extraction
- **Available N**: KCl extraction (NH4+ and NO3-)
- **pH**: 1:5 soil:water or 1:5 soil:CaCl2
- **Moisture content**: Gravimetric method

### Molecular Methods
- **DNA Extraction**: Commercial kit (e.g., DNeasy PowerSoil)
- **PCR Amplification**: ITS1/ITS2 region
  - Forward primer: ITS1F (CTTGGTCATTTAGAGGAAGTAA)
  - Reverse primer: ITS2 (GCTGCGTTCTTCATCGATGC)
- **Sequencing**: Illumina MiSeq (2 × 250 bp paired-end)
- **Quality control**: Minimum Q30 score

## Bioinformatics Pipeline

### Sequence Processing
1. **Quality filtering**: Remove low-quality reads (Q < 30)
2. **Primer trimming**: Remove primer sequences
3. **Merging**: Merge paired-end reads (minimum overlap: 20 bp)
4. **Denoising**: DADA2 or UNOISE algorithm
5. **Chimera removal**: De novo and reference-based
6. **Taxonomic assignment**: UNITE database (fungi) or SILVA (bacteria)
7. **Filtering**: Remove non-target taxa, singletons, and low-abundance OTUs

### Quality Control Thresholds
- Minimum read count per sample: 10,000 reads
- Minimum OTU abundance: 0.01% of total reads
- Samples with insufficient coverage: Remove from analysis

## Statistical Analyses

### Data Preparation
- Rarefy or use proportional abundance to account for unequal sequencing depth
- Transform environmental variables if needed (log transformation for right-skewed data)
- Check for collinearity among predictors (VIF < 5)

### Alpha Diversity
- **Metrics**: 
  - OTU richness (number of unique taxa)
  - Shannon diversity index (accounts for evenness)
  - Simpson diversity index (accounts for dominance)
- **Models**: 
  - GLMs with Gaussian or negative binomial family
  - Predictors: N, P, N:P ratio, pH, C:N
  - Random effects: Site (if hierarchical design)

### Beta Diversity
- **Distance metrics**: 
  - Bray-Curtis dissimilarity (abundance-based)
  - Jaccard distance (presence/absence)
  - UniFrac distance (phylogenetic, if applicable)
- **Ordination**: 
  - NMDS (non-metric multidimensional scaling)
  - PCoA (principal coordinates analysis)
- **Hypothesis testing**: 
  - PERMANOVA (permutational ANOVA)
  - ANOSIM (analysis of similarity)
  - Permutations: 999

### Nutrient-Diversity Relationships
- **Multivariate analyses**:
  - Redundancy analysis (RDA) for linear relationships
  - Canonical correspondence analysis (CCA) for unimodal relationships
  - Distance-based RDA (db-RDA) for non-Euclidean distances
- **Model selection**: Forward selection with permutation tests

### Indicator Species Analysis
- Identify fungal taxa associated with high/low N or P conditions
- IndVal (indicator value) analysis
- Threshold: IndVal > 0.7, p < 0.05

### Additional Analyses
- **Variance partitioning**: Separate effects of N vs. P vs. other factors
- **Functional guilds**: Analyze guild-specific responses (if data available)
  - ECM fungi
  - AM fungi
  - Saprotrophs
  - Pathogens

## Model Validation

### Assumptions
- Check residual plots for homogeneity of variance
- Test for normality (Shapiro-Wilk test, Q-Q plots)
- Assess influential points (Cook's distance)
- Multicollinearity diagnostics (VIF)

### Model Comparison
- Likelihood ratio tests
- AIC/BIC for model selection
- Cross-validation (if sample size permits)

## Visualization

### Recommended Figures
1. **Site map**: Showing sampling locations and nutrient gradients
2. **Nutrient distributions**: Box plots or violin plots
3. **Alpha diversity**: Scatter plots with regression lines
4. **Ordination plots**: NMDS with environmental vectors
5. **Heatmaps**: OTU abundance across sites
6. **Network plots**: Co-occurrence patterns (optional)

## Data Quality Checks

### Pre-analysis Checks
- [ ] All samples have sufficient sequencing depth
- [ ] No batch effects in sequencing runs
- [ ] Environmental variables have reasonable ranges
- [ ] No excessive zeros in OTU table (adjust rarefaction if needed)

### Post-analysis Checks
- [ ] Model residuals meet assumptions
- [ ] Results are robust to different distance metrics
- [ ] Outliers are identified and justified (remove or keep)
- [ ] Sensitivity analyses performed (e.g., removing rare taxa)

## Reproducibility

### Session Information
- Document R version and package versions
- Use `sessionInfo()` or `renv` for package management
- Set random seed for reproducible permutations (`set.seed()`)

### Code Documentation
- Comment complex analytical steps
- Use consistent naming conventions
- Include citations for methods in comments

## References

Key methodological papers:
- Oksanen et al. (2020) vegan: Community Ecology Package
- Callahan et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data
- Anderson (2001) A new method for non-parametric multivariate analysis of variance
- Legendre & Gallagher (2001) Ecologically meaningful transformations for ordination of species data

## Notes

- Update this document as analyses progress
- Document any deviations from planned methods
- Record decisions about data filtering or transformation
- Note any unexpected patterns or results requiring investigation
