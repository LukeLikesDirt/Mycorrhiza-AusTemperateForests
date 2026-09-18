# Mycorrhiza-AusTemperateForests

This repository holds code associated with the manuscript:

[**Soil nitrogen reduces ectomycorrhizal diversity and drives lineage-dependent arbuscular mycorrhizal responses in Australian temperate forests**]()

**Authors:**
Luke Florence<sup>1,2</sup>, John W. Morgan<sup>1</sup>, Jennifer L. Wood<sup>3</sup>, Camille Truong<sup>2,4</sup>

**Affiliations:**
1. Department of Ecological, Plant and Animal Sciences, La Trobe University, Melbourne, Victoria, Australia.
2. School of BioSciences, University of Melbourne, Melbourne, Victoria, Australia.
3. Department of Microbiology, Anatomy, Physiology and Pharmacology, La Trobe University, Melbourne, Victoria, Australia.
4. Royal Botanic Gardens Victoria, Melbourne, Victoria, Australia.

Corresponding author: Luke Florence (L.Florence@unimelb.edu.au)

## Overview

We characterised arbuscular mycorrhizal fungi (**AMF**) and ectomycorrhizal fungi (**EMF**) across a natural mineral-nitrogen gradient spanning 126 sites in the Australian Temperate Broadleaf and Mixed Forests biome, using Illumina short-read SSU metabarcoding for AMF and PacBio long-read ITS metabarcoding for EMF. Within a causal-model framework, we estimated the direct effect of mineral nitrogen on alpha diversity and relative read abundance for each guild, used threshold indicator taxa analysis (TITAN2) to identify nitrogen-sensitive taxa, and tested whether their responses were phylogenetically structured.

## Repository contents

- `generated_data/` — Contains all data required to reproduce the primary results from the manuscript.
- `bioinformatics/` — Amplicon bioinformatics pipelines that generate the AMF (Illumina SSU) and EMF (PacBio ITS) OTU tables and taxonomy used throughout the analysis. Raw sequence dataset is archived on figshare (DOI pending).
- `code/00_georef_covars/` — Builds the georeferenced covariate layer used by the diversity models: computes and predicts AM/EcM tree basal area and richness surfaces from national forest inventory data via INLA/SPDE, then extracts bioclim, aridity, soil, and those predicted surfaces at each sample site.
- `code/` — Alpha-diversity modelling, threshold indicator taxa analysis, phylogenetic signal analysis, and figure generation (see "Reproducing the analysis" below).

## Reproducing the analysis

1. **`bioinformatics/am/`, `bioinformatics/ecm/`** — Trim, denoise, classify and cluster reads into OTU tables and taxonomy (`data/amf/`, `data/emf/`).
2. **`code/00_georef_covars/01a_compute_mycorrhizal_dominance.R` → `01b_predict_mycorrhizal_dominance.R` → `01c_compute_mycorrhizal_richness.R` → `01d_predict_mycorrhizal_richness.R` → `01e_extract_georeferenced_covariates.R`** — Compute and extract AM/EcM tree basal-area and richness covariates.
3. **`code/01a_alpha_diversity_amf.R`, `01b_alpha_diversity_g_amf.R`, `01c_alpha_diversity_m_amf.R`, `01d_alpha_diversity_emf.R`, `01e_prepare_generated_data.R`** — Causal-model estimation of the mineral-nitrogen effect on alpha diversity and relative abundance for AMF and EMF.
4. **`02a_titan_amf.R`, `02b_titan_emf.R`** — Threshold indicator taxa analysis (TITAN2) identifying nitrogen-sensitive taxa for each guild.
7. **`03_phylogenetic_analysis.R`** — Tests whether nitrogen-response indicators are phylogenetically clustered (SES-MPD/SES-MNTD) within each guild.
8. **`figure_1.R`–`figure_6.R`** — Generate the manuscript's primary figures from the outputs above.