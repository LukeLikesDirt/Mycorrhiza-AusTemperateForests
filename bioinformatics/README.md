# Bioinformatics

Amplicon bioinformatics pipelines for arbuscular mycorrhizal (AM) fungi from Illumina paired-end sequencing and
ectomycorrhizal (ECM) fungi from PacBio HiFi sequencing. 

Both pipelines impliment taxon-specific cutoffs for taxonomic classification and OTU clustering.

Taxon-specific similarity cutoffs for classification and clustering were
predicted using pre-release versions of `dyna-clust-predict`:

- **AM:** [dyna-clust-predict-am](https://github.com/LukeLikesDirt/dyna-clust-predict-am) — Glomeromycota-specific cutoffs for the 18S V4 region
- **ECM:** [dyna-clust-predict-ecm](https://github.com/LukeLikesDirt/dyna-clust-predict-ecm) — ECM-specific cutoffs for the ITS region

---

## AM pipeline (`am/`)

Targets the 18S V4 region using the WANDA–AML2 primer pair. Classification and
clustering are performed against 18S V4 extracted subsets of the EUKARYOME
reference database using vsearch global alignment.

`01_prepare_reads.sh` Trim primers (cutadapt), quality filter, denoise (DADA2), and remove chimeras (vsearch)
`02_classify_asvs.sh` Classify ASVs against EUKARYOME V4 subsets with Glomeromycota-specific similarity cutoffs
`03_cluster_otus.sh` Cluster ASVs into OTUs using reference-based and de novo approaches with rank-specific cutoffs
`04_organise_data.R` Organise outputs for downstream statistical analyses

---

## ECM pipeline (`ecm/`)

Targets the ITS region (ITS2) from PacBio HiFi reads. ITS is extracted using
ITSxRust prior to denoising. Classification and clustering are performed against
rank-specific subsets of the EUKARYOME ITS reference database.

`01_trim_primers.sh` Trim primers with cutadapt
`02_extract_itsxrust.sh` Extract ITS regions using ITSxRust (HiFi preset)
`03_denoise.sh` Quality filter and denoise with DADA2
`04_remove_chimeras.sh` Remove chimeras with vsearch
`05_classify_asvs.sh` Classify ASVs against EUKARYOME ITS subsets with ECM-specific similarity cutoffs
`06_cluster_otus.sh` Cluster ASVs into OTUs using reference-based and de novo approaches with rank-specific cutoffs
`07_organise_data.R` Organise outputs for downstream statistical analyses |

---

## Environments (`envs/`)

`dyna_clust_env.yml` Main pipeline environment

`gemelli_env.yml` Environment for ordination via [gemelli][https://github.com/biocore/gemelli]