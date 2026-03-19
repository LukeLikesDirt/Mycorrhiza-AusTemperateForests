#!/bin/bash

# Script:   Taxonomic classification with vsearch --usearch_global and 
#           Glomromycota-specific similarity cutoffs
# Author:   Luke Florence
# Date:     19th March 2026
#
# Notes:
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# - Run this script from the './bioinformatics/am/' directory
#
# - vsearch --usearch_global (global alignment) is used against rank-specific
#   EUKARYOME reference subsets (species-, genus-, family-identified) and the 
#   full reference
#
# - The reference set are obtained from dyna-clust-predict-am, which extracts
#   the 18S V4 region using WANDA–AML2 (i.e. the same primers used in this study)
#   to ensure full-length global alignments and accurate similarity estimates
#
# - The 'ASVs' output is classified based on Glomromycota-specific similarity
#   cutoffs for ranks family to species predicted by dyna-clust-predict-am
#
# - dyna-clust-predict-am: https://github.com/LukeLikesDirt/dyna-clust-predict-am
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Parameters and file paths
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Input reference files
readonly REFERENCE_SEQUENCES_ALL="./data/ref_seqs/eukaryome_V4.fasta"
readonly CLASSIFICATION_FILE="./data/ref_seqs/eukaryome_V4.classification"
readonly CUTOFF_FILE="./data/ref_seqs/cutoffs_glom_V4.txt"

# Temporary rank-specific reference subsets
readonly REFERENCE_SEQUENCES_SPECIES="./tmp/eukaryome_V4_species.fasta"
readonly REFERENCE_SEQUENCES_GENUS="./tmp/eukaryome_V4_genus.fasta"
readonly REFERENCE_SEQUENCES_FAMILY="./tmp/eukaryome_V4_family.fasta"

# Temporary vsearch output files
readonly VSEARCH_SPECIES_FILE="./tmp/vsearch_species.txt"
readonly VSEARCH_GENUS_FILE="./tmp/vsearch_genus.txt"
readonly VSEARCH_FAMILY_FILE="./tmp/vsearch_family.txt"
readonly VSEARCH_ALL_FILE="./tmp/vsearch_all.txt"

# Utility scripts
readonly SUBSET_IDENTIFIED="./utils/subset_identified.R"
readonly PREPARE_CLUSTERS="./utils/prepare_glomeromycota_clusters.R"

# Input sequences
readonly IN_SEQUENCES="./data/asv_sequences.fasta"

# Output files
readonly OUT_SEQUENCES="./data/asv_sequences.fasta"
readonly OUT_CLASSIFICATION="./data/asv_classification.txt"

# Number of threads
readonly NUM_THREADS=10

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Functions
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

subset_references() {

    Rscript "$SUBSET_IDENTIFIED" \
        --input "$REFERENCE_SEQUENCES_ALL" \
        --classification "$CLASSIFICATION_FILE" \
        --rank species \
        --output "$REFERENCE_SEQUENCES_SPECIES"

    Rscript "$SUBSET_IDENTIFIED" \
        --input "$REFERENCE_SEQUENCES_ALL" \
        --classification "$CLASSIFICATION_FILE" \
        --rank genus \
        --output "$REFERENCE_SEQUENCES_GENUS"

    Rscript "$SUBSET_IDENTIFIED" \
        --input "$REFERENCE_SEQUENCES_ALL" \
        --classification "$CLASSIFICATION_FILE" \
        --rank family \
        --output "$REFERENCE_SEQUENCES_FAMILY"

}

vsearch_species() {

    vsearch \
        --usearch_global "$IN_SEQUENCES" \
        --db "$REFERENCE_SEQUENCES_SPECIES" \
        --id 0.5 \
        --strand both \
        --maxaccepts 5 \
        --maxrejects 0 \
        --userout "$VSEARCH_SPECIES_FILE" \
        --userfields "query+target+id+ql+tl+alnlen" \
        --threads "$NUM_THREADS"

}

vsearch_genus() {

    vsearch \
        --usearch_global "$IN_SEQUENCES" \
        --db "$REFERENCE_SEQUENCES_GENUS" \
        --id 0.5 \
        --strand both \
        --maxaccepts 5 \
        --maxrejects 0 \
        --userout "$VSEARCH_GENUS_FILE" \
        --userfields "query+target+id+ql+tl+alnlen" \
        --threads "$NUM_THREADS"

}

vsearch_family() {

    vsearch \
        --usearch_global "$IN_SEQUENCES" \
        --db "$REFERENCE_SEQUENCES_FAMILY" \
        --id 0.5 \
        --strand both \
        --maxaccepts 5 \
        --maxrejects 0 \
        --userout "$VSEARCH_FAMILY_FILE" \
        --userfields "query+target+id+ql+tl+alnlen" \
        --threads "$NUM_THREADS"

}

vsearch_all() {

    vsearch \
        --usearch_global "$IN_SEQUENCES" \
        --db "$REFERENCE_SEQUENCES_ALL" \
        --id 0.5 \
        --strand both \
        --maxaccepts 5 \
        --maxrejects 0 \
        --userout "$VSEARCH_ALL_FILE" \
        --userfields "query+target+id+ql+tl+alnlen" \
        --threads "$NUM_THREADS"

}

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Main script
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

echo 'Starting at:' $(date)

# Activate the conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dyna_clust_env

echo "=== CREATING RANK-SPECIFIC REFERENCE SUBSETS ==="
echo $(date)
echo ""
subset_references

echo "=== ASSIGNING TAXONOMY WITH vsearch AT RANK SPECIES ==="
echo $(date)
echo ""
vsearch_species

echo "=== ASSIGNING TAXONOMY WITH vsearch AT RANK GENUS ==="
echo $(date)
echo ""
vsearch_genus

echo "=== ASSIGNING TAXONOMY WITH vsearch AT RANK FAMILY ==="
echo $(date)
echo ""
vsearch_family

echo "=== ASSIGNING TAXONOMY WITH vsearch AGAINST FULL REFERENCE ==="
echo $(date)
echo ""
vsearch_all

echo "=== PREPARE GLOMEROMYCOTA CLUSTERS ==="
echo $(date)
echo ""
Rscript "$PREPARE_CLUSTERS" \
    --cutoffs "$CUTOFF_FILE" \
    --classification "$CLASSIFICATION_FILE" \
    --vsearch_species "$VSEARCH_SPECIES_FILE" \
    --vsearch_genus "$VSEARCH_GENUS_FILE" \
    --vsearch_family "$VSEARCH_FAMILY_FILE" \
    --vsearch_all "$VSEARCH_ALL_FILE" \
    --input "$IN_SEQUENCES" \
    --output "$OUT_CLASSIFICATION" \
    --threads "$NUM_THREADS"

conda deactivate

echo 'Finished at:' $(date)
