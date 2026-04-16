#!/bin/bash

# Script:   Dynamic OTU clustering for ECM ASVs using rank-specific
#           similarity cutoffs, reference-based and de novo approaches
# Author:   Luke Florence
# Date:     26th March 2026
#
# Notes:
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# - Run this script from the './bioinformatics/ecm/' directory
#
# - Reference-based clustering uses vsearch --usearch_global (global alignment)
#
# - De novo clustering uses vsearch --allpairs_global + union-find (exact
#   single-linkage) for ASVs with no reference match
#
# - Requires outputs from 05_classify_asvs.sh (asv_classification.txt and
#   tmp_clusters/ecm_clusters.txt)
#
# - Clustering uses ECM-specific ITS similarity cutoffs stored in
#   data/ref_seqs/cutoffs_ecm_ITS.txt
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Parameters and file paths
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Input files
readonly ASV_SEQUENCES="./data/asv_sequences.fasta"
readonly ASV_TABLE="./data/asv_table.txt"
readonly ASV_CLASSIFICATION="./data/asv_classification.txt"
readonly REF_CLASSIFICATION="./data/ref_seqs/eukaryome_ITS.classification"
readonly ECM_CLUSTERS="./tmp_clusters/ecm_clusters.txt"
readonly ECM_GENERA="./data/ref_seqs/ecm_genera.txt"
readonly FUNGAL_TRAITS="./data/ref_seqs/fungal_traits.txt"

# Utility script
readonly CLUSTER_OTUS="./utils/cluster_otus.R"

# Number of threads
readonly NUM_THREADS=10

# vsearch clustering parameters
readonly MAXACCEPTS=0
readonly MAXREJECTS=0
readonly STRAND="plus"   # "plus" or "both"

# Create output directory if it doesn't exist
mkdir -p "./output"

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Main script
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

echo 'Starting at:' $(date)

# Activate the conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dyna_clust_env

echo "=== CLUSTERING OTUs ==="
echo $(date)
echo ""

Rscript "$CLUSTER_OTUS" \
    --asv_sequences       "$ASV_SEQUENCES" \
    --asv_table           "$ASV_TABLE" \
    --asv_classification  "$ASV_CLASSIFICATION" \
    --ref_classification  "$REF_CLASSIFICATION" \
    --ecm_clusters        "$ECM_CLUSTERS" \
    --ecm_genera          "$ECM_GENERA" \
    --fungal_traits       "$FUNGAL_TRAITS" \
    --threads             "$NUM_THREADS" \
    --maxaccepts          "$MAXACCEPTS" \
    --maxrejects          "$MAXREJECTS" \
    --strand              "$STRAND"

conda deactivate

echo 'Finished at:' $(date)
