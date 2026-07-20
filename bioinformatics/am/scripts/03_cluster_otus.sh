#!/bin/bash

# Script:   Dynamic OTU clustering for the AM pool using rank-specific
#           similarity cutoffs, reference-based and de novo approaches
# Author:   Luke Florence
# Date:     20th March 2026
#
# Notes:
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# - Run this script from the './bioinformatics/am/' directory
#
# - AM fungi evaluated here span two lineages: phylum Glomeromycota, and order
#   Densosporales (Mucoromycota, Endogonomycetes) — both are AM fungi.
#
# - Reference-based clustering uses vsearch --usearch_global (global alignment)
#
# - De novo clustering uses vsearch --allpairs_global + union-find (exact
#   single-linkage) for ASVs with no reference match
#
# - Requires outputs from 02_classify_asvs.sh (asv_classification.txt and
#   tmp_clusters/am_clusters.txt)
#
# - Clustering uses Glomromycota-specific similarity cutoffs for ranks family
#   to species, predicted by dyna-clust-predict-am:
#   https://github.com/LukeLikesDirt/dyna-clust-predict-am
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Parameters and file paths
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Input files
readonly ASV_SEQUENCES="./data/asv_sequences.fasta"
readonly ASV_TABLE="./data/asv_table.txt"
readonly ASV_CLASSIFICATION="./data/asv_classification.txt"
readonly REF_CLASSIFICATION="./data/ref_seqs/eukaryome_V4_all.classification"
readonly AM_CLUSTERS="./tmp_clusters/am_clusters.txt"

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
    --am_clusters         "$AM_CLUSTERS" \
    --threads             "$NUM_THREADS" \
    --maxaccepts          "$MAXACCEPTS" \
    --maxrejects          "$MAXREJECTS" \
    --strand              "$STRAND"

conda deactivate

echo 'Finished at:' $(date)
