#!/bin/bash

# Script:   Taxonomic classification with vsearch --usearch_global and
#           ECM-specific similarity cutoffs
# Author:   Luke Florence
# Date:     26th March 2026
#
# Notes:
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# - Run this script from the './bioinformatics/ecm/' directory
#
# - vsearch --usearch_global (global alignment) is used against rank-specific
#   EUKARYOME ITS reference subsets (species-, genus-, family-identified) and
#   the full reference
#
# - ASVs are classified using ECM-specific ITS similarity cutoffs stored in
#   data/ref_seqs/cutoffs_ecm_ITS.txt
#
# - ECM vs non-ECM assignment is based on genus membership in the ECM genera
#   list (data/ref_seqs/ecm_genera.txt)
#
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Parameters and file paths
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Number of threads
readonly NUM_THREADS=12

# vsearch clustering and assignment parameters
readonly MAXACCEPTS=10000
readonly MAXREJECTS=1000
readonly MIN_SIM_SPECIES=0.99
readonly MIN_SIM_GENUS=0.95

# vsearch --id thresholds for rejecting hits
readonly MIN_SIM_SPECIES_SEARCH=0.900
readonly MIN_SIM_GENUS_SEARCH=0.700

# Input reference files
readonly REFERENCE_SEQUENCES_ALL="./data/ref_seqs/eukaryome_ITS.fasta"
readonly CLASSIFICATION_FILE="./data/ref_seqs/eukaryome_ITS.classification"
readonly CUTOFF_FILE="./data/ref_seqs/cutoffs_ecm_ITS.txt"
readonly ECM_GENERA_FILE="./data/ref_seqs/ecm_genera.txt"

# Temporary rank-specific reference subsets
readonly REFERENCE_SEQUENCES_SPECIES="./tmp/eukaryome_ITS_species.fasta"
readonly REFERENCE_SEQUENCES_GENUS="./tmp/eukaryome_ITS_genus.fasta"

# Temporary vsearch output files
readonly VSEARCH_SPECIES_FILE="./tmp/vsearch_species.txt"
readonly VSEARCH_GENUS_FILE="./tmp/vsearch_genus.txt"
readonly VSEARCH_ALL_FILE="./tmp/vsearch_all.txt"

# Filtered query FASTAs (sequences not yet classified at prior ranks)
readonly REMAINING_AFTER_SPECIES="./tmp/remaining_after_species.fasta"
readonly REMAINING_AFTER_GENUS="./tmp/remaining_after_genus.fasta"

# Utility scripts
readonly SUBSET_IDENTIFIED="./utils/subset_identified.R"
readonly CLASSIFY_ASVS="./utils/classify_ecm.R"

# Input sequences
readonly IN_SEQUENCES="./data/asv_sequences.fasta"

# Output files
readonly OUT_CLASSIFICATION="./data/asv_classification.txt"

# Create tmp directory if it doesn't exist
mkdir -p ./tmp

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

}

vsearch_species() {

    vsearch \
        --usearch_global "$IN_SEQUENCES" \
        --db "$REFERENCE_SEQUENCES_SPECIES" \
        --id "$MIN_SIM_SPECIES_SEARCH" \
        --strand plus \
        --maxaccepts "$MAXACCEPTS" \
        --maxrejects "$MAXREJECTS" \
        --userout "$VSEARCH_SPECIES_FILE" \
        --userfields "query+target+id+ql+tl+alnlen+qcov+tcov" \
        --threads "$NUM_THREADS"

}

vsearch_genus() {
    local input_fasta="$1"

    vsearch \
        --usearch_global "$input_fasta" \
        --db "$REFERENCE_SEQUENCES_GENUS" \
        --id "$MIN_SIM_GENUS_SEARCH" \
        --strand plus \
        --maxaccepts "$MAXACCEPTS" \
        --maxrejects "$MAXREJECTS" \
        --userout "$VSEARCH_GENUS_FILE" \
        --userfields "query+target+id+ql+tl+alnlen+qcov+tcov" \
        --threads "$NUM_THREADS"

}

vsearch_all() {
    local input_fasta="$1"

    vsearch \
        --usearch_global "$input_fasta" \
        --db "$REFERENCE_SEQUENCES_ALL" \
        --id 0.7 \
        --query_cov 0.7 \
        --top_hits_only \
        --strand plus \
        --maxaccepts 100 \
        --maxrejects "$MAXREJECTS" \
        --userout "$VSEARCH_ALL_FILE" \
        --userfields "query+target+id+ql+tl+alnlen+qcov+tcov" \
        --threads "$NUM_THREADS"

}

# Filter vsearch output by qcov > min_cov% (field 7, default 90) and, optionally,
# id >= min_sim (field 3, as percentage). Only the qcov filter is applied when
# min_sim is omitted.
filter_seqs() {
    local vsearch_file="$1"
    local min_sim="${2:-}"   # optional, as proportion (e.g. 0.990)
    local min_cov="${3:-90}" # minimum query coverage %, default 90

    if [[ -n "$min_sim" ]]; then
        awk -F'\t' -v min_sim="$min_sim" -v min_cov="$min_cov" \
            '$7 > min_cov && $3 >= min_sim * 100' \
            "$vsearch_file" > "${vsearch_file}.tmp"
    else
        awk -F'\t' -v min_cov="$min_cov" '$7 > min_cov' \
            "$vsearch_file" > "${vsearch_file}.tmp"
    fi
    mv "${vsearch_file}.tmp" "$vsearch_file"
}

# Remove sequences that had a hit in the previous vsearch step before the next
# search: Reads query IDs from the vsearch output from the previous step and any
# sequence listed is excluded from the downstream query.
filter_classified() {
    local vsearch_file="$1"
    local in_fasta="$2"
    local out_fasta="$3"

    awk '{print $1}' "$vsearch_file" | sort -u > "./tmp/classified_ids.txt"

    vsearch \
        --fastx_getseqs "$in_fasta" \
        --labels "./tmp/classified_ids.txt" \
        --notmatched "$out_fasta"
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
filter_seqs "$VSEARCH_SPECIES_FILE" "$MIN_SIM_SPECIES"
filter_classified "$VSEARCH_SPECIES_FILE" "$IN_SEQUENCES" "$REMAINING_AFTER_SPECIES"

echo "=== ASSIGNING TAXONOMY WITH vsearch AT RANK GENUS ==="
echo $(date)
echo ""
vsearch_genus "$REMAINING_AFTER_SPECIES"
filter_seqs "$VSEARCH_GENUS_FILE" "$MIN_SIM_GENUS"
filter_classified "$VSEARCH_GENUS_FILE" "$REMAINING_AFTER_SPECIES" "$REMAINING_AFTER_GENUS"

echo "=== ASSIGNING TAXONOMY WITH vsearch AGAINST FULL REFERENCE ==="
echo $(date)
echo ""
vsearch_all "$REMAINING_AFTER_GENUS"
filter_seqs "$VSEARCH_ALL_FILE" "" 70

echo "=== CLASSIFY ASVs AND PREPARE ECM / NON-ECM POOLS ==="
echo $(date)
echo ""
Rscript "$CLASSIFY_ASVS" \
    --cutoffs         "$CUTOFF_FILE" \
    --classification  "$CLASSIFICATION_FILE" \
    --vsearch_species "$VSEARCH_SPECIES_FILE" \
    --vsearch_genus   "$VSEARCH_GENUS_FILE" \
    --vsearch_all     "$VSEARCH_ALL_FILE" \
    --input           "$IN_SEQUENCES" \
    --output          "$OUT_CLASSIFICATION" \
    --ecm_genera      "$ECM_GENERA_FILE"

conda deactivate

echo 'Finished at:' $(date)
