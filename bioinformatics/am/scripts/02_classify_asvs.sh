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
# - The 'ASVs' output is classified based on Glomromycota-specific similarity
#   cutoffs for ranks family to species predicted by dyna-clust-predict-am,
#   available at: https://github.com/LukeLikesDirt/dyna-clust-predict-am
#
# - The reference set are obtained from dyna-clust-predict-am, which extracts
#   the 18S V4 region using WANDA–AML2 (i.e. the same primers used in this study)
#   to ensure full-length global alignments and accurate similarity estimates.
#   To be fair, full-length global alignments were needed for the prediction
#   of the Glomromycota-specific similarity cutoffs, and future work could use
#   the full-length reference sequences for classification to increase taxonomic
#   resolution for non-Glomromycota ASVs.
#
# - AM fungi evaluated here span two lineages: phylum Glomeromycota, and order
#   Densosporales (class Endogonomycetes, Mucoromycota) — both are AM fungi.
#   All of order Densosporales is included per Lutz et al. 2025 (Fungal
#   Ecology 74, 101407), co-authored by the same EUKARYOME curators
#   (Mikryukov, Tedersoo) who name the Densosporales families in Tedersoo
#   et al. 2024 (MycoKeys 107): they define putative E-AMF at order rank,
#   citing strong root-colonisation evidence for Densosporales as a whole
#   (see AM_ORDERS below).
#
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Parameters and file paths
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Number of threads
readonly NUM_THREADS=12

# vsearch clustering and assignment parameters
readonly MAXACCEPTS=0
readonly MAXREJECTS=1000
readonly MIN_SIM_SPECIES=0.995
readonly MIN_SIM_GENUS=0.980
readonly MIN_SIM_FAMILY=0.950
readonly MIN_SIM=0.750

# vsearch --id thresholds: 5% below the target cutoff to allow filter_seqs to
# apply the actual cutoff to the output (avoids missing borderline hits)
readonly MIN_SIM_SPECIES_SEARCH=0.945
readonly MIN_SIM_GENUS_SEARCH=0.930
readonly MIN_SIM_FAMILY_SEARCH=0.900

# Minimum similarity for ASVs to be classified glomeromycota (analgous to family-level classification)
readonly MIN_SIM_GLOM=0.950

# Orders outside phylum Glomeromycota to include in the AM pool as AM fungi
# (comma-separated). Densosporales (class Endogonomycetes, Mucoromycota) is
# included in full per Lutz et al. 2025 (Fungal Ecology 74, 101407), which
# defines putative E-AMF as all ASVs/OTUs assigned to order Densosporales,
# citing strong root-colonisation evidence for the order.
readonly AM_ORDERS="Densosporales"

# Input reference files
readonly REFERENCE_SEQUENCES_ALL="./data/ref_seqs/eukaryome_V4_all.fasta"
readonly CLASSIFICATION_FILE="./data/ref_seqs/eukaryome_V4_all.classification"
readonly CUTOFF_FILE="./data/ref_seqs/cutoffs_V4_am.txt"

# Temporary rank-specific reference subsets
readonly REFERENCE_SEQUENCES_SPECIES="./tmp/eukaryome_V4_species.fasta"
readonly REFERENCE_SEQUENCES_GENUS="./tmp/eukaryome_V4_genus.fasta"
readonly REFERENCE_SEQUENCES_FAMILY="./tmp/eukaryome_V4_family.fasta"

# Temporary vsearch output files
readonly VSEARCH_SPECIES_FILE="./tmp/vsearch_species.txt"
readonly VSEARCH_GENUS_FILE="./tmp/vsearch_genus.txt"
readonly VSEARCH_FAMILY_FILE="./tmp/vsearch_family.txt"
readonly VSEARCH_ALL_FILE="./tmp/vsearch_all.txt"

# Filtered query FASTAs (sequences not yet classified at prior ranks)
readonly REMAINING_AFTER_SPECIES="./tmp/remaining_after_species.fasta"
readonly REMAINING_AFTER_GENUS="./tmp/remaining_after_genus.fasta"
readonly REMAINING_AFTER_FAMILY="./tmp/remaining_after_family.fasta"

# Utility scripts
readonly SUBSET_IDENTIFIED="./utils/subset_identified.R"
readonly CLASSIFY_ASVS="./utils/classify_am.R"

# Input sequences
readonly IN_SEQUENCES="./data/asv_sequences.fasta"

# Output files
readonly OUT_SEQUENCES="./data/asv_sequences.fasta"
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

vsearch_family() {
    local input_fasta="$1"

    vsearch \
        --usearch_global "$input_fasta" \
        --db "$REFERENCE_SEQUENCES_FAMILY" \
        --id "$MIN_SIM_FAMILY_SEARCH" \
        --strand plus \
        --maxaccepts "$MAXACCEPTS" \
        --maxrejects "$MAXREJECTS" \
        --userout "$VSEARCH_FAMILY_FILE" \
        --userfields "query+target+id+ql+tl+alnlen+qcov+tcov" \
        --threads "$NUM_THREADS"

}

vsearch_all() {
    local input_fasta="$1"

    vsearch \
        --usearch_global "$input_fasta" \
        --db "$REFERENCE_SEQUENCES_ALL" \
        --id 0.8 \
        --query_cov 0.8 \
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

echo "=== ASSIGNING TAXONOMY WITH vsearch AT RANK FAMILY ==="
echo $(date)
echo ""
vsearch_family "$REMAINING_AFTER_GENUS"
filter_seqs "$VSEARCH_FAMILY_FILE" "$MIN_SIM_FAMILY"
filter_classified "$VSEARCH_FAMILY_FILE" "$REMAINING_AFTER_GENUS" "$REMAINING_AFTER_FAMILY"

echo "=== ASSIGNING TAXONOMY WITH vsearch AGAINST FULL REFERENCE ==="
echo $(date)
echo ""
vsearch_all "$REMAINING_AFTER_FAMILY"
filter_seqs "$VSEARCH_ALL_FILE" "" 80

echo "=== CLASSIFY ASVs AND PREPARE AM (GLOMEROMYCOTA + ${AM_ORDERS}) / NON-AM POOLS ==="
echo $(date)
echo ""
Rscript "$CLASSIFY_ASVS" \
    --cutoffs "$CUTOFF_FILE" \
    --classification "$CLASSIFICATION_FILE" \
    --vsearch_species "$VSEARCH_SPECIES_FILE" \
    --vsearch_genus "$VSEARCH_GENUS_FILE" \
    --vsearch_family "$VSEARCH_FAMILY_FILE" \
    --vsearch_all "$VSEARCH_ALL_FILE" \
    --input "$IN_SEQUENCES" \
    --output "$OUT_CLASSIFICATION" \
    --min_sim_glom "$MIN_SIM_GLOM" \
    --am_orders "$AM_ORDERS"

conda deactivate

echo 'Finished at:' $(date)


