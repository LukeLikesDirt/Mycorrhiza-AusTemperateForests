#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=day
#SBATCH --time=1-00:00:00
#SBATCH --output="./slurm/%x.%j.out"

# Script: Taxonomic assignment with BLASTn
# Author: Luke Florence
# Date: 5th November 2024
# Software: BLAST v2.14.1 - https://blast.ncbi.nlm.nih.gov/Blast.cgi
#
# Notes:
# --------------------------------------------------------------------------- #
# - BLAST is performed on UNITE with species information and the full UNTE 
#   with prefernce given to the species only dataset when classifying OTUs.
# - The 'ASVs' outputs corresponds to taxonomically informed dynamic 
#   clustering developed for this project. The 'OTUs' output is used to 
#   compare results between conventional 97% OTU clustering and the dynamic
#   clustering.

# Constants and subdirectories
readonly NUM_THREADS=64
readonly REFERENCE_SEQUENCES_SPECIES="./dnabarcoder/data/ref_seqs_V4_unique_species"
readonly REFERENCE_SEQUENCES_GENUS="./dnabarcoder/data/ref_seqs_V4_unique_genus"
readonly REFERENCE_SEQUENCES_FAMILY="./dnabarcoder/data/ref_seqs_V4_unique_family"
readonly REFERENCE_SEQUENCES_ALL="./dnabarcoder/data/ref_seqs_V4"
readonly ASV_SEQUENCE_FILE="../data/ASVs_filtered.fasta"
OUTPUT="../data/"

# Helper function and constants for preparing Glomeromycota clusters
readonly PREPARE_CLUSTERS="./helper_functions/prepare_glomeromycota_clusters.R"
readonly CUTOFF_FILE="./dnabarcoder/V4_cutoffs.txt"
readonly CLASSIFICATION_FILE="./dnabarcoder/data/ref_seqs_V4.classification"
readonly BLAST_SPECIES_FILE="../data/blast_species.txt"
readonly BLAST_GENUS_FILE="../data/blast_genus.txt"
readonly BLAST_FAMILY_FILE="../data/blast_family.txt"
readonly BLAST_ALL_FILE="../data/blast_all.txt"
readonly ASV_SEQUENCE_FILE="../data/ASVs_filtered.fasta"
readonly CLASSIFICATION_OUTPUT="../data/asv_classification.txt"
readonly MINLEN=400

blast_species() {

    makeblastdb \
        -in "${REFERENCE_SEQUENCES_SPECIES}.fasta" \
        -out "$REFERENCE_SEQUENCES_SPECIES" \
        -dbtype 'nucl' \
        -hash_index

    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$ASV_SEQUENCE_FILE" \
        -db "$REFERENCE_SEQUENCES_SPECIES" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$OUTPUT/blast_species.txt" \
        -num_threads "$NUM_THREADS"

}

blast_genus() {

    makeblastdb \
        -in "${REFERENCE_SEQUENCES_GENUS}.fasta" \
        -out "$REFERENCE_SEQUENCES_GENUS" \
        -dbtype 'nucl' \
        -hash_index

    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$ASV_SEQUENCE_FILE" \
        -db "$REFERENCE_SEQUENCES_GENUS" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$OUTPUT/blast_genus.txt" \
        -num_threads "$NUM_THREADS"

}

blast_family() {

    makeblastdb \
        -in "${REFERENCE_SEQUENCES_FAMILY}.fasta" \
        -out "$REFERENCE_SEQUENCES_FAMILY" \
        -dbtype 'nucl' \
        -hash_index

    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$ASV_SEQUENCE_FILE" \
        -db "$REFERENCE_SEQUENCES_FAMILY" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$OUTPUT/blast_family.txt" \
        -num_threads "$NUM_THREADS"

}

blast_all() {

    makeblastdb \
        -in "${REFERENCE_SEQUENCES_ALL}.fasta" \
        -out "$REFERENCE_SEQUENCES_ALL" \
        -dbtype 'nucl' \
        -hash_index
    
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$ASV_SEQUENCE_FILE" \
        -db "$REFERENCE_SEQUENCES_ALL" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -out "$OUTPUT/blast_all.txt" \
        -num_threads "$NUM_THREADS"

}

###############################################################################
## Main script ################################################################
###############################################################################

echo 'Starting at:' $(date)

# Activate the conda environment
source ~/.bashrc
conda activate dynamic_clustering

echo "=== ASSIGNING TAXONOMY WITH BLASTn AT RANK SPECIES ==="
echo $(date)
echo ""
blast_species

echo "=== ASSIGNING TAXONOMY WITH BLASTn AT RANK GENUS ==="
echo $(date)
echo ""
blast_genus

echo "=== ASSIGNING TAXONOMY WITH BLASTn AT RANK FAMILY ==="
echo $(date)
echo ""
blast_family

echo "=== ASSIGNING TAXONOMY WITH BLASTn AT RANK ALL ==="
echo $(date)
echo ""
blast_all

echo "=== PREPARE GLOMEROMYCOTA CLUSTERS ==="
echo $(date)
echo ""
Rscript "$PREPARE_CLUSTERS" \
    "$CUTOFF_FILE" \
    "$CLASSIFICATION_FILE" \
    "$BLAST_SPECIES_FILE" \
    "$BLAST_GENUS_FILE" \
    "$BLAST_FAMILY_FILE" \
    "$BLAST_ALL_FILE" \
    "$ASV_SEQUENCE_FILE" \
    "$CLASSIFICATION_OUTPUT" \
    "$MINLEN" \
    "$NUM_THREADS"

conda deactivate

echo 'Finished at:' $(date)