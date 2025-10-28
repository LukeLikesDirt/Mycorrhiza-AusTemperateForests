#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=slurm/%x.%j.out

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate dynamic_clustering

# INPUT AND OUTPUT SEQUENCE AND CLASSIFICATION FILES
readonly REF_SEQS_IN="../data/ref_seqs/ref_seqs_V4.fasta"
readonly REF_SEQS_OUT="./dnabarcoder/data/ref_seqs_V4.fasta"
readonly CLASSIFICATION_OUT="./dnabarcoder/data/ref_seqs_V4.classification"

# OUTPUT FILES FOR UNIQUE SEQUENCES
readonly GLOM_OUT="./dnabarcoder/data/ref_seqs_V4_unique_glom.fasta"
readonly SPECIES_OUT="./dnabarcoder/data/ref_seqs_V4_unique_species.fasta"
readonly GENUS_OUT="./dnabarcoder/data/ref_seqs_V4_unique_genus.fasta"
readonly FAMILY_OUT="./dnabarcoder/data/ref_seqs_V4_unique_family.fasta"
readonly ALL_OUT="./dnabarcoder/data/ref_seqs_V4_unique_all.fasta"

# OUTPUT FILES FOR PREVALENCE FILTERED SEQUENCES
readonly SPECIES_OUT_PREV="./dnabarcoder/data/ref_seqs_V4_unique_species_prev.fasta"
readonly GENUS_OUT_PREV="./dnabarcoder/data/ref_seqs_V4_unique_genus_prev.fasta"
readonly FAMILY_OUT_PREV="./dnabarcoder/data/ref_seqs_V4_unique_family_prev.fasta"
readonly ORDER_OUT_PREV="./dnabarcoder/data/ref_seqs_V4_unique_order_prev.fasta"
readonly CLASS_OUT_PREV="./dnabarcoder/data/ref_seqs_V4_unique_class_prev.fasta"

# PREVALENCE FILTERING CUTOFF (0.5 = 50%, 0.66 = 66%, etc.)
readonly PREVALENCE_CUTOFF="0.66"

# INPUT FOR SIMILARITY FILE
readonly SIM_GLOM="./dnabarcoder/data/ref_seqs_V4_unique_glom.sim"
readonly SIM_ALL="./dnabarcoder/data/ref_seqs_V4_unique_all.sim"

# HELPER SCRIPTS
readonly FORMAT_REF_SEQS="./helper_functions/format_fasta_classification.R"
readonly SELECT_UNIQUE_SEQS="./helper_functions/extract_unique_sequences.R"

# =============================================================================
# DNABARCODER SETUP
# =============================================================================
echo "=== SETTING UP DNABARCODER ==="
echo $(date)
echo ""

# Clone dnabarcoder from GitHub into dnabarcoder_ssu_v4
echo "Cloning dnabarcoder from GitHub..."
if [[ -d "dnabarcoder" ]]; then
    echo "Removing existing dnabarcoder directory..."
    rm -rf dnabarcoder
fi

if ! git clone https://github.com/vuthuyduong/dnabarcoder dnabarcoder; then
    echo "ERROR: Failed to clone dnabarcoder repository" >&2
    exit 1
fi

# Empty the data folder and create it if it doesn't exist
echo "Preparing dnabarcoder data folder..."
mkdir -p dnabarcoder/data
rm -rf dnabarcoder/data/*
echo "dnabarcoder setup completed!"
echo ""

# =============================================================================
# FORMAT FASTA AND CREATE CLASSIFICATION FILE
# =============================================================================
echo "=== FORMATTING FASTA AND CREATING CLASSIFICATION FILE ==="
echo $(date)
echo ""

# Check if the R script exists
if [[ ! -f "$FORMAT_REF_SEQS" ]]; then
    echo "ERROR: R script not found: $FORMAT_REF_SEQS" >&2
    echo "Please ensure the script is located at: $FORMAT_REF_SEQS" >&2
    exit 1
fi

# Check if input file exists
if [[ ! -f "$REF_SEQS_IN" ]]; then
    echo "ERROR: Input file not found: $REF_SEQS_IN" >&2
    exit 1
fi

# Execute the R script to process the FASTA and create classification file
echo "Processing FASTA file and creating classification data..."
if ! Rscript "$FORMAT_REF_SEQS" "$REF_SEQS_IN" "$REF_SEQS_OUT" "$CLASSIFICATION_OUT"; then
    echo "ERROR: R script execution failed" >&2
    exit 1
fi

echo "FASTA formatting and classification file creation completed!"
echo "Formatted FASTA saved to: $REF_SEQS_OUT"
echo "Classification file saved to: $CLASSIFICATION_OUT"
echo ""

# =============================================================================
# SELECT UNIQUE SEQUENCES FOR PREDICTION
# =============================================================================
echo "=== SELECTING UNIQUE SEQUENCES ==="
echo $(date)
echo ""

# Check if the sequence selection R script exists
if [[ ! -f "$SELECT_UNIQUE_SEQS" ]]; then
    echo "ERROR: R script not found: $SELECT_UNIQUE_SEQS" >&2
    echo "Please ensure the script is located at: $SELECT_UNIQUE_SEQS" >&2
    exit 1
fi

# Check if required input files exist
if [[ ! -f "$REF_SEQS_OUT" ]]; then
    echo "ERROR: FASTA file not found: $REF_SEQS_OUT" >&2
    exit 1
fi

if [[ ! -f "$CLASSIFICATION_OUT" ]]; then
    echo "ERROR: Classification file not found: $CLASSIFICATION_OUT" >&2
    exit 1
fi

# Execute the R script to select unique sequences
echo "Extracting unique sequences by taxonomic rank..."
if ! Rscript "$SELECT_UNIQUE_SEQS" "$REF_SEQS_OUT" "$CLASSIFICATION_OUT" "$GLOM_OUT" "$SPECIES_OUT_PREV" "$GENUS_OUT_PREV" "$FAMILY_OUT_PREV" "$ORDER_OUT_PREV" "$CLASS_OUT_PREV" "$SPECIES_OUT" "$GENUS_OUT" "$FAMILY_OUT" "$ALL_OUT" "$PREVALENCE_CUTOFF"; then
    echo "ERROR: Sequence selection script failed" >&2
    exit 1
fi

echo "Unique sequence extraction completed!"
echo "Unique sequences saved to:"
echo " - $GLOM_OUT"
echo " - $CLASS_OUT_PREV"
echo " - $ORDER_OUT_PREV"
echo " - $FAMILY_OUT_PREV"
echo " - $GENUS_OUT_PREV"
echo " - $SPECIES_OUT_PREV"
echo " - $FAMILY_OUT"
echo " - $GENUS_OUT"
echo " - $SPECIES_OUT"
echo " - $ALL_OUT"
echo ""

# =============================================================================
# PREDICT TAXONOMICALLY GUIDED CUTOFFS
# =============================================================================

echo "=== PREDICTING TAXONOMICALLY GUIDED CUTOFFS ==="
echo $(date)
echo ""

# Predict global similarity cutoffs for eukaryotes using all unique sequences
echo 'Predicting global similarity cutoffs for all taxonomic ranks in eukaryotes' $(date)
python dnabarcoder/dnabarcoder.py predict -i $ALL_OUT -c $CLASSIFICATION_OUT -st 0.85 -et 1 -s 0.001 -rank species,genus,family,order,class,phylum,kingdom -removecomplexes yes -ml 400 -sim $SIM_ALL
echo ""

# Predict global similarity cutoffs for Glomeromycota using unique and prevalence filtered Glomeromycota sequences
echo 'Predicting global similarity cutoffs for all taxonomic ranks in Glomeromycota' $(date)
echo 'Predicting class cutoff...'
python dnabarcoder/dnabarcoder.py predict -i $CLASS_OUT_PREV -c $CLASSIFICATION_OUT -st 0.85 -et 1 -s 0.001 -rank class -removecomplexes yes -ml 400 -mingroupno 3
echo ""
echo 'Predicting order cutoff...'
python dnabarcoder/dnabarcoder.py predict -i $ORDER_OUT_PREV -c $CLASSIFICATION_OUT -st 0.85 -et 1 -s 0.001 -rank order -removecomplexes yes -ml 400 -mingroupno 5
echo ""
echo 'Predicting family cutoff...'
python dnabarcoder/dnabarcoder.py predict -i $FAMILY_OUT_PREV -c $CLASSIFICATION_OUT -st 0.85 -et 1 -s 0.001 -rank family -removecomplexes yes -ml 400
echo ""
echo 'Predicting genus cutoff...'
python dnabarcoder/dnabarcoder.py predict -i $GENUS_OUT_PREV -c $CLASSIFICATION_OUT -st 0.85 -et 1 -s 0.001 -rank genus -removecomplexes yes -ml 400
echo ""
echo 'Predicting species cutoff...'
python dnabarcoder/dnabarcoder.py predict -i $SPECIES_OUT_PREV -c $CLASSIFICATION_OUT -st 0.85 -et 1 -s 0.001 -rank species -removecomplexes yes -ml 400
echo ""

# Deactivate conda environment
conda deactivate