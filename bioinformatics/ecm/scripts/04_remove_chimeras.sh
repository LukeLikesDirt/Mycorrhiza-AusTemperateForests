#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-01:00:00
#SBATCH --partition=short
#SBATCH --output=logs/%x.%j.out

export readonly IDENTITY=0.98                                        # VSEARCH: Minimum identity for preclustering during chimera removal
export readonly MAP_SCRIPT="./utils/map.pl"                          # Read mapping script (see the following link for an example pipeline that does not require a mapping file: https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline)
export readonly REF_SEQS="./data/ref_seqs/eukaryome_ITS.fasta"       # Reference sequences for chimera removal
export readonly MIN_LENGTH=250                                       # Minimum sequence length
export readonly MAX_LENGTH=1500                                      # Maximum sequence length
export CHIMERA_FILTERED_DIR="./data/chimera_filtered"
export OUTPUT="./data"
export NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# FUNCTION: Chimera removal with VSEARCH
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
chimera_filter() {
    # Remove downstream files if they exist
    local files_to_remove=("all.denovo.nonchimeras.fasta" "all.derep.fasta" "all.derep.uc" "all.nonchimeras.derep.fasta" "all.nonchimeras.fasta" "all.preclustered.fasta" "all.preclustered.uc" "all.ref.nonchimeras.fasta")

    for file in "${files_to_remove[@]}"; do
        rm -f "$CHIMERA_FILTERED_DIR/$file"
    done

    echo 'Dereplicating across samples...'
    vsearch \
        --derep_fulllength "$CHIMERA_FILTERED_DIR/all.fasta" \
        --sizein --sizeout \
        --minseqlength $MIN_LENGTH --maxseqlength $MAX_LENGTH \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.derep.uc" \
        --output "$CHIMERA_FILTERED_DIR/all.derep.fasta"

    echo 'Preclustering reads...'
    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.derep.fasta" \
        --threads $NUM_THREADS \
        --id $IDENTITY \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.preclustered.uc" \
        --centroids "$CHIMERA_FILTERED_DIR/all.preclustered.fasta"

    echo 'Starting de novo chimera detection...'
    vsearch \
        --uchime3_denovo "$CHIMERA_FILTERED_DIR/all.preclustered.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta"

    echo 'Starting reference-based chimera detection...'
    vsearch \
        --uchime_ref "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta" \
        --threads "$NUM_THREADS" \
        --db "$REF_SEQS" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta"

    echo 'Extracting all non-chimeric sequences...'
    perl "$MAP_SCRIPT" "$CHIMERA_FILTERED_DIR/all.derep.fasta" "$CHIMERA_FILTERED_DIR/all.preclustered.uc" "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta"
    perl "$MAP_SCRIPT" "$CHIMERA_FILTERED_DIR/all.fasta" "$CHIMERA_FILTERED_DIR/all.derep.uc" "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta"

    echo "Generating ASVs..."

    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta" \
        --id 1.0 \
        --threads $NUM_THREADS \
        --sizein --sizeout \
        --strand plus \
        --relabel_sha \
        --centroids "$OUTPUT/asv_sequences.fasta" \
        --otutabout "$OUTPUT/asv_table.txt"

    # Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "$OUTPUT/asv_table.txt"

    printf '\nNumber of unique sequences and ASVs\n'
    printf '    Unique non-chimeric sequence: %s\n' "$(grep -c "^>" "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta")"
    printf '    Clustered ASVs: %s\n' "$(grep -c "^>" "$OUTPUT/asv_sequences.fasta")"

}

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Main script
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

echo "Starting chimera filtering..."
echo "$(date)"
chimera_filter

echo "deactivating conda environment..."
conda deactivate

echo "=== PIPELINE COMPLETE ==="
echo "$(date)"