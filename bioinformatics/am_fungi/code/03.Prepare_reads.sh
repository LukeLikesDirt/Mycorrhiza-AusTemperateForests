#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Define directories
export INPUT="../data/raw_data"
export OUTPUT="../data/"
export PRIMERS_TRIMMED_DIR="$OUTPUT/primers_trimmed"
export CHIMERA_FILTERED_DIR="$OUTPUT/chimera_filtered"

# Make directories
mkdir -p "$OUTPUT"
mkdir -p "$PRIMERS_TRIMMED_DIR"
mkdir -p "$CHIMERA_FILTERED_DIR"

# Define constants
export readonly NUM_THREADS=32                                       # Number of threads to use
export readonly PRIMER_FWD="CAGCCGCGGTAATTCCAGCT"                    # Cutadapt: Forward primer WANDA 5'
export readonly PRIMER_REV="GAACCCAAACACTTTGGTTTCC"                  # Cutadapt: Reverse primer AML2 5' (use reverse compliments for merged reads/linked adapters AML2 3'- GGAAACCAAAGTGTTTGGGTTC -5')
export readonly MIN_OVERLAP_FWD=18                                   # Cutadapt: Overlap for the forward primer
export readonly MIN_OVERLAP_REV=20                                   # Cutadapt: Overlap for the reverse primer
export readonly MAX_ERROR_RATE=1                                     # Cutadapt & DADA2: Maximum error rate
export readonly MAX_N=0                                              # DADA2: Maximum number of ambiguous bases (N) allowed
export readonly IDENTITY=0.99                                        # VSEARCH: Minimum identity for preclustering during chimera removal
export readonly MAP_SCRIPT="./helper_functions/map.pl"               # Read mapping script (see the following link for an example pipeline that does not require a mapping file: https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline)
export readonly REF_SEQS="../data/ref_seqs/ref_seqs_V4.fasta"        # Reference sequences for chimera removal

################################################################################
# FUNCTION: Rename fastq files
################################################################################

# Rename the forward and reverse read files to have sample names as prefixes
# and read directions as suffixes (_R1 and _R2).
rename_fastq_files() {
    for file in "${INPUT}"/*.fastq.gz; do
    filename=$(basename "$file")
    
    # Extract sample name
    if [[ "$filename" == am_neg_con_* ]]; then
        base_name=$(sed 's/am_neg_con_\([0-9]*\)_FLO.*/neg\1/' <<< "$filename")
    else
        base_name=$(sed 's/am_\([^_]*\)_FLO.*/\1/' <<< "$filename")
    fi
    
    # Determine read direction
    if [[ "$filename" == *_R1_* ]]; then
        read_suffix="_R1"
    elif [[ "$filename" == *_R2_* ]]; then
        read_suffix="_R2"
    else
        echo "WARNING: Could not determine read direction for $filename" >&2
        continue
    fi
    
    # Construct new filename
    new_name="${base_name}${read_suffix}.fastq.gz"
    
    # Rename the file
    mv "$file" "${INPUT}/${new_name}"
    echo "Renamed $filename → $new_name"
    done
} 

################################################################################
# FUNCTION: Primer trimming with Cutadapt
################################################################################
trim_primers() {
    mkdir -p "$PRIMERS_TRIMMED_DIR" "$OUTPUT"
    
    echo "Trimming primers from paired reads..."
    for fwd_file in "${INPUT}"/*_R1.fastq.gz; do
        base_name=$(basename "$fwd_file" _R1.fastq.gz)
        rev_file="${INPUT}/${base_name}_R2.fastq.gz"
        
        # Check if reverse file exists
        if [ ! -f "$rev_file" ]; then
            echo "Warning: No R2 file found for ${base_name}" >> "${OUTPUT}/logfile_cutadapt.txt"
            continue
        fi
        
        echo "Processing: ${base_name}"
        
        # Trim primers from BOTH files simultaneously in paired-end mode
        cutadapt \
            -g "${PRIMER_FWD};min_overlap=${MIN_OVERLAP_FWD}" \
            -G "${PRIMER_REV};min_overlap=${MIN_OVERLAP_REV}" \
            -o "${PRIMERS_TRIMMED_DIR}/${base_name}_R1.fastq.gz" \
            -p "${PRIMERS_TRIMMED_DIR}/${base_name}_R2.fastq.gz" \
            "$fwd_file" "$rev_file" \
            -e "$MAX_ERROR_RATE" \
            --times 4 \
            --discard-untrimmed \
            --pair-filter=both \
            >> "${OUTPUT}/logfile_cutadapt.txt" 2>&1
    done
    
}

################################################################################
# FUNCTION: Quality reporting with FastQC and MultiQC
################################################################################

quality_reporting() {
    # Quality check with FastQC
    fastqc "${PRIMERS_TRIMMED_DIR}"/*.fastq.gz -o "$PRIMERS_TRIMMED_DIR" -t "$NUM_THREADS"

    # Summarise FastQC results with MultiQC
    multiqc "${PRIMERS_TRIMMED_DIR}" -o "$OUTPUT"

    # Remove intermediate files
    rm "${PRIMERS_TRIMMED_DIR}"/*fastqc.zip "${PRIMERS_TRIMMED_DIR}"/*fastqc.html
    rm -r "${OUTPUT}"/multiqc_data

    # Rename MultiQC report
    mv "${OUTPUT}/multiqc_report.html" "${OUTPUT}/logfile_multiqc_report.html"
}

################################################################################
# FUNCTION: Quality filtering and denoising with DADA2
################################################################################

dada2_processing() {

Rscript --vanilla - <<EOF

library(dada2); packageVersion("dada2")
library(Biostrings)
library(tidyverse)

path <- Sys.getenv("OUTPUT") # CHANGE ME to the directory containing the fastq files
max_ee <- as.numeric(Sys.getenv("MAX_ERROR_RATE"))
max_n <- as.numeric(Sys.getenv("MAX_N"))
# List all files
all_files <- list.files(file.path(path, "primers_trimmed"))
cat("Total files in directory:", length(all_files), "\n")

# Get forward and reverse files
fnFs <- sort(list.files(file.path(path, "primers_trimmed"), pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(file.path(path, "primers_trimmed"), pattern="_R2.fastq.gz", full.names = TRUE))

cat("Number of forward files:", length(fnFs), "\n")
cat("Number of reverse files:", length(fnRs), "\n")

# Extract sample names
sample.names.F <- sapply(basename(fnFs), function(x) sub("_R1.fastq.gz$", "", x))
sample.names.R <- sapply(basename(fnRs), function(x) sub("_R2.fastq.gz$", "", x))

# Check for mismatches
cat("\n=== Checking for unpaired files ===\n")
cat("Forward files without reverse pair:\n")
missing_R <- setdiff(sample.names.F, sample.names.R)
if(length(missing_R) > 0) {
    print(missing_R)
} else {
    cat("None\n")
}

cat("\nReverse files without forward pair:\n")
missing_F <- setdiff(sample.names.R, sample.names.F)
if(length(missing_F) > 0) {
    print(missing_F)
} else {
    cat("None\n")
}

# Only keep paired samples
sample.names <- intersect(sample.names.F, sample.names.R)
cat("\nNumber of properly paired samples:", length(sample.names), "\n")

# Filter file lists to only include paired samples
fnFs <- fnFs[sample.names.F %in% sample.names]
fnRs <- fnRs[sample.names.R %in% sample.names]

# Verify lengths match
cat("After filtering - Forward:", length(fnFs), "Reverse:", length(fnRs), "\n")

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=max_n, maxEE=c(max_ee,max_ee), rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# Learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(sample.names, out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(track) <- c("sample_id", "input", "filtered", "denoised_fwd", "denoised_rev", "merged")
data.table::fwrite(track, file.path(path, "logfile_dada2.txt"), sep="\t")

## Export the fasta file for chimera filtering in VSEARCH
#asv_seqs <- colnames(seqtab)
#asv_ids <- paste0("ASV_", seq_len(length(asv_seqs)))
#names(asv_seqs) <- asv_ids
#asv_stringset <- DNAStringSet(asv_seqs)
#writeXStringSet(asv_stringset, file.path(path, "all1.fasta"))

# Export the fasta file for chimera filtering in VSEARCH
fasta.tab <- as.data.frame(seqtab) %>%
  rownames_to_column('sample.names') %>%
  as_tibble() %>%
  pivot_longer(-sample.names, names_to = 'seq', values_to = 'size') %>%
  filter(size > 0) %>%
  ungroup() %>%
  mutate(seq.name = paste0(sample.names, '.fasta', '.', row_number(),
                           ';size=', size)) %>%
  select(seq.name, seq)

# Write and save the fasta file for the merged sequence table
seqinr::write.fasta(as.list(fasta.tab[["seq"]]), fasta.tab[["seq.name"]],
            file.path(path, 'denoised.fasta'),
                      open = 'w', nbchar = 60, as.string = FALSE)

EOF
}

################################################################################
# FUNCTION: Chimera removal with VSEARCH
################################################################################
chimera_filter() {
    # Remove downstream files if they exist
    local files_to_remove=("all.denovo.nonchimeras.fasta" "all.derep.fasta" "all.derep.uc" "all.nonchimeras.derep.fasta" "all.nonchimeras.fasta" "all.preclustered.fasta" "all.preclustered.uc" "all.ref.nonchimeras.fasta")

    for file in "${files_to_remove[@]}"; do
        rm -f "$CHIMERA_FILTERED_DIR/$file"
    done

    echo 'Dereplicating across samples...'
    vsearch \
        --derep_fulllength "$OUTPUT/denoised.fasta" \
        --sizein --sizeout \
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
    perl "$MAP_SCRIPT" "$OUTPUT/denoised.fasta" "$CHIMERA_FILTERED_DIR/all.derep.uc" "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta"

    echo "Generating ASVs..."

    vsearch \
        --cluster_unoise "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta" \
        --threads $NUM_THREADS \
        --sizein --sizeout \
        --relabel_sha \
        --uc "$OUTPUT/ASVs.uc" \
        --centroids "$OUTPUT/ASVs.fasta" \
        --biomout "$OUTPUT/ASVs.biom" \
        --otutabout "$OUTPUT/ASVs.txt"

    # Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "$OUTPUT/ASVs.txt"

    printf '\nNumber of unique sequences and ASVs\n'
    printf '    Unique non-chimeric sequence: %s\n' "$(grep -c "^>" "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta")"
    printf '    Clustered ASVs: %s\n' "$(grep -c "^>" "$OUTPUT/ASVs.fasta")"


}

################################################################################
# Main script execution
################################################################################

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate sequence_prep

echo -e "Starting fastq renaming process at $(date)"
rename_fastq_files
echo -e "Fastq renaming process completed at $(date)"
echo -e "Renamed sequences saved in: $INPUT"

echo -e "Starting primer trimming process at $(date)"
trim_primers
echo -e "Primer trimming process completed at $(date)"
echo -e "Trimmed sequences saved in: $OUTPUT"
echo -e "Cutadapt log: ${OUTPUT}/logfile_cutadapt.txt"
echo -e "Seqkit log: ${OUTPUT}/logfile_seqkit.txt"

echo -e "Starting quality reporting at $(date)"
quality_reporting
echo -e "Quality reporting completed at $(date)"
echo -e "Quality report: ${OUTPUT}/logfile_multiqc_report.html"

echo -e "Starting DADA2 processing at $(date)"
dada2_processing
echo -e "DADA2 processing completed at $(date)"
echo -e "DADA2 outputs saved in: $OUTPUT"
echo -e "DADA2 log: ${OUTPUT}/logfile_dada2.txt"
echo -e "Fasta files: ${OUTPUT}/all.fasta"

echo -e "Starting chimera filtering at $(date)"
chimera_filter
echo -e "Chimera filtering completed at $(date)"
echo -e "Chimera filtered sequences saved in: $CHIMERA_FILTERED_DIR"
echo -e "Final non-chimeric fasta: ${CHIMERA_FILTERED_DIR}/all.nonchimeras.fasta"

echo -e "All processes completed at $(date)"

# Deactivate conda environment
conda deactivate