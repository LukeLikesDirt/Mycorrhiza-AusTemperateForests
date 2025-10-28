#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# =============================================================================
# Cutadapt Primer Trimming Pipeline
# Description: Download and process EUKARYOME and MAARJAM SSU sequences in
# preparation for dnabarcoder analysis. EUKARYOME reads are trimmed with NS31
# and AML2 primers to match the MAARJAM database.
# =============================================================================

# Note:
# Apply stringent filtering to the AM-specific reverse primer (AML2) to
# EUKARYOME to ensure reads cover the target taxa. The forward primer (NG31) is
# trimmed optionally with min and max length filters to retain reads of expected
# size.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

# EUKARYOME PARAMETERS: Download URL (EUKARYOME EUK SSU v1.9.4)
readonly DOWNLOAD_URL_EUK="https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_SSU_v1.9.4.zip"
readonly DOWNLOAD_FILE_EUK="../data/ref_seqs/General_EUK_SSU_v1.9.4.zip"
readonly EXTRACTED_DIR_EUK="../data/ref_seqs/"

# MAARJAM PARAMETERS: Download URL (MAARJAM SSU v2021)
readonly DOWNLOAD_URL_MAARJAM="https://maarjam.ut.ee/resources/maarjam_database_SSU.fasta.2021.zip"
readonly DOWNLOAD_FILE_MAARJAM="../data/ref_seqs/maarjam_database_SSU.fasta.2021.zip"
readonly EXTRACTED_DIR_MAARJAM="../data/ref_seqs/"

# FILE PATHS
readonly INPUT_SEQS_EUK="../data/ref_seqs/General_EUK_SSU_v1.9.4.fasta"
readonly OUTPUT_SEQS_EUK="../data/ref_seqs/eukaryome_SSU_NS31_AML2_v1.9.4.fasta"
readonly INPUT_SEQS_MAARJAM="../data/ref_seqs/maarjam_database_SSU.fasta"
readonly OUTPUT_SEQS_MAARJAM="../data/ref_seqs/maarjam_SSU_NS31_2021.fasta"
readonly MERGED_OUTPUT="../data/ref_seqs/ref_seqs_V4.fasta"

# Create necessary directories if they don't exist
mkdir -p ../data/ref_seqs

# HELPER SCRIPTS
readonly MERGE_REFSEQS="helper_functions/reformat_merge_ref_seqs.R"

# CUTADAPT PARAMETERS
readonly LOG_FILE="../data/ref_seqs/logfile_cutadapt_eukaryome.txt"             # Log file for cutadapt
readonly NUM_THREADS=64                                                         # Number of threads for cutadapt
readonly MAX_ERROR_RATE=1                                                       # error rate
readonly PRIMER_FWD="TTGGAGGGCAAGTCTGGTGCC"                                     # NS31 5' primer
readonly PRIMER_REV="GGAAACCAAAGTGTTTGGGTTC"                                    # AML2 3' primer (reverse complement)
readonly MIN_OVERLAP_FWD=21                                                     # Minimum overlap for the forward primer
readonly MIN_OVERLAP_REV=22                                                     # Minimum overlap for the reverse primer
readonly MIN_LEN=400                                                            # Minimum length of reads to keep
readonly MAX_LEN=550                                                            # Maximum length of reads to keep

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate dynamic_clustering

# =============================================================================
# FILE DOWNLOAD: EUKARYOME
# =============================================================================

echo "=== DOWNLOADING EUKARYOME DATABASE ==="

echo "Downloading file from: $DOWNLOAD_URL_EUK"
if ! curl -o "$DOWNLOAD_FILE_EUK" "$DOWNLOAD_URL_EUK"; then
  echo "ERROR: Failed to download file from $DOWNLOAD_URL_EUK" >&2
  exit 1
fi

echo "Unzipping downloaded file..."
if ! unzip -o "$DOWNLOAD_FILE_EUK" -d "$EXTRACTED_DIR_EUK"; then
  echo "ERROR: Failed to unzip $DOWNLOAD_FILE_EUK" >&2
  exit 1
fi

# Clean up zip file
rm -f "$DOWNLOAD_FILE_EUK"

echo "Download and extraction completed successfully!"
echo ""

# =============================================================================
# PRIMER TRIMMING: EUKARYOME
# =============================================================================

# Check if input file exists
if [[ ! -f "$INPUT_SEQS_EUK" ]]; then
echo "ERROR: Input file '$INPUT_SEQS_EUK' not found!" >&2
exit 1
fi

echo "=== STARTING PRIMER TRIMMING ==="
echo "Input file: $INPUT_SEQS_EUK"
echo "Output file: $OUTPUT_SEQS_EUK"
echo "Forward primer: $PRIMER_FWD (min overlap: $MIN_OVERLAP_FWD)"
echo "Reverse primer: $PRIMER_REV (min overlap: $MIN_OVERLAP_REV)"
echo "Max error rate: $MAX_ERROR_RATE"
echo ""

# Step 1: Trim reverse primer (required - discard reads without it)
echo "Step 1: Trimming reverse primer (required)..."
cutadapt \
  -a "$PRIMER_REV;min_overlap=$MIN_OVERLAP_REV" \
  -e "$MAX_ERROR_RATE" \
  --discard-untrimmed \
  --cores "$NUM_THREADS" \
  -o "temp_fwd_trimmed.fasta" "$INPUT_SEQS_EUK" \
  >> "$LOG_FILE" 2>&1

# Step 2: Trim forward primer (optional - keeps all reads)
echo "Step 2: Trimming forward primer (optional)..."
cutadapt \
  -g "$PRIMER_FWD;min_overlap=$MIN_OVERLAP_FWD" \
  -e "$MAX_ERROR_RATE" \
  -M "$MAX_LEN" \
  -m "$MIN_LEN" \
  --cores "$NUM_THREADS" \
  -o "$OUTPUT_SEQS_EUK" "temp_fwd_trimmed.fasta" \
  >> "$LOG_FILE" 2>&1

# Clean up temporary file
rm -f temp_fwd_trimmed.fasta

# Clean up downloaded untrimmed file
rm -f "$INPUT_SEQS_EUK"

echo "Primer trimming completed successfully!"
echo "Output saved to: $OUTPUT_SEQS_EUK"

# Display summary statistics
if [[ -f "$LOG_FILE" ]]; then
echo ""
echo "=== TRIMMING SUMMARY ==="
grep -E "(Total reads processed|Reads with adapters|Reads that were too short|Reads written)" "$LOG_FILE" | tail -8
fi

echo ""

# =============================================================================
# FILE DOWNLOAD: MAARJAM
# =============================================================================

echo "=== DOWNLOADING MAARJAM DATABASE ==="

# Download the Maarjam SSU v2021 database
echo "Downloading MAARJAM file from: $DOWNLOAD_URL_MAARJAM"
if ! curl -o "$DOWNLOAD_FILE_MAARJAM" "$DOWNLOAD_URL_MAARJAM"; then
  echo "ERROR: Failed to download MAARJAM file from $DOWNLOAD_URL_MAARJAM" >&2
  exit 1
fi

# Unzip the downloaded file
echo "Unzipping MAARJAM downloaded file..."
if ! unzip -o "$DOWNLOAD_FILE_MAARJAM" -d "$EXTRACTED_DIR_MAARJAM"; then
  echo "ERROR: Failed to unzip MAARJAM file from $DOWNLOAD_FILE_MAARJAM" >&2
  exit 1
fi

# Clean up zip file
rm -f "$DOWNLOAD_FILE_MAARJAM"

echo "MAARJAM download and extraction completed successfully!"
echo ""

# =============================================================================
# LENGTH FILTERING: MAARJAM
# =============================================================================

echo "=== STARTING LENGTH FILTERING ==="

# Check if input file exists
if [[ ! -f "$INPUT_SEQS_MAARJAM" ]]; then
  echo "ERROR: Input file '$INPUT_SEQS_MAARJAM' not found!" >&2
  exit 1
fi

echo "Input file: $INPUT_SEQS_MAARJAM"
echo "Output file: $OUTPUT_SEQS_MAARJAM"
echo "Minimum length: $MIN_LEN"
echo "Maximum length: $MAX_LEN"
echo ""

# Perform length filtering using seqkit
echo "Running length filtering..."
seqkit seq -m "$MIN_LEN" -M "$MAX_LEN" -g "$INPUT_SEQS_MAARJAM" > "$OUTPUT_SEQS_MAARJAM"

# Remove original unfiltered file
rm -f "$INPUT_SEQS_MAARJAM"

echo "Length filtering completed successfully!"
echo "Output saved to: $OUTPUT_SEQS_MAARJAM"
echo ""

# =============================================================================
# REFORMAT AND MERGE FASTA FILES
# =============================================================================

echo "=== REFORMATTING AND MERGING FASTA FILES ==="

# Check if R script exists
if [[ ! -f "$MERGE_REFSEQS" ]]; then
  echo "ERROR: R script not found: $MERGE_REFSEQS" >&2
  echo "Please ensure the script is located at: $MERGE_REFSEQS" >&2
  exit 1
fi

# Check if required input files exist
if [[ ! -f "$OUTPUT_SEQS_MAARJAM" ]]; then
  echo "ERROR: MAARJAM file not found: $OUTPUT_SEQS_MAARJAM" >&2
  exit 1
fi

if [[ ! -f "$OUTPUT_SEQS_EUK" ]]; then
  echo "ERROR: EUKARYOME file not found: $OUTPUT_SEQS_EUK" >&2
  exit 1
fi

# Execute R script for reformatting and merging
echo "Executing R script for header reformatting and merging..."
Rscript "$MERGE_REFSEQS" "$OUTPUT_SEQS_MAARJAM" "$OUTPUT_SEQS_EUK" "$MERGED_OUTPUT"

# Check if R script executed successfully
if [[ $? -ne 0 ]]; then
  echo "ERROR: R script execution failed!" >&2
  exit 1
fi

echo "Reformatting and merging completed successfully!"
echo "Merged output saved to: $MERGED_OUTPUT"
echo ""

echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="