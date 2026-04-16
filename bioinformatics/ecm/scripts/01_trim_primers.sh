#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=slurm/%x.%j.out
################################################################################
# Script: 02a_extract_itsxpress.sh
#
# Purpose: Extract ITS regions from primer-trimmed single-end reads using
#          ITSxpress. Reads from data/primer_trimmed/*.fastq.gz are processed
#          per sample and extracted sequences are written to
#          data/itsxpress_extracted/<sample>.fastq.gz.
#
# Software:
#    itsxpress v2.1.4
#
# Inputs:  data/primer_trimmed/*.fastq.gz
# Outputs: data/itsxpress_extracted/<sample>.fastq.gz
# Logs:    logs/itsxpress_<JOBID_or_timestamp>.log
################################################################################

set -euo pipefail

# ------------------------------------------------------------------------------
# PARAMETER SETUP
# ------------------------------------------------------------------------------
export TAXA="All"
export REGION="ITS2"    # e.g. ITS1, ITS2, ALL — set as appropriate
export NUM_THREADS=40

# ------------------------------------------------------------------------------
# DIRECTORY SETUP
# ------------------------------------------------------------------------------
TMP_DIR="./tmp"
IN_DIR="./data/primer_trimmed"        # FIX: was missing '=' and 'data/' prefix
OUT_DIR="./data/itsxpress_extracted"
LOG_DIR="./logs"

mkdir -p "$OUT_DIR" "$TMP_DIR" "$LOG_DIR"

# ------------------------------------------------------------------------------
# LOGGING SETUP
# ------------------------------------------------------------------------------
LOG_FILE="${LOG_DIR}/itsxpress_${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "Log file: ${LOG_FILE}"

# ------------------------------------------------------------------------------
# INPUT DISCOVERY
# ------------------------------------------------------------------------------
SAMPLES=()
while IFS= read -r f; do
    [[ -n "${f}" ]] && SAMPLES+=("${f}")
done < <(ls -1 "${IN_DIR}"/*.fastq.gz 2>/dev/null | sort || true)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: No input files found. Expected ${IN_DIR}/*.fastq.gz" >&2   # FIX: broken quoting in original
    exit 1
fi

echo "Discovered ${#SAMPLES[@]} input files."

# ------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

# ------------------------------------------------------------------------------
# EXTRACT ITS USING ITSXPRESS
# ------------------------------------------------------------------------------
echo "=== START: Extracting ITS using ITSxpress ==="   # FIX: typo 'unsing'
echo "$(date)"
echo "Taxa:    ${TAXA}"
echo "Region:  ${REGION}"
echo ""

for sample_file in "${SAMPLES[@]}"; do
    if [[ ! -f "${sample_file}" ]]; then
        echo "WARNING: Skipping missing input: ${sample_file}" >&2
        continue
    fi

    sample_basename="$(basename "${sample_file}" .fastq.gz)"
    extracted_output="${OUT_DIR}/${sample_basename}.fastq.gz"

    echo "=== ${sample_basename} ==="
    echo "Input:  ${sample_file}"
    echo "Output: ${extracted_output}"
    echo "Extracting ITS (ITSxpress)..."

    itsxpress \
        --fastq "${sample_file}" \
        --single_end \
        --taxa "${TAXA}" \
        --region "${REGION}" \
        --threads "${NUM_THREADS}" \
        --outfile "${extracted_output}"

    echo "Saved ITS-extracted reads to ${extracted_output}"
    echo ""
done

echo "=== PIPELINE COMPLETE ==="
echo "$(date)"