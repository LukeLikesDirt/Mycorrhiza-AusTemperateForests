#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=logs/%x.%j.out

################################################################################
# Script: 02b_extract_itsxrust.sh
#
# Purpose: Extract ITS regions from primer-trimmed PacBio HiFi reads using
#          ITSxRust (https://github.com/ayobi/itsxrust), run in parallel with
#          02a_extract_itsxpress.sh for method comparison.
#
#          Uses --preset hifi (strict E-value 1e-10, min-anchor-score 30,
#          max-anchor-evalue 1e-8) appropriate for high-accuracy HiFi data.
#          Extracts ITS1, ITS2, and full ITS simultaneously (--region all).
#          Per-sample QC JSON and anchor TSV are written for downstream QC.
#
# Software:
#    itsxrust v0.2.2
#    hmmer 3.x (nhmmer must be on PATH via conda env)
#
# Inputs:  data/primer_trimmed/*.fastq.gz
#          data/HMMs/all.hmm  (pre-built combined all-taxa HMM profile)
# Outputs: data/itsxrust_extracted/<sample>.its1.fastq
#          data/itsxrust_extracted/<sample>.its2.fastq
#          data/itsxrust_extracted/<sample>.full.fastq   (confident + partial)
#          data/itsxrust_extracted/<sample>_full.fastq   (confident only)
#          data/itsxrust_extracted/<sample>_partial.fastq (partial only)
#          data/itsxrust_extracted/<sample>.ambiguous.fastq.gz
#          data/itsxrust_extracted/<sample>.skipped.fastq.gz
# QC:      data/itsxrust_extracted/<sample>.qc.json
#          data/itsxrust_extracted/<sample>.anchors.tsv
# Logs:    logs/itsxrust_<JOBID_or_timestamp>.log
################################################################################

# ------------------------------------------------------------------------------
# PARAMETER SETUP
# ------------------------------------------------------------------------------
HMM="./data/HMMs/all.hmm"
REGION="full"            # its1 | its2 | full | all
MIN_LEN=250              # Minimum full ITS length (adjust as needed)
MAX_LEN=1500             # Maximum full ITS length (adjust as needed)
PRESET="hifi"            # hifi preset (nc-e=1e-10, min-anchor-score=30, max-anchor-evalue=1e-8)
NUM_THREADS=40

# ------------------------------------------------------------------------------
# DIRECTORY SETUP
# ------------------------------------------------------------------------------
IN_DIR="./data/primer_trimmed"
OUT_DIR="./data/itsxrust_extracted"
LOG_DIR="./logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# ------------------------------------------------------------------------------
# LOGGING SETUP
# ------------------------------------------------------------------------------
LOG_FILE="${LOG_DIR}/itsxrust_${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "Log file: ${LOG_FILE}"

# ------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

# ------------------------------------------------------------------------------
# PREFLIGHT CHECKS
# ------------------------------------------------------------------------------
for tool in nhmmer itsxrust; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "ERROR: '${tool}' not found on PATH — ensure hmmer and itsxrust are installed in the conda env." >&2
        wait; exit 1
    fi
done

if [[ ! -f "${HMM}" ]]; then
    echo "ERROR: HMM profile not found: ${HMM}" >&2
    wait; exit 1
fi

# ------------------------------------------------------------------------------
# INPUT DISCOVERY
# ------------------------------------------------------------------------------
SAMPLES=()
while IFS= read -r f; do
    [[ -n "${f}" ]] && SAMPLES+=("${f}")
done < <(ls -1 "${IN_DIR}"/*.fastq.gz 2>/dev/null | sort || true)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: No input files found. Expected ${IN_DIR}/*.fastq.gz" >&2
    wait; exit 1
fi

echo "Discovered ${#SAMPLES[@]} input files."

# ------------------------------------------------------------------------------
# EXTRACT ITS USING ITSXRUST
# ------------------------------------------------------------------------------
echo "=== START: Extracting ITS using ITSxRust ==="
echo "$(date)"
echo "HMM profile: ${HMM}"
echo "Region:      ${REGION}"
echo "Preset:      ${PRESET}  (inc-e=1e-10, min-anchor-score=30, max-anchor-evalue=1e-8)"
echo "Threads:     ${NUM_THREADS}"
echo ""

for sample_file in "${SAMPLES[@]}"; do
    if [[ ! -f "${sample_file}" ]]; then
        echo "WARNING: Skipping missing input: ${sample_file}" >&2
        continue
    fi

    sample_basename="$(basename "${sample_file}" .fastq.gz)"
    # With --region all, --output is treated as a prefix; itsxrust appends
    # .its1, .its2, .full (+ file extension) automatically.
    output_prefix="${OUT_DIR}/${sample_basename}"

    echo "=== ${sample_basename} ==="
    echo "Input:         ${sample_file}"
    echo "Output prefix: ${output_prefix}"
    echo "Extracting ITS (ITSxRust)..."

    itsxrust extract \
        --input "${sample_file}" \
        --hmm "${HMM}" \
        --output "${output_prefix}".fastq \
        --region "${REGION}" \
        --hmmer-cpu "${NUM_THREADS}" \
        --qc-json "${output_prefix}.qc.json" \
        --anchors-tsv "${output_prefix}.anchors.tsv" \
        --preset "${PRESET}" \
        --max-full "${MAX_LEN}" \
        --min-full "${MIN_LEN}"

    # Post-filter: split the full-ITS output into confident and partial reads.
    # ITSxRust currently writes both confident and partial reads to the output
    # FASTQ. The anchors TSV contains only confident reads (all four anchors
    # found: SSU, 5.8S-start, 5.8S-end, LSU), so we use it as an allow-list.
    # FASTQ headers have the form:  @<read_id>|full:<start>-<end>
    # Anchors TSV col 1 (read_id) matches everything before "|full:".
    # This can be replaced with --confident-only once that flag is released.
    echo "Splitting ${sample_basename}.fastq into confident and partial..."

    # Confident (full-length) reads
    # Note: split() treats the delimiter as regex where | means OR, so we use
    # index() to find the literal "|full:" string instead.
    awk 'NR==FNR {
        if (FNR > 1) ids[$1] = 1
        next
    }
    FNR % 4 == 1 {
        hdr = substr($0, 2)
        pos = index(hdr, "|full:")
        keep = (pos > 0 && substr(hdr, 1, pos - 1) in ids)
    }
    keep { print }
    ' "${output_prefix}.anchors.tsv" "${output_prefix}.full.fastq" \
        > "${output_prefix}_full.fastq"

    # Partial reads (not in anchors TSV)
    awk 'NR==FNR {
        if (FNR > 1) ids[$1] = 1
        next
    }
    FNR % 4 == 1 {
        hdr = substr($0, 2)
        pos = index(hdr, "|full:")
        keep = (pos > 0 && substr(hdr, 1, pos - 1) in ids)
    }
    !keep { print }
    ' "${output_prefix}.anchors.tsv" "${output_prefix}.full.fastq" \
        > "${output_prefix}_partial.fastq"

    echo "  Confident: $(grep -c '^@' "${output_prefix}_full.fastq") reads"
    echo "  Partial:   $(grep -c '^@' "${output_prefix}_partial.fastq") reads"

    echo ""
done

# Compress all output FASTQ files in parallel
echo "Compressing output FASTQ files..."
pigz -p "$NUM_THREADS" "${OUT_DIR}"/*.fastq

echo "=== PIPELINE COMPLETE ==="
echo "$(date)"
wait
