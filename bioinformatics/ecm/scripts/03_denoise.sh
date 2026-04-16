#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=4G
#SBATCH --time=07-00:00:00
#SBATCH --partition=week
#SBATCH --output=logs/%x.%j.out

################################################################################
# Script: prepare_reads.sh
#
# Purpose: Quality filter and denoise PacBio CCS (single-end) amplicon reads
#
# Software:
#   R v4.3.3              https://www.r-project.org/
#   tidyverse v2.0.0      https://tidyverse.tidyverse.org
#   DADA2 v1.30.0         https://benjjneb.github.io/dada2/
#   VSEARCH v2.30.3       https://github.com/torognes/vsearch
#   SeqKit v2.12.0        https://bioinf.shenwei.me/seqkit/
################################################################################

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# PARAMETER SETUP
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Quality filtering:
export MAXEE_RATE=0.001 # Maximum expected error rate (0.001 = Q30)
export MAXN=0           # Discard reads containing ambiguous (N) bases
export QMAX=93          # Maximum quality score to accept: 93 for PacBio CCS
export MINLEN=250

# Threading
export NUM_THREADS=64

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# DIRECTORY SETUP
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
export IN_DIR="./data/itsxrust_extracted"
export QUAL_FILTERED_DIR="./tmp"
export DENOISED_DIR="./data/denoised"
export CHIMERA_FILTERED_DIR="./data/chimera_filtered"

mkdir -p $IN_DIR $QUAL_FILTERED_DIR $DENOISED_DIR $CHIMERA_FILTERED_DIR

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# FUNCTIONS
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Quality filter trimmed reads using VSEARCH:
#   MAXEE_RATE — maximum expected errors per base (proportional, scales with
#                read length)
#   MAXN       — maximum number of ambiguous (N) bases allowed per read
#   QMAX       — maximum quality score to accept; 93 prevents VSEARCH from
#                rejecting high-quality PacBio CCS scores
qual_filter() {
    mkdir -p "$QUAL_FILTERED_DIR"

    echo "Quality filtering reads..."
    for file in "${IN_DIR}"/*[0-9].fastq.gz; do
        base_name=$(basename "$file" .fastq.gz)
        vsearch \
            --fastq_filter "$file" \
            --fastq_maxee_rate "$MAXEE_RATE" \
            --fastq_maxns "$MAXN" \
            --fastq_qmax "$QMAX" \
            --fastq_minlen "$MINLEN" \
            --fastqout "${QUAL_FILTERED_DIR}/${base_name}.fastq" 2>> "logs/logfile_vsearch.txt"
    done
    pigz -f -p "$NUM_THREADS" "${QUAL_FILTERED_DIR}"/*.fastq
}
# Denoise reads using DADA2 in PacBio CCS single-end mode:
#   dada2:::PacBioErrfun — CCS-specific error estimation function; replaces the
#                          default loessErrfun which was designed for Illumina
#                          (Callahan et al. 2019, Nucleic Acids Research)
#   BAND_SIZE=32         — wider alignment band to accommodate PacBio-length indels
#   pool=TRUE            — globally pools all samples before denoising; maximises
#                          sensitivity to rare variants by sharing read information
#                          across samples (slower but more sensitive than pseudo)
#   MIN_ABUNDANCE=1      — retains singletons
#
# Outputs: seqtab.rds, logfile_dada2.txt, seqtab.txt, all.fasta
denoise() {
    mkdir -p "$DENOISED_DIR" "$CHIMERA_FILTERED_DIR"
    echo "Denoising reads using DADA2 (PacBio CCS single-end)..."

    Rscript - <<EOF
        suppressPackageStartupMessages(require(dada2))
        suppressPackageStartupMessages(require(tidyverse))

        filtFs <- sort(list.files("$QUAL_FILTERED_DIR", pattern=".fastq.gz", full.names = TRUE))
        # Drop empty files (samples where all reads were filtered out)
        filtFs <- filtFs[file.size(filtFs) > 100]
        if (length(filtFs) == 0) stop("No reads remain after quality filtering.")
        sample.names <- gsub(".fastq.gz", "", basename(filtFs))
        names(filtFs) <- sample.names

        errF <- learnErrors(filtFs,
                            multithread             = $NUM_THREADS,
                            errorEstimationFunction = dada2::PacBioErrfun,
                            BAND_SIZE               = 32)

        dadaFs <- dada(filtFs,
                       err          = errF,
                       multithread  = $NUM_THREADS,
                       pool         = TRUE,
                       MIN_ABUNDANCE = 1,
                       BAND_SIZE    = 32)

        seqtab <- makeSequenceTable(dadaFs)

        getN <- function(x) sum(getUniques(x))
        track <- cbind(sapply(dadaFs, getN), rowSums(seqtab))
        colnames(track) <- c("denoised", "tabled")

        saveRDS(seqtab, file.path("$DENOISED_DIR", "seqtab.rds"))
        write.table(track, "logs/logfile_dada2.txt",
                    sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

        # Convert seqtab to tab-separated format with VSEARCH-style ;size= annotations
        # Header format: SAMPLENAME.N;size=N
        as.data.frame(seqtab) %>%
            rownames_to_column('sample.names') %>%
            pivot_longer(-sample.names, names_to='seq', values_to='size') %>%
            filter(size > 0) %>%
            ungroup() %>%
            mutate(seq.name=paste0(sample.names, ".", row_number(), ";size=", size)) %>%
            select(seq.name, seq) %>%
            write.table(., file=file.path("$DENOISED_DIR", "seqtab.txt"),
                        sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
EOF
    if [ $? -ne 0 ]; then
        echo "R script failed"
        exit 1
    fi

    # Convert tab-separated output to FASTA for downstream VSEARCH steps
    seqkit tab2fx "$DENOISED_DIR/seqtab.txt" -o "$CHIMERA_FILTERED_DIR/all.fasta"
}

# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# MAIN
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

echo "Starting denoising pipeline..."
echo "Parameters:"
echo "MAXEE_RATE: ${MAXEE_RATE}"
echo "MAXN:       ${MAXN}"
echo "QMAX:       ${QMAX}"
echo "MINLEN:     ${MINLEN}"
echo "NUM_THREADS: ${NUM_THREADS}"
echo "" 

echo "Input directory: ${IN_DIR}"
echo "Quality-filtered output directory: ${QUAL_FILTERED_DIR}"
echo "Denoised output directory: ${DENOISED_DIR}"
echo "Chimera-filtered output directory: ${CHIMERA_FILTERED_DIR}"
echo ""

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

echo "=== START: Quality filtering ==="
echo "$(date)"
qual_filter
echo "=== END: Quality filtering ==="
echo ""

echo "=== START: Denoising ==="
echo "$(date)"
denoise
echo "=== END: Denoising ==="
echo ""

echo "deactivating conda environment..."
conda deactivate

# Remove the tmp dir
rm -rf "$QUAL_FILTERED_DIR"

echo "=== PIPELINE COMPLETE ==="
echo "$(date)"