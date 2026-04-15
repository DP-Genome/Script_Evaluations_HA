#!/bin/bash
# align_sort_rmdup_robust.sh

set -euo pipefail

# -------------------------------
# Usage check
# -------------------------------
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <R1_fastq_directory>"
    exit 1
fi

INPUT_DIR="$1"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "[ERROR] Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# -------------------------------
# Configuration
# -------------------------------
REF="/mnt/server_share/Hanan/human_pipeline/References/human_references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

THREADS=8
MAX_JOBS=1
SORT_MEM="1G"

OUTDIR="./mapped_output"
Q30DIR="./mapped_output_q30"
LOGDIR="./logs"
TMPDIR="./tmp"

mkdir -p "$OUTDIR" "$Q30DIR" "$LOGDIR" "$TMPDIR"

PIPELINE_LOG="$LOGDIR/pipeline.log"

# -------------------------------
# Logging
# -------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$PIPELINE_LOG"
}

# -------------------------------
# Function: process one sample
# -------------------------------
map_sample() {
    R1="$1"
    SAMPLE=$(basename "$R1" _R1.fastq.gz)

    # Adjust this if your R2 path differs
    R2="../1_control_samples/R2_fastqs/${SAMPLE}_R2.fastq.gz"

    SAMPLE_LOG="$LOGDIR/${SAMPLE}.log"
    exec > >(tee -a "$SAMPLE_LOG") 2>&1

    echo "=============================="
    echo "Sample: $SAMPLE"
    echo "Start: $(date)"
    echo "=============================="

    if [[ ! -f "$R2" ]]; then
        echo "[ERROR] Missing R2: $R2"
        return 1
    fi

    SAM="$TMPDIR/${SAMPLE}.sam"
    NAME_SORT="$TMPDIR/${SAMPLE}.namesort.bam"
    FIXMATE="$TMPDIR/${SAMPLE}.fixmate.bam"
    COORD_SORT="$TMPDIR/${SAMPLE}.coordsort.bam"

    DEDUP="$OUTDIR/${SAMPLE}.dedup.bam"
    FILTERED="$Q30DIR/${SAMPLE}.bam"

    cleanup() {
        rm -f "$SAM" "$NAME_SORT" "$FIXMATE" "$COORD_SORT"
    }
    trap cleanup EXIT

    echo "[INFO] Running BWA-MEM"
    ~/bwa-0.7.19/bwa mem -t "$THREADS" "$REF" "$R1" "$R2" > "$SAM"

    echo "[INFO] Name sorting"
    samtools sort -n -@ "$THREADS" -m "$SORT_MEM" \
        -T "$TMPDIR/${SAMPLE}_namesort" \
        -o "$NAME_SORT" "$SAM"

    echo "[INFO] Running fixmate"
    samtools fixmate -m "$NAME_SORT" "$FIXMATE"

    echo "[INFO] Coordinate sorting"
    samtools sort -@ "$THREADS" -m "$SORT_MEM" \
        -T "$TMPDIR/${SAMPLE}_coordsort" \
        -o "$COORD_SORT" "$FIXMATE"

    echo "[INFO] Marking duplicates"
    samtools markdup -r -@ "$THREADS" "$COORD_SORT" "$DEDUP"

    echo "[INFO] Indexing deduplicated BAM"
    samtools index -@ "$THREADS" "$DEDUP"

    echo "[INFO] Filtering MAPQ ≥ 30"
    samtools view -@ "$THREADS" -b -q 30 "$DEDUP" > "$FILTERED"

    echo "[INFO] Indexing filtered BAM"
    samtools index -@ "$THREADS" "$FILTERED"

    echo "[INFO] Calculating coverage"

    cov_dedup=$(samtools depth -@ "$THREADS" -a "$DEDUP" | \
        awk '{sum+=$3} END {if (NR>0) printf "%.6f", sum/NR; else print 0}')

    cov_q30=$(samtools depth -@ "$THREADS" -a "$FILTERED" | \
        awk '{sum+=$3} END {if (NR>0) printf "%.6f", sum/NR; else print 0}')

    echo -e "${SAMPLE}\t${cov_dedup}\t${cov_q30}" > "coverage.${SAMPLE}.tmp"

    echo "[DONE] Finished $SAMPLE at $(date)"
}

export -f map_sample
export REF THREADS SORT_MEM TMPDIR OUTDIR Q30DIR LOGDIR

# -------------------------------
# Find input files
# -------------------------------
R1_FILES=$(find "$INPUT_DIR" -type f -name "*_R1.fastq.gz")

if [[ -z "$R1_FILES" ]]; then
    echo "[ERROR] No R1 files found in $INPUT_DIR"
    exit 1
fi

log "Found R1 files:"
echo "$R1_FILES" | tee -a "$PIPELINE_LOG"

# -------------------------------
# Run pipeline
# -------------------------------
log "Starting parallel processing"

parallel --bar -j "$MAX_JOBS" map_sample ::: $R1_FILES

log "All samples completed"

# -------------------------------
# Merge coverage
# -------------------------------
log "Merging coverage results"

echo -e "Sample\tDedup_Coverage\tMAPQ30_Coverage" > coverage_summary.tsv
cat coverage.*.tmp | sort -u >> coverage_summary.tsv
rm coverage.*.tmp

log "Pipeline finished successfully"
