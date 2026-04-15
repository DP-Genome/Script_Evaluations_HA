#!/bin/bash
# align_sort_rmdup_MAsh

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <R1_fastq_directory>"
    exit 1
fi

INPUT_DIR="$1"
REF="/mnt/server_share/Hanan/human_pipeline/References/human_references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# ---------------------------------------------------------
# HARDWARE TUNING (64 Cores / 256GB RAM)
# ---------------------------------------------------------
MAX_JOBS=4          # Process 4 human genomes simultaneously
THREADS=15          # 15 threads per genome (60 cores total, 4 for OS)
SORT_MEM="3G"       # 3GB per thread = 45GB per job = 180GB total RAM used.

OUTDIR="./mapped_output"
LOGDIR="./logs"
TMPDIR="./tmp"
mkdir -p "$OUTDIR" "$LOGDIR" "$TMPDIR"

# ---------------------------------------------------------
# Pipeline Function
# ---------------------------------------------------------
map_sample() {
    R1="$1"
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    
    # Dynamically find R2 in the same directory structure
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

    SAMPLE_LOG="$LOGDIR/${SAMPLE}.log"
    exec > >(tee -a "$SAMPLE_LOG") 2>&1

    echo "[$(date)] Starting $SAMPLE"

    # Define paths
    COORD_SORT="$TMPDIR/${SAMPLE}.coordsort.bam"
    DEDUP="$OUTDIR/${SAMPLE}.dedup.bam"
    RG_STR="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"

    # STEP 1: Align & Coordinate Sort in one streamed command
    # This keeps data in RAM and writes directly to the sorted BAM
    ~/bwa-0.7.19/bwa mem -t "$THREADS" -R "$RG_STR" "$REF" "$R1" "$R2" | \
        samtools sort -@ "$THREADS" -m "$SORT_MEM" -T "$TMPDIR/${SAMPLE}_sort" -o "$COORD_SORT" -

    # STEP 2: Mark Duplicates with Sambamba (faster, no name-sort required)
    sambamba markdup -r -t "$THREADS" "$COORD_SORT" "$DEDUP"

    # STEP 3: Index the final BAM
    samtools index -@ "$THREADS" "$DEDUP"

    # Cleanup intermediate BAM to save disk space
    rm -f "$COORD_SORT"

    echo "[$(date)] Finished $SAMPLE"
}
export -f map_sample
export REF THREADS SORT_MEM TMPDIR OUTDIR LOGDIR

# ---------------------------------------------------------
# Execution
# ---------------------------------------------------------
# Find files and pass safely to parallel
find "$INPUT_DIR" -type f -name "*_R1.fastq.gz" -print0 | \
    parallel -0 --bar -j "$MAX_JOBS" map_sample

echo "All samples completed."
