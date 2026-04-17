#!/bin/bash
# align_sort_rmdup_MA_new.sh

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <R1_fastq_directory>"
    exit 1
fi

INPUT_DIR="$1"
REF="/mnt/server_share/Hanan/human_pipeline/References/human_references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# ---------------------------------------------------------
# HARDWARE TUNING (64 Cores / 128 Threads / 256GB RAM)
# ---------------------------------------------------------
# Processing 18 samples. 3 rounds of 6 jobs.
MAX_JOBS=6          # Process 6 human genomes simultaneously
THREADS=8           # 8 threads per tool in the pipe (3x8=24 threads per job)
                    # 6 jobs * 24 threads = 144 threads total (~1.1x oversubscription of 128 threads)
SORT_MEM="4G"       # 4GB per thread = 32GB per job for sorting.
                    # 6 jobs * (32GB sort + 6GB BWA index) = 228GB total RAM used.

OUTDIR="./mapped_output"
LOGDIR="./logs"
TMPDIR="./tmp"
mkdir -p "$OUTDIR" "$LOGDIR" "$TMPDIR"

# Tool paths (absolute paths for robustness)
SAMBAMBA="/mnt/server_share/Hanan/human_pipeline/References/softwares/sambamba-1.0.1"
MOSDEPTH="/mnt/server_share/Hanan/human_pipeline/References/softwares/mosdepth"

# Summary file
SUMMARY_FILE="coverage_summary_q30.txt"
echo -e "Sample\tCoverage_Before_Q30\tCoverage_After_Q30" > "$SUMMARY_FILE"

# ---------------------------------------------------------
# Pipeline Function
# ---------------------------------------------------------
map_sample() {
    R1="$1"
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

    SAMPLE_LOG="$LOGDIR/${SAMPLE}.log"
    exec > >(tee -a "$SAMPLE_LOG") 2>&1

    echo "[$(date)] Starting $SAMPLE"

    # Define paths
    DEDUP_BAM="$OUTDIR/${SAMPLE}.dedup.bam"
    FINAL_BAM="$OUTDIR/${SAMPLE}.q30.bam"
    RG_STR="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"

    # STEP 1: Align, Sort, and MarkDuplicates
    # We use an intermediate sorted BAM for sambamba as it requires seeking.
    # We chain these with && to maintain logical flow in a single "pipeline" operation.
    SORTED_BAM="$TMPDIR/${SAMPLE}.sorted.bam"

    ~/bwa-0.7.19/bwa mem -t "$THREADS" -R "$RG_STR" "$REF" "$R1" "$R2" | \
        samtools sort -@ "$THREADS" -m "$SORT_MEM" -T "$TMPDIR/${SAMPLE}_sort" -o "$SORTED_BAM" - && \
        "$SAMBAMBA" markdup -t "$THREADS" "$SORTED_BAM" "$DEDUP_BAM" && \
        samtools index -@ "$THREADS" "$DEDUP_BAM" && \
        rm -f "$SORTED_BAM"

    echo "[$(date)] Calculating coverage before Q30 for $SAMPLE"
    # STEP 2: Calculate coverage before Q30
    "$MOSDEPTH" -t "$THREADS" -x -n "$OUTDIR/${SAMPLE}_before" "$DEDUP_BAM"
    COV_BEFORE=$(grep "total" "$OUTDIR/${SAMPLE}_before.mosdepth.summary.txt" | awk '{print $4}')

    echo "[$(date)] Calculating coverage after Q30 for $SAMPLE"
    # STEP 3: Calculate coverage after Q30 (using -Q 30 flag in mosdepth)
    "$MOSDEPTH" -t "$THREADS" -x -n -Q 30 "$OUTDIR/${SAMPLE}_after" "$DEDUP_BAM"
    COV_AFTER=$(grep "total" "$OUTDIR/${SAMPLE}_after.mosdepth.summary.txt" | awk '{print $4}')

    echo "[$(date)] Filtering for Q30 and indexing for $SAMPLE"
    # STEP 4: Filter out reads below Q30 for the final output
    samtools view -@ "$THREADS" -b -q 30 -o "$FINAL_BAM" "$DEDUP_BAM"
    samtools index -@ "$THREADS" "$FINAL_BAM"

    # Save results to a temporary file (to be merged later to avoid race conditions)
    echo -e "${SAMPLE}\t${COV_BEFORE}\t${COV_AFTER}" > "$OUTDIR/${SAMPLE}.cov_stats"

    # Cleanup intermediate BAM if desired (keeping the q30.bam as the final)
    rm -f "$DEDUP_BAM" "$DEDUP_BAM.bai"
    rm -f "$OUTDIR/${SAMPLE}_before"* "$OUTDIR/${SAMPLE}_after"*

    echo "[$(date)] Finished $SAMPLE"
}

export -f map_sample
export REF THREADS SORT_MEM TMPDIR OUTDIR LOGDIR SAMBAMBA MOSDEPTH

# ---------------------------------------------------------
# Execution
# ---------------------------------------------------------
# Find files and pass safely to parallel
find "$INPUT_DIR" -type f -name "*_R1.fastq.gz" -print0 | \
    parallel -0 --bar -j "$MAX_JOBS" map_sample

# Combine all coverage stats into the final summary file
cat "$OUTDIR"/*.cov_stats >> "$SUMMARY_FILE"
rm -f "$OUTDIR"/*.cov_stats

echo "All samples completed. Results in $SUMMARY_FILE"
