#!/bin/bash
# ==============================================================================
# extract_metrics.sh
# Bash Orchestration Script for Low-Coverage Genotype Imputation Benchmarking
# ==============================================================================

set -euo pipefail

# --- Configuration & Paths ---
PROJECT_DIR="/mnt/server_share/Hanan/human_pipeline/middle_eastern_study"
REF_DIR="/mnt/server_share/Hanan/human_pipeline/References"
REF_FASTA="${REF_DIR}/human_references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
BAM_DIR="${PROJECT_DIR}/control_bams"
OUTPUT_DIR="${PROJECT_DIR}/validation_output_2"
TMP_DIR="${OUTPUT_DIR}/tmp"

# Reference Panels and Coverage Depths to evaluate
PANELS=("1kgp" "hgdp1kgp")
DEPTHS=("0.1x" "0.25x" "0.5x" "1x" "2x" "4x")
GP_FILTERS=("raw" "gp90" "gp99")

# Hardware & Multi-threading parameters (Tuned for 128 virtual CPUs / 64 Cores)
# We use 45 parallel jobs of single-threaded mpileup pipes to utilize 45 cores total.
MAX_PARALLEL_JOBS=45

# --- Software Discovery ---
echo "[$(date)] Discovering software tools..."
BCFTOOLS=$(which bcftools 2>/dev/null || find "${REF_DIR}/softwares" -name "bcftools" -type f -executable | head -n 1)
PARALLEL=$(which parallel 2>/dev/null || find "${REF_DIR}/softwares" -name "parallel" -type f -executable | head -n 1)

if [[ -z "$BCFTOOLS" ]]; then
    echo "ERROR: bcftools not found!" >&2
    exit 1
fi
if [[ -z "$PARALLEL" ]]; then
    echo "ERROR: GNU parallel not found!" >&2
    exit 1
fi

echo "Using BCFtools: $BCFTOOLS"
echo "Using GNU Parallel: $PARALLEL"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR"
mkdir -p "${PROJECT_DIR}/control_vcfs"
mkdir -p "${PROJECT_DIR}/imputed_GP_filtered_vcfs"

# --- Step 1: Consolidate BAM files ---
echo "[$(date)] Locating ground truth BAM files..."
BAM_LIST="${OUTPUT_DIR}/ground_truth_bams.txt"
find "$BAM_DIR" -name "*.bam" > "$BAM_LIST"
NUM_BAMS=$(wc -l < "$BAM_LIST")
echo "Found $NUM_BAMS ground truth BAM files."

if [[ $NUM_BAMS -eq 0 ]]; then
    echo "ERROR: No BAM files found in $BAM_DIR" >&2
    exit 1
fi

# --- Step 2: Extract sites to benchmark ---
echo "[$(date)] Extracting target variant sites for genotyping (union of both reference panels)..."
TARGET_HGDP="${PROJECT_DIR}/imputed_vcfs_hgdp1kgp/impute_1x_work/imputation_output/imputed_samples/imputed_middle_eastern_1x_biallelic.bcf"
TARGET_1KGP="${PROJECT_DIR}/imputed_vcfs_1kgp/impute_1x_work/imputation_output/imputed_samples/imputed_middle_eastern_1x_biallelic.bcf"
SITES_BED="${OUTPUT_DIR}/target_sites_union.bed"

if [[ -f "$TARGET_HGDP" && -f "$TARGET_1KGP" ]]; then
    echo "Extracting and merging sites from HGDP+1KGP and 1KGP panels..."
    $BCFTOOLS query -f '%CHROM\t%POS\t%POS\n' "$TARGET_HGDP" > "${OUTPUT_DIR}/hgdp_sites.tmp"
    $BCFTOOLS query -f '%CHROM\t%POS\t%POS\n' "$TARGET_1KGP" > "${OUTPUT_DIR}/1kgp_sites.tmp"
    # Sort genomic-wise by chromosome (version/numerical) and position (numerical) to match BCF order
    cat "${OUTPUT_DIR}/hgdp_sites.tmp" "${OUTPUT_DIR}/1kgp_sites.tmp" | sort -k1,1V -k2,2n -u > "$SITES_BED"
    rm -f "${OUTPUT_DIR}/hgdp_sites.tmp" "${OUTPUT_DIR}/1kgp_sites.tmp"
else
    echo "WARNING: Target imputed BCFs not found at standard locations. Finding any biallelic BCF to extract..."
    TARGET_BCF=$(find "${PROJECT_DIR}" -name "*_biallelic.bcf" | head -n 1)
    if [[ -n "$TARGET_BCF" ]]; then
        echo "Extracting sites from: $TARGET_BCF"
        $BCFTOOLS query -f '%CHROM\t%POS\t%POS\n' "$TARGET_BCF" | sort -k1,1V -k2,2n -u > "$SITES_BED"
    else
        echo "ERROR: Could not find any biallelic BCF to extract sites." >&2
        exit 1
    fi
fi
echo "Extracted $(wc -l < "$SITES_BED") target sites in union."

# --- Step 3: Genotype Ground Truth BAMs at Target Sites ---
export BCFTOOLS REF_FASTA BAM_LIST TMP_DIR

# Split target sites into chunks of 200,000 lines
echo "[$(date)] Splitting target sites into chunks of 200,000 lines..."
CHUNK_PREFIX="${TMP_DIR}/sites_chunk_"
split -l 200000 -d -a 4 "$SITES_BED" "$CHUNK_PREFIX"

# Get list of chunk files
CHUNK_FILES=( $(ls ${TMP_DIR}/sites_chunk_*) )
NUM_CHUNKS=${#CHUNK_FILES[@]}
echo "Split sites into $NUM_CHUNKS chunks."

call_variants_by_chunk() {
    local chunk_file="$1"
    local chunk_id=$(basename "$chunk_file" | sed 's/sites_chunk_//')
    local chunk_raw_bcf="${TMP_DIR}/truth_chunk_${chunk_id}.bcf"

    # Skip indels (-I) to speed up variant calling significantly, as they are filtered out in step 4 anyway.
    $BCFTOOLS mpileup -I -Ou -f "$REF_FASTA" -R "$chunk_file" -q 30 -Q 30 -a FORMAT/DP,FORMAT/AD -b "$BAM_LIST" | \
        $BCFTOOLS call -m -G - -Ob -o "$chunk_raw_bcf"
}

export -f call_variants_by_chunk

echo "[$(date)] Calling genotypes on ground truth BAMs in parallel across chunks..."
printf "%s\n" "${CHUNK_FILES[@]}" | $PARALLEL -j "$MAX_PARALLEL_JOBS" call_variants_by_chunk {}

# Concatenate chromosomal truth files
echo "[$(date)] Concatenating chunk truth files..."
TRUTH_RAW_BCF="${PROJECT_DIR}/control_vcfs/ground_truth_raw.bcf"
$BCFTOOLS concat -Ob -o "$TRUTH_RAW_BCF" "${TMP_DIR}"/truth_chunk_*.bcf
$BCFTOOLS index -f "$TRUTH_RAW_BCF"

# Filter to retain only biallelic SNPs and compute MAF using fill-tags plugin
echo "[$(date)] Filtering truth calls and calculating MAF..."
TRUTH_FILT_BCF="${PROJECT_DIR}/control_vcfs/ground_truth_biallelic.bcf"
$BCFTOOLS view -v snps -m2 -M2 "$TRUTH_RAW_BCF" -Ou | \
    $BCFTOOLS +fill-tags -Ou -- -t MAF | \
    $BCFTOOLS view -Ob -o "$TRUTH_FILT_BCF"
$BCFTOOLS index -f "$TRUTH_FILT_BCF"

# --- Step 4: Extract Genotypes from Ground Truth BCF ---
echo "[$(date)] Extracting tidy genotypes from ground truth..."
TRUTH_TSV="${OUTPUT_DIR}/truth_extracted.tsv.gz"
$BCFTOOLS query -f '[%CHROM:%POS:%REF:%ALT\t%SAMPLE\t%GT\t%INFO/MAF\n]' "$TRUTH_FILT_BCF" | gzip > "$TRUTH_TSV"
echo "Truth genotypes exported to $TRUTH_TSV"

# --- Step 5: Process Imputed BCFs and Apply GP Filters ---
extract_imputed_task() {
    local panel="$1"
    local depth="$2"
    local gp_filter="$3"
    local output_dir="$4"
    local project_dir="$5"
    local bcftools_bin="$6"

    local input_bcf="${project_dir}/imputed_vcfs_${panel}/impute_${depth}_work/imputation_output/imputed_samples/imputed_middle_eastern_${depth}_biallelic.bcf"
    local out_tsv="${output_dir}/imputed_${panel}_${depth}_${gp_filter}.tsv.gz"
    local out_bcf="${project_dir}/imputed_GP_filtered_vcfs/imputed_${panel}_${depth}_${gp_filter}.bcf"

    if [[ ! -f "$input_bcf" ]]; then
        echo "WARNING: File not found: $input_bcf. Skipping."
        return
    fi

    echo "Processing panel: $panel | depth: $depth | filter: $gp_filter..."

    if [[ "$gp_filter" == "raw" ]]; then
        $bcftools_bin view -Ob -o "$out_bcf" "$input_bcf"
    elif [[ "$gp_filter" == "gp90" ]]; then
        $bcftools_bin +setGT "$input_bcf" -Ob -o "$out_bcf" -- -t q -n . -i 'MAX(FMT/GP)<=0.90'
    elif [[ "$gp_filter" == "gp99" ]]; then
        $bcftools_bin +setGT "$input_bcf" -Ob -o "$out_bcf" -- -t q -n . -i 'MAX(FMT/GP)<=0.99'
    fi

    $bcftools_bin index -f "$out_bcf"
    $bcftools_bin query -f '[%CHROM:%POS:%REF:%ALT\t%SAMPLE\t%GT\t%GP\t%DS\n]' "$out_bcf" | gzip > "$out_tsv"
}

export -f extract_imputed_task

# Generate parameter task grid
TASK_FILE="${OUTPUT_DIR}/task_grid.txt"
> "$TASK_FILE"

for panel in "${PANELS[@]}"; do
    for depth in "${DEPTHS[@]}"; do
        for gp_filter in "${GP_FILTERS[@]}"; do
            echo "$panel $depth $gp_filter" >> "$TASK_FILE"
        done
    done
done

echo "[$(date)] Extracting imputed data in parallel..."
$PARALLEL --colsep ' ' -j "$MAX_PARALLEL_JOBS" \
    extract_imputed_task {1} {2} {3} "$OUTPUT_DIR" "$PROJECT_DIR" "$BCFTOOLS" < "$TASK_FILE"

# Cleanup intermediate files
rm -rf "$TMP_DIR"
rm -f "$TASK_FILE"

echo "=============================================================================="
echo " Bash extraction pipeline completed successfully!"
echo " Raw outputs located in: $OUTPUT_DIR"
echo " Filtered BCFs located in: ${PROJECT_DIR}/imputed_GP_filtered_vcfs"
echo "=============================================================================="
