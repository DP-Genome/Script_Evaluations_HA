#!/bin/bash
# ==============================================================================
# extract_metrics.sh
# Bash Orchestration Script for Low-Coverage Genotype Imputation Benchmarking
# ==============================================================================

set -euo pipefail

# Configuration & Paths
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
MAX_PARALLEL_JOBS=24

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



# --- Step 3: Genotype Ground Truth BAMs at Target Sites ---
HGDP1K_VCF_DIR="${REF_DIR}/glimpse2_references/hgdp1kgp/vcf_sites"
HGDP1K_TSV_DIR="${REF_DIR}/glimpse2_references/hgdp1kgp/tsv_sites"
KG_VCF_DIR="${REF_DIR}/glimpse2_references/1kgp/vcf_sites"
KG_TSV_DIR="${REF_DIR}/glimpse2_references/1kgp/tsv_sites"

call_truth_by_chr_panel() {
    local panel="$1"   # hgdp1kgp OR 1kgp
    local chr="$2"

    local vcf_dir tsv_dir

    if [[ "$panel" == "hgdp1kgp" ]]; then
        vcf_dir="$HGDP1K_VCF_DIR"
        tsv_dir="$HGDP1K_TSV_DIR"
    else
        vcf_dir="$KG_VCF_DIR"
        tsv_dir="$KG_TSV_DIR"
    fi

    local sites_vcf="${vcf_dir}/${panel}_chr${chr}.sites.vcf.gz"
    local alleles_tsv="${tsv_dir}/${panel}_chr${chr}.sites.tsv.gz"

    local out_bcf="${TMP_DIR}/truth_${panel}_chr${chr}.bcf"

    echo "[$(date)] Processing ${panel} chr${chr}"

    $BCFTOOLS mpileup -f "$REF_FASTA" -I -E -r "chr${chr}" -T "$sites_vcf" -q 30 -Q 30 -a FORMAT/DP,FORMAT/AD -Ou -b "$BAM_LIST" | \
    $BCFTOOLS call -Aim -C alleles -T "$alleles_tsv" -Ob -o "$out_bcf"
}

export -f call_truth_by_chr_panel
export BCFTOOLS REF_FASTA BAM_LIST TMP_DIR

# run variant calling for the truth files for each panel
echo "[$(date)] Running HGDP+1KGP truth calling..."
printf "%s\n" {1..22} | parallel -j "$MAX_PARALLEL_JOBS" call_truth_by_chr_panel hgdp1kgp {}
echo "[$(date)] Running 1KGP truth calling..."
printf "%s\n" {1..22} | parallel -j "$MAX_PARALLEL_JOBS" call_truth_by_chr_panel 1kgp {}

# Concatenate for each panel
TRUTH_HGDP1K="${PROJECT_DIR}/control_vcfs/truth_hgdp1kgp.bcf"
files_hgdp=()
for chr in {1..22}; do
    files_hgdp+=("${TMP_DIR}/truth_hgdp1kgp_chr${chr}.bcf")
done
$BCFTOOLS concat -Ob -o "$TRUTH_HGDP1K" "${files[@]}"
$BCFTOOLS index -f "$TRUTH_HGDP1K"

TRUTH_1KG="${PROJECT_DIR}/control_vcfs/truth_1kgp.bcf"
files_1kg=()
for chr in {1..22}; do
    files_1kg+=("${TMP_DIR}/truth_1kgp_chr${chr}.bcf")
done
$BCFTOOLS concat -Ob -o "$TRUTH_1KG" "${files[@]}"
$BCFTOOLS index -f "$TRUTH_1KG"


# Filter to retain only biallelic SNPs and compute MAF using fill-tags plugin
echo "[$(date)] Filtering truth calls and calculating MAF for the hgdp1kgp truth set..."
$BCFTOOLS norm -m+ "$TRUTH_HGDP1K" -Ou | \
    $BCFTOOLS view -v snps -m2 -M2 -Ou | \
    $BCFTOOLS +fill-tags -Ou -- -t MAF | \
    $BCFTOOLS view -Ob -o "${PROJECT_DIR}/control_vcfs/truth_hgdp1kgp_biallelic.bcf"
$BCFTOOLS index -f "${PROJECT_DIR}/control_vcfs/truth_hgdp1kgp_biallelic.bcf"


echo "[$(date)] Filtering truth calls and calculating MAF for the 1kgp truth set..."
$BCFTOOLS norm -m+ "$TRUTH_1KG" -Ou | \
    $BCFTOOLS view -v snps -m2 -M2 -Ou | \
    $BCFTOOLS +fill-tags -Ou -- -t MAF | \
    $BCFTOOLS view -Ob -o "${PROJECT_DIR}/control_vcfs/truth_1kgp_biallelic.bcf"
$BCFTOOLS index -f "${PROJECT_DIR}/control_vcfs/truth_1kgp_biallelic.bcf"


# --- Step 4: Extract Genotypes from Ground Truth BCF --- CHECK
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

# export -f extract_imputed_task

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

#echo "[$(date)] Extracting imputed data in parallel..."
#$PARALLEL --colsep ' ' -j "$MAX_PARALLEL_JOBS" \
#    extract_imputed_task {1} {2} {3} "$OUTPUT_DIR" "$PROJECT_DIR" "$BCFTOOLS" < "$TASK_FILE"

# Cleanup intermediate files
rm -rf "$TMP_DIR"
#rm -f "$TASK_FILE"

echo "=============================================================================="
echo " Bash extraction pipeline completed successfully!"
echo " Raw outputs located in: $OUTPUT_DIR"
echo " Filtered BCFs located in: ${PROJECT_DIR}/imputed_GP_filtered_vcfs"
echo "=============================================================================="
