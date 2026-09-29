#!/bin/bash
# Version: 1.1.1
#   1.1.1 - ldprio v1.5 compatible (see pca_admixture.sh).
#   1.1.0 - new --transversions true|false: restrict PCA/Admixture to transversion SNPs.
#   1.0.1 - fixed admixture/PCA prefix bugs with --prune ldprio|indep, single seed shared by all
#           pileupCaller outputs, --cores option (replaces hard-coded threads), safe config
#           loader (no eval), input validation, exported log_msg in all variants.
#   1.0.0 - initial tagged release (see OLD_VERSION_SCRIPTS/)
set -euo pipefail
VERSION="1.1.1"

# ======== Logger setup ========
LOG_FILE="pipeline_$(date +%Y%m%d_%H%M%S).log"
export LOG_FILE
log_msg() {
    local msg="$1"
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo "[$timestamp] $msg" | tee -a "$LOG_FILE"
}
export -f log_msg
# ======== Set fixed paths ========
REF_DIR="/mnt/server_share/Hanan/human_pipeline/References"

# ======== Software discovery ========
log_msg "pileupCaller_glimpse.sh version ${VERSION}"
log_msg "Discovering software tools..."
SOFTWARES_DIR="${REF_DIR}/softwares"

# Helper to find binary
find_tool() {
    local name="$1"
    local path=$(find "$SOFTWARES_DIR" -name "$name" -type f -executable | head -n 1)
    if [[ -z "$path" ]]; then
        # Fallback to system path
        if command -v "$name" >/dev/null 2>&1; then
            echo "$name"
        else
            log_msg "ERROR: Tool '$name' not found in $SOFTWARES_DIR or system PATH." >&2
            exit 1
        fi
    else
        echo "$path"
    fi
}

SAMTOOLS=$(find_tool "samtools")
BCFTOOLS=$(find_tool "bcftools")
BGZIP=$(find_tool "bgzip")
TABIX=$(find_tool "tabix")
PARALLEL=$(find_tool "parallel")
PLINK_BIN=$(find_tool "plink")
PILEUPCALLER_BIN=$(find_tool "pileupCaller-linux")
ADMIXTURE_BIN=$(find_tool "admixture")

# GLIMPSE tools (prefer static binaries)
GLIMPSE2_PHASE=$(find "$SOFTWARES_DIR"/GLIMPSE2 -name "GLIMPSE2_phase_static" -type f -executable | head -n 1)
GLIMPSE2_LIGATE=$(find "$SOFTWARES_DIR"/GLIMPSE2 -name "GLIMPSE2_ligate_static" -type f -executable | head -n 1)
GLIMPSE1_PHASE=$(find "$SOFTWARES_DIR"/GLIMPSE1 -name "GLIMPSE_phase_static" -type f -executable | head -n 1)
GLIMPSE1_LIGATE=$(find "$SOFTWARES_DIR"/GLIMPSE1 -name "GLIMPSE_ligate_static" -type f -executable | head -n 1)

# EIGENSTRAT
CONVERTF_BIN=$(find "$SOFTWARES_DIR"/EIG-8.0.0 -name "convertf" -type f -executable | head -n 1)
SMARTPCA_BIN=$(find "$SOFTWARES_DIR"/EIG-8.0.0 -name "smartpca" -type f -executable | head -n 1)

for tool_var in GLIMPSE2_PHASE GLIMPSE2_LIGATE GLIMPSE1_PHASE GLIMPSE1_LIGATE CONVERTF_BIN SMARTPCA_BIN; do
    if [[ -z "${!tool_var}" ]]; then
        log_msg "ERROR: could not find $tool_var under $SOFTWARES_DIR."
        exit 1
    fi
done

log_msg "Software discovery completed."

# ======= Process arguments =======

# STEP1: Parse long options - set defaults
INPUT_FILE=""
PSEUDOHAP="true"
PANEL="wgs"
IMPUTE="false"
PLOT="false"
XCHR="true"
THREADS=2      # threads per GLIMPSE chunk job
CORES=24       # total cores for samtools, bcftools, smartpca and admixture
PRUNE="ldprio"
SEED=12345     # one seed shared by all pileupCaller runs so EIGENSTRAT/PLINK/VCF agree
TRANSVERSIONS="false"
export CORES PRUNE

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)     INPUT_FILE="$2"; shift 2 ;;
        --pseudoHap) PSEUDOHAP="$2";  shift 2 ;;
        --panel)     PANEL="$2";      shift 2 ;;
        --impute)    IMPUTE="$2";     shift 2 ;;
        --plot)      PLOT="$2";       shift 2 ;;
        --xchr)      XCHR="$2";       shift 2 ;;
        --threads)   THREADS="$2";    shift 2 ;;
        --cores)     CORES="$2";      shift 2 ;;
        --prune)     PRUNE="$2";      shift 2 ;;
        --seed)      SEED="$2";       shift 2 ;;
        --transversions) TRANSVERSIONS="$2"; shift 2 ;;
        --version)   echo "$VERSION"; exit 0 ;;
        --help)
            echo "Usage: $0 --input file.txt --pseudoHap true|false --panel twist|wgs --impute true|false --plot false|pca|admixture|both --xchr true|false --prune indep|ldprio --threads 2 --cores 24 --seed 12345 --transversions true|false [--version]"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Validate option values
case "$PLOT"  in false|pca|admixture|both) ;; *) log_msg "Error: --plot must be false|pca|admixture|both"; exit 1 ;; esac
case "$PRUNE" in indep|ldprio) ;;             *) log_msg "Error: --prune must be indep|ldprio"; exit 1 ;; esac
case "$TRANSVERSIONS" in true|false) ;; *) log_msg "Error: --transversions must be true|false"; exit 1 ;; esac
case "$PANEL" in wgs|twist) ;;                *) log_msg "Error: --panel must be twist|wgs"; exit 1 ;; esac

# STEP 2: Validate input/config file
if [[ -z "$INPUT_FILE" ]]; then
    log_msg "Error: --input file is required. For more information see README file."
    exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    log_msg "Error: '$INPUT_FILE' not found."
    exit 1
fi

# STEP 3: Load input/config file variables
log_msg "Loading parameters from $INPUT_FILE..."
while IFS='=' read -r key value || [[ -n "$key" ]]; do
    [[ $key =~ ^# ]] && continue   # skip comments
    [[ -z $key ]] && continue      # skip empty lines
    # Trim surrounding whitespace and carriage returns (no eval, no xargs)
    key="${key//$'\r'/}";     key="${key#"${key%%[![:space:]]*}"}";     key="${key%"${key##*[![:space:]]}"}"
    value="${value//$'\r'/}"; value="${value#"${value%%[![:space:]]*}"}"; value="${value%"${value##*[![:space:]]}"}"
    [[ -z "$key" ]] && continue
    if [[ ! "$key" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]; then
        log_msg "Error: invalid variable name '$key' in $INPUT_FILE"
        exit 1
    fi
    printf -v "$key" '%s' "$value"
done < "$INPUT_FILE"

for req in INPUT_BAM SAMPLE_NAME SAMPLE_PREFIX; do
    if [[ -z "${!req:-}" ]]; then
        log_msg "Error: $req is missing from $INPUT_FILE"
        exit 1
    fi
done

# STEP 4: Clean the input bam and save in list
log_msg "Filtering BAM files..."
BAM_DIR="clean_bams"
mkdir -p "$BAM_DIR"

BAM_LIST="bam_list.txt"
> "$BAM_LIST" # clear existing file

IFS=',' read -ra ADDR <<< "$INPUT_BAM"

# Loop through the resulting array
declare -A SEEN_BAM_NAMES=()
for input_bam in "${ADDR[@]}"; do
    filename=$(basename "$input_bam")
    output_bam="$BAM_DIR/${filename%.bam}_filtered.bam"

    if [[ -n "${SEEN_BAM_NAMES[$filename]:-}" ]]; then
        log_msg "ERROR: duplicate BAM file name '$filename' would overwrite its filtered output in $BAM_DIR."
        exit 1
    fi
    SEEN_BAM_NAMES[$filename]=1

    log_msg "Processing $input_bam -> $output_bam"
    $SAMTOOLS view -@ "$CORES" -b -h -F 0x400 -q 30 -o "$output_bam" "$input_bam" &>> "$LOG_FILE"
    echo "./$output_bam" >> "$BAM_LIST"
    $SAMTOOLS index "$output_bam" &>> "$LOG_FILE"
done
log_msg "Filtering completed. Filtered BAMs are in: $BAM_DIR"

# ======== PART 1: CREATING A PSEUDOHAPLOID ========

if [[ "$PSEUDOHAP" == "true" ]]; then

    log_msg "Starting PseudoHaploid creation..."

    # STEP 1: Choose BED/SNP files for mpileup/PileupCaller
    if [[ "$PANEL" == "wgs" ]]; then
        log_msg "Panel: WGS"
       BED_FILE="${REF_DIR}/hgdp_references/4M_hgdp_all_chr_withX.bed"
       SNP_FILE="${REF_DIR}/hgdp_references/eigenstrat_files/hgdp_all_chr_withX.snp"
    elif [[ "$PANEL" == "twist" ]]; then
        log_msg "Panel: TWIST 1240k"
        BED_FILE="${REF_DIR}/twist_references/twist_1240k_hg38.bed"
        SNP_FILE="${REF_DIR}/twist_references/twist_1240k_hg38.snp"
    else
        log_msg "Error: --panel 'twist' or 'wgs' required."
        exit 1
    fi

    # STEP 2: Run Pileupcaller
    FASTA_REF="${REF_DIR}/human_references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
    PILEUP_CALLER_DIR="./pileupCaller_output"
    mkdir -p "$PILEUP_CALLER_DIR"

    MPILEUP_CMD="$SAMTOOLS mpileup -B -R -q 25 -Q 20 -f $FASTA_REF -b $BAM_LIST -l $BED_FILE"

    # The same --seed is used for every output format so all three contain identical calls
    log_msg "Running PileupCaller for EigenStrat output (seed $SEED)..."
    ($MPILEUP_CMD | grep -v "^chrM" | $PILEUPCALLER_BIN --randomHaploid --seed "$SEED" --sampleNames "$SAMPLE_NAME" \
        --samplePopName pseudo --snpFile "$SNP_FILE" --eigenstratOut "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}") &>> "$LOG_FILE"

    # Rename chromosomes and positions in SNP file
    awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $2); $1=$2"_"$4; print}' \
        "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.snp" > tmp && mv tmp "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.snp"

    log_msg "Running PileupCaller for PLINK output (seed $SEED)..."
    ($MPILEUP_CMD | grep -v "^chrM" | $PILEUPCALLER_BIN --randomHaploid --seed "$SEED" --sampleNames "$SAMPLE_NAME" \
        --samplePopName pseudo --snpFile "$SNP_FILE" --plinkOut "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}") &>> "$LOG_FILE"

    # Rename chromosomes and positions in BIM file
    awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $1); $2=$1"_"$4; print}' \
        "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.bim" > tmp && mv tmp "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.bim"

    # create an autosomal version for plotting
    $PLINK_BIN --bfile "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}" --not-chr X --make-bed --out "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}_noX" &>> "$LOG_FILE"

    log_msg "Running PileupCaller for VCF output (seed $SEED)..."
    ($MPILEUP_CMD | grep -v "^chrM" | $PILEUPCALLER_BIN --randomHaploid --seed "$SEED" --sampleNames "$SAMPLE_NAME" --samplePopName pseudoHaploid \
    --snpFile "$SNP_FILE" --vcf | $BGZIP > "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.vcf.gz") &>> "$LOG_FILE"
    $TABIX -p vcf "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.vcf.gz" &>> "$LOG_FILE"

    log_msg "✅ PseudoHaploid files created for $SAMPLE_PREFIX"
fi

# ======== PART 2: IMPUTAION using GLIMPSE ========

if [[ "$IMPUTE" == "true" ]]; then
    log_msg "Starting GLIMPSE imputation..."

    # Number of concurrent chromosome and chunk jobs
    MAX_CHR_PARALLEL=2
    MAX_CHUNK_PARALLEL=8

    # set paths
    IMPUTE_DIR="./imputation_output"
    SPLIT_DIR="${REF_DIR}/glimpse2_references/hgdp1kgp/ref_panel_split"
    IMPUTED_CHUNKS_DIR="${IMPUTE_DIR}/imputed_chunks"
    CHUNKS_DIR="${REF_DIR}/glimpse2_references/hgdp1kgp/chunks"
    mkdir -p "$IMPUTE_DIR" "$IMPUTED_CHUNKS_DIR"

    # Add sample name column in the BAM_LIST file
    IFS=',' read -r -a SAMPLE_ARRAY <<< "$SAMPLE_NAME"
    NAMES_FILE="${IMPUTE_DIR}/sample_names.txt"
    > "$NAMES_FILE"

    i=0
    while read -r bam; do
        sample="${SAMPLE_ARRAY[i]}"
        echo -e "$bam\t$sample" >> "$NAMES_FILE"
        ((i++)) || true
    done < "$BAM_LIST"

    export NAMES_FILE SPLIT_DIR IMPUTED_CHUNKS_DIR CHUNKS_DIR MAX_CHUNK_PARALLEL GLIMPSE2_PHASE THREADS LOG_FILE

    # function to impute each autosomal chunk
    process_chunk_autosome() {
      LINE="$1"
      CHR=$(echo "$LINE" | awk '{print $2}')
      IRG=$(echo "$LINE" | awk '{print $3}')
      ORG=$(echo "$LINE" | awk '{print $4}')
      REGS=$(echo "$IRG" | cut -d: -f2 | cut -d- -f1)
      REGE=$(echo "$IRG" | cut -d: -f2 | cut -d- -f2)

      REF_BIN="${SPLIT_DIR}/hgdp1kgp_${CHR}_${REGS}_${REGE}.bin"
      OUT_BCF="${IMPUTED_CHUNKS_DIR}/imputed_${CHR}_${REGS}_${REGE}.bcf"

      $GLIMPSE2_PHASE \
        --bam-list "$NAMES_FILE" \
        --reference "$REF_BIN" \
        --output "$OUT_BCF" \
        --thread "$THREADS" &>> "$LOG_FILE"
    }
    export -f process_chunk_autosome

    process_chr() {
      CHR="$1"
      CHUNKS_FILE="${CHUNKS_DIR}/chunks.chr${CHR}.txt"
      log_msg "    Running GLIMPSE2 imputation for chromosome $CHR..."
      cat "$CHUNKS_FILE" | parallel -j "$MAX_CHUNK_PARALLEL" process_chunk_autosome
    }
    export -f process_chr

    log_msg "STEP 1/4: Running imputation for the autosomes..."
    $PARALLEL -j "$MAX_CHR_PARALLEL" process_chr ::: $(seq 1 22)
    log_msg "✅ Per-chunk autosomal imputation completed."

    # Ligate each chromosome
    log_msg "STEP 2/4: Ligating imputed chunks per autosome..."
    CHR_LIG_PARALLEL=22
    IMPUTED_CHR_DIR="${IMPUTE_DIR}/imputed_chromosomes"
    LISTS_DIR="${IMPUTED_CHR_DIR}/list_files"
    mkdir -p "$IMPUTED_CHR_DIR" "$LISTS_DIR"

    export IMPUTED_CHR_DIR LISTS_DIR BCFTOOLS GLIMPSE2_LIGATE

    ligate_chr() {
      local CHR="$1"
      echo "    Ligating chromosome $CHR..."
      LIST_FILE="${LISTS_DIR}/list.chr${CHR}.txt"
      mapfile -t IMPUTED_FILES < <(ls -1v "${IMPUTED_CHUNKS_DIR}/imputed_chr${CHR}_"*.bcf 2>/dev/null || true)
      printf "%s\n" "${IMPUTED_FILES[@]}" > "$LIST_FILE"

      $GLIMPSE2_LIGATE --input "$LIST_FILE" --output "${IMPUTED_CHR_DIR}/imputed_chr${CHR}_ligated.bcf" &>> "$LOG_FILE"
      $BCFTOOLS index -f "${IMPUTED_CHR_DIR}/imputed_chr${CHR}_ligated.bcf" &>> "$LOG_FILE"
    }
    export -f ligate_chr
    $PARALLEL -j "$CHR_LIG_PARALLEL" ligate_chr ::: $(seq 1 22)
    log_msg "✅ Per-chromosome autosomal ligation completed."

    # Integrated X-chromosome imputation logic
    if [[ "$XCHR" == "true" ]]; then
        log_msg "STEP 3/4: Starting imputation for the X chromosome (GLIMPSE1)..."
        FASTA_REF="${REF_DIR}/human_references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

        IFS=',' read -r -a NAME_ARRAY <<< "$SAMPLE_NAME"
        IFS=',' read -r -a SEX_ARRAY <<< "$SAMPLE_SEX"

        GL_PLOIDY="${IMPUTE_DIR}/GL_ploidy.txt"
        > "$GL_PLOIDY"

        for i in "${!NAME_ARRAY[@]}"; do
            NAME="${NAME_ARRAY[i]}"
            SEX="${SEX_ARRAY[i]}"
            [[ "$SEX" == "M" ]] && PLOIDY=1 || PLOIDY=2
            echo -e "chrX\t2781514\t155700519\t${SEX}\t${PLOIDY}" >> "$GL_PLOIDY"
        done

        log_msg "Computing genotype likelihoods for chrX..."
        SITES_DIR="${REF_DIR}/glimpse1_references/sites_files"
        OUT_VCF_GL="${IMPUTE_DIR}/${SAMPLE_PREFIX}_GL_chrX.vcf.gz"

        ($BCFTOOLS mpileup -f "$FASTA_REF" -b "$BAM_LIST" --ignore-RG -I -E -a FORMAT/DP -T "${SITES_DIR}/hgdp1kgp_chrX.sites.vcf.gz" -r chrX -Ou | \
            $BCFTOOLS call --ploidy GRCh38 --ploidy-file "$GL_PLOIDY" -Aim -C alleles -T "${SITES_DIR}/hgdp1kgp_chrX.sites.tsv.gz" -Oz -o "$OUT_VCF_GL") &>> "$LOG_FILE"

        $BCFTOOLS reheader -s "$NAMES_FILE" "$OUT_VCF_GL" -o "$OUT_VCF_GL.tmp" &>> "$LOG_FILE" && mv "$OUT_VCF_GL.tmp" "$OUT_VCF_GL"
        $BCFTOOLS index -f "$OUT_VCF_GL" &>> "$LOG_FILE"

        SAMPLE_PLOIDY="${IMPUTE_DIR}/sample_ploidy.txt"
        > "$SAMPLE_PLOIDY"
        for i in "${!NAME_ARRAY[@]}"; do
            [[ "${SEX_ARRAY[$i]}" == "M" || "${SEX_ARRAY[$i]}" == "m" ]] && PL=1 || PL=2
            printf "%s\t%d\n" "${NAME_ARRAY[$i]}" "$PL" >> "$SAMPLE_PLOIDY"
        done

        log_msg "Imputing X chromosome chunks..."
        export GLIMPSE1_PHASE REF_DIR OUT_VCF_GL SAMPLE_PLOIDY IMPUTED_CHUNKS_DIR BCFTOOLS THREADS LOG_FILE
        process_chunk_chrX() {
          LINE="$1"
          CHR=$(echo "$LINE" | awk '{print $2}')
          IRG=$(echo "$LINE" | awk '{print $3}')
          ORG=$(echo "$LINE" | awk '{print $4}')
          REGS=$(echo "$IRG" | cut -d: -f2 | cut -d- -f1)
          REGE=$(echo "$IRG" | cut -d: -f2 | cut -d- -f2)

          REF_X_BCF="${REF_DIR}/glimpse1_references/ref_panel_bcf/hgdp1kgp_chrX.bcf"
          MAP_X="${REF_DIR}/glimpse2_references/map_files_hg38/chrX.b38.gmap.gz"
          OUT_BCF="${IMPUTED_CHUNKS_DIR}/imputed_${CHR}_${REGS}_${REGE}.bcf"

          $GLIMPSE1_PHASE --input "$OUT_VCF_GL" --reference "$REF_X_BCF" --map "$MAP_X" \
              --input-region "${IRG}" --output-region "${ORG}" --samples-file "$SAMPLE_PLOIDY" \
              --output "$OUT_BCF" --thread "$THREADS" &>> "$LOG_FILE"
          $BCFTOOLS index -f "$OUT_BCF" &>> "$LOG_FILE"
        }
        export -f process_chunk_chrX
        CHUNKS_DIR_X="${REF_DIR}/glimpse2_references/hgdp1kgp/chunks"
        [[ -f "${CHUNKS_DIR_X}/chunks.chrX.txt" ]] || echo "Missing chrX chunks"

        cat "${CHUNKS_DIR_X}/chunks.chrX.txt" | $PARALLEL -j "$MAX_CHUNK_PARALLEL" process_chunk_chrX

        LST_X="${LISTS_DIR}/list.chrX.txt"
        ls "${IMPUTED_CHUNKS_DIR}"/imputed_chrX_*.bcf > "${LST_X}"
        $GLIMPSE1_LIGATE --input "${LST_X}" --output "${IMPUTED_CHR_DIR}/imputed_chrX_ligated.bcf" &>> "$LOG_FILE"
        $BCFTOOLS index -f "${IMPUTED_CHR_DIR}/imputed_chrX_ligated.bcf" &>> "$LOG_FILE"
        log_msg "✅ X chromosome imputation completed."
    fi

    # concatenate imputed chromosomes into a single file
    log_msg "STEP 4/4: Concatenating and filtering imputed chromosomes..."
    IMPUTED_SAMPLE_DIR="${IMPUTE_DIR}/imputed_samples"
    mkdir -p "$IMPUTED_SAMPLE_DIR"

    OUT_CONCAT_BCF="${IMPUTED_SAMPLE_DIR}/imputed_${SAMPLE_PREFIX}.bcf"
    CHR_LIG_BCF=$(ls -1v "${IMPUTED_CHR_DIR}/imputed_chr"*"_ligated.bcf")

    $BCFTOOLS concat -Ob -o "$OUT_CONCAT_BCF" $CHR_LIG_BCF &>> "$LOG_FILE"
    $BCFTOOLS index -f "$OUT_CONCAT_BCF" &>> "$LOG_FILE"

    # apply mask
    MASK_BED="${REF_DIR}/hgdp_references/hgdp_wgs.20190516.mask.bed"
    MASKED_BCF="${IMPUTED_SAMPLE_DIR}/imputed_${SAMPLE_PREFIX}_masked.bcf"
    $BCFTOOLS view --regions-file "$MASK_BED" --threads "$CORES" -Ob -o "$MASKED_BCF" "$OUT_CONCAT_BCF" &>> "$LOG_FILE"
    $BCFTOOLS index "$MASKED_BCF" &>> "$LOG_FILE"

    # filter for biallelic snps only
    BIALLELIC_BCF="${IMPUTED_SAMPLE_DIR}/imputed_${SAMPLE_PREFIX}_biallelic.bcf"
    ($BCFTOOLS norm -m+ "$MASKED_BCF" -Ou | $BCFTOOLS view -v snps -m2 -M2 -Ob -o "$BIALLELIC_BCF") &>> "$LOG_FILE"
    $BCFTOOLS index "$BIALLELIC_BCF" &>> "$LOG_FILE"

    # set GP<0.99 as missing and filter out MAF<0.1
    FILTERED_BCF="${IMPUTED_SAMPLE_DIR}/imputed_${SAMPLE_PREFIX}_filtered.bcf"
    ($BCFTOOLS +setGT "$BIALLELIC_BCF" -- -t q -n . -i 'MAX(FMT/GP)<0.99' | $BCFTOOLS view -i 'RAF>0.1 & RAF<0.9' -Ob -o "$FILTERED_BCF") &>> "$LOG_FILE"
    $BCFTOOLS index "$FILTERED_BCF" &>> "$LOG_FILE"

    # Generate plink files and rename chromosomes and positions in BIM file
    IMPUTED_PLINK_DIR="${IMPUTE_DIR}/plink_files"
    mkdir -p "$IMPUTED_PLINK_DIR"
    $PLINK_BIN --bcf "$FILTERED_BCF" --const-fid imputed --make-bed --out "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered" &>> "$LOG_FILE"
    awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $1); $2=$1"_"$4; print}' \
        "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered.bim" > tmp && mv tmp "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered.bim"
    # create an autosomal version for plotting
    $PLINK_BIN --bfile "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered" --not-chr X --make-bed --out "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered_noX" &>> "$LOG_FILE"

    log_msg "✅ Imputation and filtering completed successfully."
    fi

# ======== PART 3: Plotting ========
if [[ "$PLOT" != "false" ]]; then
    log_msg "Launching PCA/Admixture analysis..."
    # Export discovered bins for sub-script
    export PLINK_BIN CONVERTF_BIN SMARTPCA_BIN ADMIXTURE_BIN CORES PRUNE TRANSVERSIONS
    source "${REF_DIR}/../Pipelines/pca_admixture.sh"
fi

log_msg "Pipeline execution finished."
