#!/bin/bash
set -euo pipefail

# ======== Set fixed paths ========

REF_DIR="/mnt/server_share/Hanan/human_pipeline/References"

# ======= Process arguments =======

# STEP1: Parse long options - set defaults

INPUT_FILE=""
PSEUDOHAP="true"
PANEL="wgs"
IMPUTE="false"
PLOT="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)     INPUT_FILE="$2"; shift 2 ;;
        --pseudoHap) PSEUDOHAP="$2";  shift 2 ;;
        --panel)     PANEL="$2";      shift 2 ;;
        --impute)    IMPUTE="$2";     shift 2 ;;
        --plot)      PLOT="$2";       shift 2 ;;
        --help)
            echo "Usage: $0 --input file.txt --pseudoHap true|false --panel twist|wgs --impute true|false --plot false|pca|admixture|both"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# STEP 2: Validate input/config file

if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: --input file is required. For more information see README file."
    exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: '$INPUT_FILE' not found."
    exit 1
fi

# STEP 3: Load input/config file variables
echo
echo "Parameters provided:"
echo "--------------------"
while IFS='=' read -r key value; do
    [[ $key =~ ^# ]] && continue   # skip comments
    [[ -z $key ]] && continue      # skip empty lines
    eval "$key=\"$value\""
    echo "  $key = '$value'"
done < "$INPUT_FILE"

# STEP 4: Clean the input bam and save in list

BAM_DIR="clean_bams"
mkdir -p "$BAM_DIR"

BAM_LIST="bam_list.txt"
> "$BAM_LIST" # clear existing file

IFS=',' read -ra ADDR <<< "$INPUT_BAM"
echo
echo "Removing marked duplicates and setting quality thresholds..."
echo

# Loop through the resulting array
for input_bam in "${ADDR[@]}"; do
    filename=$(basename "$input_bam")
    output_bam="$BAM_DIR/${filename%.bam}_filtered.bam"

    samtools view -@ 48 -b -h -F 0x400 -q 30 -Q 30 -o "$output_bam" "$input_bam"
    echo "./$output_bam" >> "$BAM_LIST"
    samtools index "$output_bam"
done

echo "Done. All filtered BAMs are in: $BAM_DIR"

# ======== PART 1: CREATING A PSEUDOHAPLOID ========

if [[ "$PSEUDOHAP" == "true" ]]; then

    # STEP 1: Choose BED/SNP files for mpileup/PileupCaller

    if [[ "$PANEL" == "wgs" ]]; then
	echo "Panel used: none. Running for WGS..."
        BED_FILE="${REF_DIR}/hgdp_references/4M_hgdp_all_chr_withX.bed"
        SNP_FILE="${REF_DIR}/hgdp_references/eigenstrat_files/hgdp_all_chr_withX.snp"

    elif [[ "$PANEL" == "twist" ]]; then
	echo "Panel used: TWSIT 1240k."
        BED_FILE="${REF_DIR}/twist_references/twist_1240k_hg38.bed"
        SNP_FILE="${REF_DIR}/twist_references/twist_1240k_hg38.snp"

    else
        echo "Error: --panel 'twist' or 'wgs' required."
        exit 1
    fi

    # STEP 2: Run Pileupcaller

    FASTA_REF="${REF_DIR}/human_references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    PILEUP_CALLER_DIR="./pileupCaller_output"
    mkdir -p "$PILEUP_CALLER_DIR"

    MPILEUP_CMD="samtools mpileup -B -R -q 25 -Q 20 -f $FASTA_REF -b $BAM_LIST -l $BED_FILE"

    echo
    echo "Running PileupCaller for EigenStrat output..."
    echo "---------------------------------------------"
    $MPILEUP_CMD | grep -v "^chrM" | $PILEUPCALLER_BIN --randomHaploid --sampleNames "$SAMPLE_NAME" \
        --samplePopName pseudoHaploid --snpFile "$SNP_FILE" --eigenstratOut "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}"

    # Rename chromosomes and positions in SNP file
    awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $2); $1=$2"_"$4; print}' \
        "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.snp" > tmp && mv tmp "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.snp"

    echo
    echo "Running PileupCaller for PLINK output..."
    echo "----------------------------------------"
    $MPILEUP_CMD | grep -v "^chrM" | $PILEUPCALLER_BIN --randomHaploid --sampleNames "$SAMPLE_NAME" \
        --samplePopName pseudoHaploid --snpFile "$SNP_FILE" --plinkOut "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}"

    # Rename chromosomes and positions in BIM file
    awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $1); $2=$1"_"$4; print}' \
        "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.bim" > tmp && mv tmp "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.bim"

    # create an autosomal version for plotting
    $PLINK_BIN --bfile "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}" --not-chr X --make-bed --out "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}_noX"

    echo
    echo "Running PileupCaller for VCF output..."
    echo "---------------------------------------------"
    $MPILEUP_CMD | grep -v "^chrM" | $PILEUPCALLER_BIN --randomHaploid --sampleNames "$SAMPLE_NAME" --samplePopName pseudoHaploid \
    --snpFile "$SNP_FILE" --vcf | bgzip > "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.vcf.gz"
    tabix -p vcf "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.vcf.gz"

    echo
    echo "✅ PseudoHaploid files created for $SAMPLE_PREFIX"
fi

# ======== PART 2: IMPUTAION using GLIMPSE 2 ========

if [[ "$IMPUTE" == "true" ]]; then

    echo
    echo "Running GLIMPSE2 imputation..."
    echo
    # STEP 1: Impute all chunks for each BAM input

    echo "STEP 1/4: Running imputation for the input BAMs..."

    # Number of concurrent chromosome and chunk jobs
    MAX_CHR_PARALLEL=8
    MAX_CHUNK_PARALLEL=15

    # set paths
    IMPUTE_DIR="./imputation_output"
    SPLIT_DIR="${REF_DIR}/glimpse2_references/ref_panel_split/"
    IMPUTED_CHUNKS_DIR="${IMPUTE_DIR}/imputed_chunks"
    CHUNKS_DIR="${REF_DIR}/glimpse2_references/chunks"
    mkdir -p "$IMPUTE_DIR" "$IMPUTED_CHUNKS_DIR"

    # Add sample name column in the BAM_LIST file to avoid (name too long error in smartPCA)
    IFS=',' read -r -a SAMPLE_ARRAY <<< "$SAMPLE_NAME"

    NAMES_FILE="${IMPUTE_DIR}/sample_names.txt"
    rm -f NAMES_FILE

    i=0
    while read -r bam; do
        sample="${SAMPLE_ARRAY[i]}"
        echo -e "$bam\t$sample" >> "$NAMES_FILE"
        ((i++)) || true
    done < "$BAM_LIST"

    echo "Using sample names:"
    cat "$NAMES_FILE"

    export GLIMPSE2_BIN NAMES_FILE SPLIT_DIR IMPUTED_CHUNKS_DIR CHUNKS_DIR MAX_CHUNK_PARALLEL

    # function to impute each autosomal chunk
    process_chunk_autosome() {
      LINE="$1"
      CHR=$(echo "$LINE" | awk '{print $2}')
      IRG=$(echo "$LINE" | awk '{print $3}')
      ORG=$(echo "$LINE" | awk '{print $4}')
      REGS=$(echo "$IRG" | cut -d: -f2 | cut -d- -f1)
      REGE=$(echo "$IRG" | cut -d: -f2 | cut -d- -f2)

      REF_BIN="${SPLIT_DIR}/${CHR}_${CHR}_${REGS}_${REGE}.bin"
      OUT_BCF="${IMPUTED_CHUNKS_DIR}/imputed_${CHR}_${REGS}_${REGE}.bcf"

      $GLIMPSE2_BIN/GLIMPSE2_phase_static \
        --bam-list "$NAMES_FILE" \
        --reference "$REF_BIN" \
        --output "$OUT_BCF"
    }

    export -f process_chunk_autosome

    # funtion to process chromosome
    process_chr() {
      CHR="$1"
      CHUNKS_FILE="${CHUNKS_DIR}/chunks.chr${CHR}.txt"
      echo "    Running imputation for chromosome $CHR"
      cat "$CHUNKS_FILE" | parallel -j "$MAX_CHUNK_PARALLEL" process_chunk_autosome
    }

    export -f process_chr

    # run autosomal chromosomes in parallel
    parallel -j "$MAX_CHR_PARALLEL" process_chr ::: $(seq 1 22)

    echo "✅ Per-chunk imputation completed."
    echo

    # STEP 2: Ligate imputed chunks for each chromosome

    echo "STEP 2/4: Ligating imputed chunks per autosome..."

    # number of parallel jobs
    CHR_LIG_PARALLEL=22

    IMPUTED_CHR_DIR="${IMPUTE_DIR}/imputed_chromosomes"
    LISTS_DIR="${IMPUTED_CHR_DIR}/list_files"
    mkdir -p "$IMPUTED_CHR_DIR" "$LISTS_DIR"

    export IMPUTED_CHR_DIR LISTS_DIR

    # function to ligate imputed chromosomes
    ligate_chr() {

      local CHR="$1"
      echo "    Ligating chromosome $CHR..."
      # Create the list file
      LIST_FILE="${LISTS_DIR}/list.chr${CHR}.txt"

      # Populate array of BCF chunk files
      mapfile -t IMPUTED_FILES < <(ls -1v "${IMPUTED_CHUNKS_DIR}/imputed_chr${CHR}_"*.bcf 2>/dev/null || true)
      printf "%s\n" "${IMPUTED_FILES[@]}" > "$LIST_FILE"

      # Ligate
      $GLIMPSE2_BIN/GLIMPSE2_ligate_static \
        --input "$LIST_FILE" \
        --output "${IMPUTED_CHR_DIR}/imputed_chr${CHR}_ligated.bcf"

      bcftools index -f "${IMPUTED_CHR_DIR}/imputed_chr${CHR}_ligated.bcf"
      echo "    ✅ Done with chr${CHR}"
    }

    export -f ligate_chr
    # run ligation function
    parallel -j "$CHR_LIG_PARALLEL" ligate_chr ::: $(seq 1 22)

    echo "✅ Per-chromosome ligatation completed."


    # STEP 3: Impute the X chromosome
    echo
    echo " Starting imputation for the X chromosome..."
    source "${REF_DIR}/../Pipelines/impute_chrX.sh"


    # STEP 4: Concatenation
    echo
    echo "STEP 4/4: Concatenating imputed chromosomes..."
    echo
    IMPUTED_SAMPLE_DIR="${IMPUTE_DIR}/imputed_samples"
    mkdir -p "$IMPUTED_SAMPLE_DIR"

    # concatenation
    OUT_CONCAT_BCF="${IMPUTED_SAMPLE_DIR}/imputed_${SAMPLE_PREFIX}.bcf"
    CHR_LIG_BCF=$(ls -1v "${IMPUTED_CHR_DIR}/imputed_chr"*"_ligated.bcf")

    bcftools concat -Ob -o "$OUT_CONCAT_BCF" $CHR_LIG_BCF
    bcftools index -f "$OUT_CONCAT_BCF"

    echo "✅ Concatenated all chromosomes for ${SAMPLE_PREFIX}"
    echo "END OF IMPUTATION."

    # STEP 4: Filtering the imputed sample
    echo
    echo "Masking the concatenated BCF..."
    echo
    MASK_BED="${REF_DIR}/hgdp_references/hgdp_wgs.20190516.mask.bed"
    MASKED_BCF="${IMPUTED_SAMPLE_DIR}/imputed_${SAMPLE_PREFIX}_masked.bcf"

    bcftools view --regions-file "$MASK_BED" --threads 24 -Ob -o "$MASKED_BCF" "$OUT_CONCAT_BCF"
    bcftools index "$MASKED_BCF"

    echo
    echo "Collapsing and filtering multiallelic sites..."
    echo
    BIALLELIC_BCF="${IMPUTED_SAMPLE_DIR}/imputed_${SAMPLE_PREFIX}_biallelic.bcf"

    bcftools norm -m+ "$MASKED_BCF" -Ou | bcftools view -v snps -m2 -M2 -Ob -o "$BIALLELIC_BCF"
    bcftools index "$BIALLELIC_BCF"

    echo
    echo "Setting GP<0.99 as '.' and filtering for RAF [0.1-0.9]..."
    echo
    FILTERED_BCF="${IMPUTED_SAMPLE_DIR}/imputed_${SAMPLE_PREFIX}_filtered.bcf"

    bcftools +setGT "$BIALLELIC_BCF" -- -t q -n . -i 'MAX(FMT/GP)<0.99' | bcftools view -i 'RAF>0.1 & RAF<0.9' -Ob -o "$FILTERED_BCF"
    bcftools index "$FILTERED_BCF"

    echo "✅ Done filtering imputed BCFs."

    # STEP 5: Creating plink files from the filtered BCF
    echo
    echo "Creating plink files from '$FILTERED_BCF'"
    echo
    IMPUTED_PLINK_DIR="${IMPUTE_DIR}/plink_files"
    mkdir -p "$IMPUTED_PLINK_DIR"

    $PLINK_BIN --bcf "$FILTERED_BCF" --double-id --make-bed --out "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered"

    # Rename chromosomes and positions in BIM file
    awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $1); $2=$1"_"$4; print}' \
        "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered.bim" > tmp && mv tmp "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered.bim"
    # create an autosomal version for plotting
    $PLINK_BIN --bfile "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered" --not-chr X --make-bed --out "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered_noX"


    echo echo "✅ Created plink files."

fi # end of [if impute==true]


# ======== PART 3: Run the PCA/Admixture plotting script ========

if [[ "$PLOT" == "pca" || "$PLOT" == "admixture" || "$PLOT" == "both" ]]; then
    source "${REF_DIR}/../Pipelines/pca_admixture.sh"
fi


