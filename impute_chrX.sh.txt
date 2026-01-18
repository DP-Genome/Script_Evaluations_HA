#!/usr/bin/env bash

set -euo pipefail

# ======= Imputation of the X chromosome =======

# STEP 1: Computing Genotype Likelihoods (GL)

# Create the ploidy file required for bcftools call

IFS=',' read -r -a NAME_ARRAY <<< "$SAMPLE_NAME"
IFS=',' read -r -a SEX_ARRAY <<< "$SAMPLE_SEX"

GL_PLOIDY="${IMPUTE_DIR}/GL_ploidy.txt"
rm -f "$GL_PLOIDY" # Clear file if it exists

for i in "${!NAME_ARRAY[@]}"; do
    NAME="${NAME_ARRAY[i]}"
    SEX="${SEX_ARRAY[i]}"

    if [[ "$SEX" == "M" ]]; then
        PLOIDY=1
    elif [[ "$SEX" == "F" ]]; then
        PLOIDY=2
    else
        PLOIDY="unknown"
    fi

    echo -e "chrX\t2781514\t155700519\t${SEX}\t${PLOIDY}" >> "$GL_PLOIDY"
done

echo "Created $GL_PLOIDY:"
cat "$GL_PLOIDY"

# Calculate the genotype likelihoods
echo "Computing genotype likelihoods..."

SITES_DIR="${REF_DIR}/glimpse1_references/sites_files"
OUT_VCF_GL="${IMPUTE_DIR}/${SAMPLE_PREFIX}_GL_chrX.vcf.gz"

bcftools mpileup -f "$FASTA_REF" -b "$BAM_LIST" --ignore-RG -I -E -a FORMAT/DP -T "${SITES_DIR}/hgdp1kgp_chrX.sites.vcf.gz" -r chrX -Ou | \
    bcftools call --ploidy GRCh38 --ploidy-file "$GL_PLOIDY" -Aim -C alleles -T "${SITES_DIR}/hgdp1kgp_chrX.sites.tsv.gz" -Oz -o "$OUT_VCF_GL"

# rename the samples to match the imputed autosomes
bcftools reheader -s "$NAMES_FILE" "$OUT_VCF_GL" -o "$OUT_VCF_GL.tmp" && mv "$OUT_VCF_GL.tmp" "$OUT_VCF_GL"
bcftools index -f "$OUT_VCF_GL"


# STEP 2: Impute the X chromosomse chunks

# First, create the sample ploidy information file
SAMPLE_PLOIDY="${IMPUTE_DIR}/sample_ploidy.txt"
> "$SAMPLE_PLOIDY"

IFS=',' read -ra SAMPLE_NAMES <<< "$SAMPLE_NAME"
IFS=',' read -ra SAMPLE_SEXES <<< "$SAMPLE_SEX"

if [[ ${#SAMPLE_NAMES[@]} -ne ${#SAMPLE_SEXES[@]} ]]; then
    echo "ERROR: SAMPLE_NAME and SAMPLE_SEX length mismatch" >&2
    exit 1
fi

for i in "${!SAMPLE_NAMES[@]}"; do
    case "${SAMPLE_SEXES[$i]}" in
        M|m) PLOIDY=1 ;;
        F|f) PLOIDY=2 ;;
        *)
            echo "ERROR: Invalid sex '${SAMPLE_SEXES[$i]}' for sample ${SAMPLE_NAMES[$i]}" >&2
            exit 1
            ;;
    esac
    printf "%s\t%d\n" "${SAMPLE_NAMES[$i]}" "$PLOIDY" >> "$SAMPLE_PLOIDY"
done
echo
echo "Sample ploidy file created:"
cat "$SAMPLE_PLOIDY"
echo

# Define chunks imputation function
echo "Imputing X chromosome chunks..."

export GLIMPSE1_BIN REF_DIR OUT_VCF_GL SAMPLE_PLOIDY

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

  $GLIMPSE1_BIN/GLIMPSE_phase_static \
      --input "$OUT_VCF_GL" \
      --reference "$REF_X_BCF" \
      --map "$MAP_X" \
      --input-region "${IRG}" \
      --output-region "${ORG}" \
      --samples-file "$SAMPLE_PLOIDY" \
      --output "$OUT_BCF"

  bcftools index -f "$OUT_BCF"
}

export -f process_chunk_chrX

# Run the imputation function
cat "${CHUNKS_DIR}/chunks.chrX.txt" | parallel -j 25 process_chunk_chrX

# STEP 3: Ligate the X chromosomse chunks
LST="${LISTS_DIR}/list.chrX.txt"
ls "${IMPUTED_CHUNKS_DIR}"/imputed_chrX_*.bcf > "${LST}"

$GLIMPSE1_BIN/GLIMPSE_ligate_static --input "${LST}" --output "${IMPUTED_CHR_DIR}/imputed_chrX_ligated.bcf"
