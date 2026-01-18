#!/bin/bash
set -euo pipefail

# ======== Set fixed paths ========

REF_DIR="/mnt/server_share/Hanan/human_pipeline/References"
HG38_REF="${REF_DIR}/human_references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
HG19_REF="${REF_DIR}/human_references/other_files/hs37d5_noChr.fa"
CHAIN_FILE="${REF_DIR}/human_references/other_files/hg38ToHg19.over.chain.gz"

REF_RSID="${REF_DIR}/twist_references/twist_1240k_annotation_file.txt.gz"

# ======= Process Input =======

# Check that exactly 1 argument was provided
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input_bcf>"
    exit 1
fi

# Read the input BCF
INPUT_HG38_BCF="$1"

# ======= STEP 1: liftOver hg38 →  hg19 =======
echo
echo "Running liftover hg38 →  hg19..."
echo

OUTPUT_DIR="./GEDmatch_input"
BCF_RSID="${OUTPUT_DIR}/$(basename "${INPUT_HG38_BCF%.bcf}")_withRSID.bcf"
OUTPUT_HG19_BCF="${OUTPUT_DIR}/$(basename "${INPUT_HG38_BCF%.bcf}")_hg19.bcf"
mkdir -p "$OUTPUT_DIR"

# remove DS values to avoid liftover error
# annotate the bcf with RSIDs - needed for STEP 2
bcftools annotate -x FORMAT/DS -a "$REF_RSID" -h <(echo '##annotate=file,ID=ID,Number=1,Type=String,Description="RSID annotation"') \
        -c CHROM,POS,ID -Ob -o "$BCF_RSID" "$INPUT_HG38_BCF"
bcftools index --threads 24 "$BCF_RSID"

# liftover
bcftools +liftover "$BCF_RSID" -- -s "$HG38_REF" -f "$HG19_REF" \
        -c "$CHAIN_FILE" | bcftools sort -Ob -o "$OUTPUT_HG19_BCF"

# index and convert to vcf
OUTPUT_HG19_VCF="${OUTPUT_HG19_BCF%.bcf}.vcf"

bcftools index --threads 24 "$OUTPUT_HG19_BCF"
bcftools view "$OUTPUT_HG19_BCF" -Ov -o "$OUTPUT_HG19_VCF"
bgzip -f "$OUTPUT_HG19_VCF"
bcftools index --threads 24 "${OUTPUT_HG19_VCF}.gz"

echo
echo "✅ LiftOver complete. Output files:"
echo " BCF: $OUTPUT_HG19_BCF"
echo " VCF: ${OUTPUT_HG19_VCF}.gz"
echo
# ======= STEP 2: create GEDmatch input file =======

# create one GEDmatch input file per sample

bcftools query -l "${OUTPUT_HG19_VCF}.gz" > "${OUTPUT_DIR}/samples.txt"

for S in $(cat "${OUTPUT_DIR}/samples.txt"); do
    echo "Generating GEDmatch input file for sample $S ..."

    bcftools query -s "$S" -f '%ID\t%CHROM\t%POS\t[%TGT]\n' "${OUTPUT_HG19_VCF}.gz" \
    | awk 'BEGIN{OFS="\t"}{
        g=$4
        g=gensub(/\|/,"","g",g)
        g=gensub(/\//,"","g",g)
        g=gensub(/\./,"N","g",g)
        print $1,$2,$3,g
    }' > "${OUTPUT_DIR}/GEDmatch_input_${S}.txt"

    # add header
    sed -i '1irsid\tchromosome\tposition\tgenotype' "${OUTPUT_DIR}/GEDmatch_input_${S}.txt"
    # remove positions with missing genotypes and positions
    awk 'NR==1 || $4 != "NN"' "${OUTPUT_DIR}/GEDmatch_input_${S}.txt" | \
    awk 'NR==1 || $1 != "\."' \
    > tmp && mv tmp "${OUTPUT_DIR}/GEDmatch_input_${S}.txt"


done

echo
echo "✅ GEDmatch input file created for each sample."
echo
