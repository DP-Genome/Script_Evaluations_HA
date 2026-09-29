#!/bin/bash
set -euo pipefail

# ======== Logger setup ========
LOG_FILE="gedmatch_export_$(date +%Y%m%d_%H%M%S).log"
log_msg() {
    local msg="$1"
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo "[$timestamp] $msg" | tee -a "$LOG_FILE"
}

# ======== Software discovery ========
log_msg "Discovering software tools..."
SOFTWARES_DIR="../References/softwares"
find_tool() {
    local name="$1"
    local path=$(find "$SOFTWARES_DIR" -name "$name" -type f -executable | head -n 1)
    if [[ -z "$path" ]]; then
        if command -v "$name" >/dev/null 2>&1; then echo "$name"; else log_msg "ERROR: '$name' not found."; exit 1; fi
    else
        echo "$path"
    fi
}

# ======== Set fixed paths ========
REF_DIR="../References"

HG38_REF="${REF_DIR}/human_references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
HG19_REF="${REF_DIR}/human_references/other_files/hs37d5_noChr.fa"
CHAIN_FILE="${REF_DIR}/human_references/other_files/hg38ToHg19.over.chain.gz"

# RSID Annotation sources (hg38 based)
RSID_BIM="./23andMe_variants.txt"
# Subset references
FIGG_SNPS_HG37="${REF_DIR}/GEDmatch_references/snps_23andME_v5.txt"

# ======= Process Input =======
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input_bcf>"
    exit 1
fi
INPUT_HG38_BCF="$1"
OUTPUT_DIR="./GEDmatch_input"
mkdir -p "$OUTPUT_DIR"

# ======= STEP 0: Prepare RSID Annotation =======
# Convert BIM to a format bcftools annotate likes (CHROM POS RSID)
ANNOT_FILE="${OUTPUT_DIR}/rsid_annotation.txt.gz"
if [[ ! -f "$ANNOT_FILE" || ! -f "$ANNOT_FILE.tbi" ]]; then
    log_msg "Preparing RSID annotation file..."
    # Ensure CHROM names match input (chr1, chr2...)
    awk 'BEGIN{OFS="\t"}{print $1, $4, $2}' "$RSID_BIM" | sort -k1,1 -k2,2n | bgzip -c > "$ANNOT_FILE"
    tabix -s1 -b2 -e2 "$ANNOT_FILE"
fi

# ======= STEP 2: LiftOver hg38 → hg19 =======
log_msg "STEP 2: Running liftover..."
BCF_NAME=$(basename "${INPUT_HG38_BCF%.bcf}")

liftover() {
    local in_bcf="$1"
    local out_bcf="$2"

    bcftools annotate -x FORMAT/DS -Ou "$in_bcf" | \
    bcftools +liftover -Ou -- -s "$HG38_REF" -f "$HG19_REF" -c "$CHAIN_FILE" | \
    bcftools sort -Ob -o "$out_bcf" &>> "$LOG_FILE"

    bcftools index --threads 24 "$out_bcf" &>> "$LOG_FILE"
}

BCF_WG_HG19="${OUTPUT_DIR}/${BCF_NAME}_WG_hg19.bcf"

log_msg "  Lifting over Whole Genome..."
liftover "$INPUT_HG38_BCF" "$BCF_WG_HG19"

# ======= STEP 1: Annotate RSIDs and Subset =======
log_msg "STEP 1: Annotating and subsetting BCF..."
BCF_ANNOT="${OUTPUT_DIR}/${BCF_NAME}_annotated.bcf"

# 1.1 Annotate with comprehensive RSIDs
log_msg "  Annotating with HGDP RSIDs..."
bcftools annotate -a "$ANNOT_FILE" -c CHROM,POS,ID -Ob -o "$BCF_ANNOT" "$BCF_WG_HG19" &>> "$LOG_FILE"


# ======= STEP 3: Subset FIGG from WG hg19 =======
log_msg "STEP 3: Subsetting for FIGG SNPs (hg19)..."
BCF_FIGG_HG19="${OUTPUT_DIR}/${BCF_NAME}_FIGG_hg19.bcf"
# Create FIGG positions file (CHROM POS)
FIGG_POS="${OUTPUT_DIR}/FIGG_positions.txt"
awk '{print $2"\t"$3}' "$FIGG_SNPS_HG37" > "$FIGG_POS"
bcftools view -T "$FIGG_POS" -Ob -o "$BCF_FIGG_HG19" "$BCF_ANNOT" &>> "$LOG_FILE"
bcftools index --threads 24 "$BCF_FIGG_HG19" &>> "$LOG_FILE"



# ======= STEP 4: Create GEDmatch output files =======
log_msg "STEP 4: Creating GEDmatch input files..."

create_gedmatch_txt() {
    local in_bcf="$1"
    local sample="$2"
    local type_label="$3"
    local out_file="${OUTPUT_DIR}/GEDmatch_input_${sample}_${type_label}.txt"

    log_msg "  Generating $type_label file for sample $sample..."

    # Query data and convert format
    bcftools query -s "$sample" -f '%ID\t%CHROM\t%POS\t[%TGT]\n' "$in_bcf" \
    | awk 'BEGIN{OFS="\t"}{
        g=$4
        gsub(/[|\/]/, "", g)
        gsub(/\./, "N", g)
        print $1, $2, $3, g
    }' > "$out_file.tmp"

    # Add header and filter missing alleles
    (echo -e "rsid\tchromosome\tposition\tgenotype" && \
     awk 'BEGIN{OFS="\t"}{
         if ($4 ~ /N/)
             $4="--"
         print
     }' "$out_file.tmp") > "$out_file"

    rm "$out_file.tmp"
}

SAMPLES=$(bcftools query -l "$BCF_ANNOT")
for S in $SAMPLES; do
    create_gedmatch_txt "$BCF_ANNOT" "$S" "WG"
    create_gedmatch_txt "$BCF_FIGG_HG19" "$S" "FIGG"
done

log_msg "✅ GEDmatch export completed for all samples and types."
