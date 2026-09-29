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

BCFTOOLS=$(find_tool "bcftools")
BGZIP=$(find_tool "bgzip")
TABIX=$(command -v tabix || find_tool "tabix")

# ======== Set fixed paths ========
REF_DIR="../References"

HG38_REF="${REF_DIR}/human_references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
HG19_REF="${REF_DIR}/human_references/other_files/hs37d5_noChr.fa"
CHAIN_FILE="${REF_DIR}/human_references/other_files/hg38ToHg19.over.chain.gz"

# RSID Annotation sources (hg38 based)
# Primary: Comprehensive HGDP-based RSID file
RSID_BIM="${REF_DIR}/hgdp_references/plink_files/hgdp_all_chr_withRSID.bim"
# Secondary: Twist panel annotation file
TWIST_RSID_REF="${REF_DIR}/twist_references/twist_1240k_annotation_file.txt.gz"

# Subset references
TWIST_BED_HG38="${REF_DIR}/twist_references/twist_1240k_hg38_chr.bed"
FIGG_SNPS_HG37="${REF_DIR}/GEDmatch_references/FIGG_snps_hg37.txt"

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
    awk 'BEGIN{OFS="\t"}{print "chr"$1, $4, $2}' "$RSID_BIM" | sort -k1,1 -k2,2n | $BGZIP -c > "$ANNOT_FILE"
    $TABIX -s1 -b2 -e2 "$ANNOT_FILE"
fi

# ======= STEP 1: Annotate RSIDs and Subset =======
log_msg "STEP 1: Annotating and subsetting BCF..."
BCF_NAME=$(basename "${INPUT_HG38_BCF%.bcf}")
BCF_ANNOT="${OUTPUT_DIR}/${BCF_NAME}_annotated.bcf"

# 1.1 Annotate with comprehensive RSIDs
log_msg "  Annotating with HGDP RSIDs..."
$BCFTOOLS annotate -a "$ANNOT_FILE" -c CHROM,POS,ID -Ob -o "$BCF_ANNOT.tmp" "$INPUT_HG38_BCF" &>> "$LOG_FILE"

# 1.2 Fill in missing RSIDs with Twist panel RSIDs
log_msg "  Filing in missing RSIDs with Twist panel..."
$BCFTOOLS annotate -a "$TWIST_RSID_REF" -c CHROM,POS,ID --collapse none -Ob -o "$BCF_ANNOT" "$BCF_ANNOT.tmp" &>> "$LOG_FILE"
rm "$BCF_ANNOT.tmp"
$BCFTOOLS index --threads 24 "$BCF_ANNOT" &>> "$LOG_FILE"

# 1.3 Create Twist hg38 subset
BCF_TWIST_HG38="${OUTPUT_DIR}/${BCF_NAME}_twist_hg38.bcf"
log_msg "  Subsetting for Twist panel (hg38)..."
$BCFTOOLS view -T "$TWIST_BED_HG38" -Ob -o "$BCF_TWIST_HG38" "$BCF_ANNOT" &>> "$LOG_FILE"
$BCFTOOLS index --threads 24 "$BCF_TWIST_HG38" &>> "$LOG_FILE"

# ======= STEP 2: LiftOver hg38 → hg19 =======
log_msg "STEP 2: Running liftover..."

lift_over() {
    local in_bcf="$1"
    local out_bcf="$2"
    $BCFTOOLS +liftover "$in_bcf" -- -s "$HG38_REF" -f "$HG19_REF" -c "$CHAIN_FILE" \
    | $BCFTOOLS sort -Ob -o "$out_bcf" &>> "$LOG_FILE"
    $BCFTOOLS index --threads 24 "$out_bcf" &>> "$LOG_FILE"
}

BCF_WG_HG19="${OUTPUT_DIR}/${BCF_NAME}_WG_hg19.bcf"
BCF_TWIST_HG19="${OUTPUT_DIR}/${BCF_NAME}_twist_hg19.bcf"

log_msg "  Lifting over Whole Genome..."
lift_over "$BCF_ANNOT" "$BCF_WG_HG19"
log_msg "  Lifting over Twist subset..."
lift_over "$BCF_TWIST_HG38" "$BCF_TWIST_HG19"

# ======= STEP 3: Subset FIGG from WG hg19 =======
log_msg "STEP 3: Subsetting for FIGG SNPs (hg19)..."
BCF_FIGG_HG19="${OUTPUT_DIR}/${BCF_NAME}_FIGG_hg19.bcf"
# Create FIGG positions file (CHROM POS)
FIGG_POS="${OUTPUT_DIR}/FIGG_positions.txt"
awk '{print $2"\t"$3}' "$FIGG_SNPS_HG37" > "$FIGG_POS"
$BCFTOOLS view -T "$FIGG_POS" -Ob -o "$BCF_FIGG_HG19" "$BCF_WG_HG19" &>> "$LOG_FILE"
$BCFTOOLS index --threads 24 "$BCF_FIGG_HG19" &>> "$LOG_FILE"

# ======= STEP 4: Create GEDmatch output files =======
log_msg "STEP 4: Creating GEDmatch input files..."

create_gedmatch_txt() {
    local in_bcf="$1"
    local sample="$2"
    local type_label="$3"
    local out_file="${OUTPUT_DIR}/GEDmatch_input_${sample}_${type_label}.txt"

    log_msg "  Generating $type_label file for sample $sample..."
    
    # Query data and convert format
    $BCFTOOLS query -s "$sample" -f '%ID\t%CHROM\t%POS\t[%TGT]\n' "$in_bcf" \
    | awk 'BEGIN{OFS="\t"}{
        g=$4
        gsub(/[|\/]/, "", g)
        gsub(/\./, "N", g)
        print $1, $2, $3, g
    }' > "$out_file.tmp"

    # Add header and filter missing alleles
    (echo -e "rsid\tchromosome\tposition\tgenotype" && \
     awk '$4 !~ /N/' "$out_file.tmp") > "$out_file"

    rm "$out_file.tmp"
}

SAMPLES=$($BCFTOOLS query -l "$BCF_WG_HG19")
for S in $SAMPLES; do
    create_gedmatch_txt "$BCF_WG_HG19" "$S" "WG"
    create_gedmatch_txt "$BCF_TWIST_HG19" "$S" "Twist"
    create_gedmatch_txt "$BCF_FIGG_HG19" "$S" "FIGG"
done

log_msg "✅ GEDmatch export completed for all samples and types."
