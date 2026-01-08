#|!/bin/bash

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Check that exactly 2 arguments were provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <BAM1> <BAM2>"
    exit 1
fi

# Read inputs
BAM1="$1"
BAM2="$2"

# To mask
MASK_FILE="../../References/hgdp_references/hgdp_wgs.20190516.mask.bed"
BAM_MASKED1="${BAM1%.bam}.masked.bam"
BAM_MASKED2="${BAM2%.bam}.masked.bam"

# For mpileup
BED_FILE="/mnt/server_share/Hanan/human_pipeline/ancient_DNA/contamination_estimation/hapCon_hg38_trial/bed_liftover/output_1kghgdp_hg38_chrX.bed"

# Chromosome X Length (hg38 standard).
# Adjust this if you are using a specific reference build where lengths differ.
CHROM_X_LEN=156040895

# Ratios to process (S1 fraction). The S2 fraction is calculated as (1 - S1).
RATIOS=(0.95 0.90 0.85 0.80)

# Coverages to generate
COVERAGES=(0.1 0.25 0.5 1.0)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Function to get total mapped reads in ChrX
get_read_count() {
    local bam=$1
    # Uses idxstats to get count specifically for 'chrX' or 'X'
    # Sums up counts for X to handle potential chr prefix differences
    samtools idxstats "$bam" | grep -E "^(chr)?X\b" | cut -f3 | awk '{s+=$1} END {print s}'
}

# Function to estimate read length (checks first 100 reads)
get_read_length() {
    local bam=$1
    samtools view "$bam" | head -n 100 | awk '{sum+=length($10)} END {print int(sum/NR)}'
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

echo "Starting Analysis..."
echo "Input 1: $BAM1"
echo "Input 2: $BAM2"

# 1. Mask the Input Files
echo "Masking ${BAM1}..."
samtools view -@ 96 -L "$BED_FILE" -b "$BAM1" > "$BAM_MASKED1"
samtools index -@ 24 "$BAM_MASKED1"

echo "Masking ${BAM2}..."
samtools view -@ 96 -L "$BED_FILE" -b "$BAM2" > "$BAM_MASKED2"
samtools index -@ 24 "$BAM_MASKED2"

# 2. Gather Input Statistics
echo "Gathering input statistics..."
TOTAL_READS_1=$(get_read_count "$BAM_MASKED1")
TOTAL_READS_2=$(get_read_count "$BAM_MASKED2")
READ_LEN=$(get_read_length "$BAM_MASKED1")

echo "Reads in BAM1 (ChrX): $TOTAL_READS_1"
echo "Reads in BAM2 (ChrX): $TOTAL_READS_2"
echo "Estimated Read Length: $READ_LEN bp"
echo "Target Genome Size (ChrX): $CHROM_X_LEN bp"
echo "------------------------------------------------"

# Loop through Coverages
for COV in "${COVERAGES[@]}"; do

    # Calculate Total Reads needed for this specific coverage
    # Formula: (Coverage * GenomeLength) / ReadLength
    REQUIRED_TOTAL_READS=$(echo "$COV * $CHROM_X_LEN / $READ_LEN" | bc)
    echo "Processing Target Coverage: ${COV}x (Requires approx $REQUIRED_TOTAL_READS reads)"

    # Loop through Ratios
    for R1 in "${RATIOS[@]}"; do

        # Calculate R2 (Sample 2 ratio) using bc for float subtraction
        R2=$(echo "1 - $R1" | bc)

        # Construct Output Filename
        OUT_BASE="merged_r${R1}_cov${COV}"
        OUT_BAM="${OUT_BASE}.bam"

        echo "  > Ratio: S1=$R1 / S2=$R2"

        # Calculate exact reads needed from each sample
        NEEDED_S1=$(echo "$REQUIRED_TOTAL_READS * $R1" | bc | awk '{print int($1)}')
        NEEDED_S2=$(echo "$REQUIRED_TOTAL_READS * $R2" | bc | awk '{print int($1)}')

        # Calculate Sampling Fraction for samtools view -s
        # We assume a fixed seed (e.g., 42) + the fraction
        # If needed reads > available reads, set fraction to 1 (take all)

        if (( $(echo "$NEEDED_S1 > $TOTAL_READS_1" | bc -l) )); then
            echo "    WARNING: Not enough reads in BAM_MASKED1 for requested coverage. Taking all."
            FRAC_1=1
        else
            FRAC_1=$(echo "$NEEDED_S1 / $TOTAL_READS_1" | bc -l)
        fi

        if (( $(echo "$NEEDED_S2 > $TOTAL_READS_2" | bc -l) )); then
            echo "    WARNING: Not enough reads in BAM_MASKED2 for requested coverage. Taking all."
            FRAC_2=1
        else
            FRAC_2=$(echo "$NEEDED_S2 / $TOTAL_READS_2" | bc -l)
        fi

        # Temporary files for subsampled chunks
        TMP_S1="tmp_s1_${OUT_BASE}.bam"
        TMP_S2="tmp_s2_${OUT_BASE}.bam"

        # Subsample Sample 1
        # Seed 42 is added to the integer part
        samtools view -b -s 42$FRAC_1 "$BAM_MASKED1" > "$TMP_S1" &
        PID1=$!

        # Subsample Sample 2
        # Use different seed (e.g. 43) just to be safe, though inputs differ
        samtools view -b -s 43$FRAC_2 "$BAM_MASKED2" > "$TMP_S2" &
        PID2=$!

        wait $PID1 $PID2

        # Merge the subsampled files
        # -f forces overwrite, -c attaches RG tags if present (optional but good practice)
        samtools merge -f "$OUT_BAM" "$TMP_S1" "$TMP_S2"
        samtools index "$OUT_BAM"

        rm "$TMP_S1" "$TMP_S2"

        # ----------------------------------------------------------------------
        # Run mpileup
        # ----------------------------------------------------------------------
        echo "    Running mpileup for $OUT_BAM..."
        samtools mpileup --positions "$BED_FILE" -r chrX -q 20 -Q 30 -o "${OUT_BASE}.mpileup" "$OUT_BAM"

    done
done

echo "------------------------------------------------"
echo "All 16 combinations processed successfully."
