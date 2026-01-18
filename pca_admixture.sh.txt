#!/usr/bin/env bash

set -euo pipefail

# ======= PART 1: Create the merged plink files =======

# STEP 1: set file paths

# Reference plink prefix
REF_PLINK_PREFIX="${REF_DIR}/hgdp_references/plink_files/hgdp_all_chr"
# output dir
mkdir -p "merged_plink"
# PLINK sample/reference merged output
PLINK_MERGED_PREFIX="./merged_plink/merged_ref"


# STEP 2: Create list of plink files to be merged with reference
echo
echo "Creating merge_list.txt..."
echo
MERGE_LIST="merge_list.txt"
> "$MERGE_LIST"

if [[ "$PSEUDOHAP" == "true" ]]; then
    echo "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}_noX" >> "$MERGE_LIST"
fi

if [[ "$IMPUTE" == "true" ]]; then
    echo "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered_noX" >> "$MERGE_LIST"
fi

# check list
if [[ ! -s "$MERGE_LIST" ]]; then
    echo "ERROR: merge list is empty. Nothing to merge." >&2
    exit 1
fi
echo
echo "Merging the following datasets:"
cat "$MERGE_LIST"

# STEP 3: Merge plinks
echo
echo "Merging step 1/2: Merging sample with reference plink files..."
echo "--------------------------------------------------------------"
$PLINK_BIN --bfile "$REF_PLINK_PREFIX" --merge-list "$MERGE_LIST" --make-bed --out "${PLINK_MERGED_PREFIX}"
echo
echo "Merging step 2/2: Performing LD pruning"
echo "-------------------------------"
echo "Pruning..."
$PLINK_BIN --bfile "${PLINK_MERGED_PREFIX}" --indep-pairwise 50 5 0.2 --out "${PLINK_MERGED_PREFIX}_pruned"
$PLINK_BIN --bfile "${PLINK_MERGED_PREFIX}" --extract "${PLINK_MERGED_PREFIX}_pruned.prune.in" --make-bed --out "${PLINK_MERGED_PREFIX}_pruned"

# ======== PART 2: Plot PCA ========

if [[ "$PLOT" == "pca" || "$PLOT" == "both" ]]; then
    echo
    echo "Runnin smart PCA..."
    echo

    # Population list to project on
    POPULATION_LIST="${REF_DIR}/hgdp_references/hgdp_population/populationList.txt"
    SUBPOPULATION_LIST="${REF_DIR}/hgdp_references/hgdp_population/subpopulationList.txt"

    # Define output paths
    mkdir -p "smart_pca"
    PCA_MERGED_PREFIX="./smart_pca/merged_ref"

    # Define parameter files to be generated
    PAR_PCA_CONVERTF="par_pca.convertf"
    PAR_SMARTPCA="par.smartpcaProject"

    echo "=== Generating convertf parameter file: $PAR_PCA_CONVERTF ==="
    cat <<EOF > $PAR_PCA_CONVERTF
    genotypename:       ${PLINK_MERGED_PREFIX}_pruned.bed
    snpname:            ${PLINK_MERGED_PREFIX}_pruned.bim
    indivname:          ${PLINK_MERGED_PREFIX}_pruned.fam
    outputformat:       EIGENSTRAT
    genotypeoutname:    ${PCA_MERGED_PREFIX}.geno
    snpoutname:         ${PCA_MERGED_PREFIX}.snp
    indivoutname:       ${PCA_MERGED_PREFIX}.ind
EOF

    echo "=== Generating smartpca parameter file: $PAR_SMARTPCA ==="
    cat <<EOF > $PAR_SMARTPCA
    genotypename: ${PCA_MERGED_PREFIX}.geno
    snpname:      ${PCA_MERGED_PREFIX}.snp
    indivname:    ${PCA_MERGED_PREFIX}_withPop.ind
    poplistname:  ${POPULATION_LIST}
    evecoutname:  ${PCA_MERGED_PREFIX}.evec
    evaloutname:  ${PCA_MERGED_PREFIX}.eval
    lsqproject: YES
    shrink mode: YES
    numoutlieriter: 0
    numoutevec: 10
    numthreads: 24
EOF

    #  Run smart PCA
    echo
    echo "PCA step 1/2: Converting to EIGENSTRAT geno/snp/ind format..."
    echo "-------------------------------------------------------------"
    $CONVERTF_BIN -p "$PAR_PCA_CONVERTF"
    echo
    echo "Assigning population labels..."
    awk 'NR==FNR { pop[$1]=$2; next }          # first file: build lookup table
        $3=="???" && ($1 in pop) { $3=pop[$1] }       # replace ??? with population if known
        $3=="???" && !($1 in pop) { $3=$1 }  # replace ??? with sample name if unknown
        { print }
    ' "${REF_DIR}/hgdp_references/population_labels.tsv" "${PCA_MERGED_PREFIX}.ind" > "${PCA_MERGED_PREFIX}_withPop.ind"

    echo
    echo "PCA step 2/2: Running smartPCA..."
    echo "---------------------------------"
    $SMARTPCA_BIN -p "$PAR_SMARTPCA"

    # add the subpopultion hgdp names to the evec file
    paste "${PCA_MERGED_PREFIX}.evec" <(awk '{print $2}' "$SUBPOPULATION_LIST") > tmp.evec
    # same for our samples
    awk '{
      if ($NF ~ /:/) {
        split($NF, a, ":");
        print $0 "\t" a[2];
      } else {
        print;
      }
    }' tmp.evec > "${PCA_MERGED_PREFIX}_labeled.evec"

    rm tmp.evec


    # print number of SNPs in each step
    echo
    echo "Counting all SNPS"
    echo "-----------------"
    echo "  --> HGDP reference input:" $(wc -l < "${REF_PLINK_PREFIX}.bim")
    echo "  --> SNPs in the Pseudohaplotype:" $(wc -l < "${PILEUP_CALLER_DIR}/${SAMPLE_PREFIX}.bim")
    echo "  --> SNPs in the imputed sample:" $(wc -l < "${IMPUTED_PLINK_DIR}/imputed_${SAMPLE_PREFIX}_filtered.bim")
    echo "  --> Merged PLINK files:" $(wc -l < "${PLINK_MERGED_PREFIX}.bim")
    echo "  --> After pruning (PCA input):" $(wc -l < "${PLINK_MERGED_PREFIX}_pruned.bim")
    echo
    echo "=== PCA FINISHED SUCCESSFULLY ==="
    echo

fi # end of if [[ PLOT == pca || PLOT == both ]]

# ======== PART 3: Plot Admixture ========

if [[ "$PLOT" == "admixture" || "$PLOT" == "both" ]]; then

    # output dir
    ADMIXTURE_DIR="./admixture"
    mkdir -p "$ADMIXTURE_DIR"
    # Admixture metadata file for annotation
    ADMIXTURE_LABELS="${REF_DIR}/hgdp_references/hgdp_population/admixture_labels.txt"
    # Define output prefix
    ADMIXTURE_PREFIX="./$ADMIXTURE_DIR/merged_ref"

    echo
    echo "Running Admixture..."
    echo "--------------------"
    echo " Converting input to PLINK..."
    for k in {5..14}; do $ADMIXTURE_BIN -j120 --cv ${PLINK_MERGED_PREFIX}_pruned.bed $k 2>&1 | tee admixture_K${k}.log; done
    grep '^CV error' admixture_K*.log > cv_comparison.txt

    # add population labels
    cp "$ADMIXTURE_LABELS" "${ADMIXTURE_DIR}/"

    # append sample names to the labels file
    IFS=',' read -ra SAMPLES <<< "$SAMPLE_NAME"
    for s in "${SAMPLES[@]}"; do
        echo -e "$s\t$s\t$s\tM" >> "${ADMIXTURE_DIR}/admixture_labels.txt"
    done

    # paste to each Q file
    for Qfile in *.Q; do paste "${ADMIXTURE_DIR}/admixture_labels.txt" "$Qfile" > "${Qfile}.withPop"; done

    mkdir -p "${ADMIXTURE_DIR}/Qfiles_withpop"
    mv *withPop "${ADMIXTURE_DIR}/Qfiles_withpop/"
    mv *.{P,Q,log} "$ADMIXTURE_DIR" && mv cv_comparison.txt "$ADMIXTURE_DIR"

    echo
    echo "=== ADMIXTURE FINISHED SUCCESSFULLY ==="


fi # end of if [[ PLOT == admixture || PLOT == both ]]


