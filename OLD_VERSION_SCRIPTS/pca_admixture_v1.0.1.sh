#!/usr/bin/env bash
# Version: 1.0.1
#   1.0.1 - PCA/Admixture use the pruned set that matches --prune (ldprio or indep);
#           new-sample list no longer hard-coded to line 930; admixture runs in its own
#           directory (no more moving *.log); threads come from $CORES; PRUNE/PLOT defaults.
#   1.0.0 - initial tagged release (see OLD_VERSION_SCRIPTS/)

set -euo pipefail
: "${PRUNE:=ldprio}"
: "${PLOT:=both}"
: "${CORES:=24}"

# ======== Logger and Software Discovery ========
# These variables should be exported from the calling script (pileupCaller_glimpse.sh)
# If not, we attempt to find them here for standalone execution.

if ! command -v log_msg >/dev/null 2>&1; then
    LOG_FILE="${LOG_FILE:-pca_admixture_$(date +%Y%m%d_%H%M%S).log}"
    log_msg() {
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] $1" | tee -a "$LOG_FILE"
    }
fi

SOFTWARES_DIR="../References/softwares"
find_tool() {
    local name="$1"
    local path=$(find "$SOFTWARES_DIR" -name "$name" -type f -executable | head -n 1)
    [[ -z "$path" ]] && echo "$name" || echo "$path"
}

: "${PLINK_BIN:=$(find_tool "plink")}"
: "${CONVERTF_BIN:=$(find_tool "convertf")}"
: "${SMARTPCA_BIN:=$(find_tool "smartpca")}"
: "${ADMIXTURE_BIN:=$(find_tool "admixture")}"

# ======= PART 1: Create the merged plink files =======

log_msg "PART 1: Creating merged PLINK files..."

# Reference plink prefix
REF_PLINK_PREFIX="${REF_DIR}/hgdp_references/plink_files/hgdp_all_chr"

# PLINK sample/reference merged output
mkdir -p "merged_plink"
PLINK_MERGED_PREFIX="./merged_plink/merged_ref"
PLINK_MERGED_PREFIX_LDPRIO="./merged_plink_ldprio/merged_ref"

# Prefix of the pruned fileset used by PCA and Admixture (depends on --prune)
case "$PRUNE" in
    indep)  PAR_PREFIX="${PLINK_MERGED_PREFIX}_pruned" ;;
    ldprio) PAR_PREFIX="${PLINK_MERGED_PREFIX_LDPRIO}_ldprioPruned" ;;
    *) log_msg "ERROR: PRUNE must be indep or ldprio (got '$PRUNE')"; exit 1 ;;
esac

# STEP 2: Create list of plink files to be merged with reference
log_msg "Creating merge_list.txt..."
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
    log_msg "ERROR: merge list is empty. Nothing to merge." >&2
    exit 1
fi
log_msg "Merging sample datasets with HGDP reference..."
$PLINK_BIN --bfile "$REF_PLINK_PREFIX" --merge-list "$MERGE_LIST" --make-bed --out "${PLINK_MERGED_PREFIX}" &>> "$LOG_FILE"

if [[ "$PRUNE" == "indep" ]]; then
    log_msg "Performing LD pruning (r2 > 0.2)..."
    # save pruned output in the same file as the plink files
    $PLINK_BIN --bfile "${PLINK_MERGED_PREFIX}" --maf 0.01 --geno 0.05 --indep-pairwise 50 5 0.2 --out "${PLINK_MERGED_PREFIX}_pruned" &>> "$LOG_FILE"
    $PLINK_BIN --bfile "${PLINK_MERGED_PREFIX}" --extract "${PLINK_MERGED_PREFIX}_pruned.prune.in" --make-bed --out "${PLINK_MERGED_PREFIX}_pruned" &>> "$LOG_FILE"
fi

# Samples added on top of the reference panel (= in the merged .fam but not in the reference .fam)
mkdir -p "merged_plink_ldprio"
LDPRIO_SAMPLES="./merged_plink_ldprio/ldprioSamples.txt"
awk 'NR==FNR { ref[$1 "\t" $2]; next } !(($1 "\t" $2) in ref) { print $1 "\t" $2 }' \
    "${REF_PLINK_PREFIX}.fam" "${PLINK_MERGED_PREFIX}.fam" > "$LDPRIO_SAMPLES"
if [[ ! -s "$LDPRIO_SAMPLES" ]]; then
    log_msg "ERROR: no new samples found after merging with the reference."
    exit 1
fi
log_msg "New (non-reference) samples: $(wc -l < "$LDPRIO_SAMPLES")"

if [[ "$PRUNE" == "ldprio" ]]; then

    log_msg "Performing ldprio pruning..."

    python3 "${REF_DIR}/softwares/ldprio.py" --bfile "${PLINK_MERGED_PREFIX}" --priority-samples "${LDPRIO_SAMPLES}" \
        --ld-samples "${REF_DIR}/hgdp_references/modernSamples.txt" --autosomes-only --make-bed --out "${PLINK_MERGED_PREFIX_LDPRIO}_ldprioPruned"
fi
# ======== PART 2: Plot PCA ========

if [[ "$PLOT" == "pca" || "$PLOT" == "both" ]]; then
    log_msg "PART 2: Running SmartPCA..."

    # Population list to project on
    POP_LIST_REF="${REF_DIR}/hgdp_references/hgdp_population/populationList.txt"
    SUBPOP_LIST_REF="${REF_DIR}/hgdp_references/hgdp_population/subpopulationList.txt"
    LABELS_REF="${REF_DIR}/hgdp_references/hgdp_population/population_labels.tsv"

    # Define output paths
    mkdir -p "smart_pca" "smart_pca_CI"
    PCA_MERGED_PREFIX="./smart_pca/merged_ref"
    PCA_MERGED_PREFIX_CI="./smart_pca_CI/merged_ref"

    # Define parameter files to be generated
    PAR_PCA_CONVERTF="par_pca.convertf"
    PAR_SMARTPCA="par.smartpcaProject"
    PAR_PCA_CONVERTF_CI="par_pca_CI.convertf"
    PAR_SMARTPCA_CI="par_CI.smartpcaProject"

    # for smart pca without confidence ellipses
    log_msg "Generating parameter files..."
    cat <<EOF > $PAR_PCA_CONVERTF
    genotypename:       ${PAR_PREFIX}.bed
    snpname:            ${PAR_PREFIX}.bim
    indivname:          ${PAR_PREFIX}.fam
    outputformat:       EIGENSTRAT
    genotypeoutname:    ${PCA_MERGED_PREFIX}.geno
    snpoutname:         ${PCA_MERGED_PREFIX}.snp
    indivoutname:       ${PCA_MERGED_PREFIX}.ind
EOF

    cat <<EOF > $PAR_SMARTPCA
    genotypename: ${PCA_MERGED_PREFIX}.geno
    snpname:      ${PCA_MERGED_PREFIX}.snp
    indivname:    ${PCA_MERGED_PREFIX}_withPop.ind
    poplistname:  ${POP_LIST_REF}
    evecoutname:  ${PCA_MERGED_PREFIX}.evec
    evaloutname:  ${PCA_MERGED_PREFIX}.eval
    lsqproject: YES
    newshrink: YES
    numoutlieriter: 0
    numoutevec: 10
    numthreads: ${CORES}
EOF

    # for smart pca with confidence ellipses
    log_msg "Generating parameter files..."
    cat <<EOF > $PAR_PCA_CONVERTF_CI
    genotypename:       ${PAR_PREFIX}.bed
    snpname:            ${PAR_PREFIX}.bim
    indivname:          ${PAR_PREFIX}.fam
    outputformat:       EIGENSTRAT
    genotypeoutname:    ${PCA_MERGED_PREFIX_CI}.geno
    snpoutname:         ${PCA_MERGED_PREFIX_CI}.snp
    indivoutname:       ${PCA_MERGED_PREFIX_CI}.ind
EOF

    awk '{print $1 ":" $2}' "$LDPRIO_SAMPLES" > ./smart_pca_CI/sampleList_CI.txt
    cat <<EOF > $PAR_SMARTPCA_CI
    genotypename: ${PCA_MERGED_PREFIX_CI}.geno
    snpname:      ${PCA_MERGED_PREFIX_CI}.snp
    indivname:    ${PCA_MERGED_PREFIX_CI}_withPop.ind
    poplistname:  ${POP_LIST_REF}
    evecoutname:  ${PCA_MERGED_PREFIX_CI}_CI.evec
    evaloutname:  ${PCA_MERGED_PREFIX_CI}_CI.eval
    lsqproject: YES
    newshrink: YES
    numoutlieriter: 0
    numoutevec: 2
    elllistname:  ./smart_pca_CI/sampleList_CI.txt
    elloutname:  confEllipses
    ellconf:  0.95
    numthreads: ${CORES}
EOF

    # convert
    log_msg "Converting to EIGENSTRAT format..."
    $CONVERTF_BIN -p "$PAR_PCA_CONVERTF" &>> "$LOG_FILE"
    $CONVERTF_BIN -p "$PAR_PCA_CONVERTF_CI" &>> "$LOG_FILE"

    log_msg "Assigning population labels..."
    for PREFIX in "$PCA_MERGED_PREFIX" "$PCA_MERGED_PREFIX_CI"; do
        awk '
        NR==FNR { pop[$1]=$2; next }

        $3=="???" {
            if ($1 in pop) {
                $3 = pop[$1]
            } else if (index($1, ":") > 0) {
                split($1, a, ":")
                $3 = a[1]
            } else {
                $3 = $1
            }
        }
        { print }
        ' "$LABELS_REF" "${PREFIX}.ind" > "${PREFIX}_withPop.ind"
    done

    # run smart PCA
    log_msg "Executing SmartPCA..."
    $SMARTPCA_BIN -p "$PAR_SMARTPCA" &>> "$LOG_FILE"
    $SMARTPCA_BIN -p "$PAR_SMARTPCA_CI" &>> "$LOG_FILE"

    log_msg "Fixing PCA labels in .evec files..."

    for PREFIX in "$PCA_MERGED_PREFIX" "$PCA_MERGED_PREFIX_CI"; do
        # Robustly join evec with subpopulation labels using ID lookup
        awk '
        NR==FNR { subpop[$1]=$2; next }
        {
            if ($1 ~ /^#/) { print; next }

            # Split ID (e.g., imputed:firstSample)
            if (index($1, ":") > 0) {
                split($1, a, ":")
                id = a[2]
                sp = (id in subpop) ? subpop[id] : id
            } else {
                id = $1
                sp = (id in subpop) ? subpop[id] : $NF
            }

            print $0 "\t" sp
        }
        ' "$SUBPOP_LIST_REF" "${PREFIX}.evec" > "${PREFIX}_labeled.evec"

        # Reorder evec: HGDP first, then targets
        (
            grep "^[[:space:]]*#" "${PREFIX}_labeled.evec"
            grep "[[:space:]]HGDP" "${PREFIX}_labeled.evec" || true
            grep -v "^[[:space:]]*#" "${PREFIX}_labeled.evec" \
                | grep -v "[[:space:]]HGDP" || true
        ) > "${PREFIX}_labeled_ordered.evec"

        mv "${PREFIX}_labeled_ordered.evec" "${PREFIX}_labeled.evec"

        log_msg "Final evec: ${PREFIX}_labeled.evec"
    done

    log_msg "PCA label fixing completed."
fi

# ======== PART 3: Plot Admixture ========

if [[ "$PLOT" == "admixture" || "$PLOT" == "both" ]]; then
    log_msg "PART 3: Running Admixture on ${PAR_PREFIX}..."

    ADMIXTURE_DIR="./admixture"
    mkdir -p "$ADMIXTURE_DIR"
    ADMIXTURE_LABELS_REF="${REF_DIR}/hgdp_references/hgdp_population/admixture_labels.txt"

    # Admixture writes its .Q/.P files into the current directory, so run it inside ADMIXTURE_DIR
    ADMIX_INPUT="$(realpath "${PAR_PREFIX}.bed")"
    ADMIX_FAM="$(realpath "${PAR_PREFIX}.fam")"

    log_msg "Running Admixture for K=5 to 14..."
    : > "${ADMIXTURE_DIR}/cv_comparison.txt"
    for k in {5..14}; do
        log_msg "  Processing K=$k"
        (cd "$ADMIXTURE_DIR" && "$ADMIXTURE_BIN" -j"${CORES}" --cv "$ADMIX_INPUT" "$k") > "${ADMIXTURE_DIR}/admixture_K${k}.log" 2>&1
        grep '^CV error' "${ADMIXTURE_DIR}/admixture_K${k}.log" >> "${ADMIXTURE_DIR}/cv_comparison.txt" || true
    done

    log_msg "Fixing Admixture labels..."
    # Generate the label file in the EXACT order of the .fam file used for Admixture
    RUN_LABELS="${ADMIXTURE_DIR}/run_order_labels.txt"
    awk '
    NR==FNR {
        subpop[$1]=$2; pop[$1]=$3; gender[$1]=$4; next
    }
    {
        iid = $2
        fid = $1
        # Label as FID:IID for target samples (e.g., imputed:firstSample)
        label = (fid == iid) ? iid : fid ":" iid
        if (fid == "pseudo" || fid == "imputed") {
            sp = iid
            p  = fid
        } else {
            sp = (iid in subpop) ? subpop[iid] : iid
            p  = (iid in pop)    ? pop[iid]    : iid
        }
        g  = (iid in gender) ? gender[iid] : "U"
        print label "\t" sp "\t" p "\t" g
    }
    ' "$ADMIXTURE_LABELS_REF" "$ADMIX_FAM" > "$RUN_LABELS"

    log_msg "Annotating .Q files..."
    mkdir -p "${ADMIXTURE_DIR}/Qfiles_withpop"
    for Qfile in "${ADMIXTURE_DIR}"/*.Q; do
        out="${ADMIXTURE_DIR}/Qfiles_withpop/$(basename "$Qfile").withPop"
        paste "$RUN_LABELS" "$Qfile" > "${out}.tmp"
        # Reorder Admixture output: HGDP first, then targets
        (grep "HGDP" "${out}.tmp" || true; grep -v "HGDP" "${out}.tmp" || true) > "$out"
        rm -f "${out}.tmp"
    done

    log_msg "Admixture analysis completed."
fi
