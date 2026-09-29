#!/usr/bin/env bash
# Version: 1.2.0
#   1.2.0 - options: PCA_EVECS (--pca-evecs, eigenvectors written by smartpca, default 10) and RUN_CI
#           (--ci true|false, skip the confidence-ellipse smartpca run). Fixed CI evec labelling, which
#           read a non-existent file and aborted the script; it now writes <prefix>_CI_labeled.evec.
#   1.1.1 - ldprio v1.5 compatible (--autosomes-only removed in ldprio; --plink passed explicitly;
#           ldprio output goes to the pipeline log).
#   1.1.0 - optional transversion-only filter (TRANSVERSIONS=true / --transversions true) applied
#           to the merged reference+sample set before LD pruning, PCA and Admixture; outputs get a _tv suffix.
#   1.0.1 - PCA/Admixture use the pruned set that matches --prune (ldprio or indep);
#           new-sample list no longer hard-coded to line 930; admixture runs in its own
#           directory (no more moving *.log); threads come from $CORES; PRUNE/PLOT defaults.
#   1.0.0 - initial tagged release (see OLD_VERSION_SCRIPTS/)

set -euo pipefail
: "${PRUNE:=ldprio}"
: "${PLOT:=both}"
: "${CORES:=24}"
: "${TRANSVERSIONS:=false}"
: "${PCA_EVECS:=10}"
: "${RUN_CI:=true}"
TV_SUFFIX=""

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

# ======== Optional transversion-only filter ========
# Keeps only A/C, A/T, G/C, G/T SNPs. Transitions (A/G, C/T) and non-ACGT alleles are removed,
# which avoids post-mortem damage (C>T / G>A) affecting PCA and Admixture.
if [[ "$TRANSVERSIONS" == "true" ]]; then
    log_msg "Filtering merged dataset to transversions only..."
    TV_IDS="${PLINK_MERGED_PREFIX}_transversion_ids.txt"
    awk '{ a=toupper($5); b=toupper($6)
           if (a !~ /^[ACGT]$/ || b !~ /^[ACGT]$/ || a == b) next
           p = (a < b) ? a b : b a
           if (p == "AG" || p == "CT") next
           print $2 }' "${PLINK_MERGED_PREFIX}.bim" > "$TV_IDS"
    log_msg "Transversions kept: $(wc -l < "$TV_IDS") of $(wc -l < "${PLINK_MERGED_PREFIX}.bim") variants"
    if [[ ! -s "$TV_IDS" ]]; then
        log_msg "ERROR: no transversions found in merged dataset."
        exit 1
    fi
    $PLINK_BIN --bfile "${PLINK_MERGED_PREFIX}" --extract "$TV_IDS" --make-bed --out "${PLINK_MERGED_PREFIX}_tv" &>> "$LOG_FILE"
    PLINK_MERGED_PREFIX="${PLINK_MERGED_PREFIX}_tv"
    PLINK_MERGED_PREFIX_LDPRIO="${PLINK_MERGED_PREFIX_LDPRIO}_tv"
    TV_SUFFIX="_tv"
fi

# Prefix of the pruned fileset used by PCA and Admixture (depends on --prune)
case "$PRUNE" in
    indep)  PAR_PREFIX="${PLINK_MERGED_PREFIX}_pruned" ;;
    ldprio) PAR_PREFIX="${PLINK_MERGED_PREFIX_LDPRIO}_ldprioPruned" ;;
    *) log_msg "ERROR: PRUNE must be indep or ldprio (got '$PRUNE')"; exit 1 ;;
esac

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
        --ld-samples "${REF_DIR}/hgdp_references/modernSamples.txt" --plink "$PLINK_BIN" --make-bed --out "${PLINK_MERGED_PREFIX_LDPRIO}_ldprioPruned" &>> "$LOG_FILE"
fi
# ======== PART 2: Plot PCA ========

if [[ "$PLOT" == "pca" || "$PLOT" == "both" ]]; then
    log_msg "PART 2: Running SmartPCA..."

    # Population list to project on
    POP_LIST_REF="${REF_DIR}/hgdp_references/hgdp_population/populationList.txt"
    SUBPOP_LIST_REF="${REF_DIR}/hgdp_references/hgdp_population/subpopulationList.txt"
    LABELS_REF="${REF_DIR}/hgdp_references/hgdp_population/population_labels.tsv"

    # Define output paths
    mkdir -p "smart_pca"
    [[ "$RUN_CI" == "true" ]] && mkdir -p "smart_pca_CI"
    PCA_MERGED_PREFIX="./smart_pca/merged_ref${TV_SUFFIX}"
    PCA_MERGED_PREFIX_CI="./smart_pca_CI/merged_ref${TV_SUFFIX}"

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
    numoutevec: ${PCA_EVECS}
    numthreads: ${CORES}
EOF

    if [[ "$RUN_CI" == "true" ]]; then
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
    fi

    # convert
    log_msg "Converting to EIGENSTRAT format..."
    $CONVERTF_BIN -p "$PAR_PCA_CONVERTF" &>> "$LOG_FILE"
    [[ "$RUN_CI" == "true" ]] && $CONVERTF_BIN -p "$PAR_PCA_CONVERTF_CI" &>> "$LOG_FILE"

    log_msg "Assigning population labels..."
    PCA_PREFIXES=("$PCA_MERGED_PREFIX")
    [[ "$RUN_CI" == "true" ]] && PCA_PREFIXES+=("$PCA_MERGED_PREFIX_CI")
    for PREFIX in "${PCA_PREFIXES[@]}"; do
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
    [[ "$RUN_CI" == "true" ]] && $SMARTPCA_BIN -p "$PAR_SMARTPCA_CI" &>> "$LOG_FILE"

    log_msg "Fixing PCA labels in .evec files..."

    EVEC_IN=("${PCA_MERGED_PREFIX}.evec");            EVEC_OUT=("${PCA_MERGED_PREFIX}_labeled.evec")
    if [[ "$RUN_CI" == "true" ]]; then
        EVEC_IN+=("${PCA_MERGED_PREFIX_CI}_CI.evec"); EVEC_OUT+=("${PCA_MERGED_PREFIX_CI}_CI_labeled.evec")
    fi
    for i in "${!EVEC_IN[@]}"; do
        EVEC_FILE="${EVEC_IN[i]}"; LABELED="${EVEC_OUT[i]}"
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
        ' "$SUBPOP_LIST_REF" "$EVEC_FILE" > "$LABELED"

        # Reorder evec: HGDP first, then targets
        (
            grep "^[[:space:]]*#" "$LABELED" || true
            grep "[[:space:]]HGDP" "$LABELED" || true
            grep -v "^[[:space:]]*#" "$LABELED" \
                | grep -v "[[:space:]]HGDP" || true
        ) > "${LABELED}.tmp"

        mv "${LABELED}.tmp" "$LABELED"

        log_msg "Final evec: $LABELED"
    done

    log_msg "PCA label fixing completed."
fi

# ======== PART 3: Plot Admixture ========

if [[ "$PLOT" == "admixture" || "$PLOT" == "both" ]]; then
    log_msg "PART 3: Running Admixture on ${PAR_PREFIX}..."

    ADMIXTURE_DIR="./admixture${TV_SUFFIX}"
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
