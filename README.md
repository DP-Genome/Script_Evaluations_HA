# Human Genomics Pipeline

A comprehensive, modular pipeline for processing human genomic data from BAM files to downstream population genetics analyses. It includes pseudohaploid calling, genotype imputation (autosomal and X-chromosome), PCA, Admixture analysis, and GEDmatch export.

## Architecture

The pipeline follows a modular, shell-based orchestration pattern designed for genomic data processing.

1.  **Orchestration Layer (`pileupCaller_glimpse.sh`)**: Acts as the main controller, handling argument parsing, configuration loading, and sequential execution of pipeline stages.
2.  **Preprocessing & Quality Control**: Employs `samtools` for BAM filtering (quality and duplicate removal).
3.  **Variant Representation Stage**:
    *   **Pseudohaploid Calling**: Uses `pileupCaller` to generate random haploid genotypes, a common technique for low-coverage sequencing data.
    *   **Imputation Engine**: A hybrid approach using `GLIMPSE2` for autosomes and `GLIMPSE1` for the X-chromosome (`impute_chrX.sh`), optimized for low-coverage sequencing.
4.  **Analytical Layer (`pca_admixture.sh`)**:
    *   **Dimensionality Reduction**: Implements `smartpca` (EIGENSTRAT) for PCA with projection capability onto reference panels.
    *   **Ancestry Estimation**: Executes `admixture` across a range of K-values (5-14) with cross-validation.
5.  **Visualization Layer (R scripts)**:
    *   `PCAplot.R`: Generates interactive (plotly) and static (ggplot2) PCA and UMAP visualizations.
    *   `admixturePlot.R`: Produces ordered ancestry barplots across multiple K-values.
6.  **Integration & Export (`GEDmatch_export.sh`)**: Handles genomic coordinate liftover (hg38 to hg19) and formats data for GEDmatch upload.

## Pipeline Structure

### Core Scripts
- **`pileupCaller_glimpse.sh`**: The main driver script. Orchestrates the entire workflow from BAM filtering to imputation and population analysis.
- **`impute_chrX.sh`**: Specialized script for X-chromosome imputation using GLIMPSE1, integrated into the main pipeline.
- **`pca_admixture.sh`**: Performs PCA and Admixture runs, including merging with reference datasets and LD pruning.
- **`GEDmatch_export.sh`**: Converts hg38 BCF files to hg19 and generates text files ready for GEDmatch.

### Visualization Scripts
- **`PCAplot.R`**: R script for generating static and interactive PCA and UMAP plots.
- **`admixturePlot.R`**: R script for generating Admixture barplots across different K values.

### Configuration
- **`input_template.txt`**: Template for the required configuration file containing sample metadata and software paths.

## Prerequisites

The pipeline requires the following software to be installed and paths configured in the input file:
- `samtools` & `bcftools`
- `parallel` (GNU Parallel)
- `pileupCaller`
- `GLIMPSE1` & `GLIMPSE2`
- `plink` (1.9)
- `EIGENSTRAT` (convertf, smartpca)
- `admixture`
- `R` (with `ggplot2`, `plotly`, `uwot`, `RColorBrewer`)

## Usage

### 1. Prepare Input Configuration
Copy `input_template.txt` and fill in your sample details and software paths.

```bash
# Example configuration
INPUT_BAM=/path/to/sample1.bam,/path/to/sample2.bam
SAMPLE_NAME=Sample1,Sample2
SAMPLE_SEX=M,F
SAMPLE_PREFIX=my_project
...
```

### 2. Run the Pipeline
The main script `pileupCaller_glimpse.sh` accepts several flags to control the workflow:

```bash
./pileupCaller_glimpse.sh --input config.txt --pseudoHap true --panel wgs --impute true --plot both
```

**Arguments:**
- `--input`: Path to your configuration file (Required).
- `--pseudoHap`: Set to `true` to generate pseudohaploid calls.
- `--panel`: Choose `wgs` or `twist` (for 1240k panel).
- `--impute`: Set to `true` to run GLIMPSE imputation.
- `--plot`: Choose `pca`, `admixture`, `both`, or `false`.

### 3. Export to GEDmatch (Optional)
If you have an imputed BCF file, you can generate GEDmatch-compatible files:

```bash
./GEDmatch_export.sh imputation_output/imputed_samples/my_project_filtered.bcf
```

## Directory Outputs

- `clean_bams/`: Filtered and indexed BAM files.
- `pileupCaller_output/`: Pseudohaploid genotypes in Eigenstrat, PLINK, and VCF formats.
- `imputation_output/`: 
    *   `imputed_chunks/`: Intermediate BCF chunks.
    *   `imputed_chromosomes/`: Ligated per-chromosome BCFs.
    *   `imputed_samples/`: Final concatenated and filtered BCFs.
    *   `plink_files/`: Imputed data converted to PLINK format.
- `merged_plink/`: Samples merged with reference panels for analysis.
- `smart_pca/`: PCA results (`.evec`, `.eval`, logs).
- `admixture/`: Admixture Q-files and CV logs.
- `GEDmatch_input/`: Formatted text files and hg19 BCFs for GEDmatch.

## Reference Files
The pipeline expects a structured `References` directory containing:
- HGDP reference PLINK files and population labels.
- GLIMPSE2 reference panels (split into chunks).
- Human reference genomes (hg38 and hg19) and liftover chains.
- SNP lists for Twist 1240k and HGDP 4M sites.
