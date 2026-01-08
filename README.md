# hapCon_hg38

Below are the steps implemented to convert hapCon from an hg19-based software into a hg38-based one.
This enables the contamination estimation of BAM files that were generated using the GRCh38 reference.

----------------------------------------------------------

### 1. Download and preprocess the X chromosome reference
#### 1.1. Download non-PAR region only (bcf + index) from the HGDP+1kGP - https://console.cloud.google.com/storage/browser/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2;tab=objects?inv=1&invt=Ab4IFA&prefix=&forceOnObjectsSortingFiltering=false
#### 1.2. Filter out unwanted regions using the masking regions (.bed file from the PCA project)
`grep chrX ../../../PCA_admixture/HGDP_samples/hgdp_wgs.20190516.mask.bed > hgdp_wgs_chrX.20190516.mask.bed`<br><br>
`bcftools view --regions-file ../hgdp_wgs_chrX.20190516.mask.bed --threads 36 -O z -o hgdp_chrX_nonPar_masked.vcf.gz hgdp1kg_phased_haplotypes_chrX_non_par.bcf`<br><br>
`tabix -p vcf hgdp_chrX_nonPar_masked.vcf.gz`

#### 1.3. Filter out the rare alleles (MAF < 0.05)
`bcftools view -m2 -M2 -v snps hgdp_chrX_nonPar_masked.vcf.gz | bcftools filter -i 'MAF>0.05' -Oz -o hgdp_chrX_nonPar_masked_filtered.vcf.gz`<br><br>
`tabix -p vcf hgdp_chrX_nonPar_masked_filtered.vcf.gz`

#### 1.4. check the number of filtered records
`bcftools index --nrecords hgdp1kg_phased_haplotypes_chrX_non_par.bcf`,<br>
`output: 3389882`<br><br>
`bcftools index --nrecords hgdp1kg_chrX_nonPar_masked.vcf.gz`<br>
`output: 2493425` (896,457 sites filtered out - 73.5% retained)<br><br>
`bcftools index --nrecords hgdp1kg_chrX_nonPar_masked_filtered.vcf.gz`<br>
`output: 137155` (2,356,270 sites filtered out - 5.5% retained)

#### 1.5. from the 4091 samples in the vcf, subset the samples mentioned in samples_1kgp_hgdp_hg38_final.csv
`bcftools view -S creating_1kgp_hgdp_metadata/samples_1kgp_hgdp_hg38_final.csv -Oz -o ref_chrX/hgdp1kg_chrX_final_subset.vcf.gz ref_chrX/hgdp1kg_chrX_nonPar_masked_filtered.vcf.gz`<br><br>
`bcftools index ref_chrX/hgdp1kg_chrX_final_subset.vcf.gz`

#### 1.6. rename 'chrX' to 'X'
`bcftools annotate --rename-chrs chrX_rename.txt hgdp1kg_chrX_final_subset.vcf.gz | bgzip > hgdp1kg_X_final_subset.vcf.gz`<br><br>
`bcftools index hgdp1kg_X_final_subset.vcf.gz`

----------------------------------------------------------

### 2. Prepare the hdf5 file required for the hapCON

*Notes:*
- the hdf5 file provided with the hapCON script is hg37-based; we are creating another one that is hg38-based
- the vcf-to-hdf5 format conversion script was requested from the author and was provided in the hapCON github issues section
- the script was then amended by removing the filtering steps (vcf is already filtered from the steps above). Option for arguments was also removed; here we are only using one chromosome (X) not looping from 1-22 as in their script
- the script can only be run locally, not on the server
- the 'hapsburg' directory was obtained from the hapCON release downloaded from github and is needed and called by the script
  
#### 2.1. Prepare the `LD_map.snp` file for the LD (linkage disequilibrium) merge funcion
- This file has to have 6 tab-delimited columns as follows: snpID, chr, map (genetic position in cM), pos (physical base-pair position), ref, alt.
- The file was created from two starting files:<br>
  a) GLIMPSE `chrX.b38.gmap` file obtaied from their github page. It contains 3 columns: pos, chr, cM.<br>
  b) `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' hgdp1kg_X_final_subset.vcf.gz > vcf_alleles.tsv` which results in the following 4 columns: chr, pos, ref, alt.
- `chrX.b38.gmap` has 89,237 positions while `vcf_alleles.tsv` has 137,155 SNPs.
- Manually remove column headers from `chrX.b38.gmap` using `vi` then sort both files by positions:<br>
`sort -k1,1n chrX.b38.gmap > map.sorted && sort -k2,2n vcf_alleles.tsv > vcf.sorted`
- Paste columns taking into consideration the required column order:
```
awk '
# First file: map.sorted
NR==FNR && FNR>1 {
    map[$1] = $3;   # store genetic position (cM) indexed by physical position
    chr[$1] = $2;   # store chromosome name (if needed)
    next
}

# Second file: vcf.sorted
FNR>1 {
    pos = $2
    if (pos in map) {
        snpID = chr[pos] "_" pos      # e.g., X_2785350
        print snpID, chr[pos], map[pos], pos, $3, $4
    }
}
' map.sorted vcf.sorted > LD_map.snp
```
- There are 59,438 (overlapping) positions in the output file.


#### 2.2. Run the vcf-to-hdf5 conversion script (locally)
- Input and output paths are specified within the script
- Line 221 of the h5_functions.py file (in the hapsburg directory) was amended by adding `+ ["X"]` to `chs=list(range(1,23))`
- Run the script: `python3 prep_h5_ref_panel_from_vcf.py`
- Output:
```
...
Loaded 59438 Chr.X 1240K SNPs.
Intersection 59421 out of 137155 HDF5 SNPs
Interpolating 77734 variants.
Finished Chromosome X.
Adding map to HDF5...
We did it. Finished.
Run Finished Successfully!
```
- The output hdf5 file has a number groups/datasets (see `h5dump -n  hgdp1kg_X.hdf5`)

----------------------------------------------------------

### 3. Create the hg38 version of the DA43 bam input

#### 3.1. Create the bam file from fastq<br> 
`conda activate bwa_env`<br><br>
`REF="../../../../HG_analysis_pipeline/PCA_admixture/HGDP_samples/GRCh38_full_analysis_set_plus_decoy_hla.fa" && FASTQ="DA43.fastq.gz"`<br>
`bwa index "$REF"`<br>
`bwa aln -l 1024 -n 0.01 -o 2 "$REF" "$FASTQ" > "./DA43.sai"`<br>
`bwa samse "$REF" "DA43_hg38.sai" "$FASTQ" > "DA43_hg38.sam"`<br><br>
Convert to bam and sort: `samtools view -Sb DA43_hg38.sam | samtools sort -o DA43_hg38_sorted.bam && samtools index DA43_hg38_sorted.bam`<br><br>
Remove duplicates (none here): `samtools rmdup -s DA43_hg38_sorted.bam DA43_hg38_rmdup.bam && samtools index DA43_hg38_rmdup.bam`<br><br>

#### 3.2. Rename the Chr
`samtools view -H DA43_hg38_sorted.bam > header_chrX.txt`<br>
`sed 's/SN:chrX/SN:X/' header_chrX.txt > header_X.txt`<br>
`samtools reheader header_X.txt DA43_hg38_sorted.bam > DA43_hg38_X.bam`<br>
`samtools sort -o DA43_hg38_X.sorted.bam DA43_hg38_X.bam && samtools index DA43_hg38_X.sorted.bam`<br>

----------------------------------------------------------

### 4. Create the metadata file needed for hapCon
- This step was done manually<br>
final file: `g38_metadata.csv`

----------------------------------------------------------

### 5. Test run hapCON with the newly created hdf5 file

- A tutorial for running hapCON along with a test sample is found here: https://github.com/hringbauer/hapROH/blob/master/Notebooks/Vignettes/hapCon_vignette.ipynb (the "hapCon with 1000G Reference Panel" section was followed)
- hapConX doesn't run on the server, move to  local directory and execute

#### 5.1. Create a hg37 version of the hapCon bedfile using `liftover`
- Since the chain file used by liftover uses 'chrX' instead of 'X', rename the chr in the input bed file using 'sed':<br>
`sed 's/^X/chrX/' ../hapCon/reference/chrX_1000G.bed > chrX_1000G_chrRenamed.bed`<br><br>
`./liftOver chrX_1000G_chrRenamed.bed hg19ToHg38.over.chain output_1kghgdp_hg38_chrX.bed unmapped.bed`<br>
- The `unmapped.bed` is where unmapped sites are listed<br>
- Name it back to 'X': `sed 's/^chrX/X/' output_1kghgd_hg38_chrX.bed > output_1kghgdp_hg38_X.bed`

#### 5.2. Run HapCon
`source ~/venv/bin/activate`<br><br>
`samtools mpileup --positions ../bed_liftover/output_1kghgdp_hg38_X.bed -r X -q 30 -Q 30 -o DA43_hg38_X.mpileup DA43_hg38_X.sorted.bam`<br><br>
`~/venv/bin/hapConX -m DA43_hg38_X.mpileup -r hgdp1kg_X_clean.hdf5 --meta hg38_metadata.csv`<br>

----------------------------------------------------------










