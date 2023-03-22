#!/bin/bash

# Compare candidate SVs set with genotyped SVs. Since vg does not preserve info from the candidates VCF, we merge the genotyped set with the candidate set in R in order to assign a SVTYPE, an END and other info to the genotyped SVs (when possible).

# Because END and SVLEN were lost for many INVs when merging SR and LR SVs, we retrieve this info from the original SR and LR VCFs. Add these files to 05_candidates/raw/ and correct RAW_SR_VCF and RAW_LR_VCF variables if needed.
# These are the same VCFs files used in the script summarize_plot.sh from the SVs_long_reads and SVs_short_reads pipelines (https://github.com/LaurieLecomte/SVs_short_reads/blob/main/01_scripts/utils/summarize_plot.sh and https://github.com/LaurieLecomte/SVs_long_reads/blob/main/01_scripts/utils/summarize_plot.sh)


# valeria
# srun -c 1 -p ibis_small --time=1-00:00:00 --mem=50G -J compare_summarize_plot -o log/compare_summarize_plot_%j.log /bin/sh 01_scripts/utils/compare_summarize_plot.sh &

# manitou
# srun -c 1 -p small --time=1-00:00:00 --mem=50G -J compare_summarize_plot -o log/compare_summarize_plot_%j.log /bin/sh 01_scripts/utils/compare_summarize_plot.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
FASTQ_DIR="04_reads"

VCF_DIR="05_candidates"
INPUT_VCF="$VCF_DIR/raw/merged_SUPP2.ready.vcf"

GRAPH_DIR="06_graph"
INDEX_DIR="$GRAPH_DIR/index"
ALIGNED_DIR="$GRAPH_DIR/aligned"
SNARLS_DIR="$GRAPH_DIR/snarls"
PACKS_DIR="$GRAPH_DIR/packs"

CALLS_DIR="07_calls"
MERGED_DIR="08_merged"
FILT_DIR="09_filtered"

CANDIDATES_VCF="$VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.vcf.gz"

TMP_DIR="tmp"

REGIONS_EX="02_infos/excl_chrs.txt"

MIN_GQ=5
MAX_GQ=256 # deduct max allowed GQ from VCF, e.g. bcftools query -f "[%GQ\n]" $CALLS_DIR/raw/"$SAMPLE".vcf | sort -n | uniq | tail -n1

MIN_DP=4
MAX_DP=80 # arbitrary threshold of 5 * anticipated SR sequencing coverage, 5 * 16 

MIN_MAF=0.05
MAX_MISS=0.5

RAW_SR_VCF="$VCF_DIR/raw/merged_delly_manta_smoove.sorted.vcf"
RAW_LR_VCF="$VCF_DIR/raw/merged_sniffles_svim_nanovar.sorted.vcf"
GENO_VCF="$MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.vcf.gz"
FILT_GENO_VCF="$FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz"

# LOAD REQUIRED MODULES
module load bcftools/1.13

# 1. Convert inputs to tables
bcftools query -f '%CHROM\t%POS\t%ID\t%SVTYPE\t%SVLEN\t%END\t%SUPP\t%SUPP_VEC\n' $RAW_SR_VCF > $VCF_DIR/raw/"$(basename -s .sorted.vcf $RAW_SR_VCF)".table

bcftools query -f '%CHROM\t%POS\t%ID\t%SVTYPE\t%SVLEN\t%END\t%SUPP\t%SUPP_VEC\n' $RAW_LR_VCF > $VCF_DIR/raw/"$(basename -s .sorted.vcf $RAW_LR_VCF)".table
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SVTYPE\t%SVLEN\t%END\t%SUPP\t%SUPP_VEC\n' $CANDIDATES_VCF > $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".simpl.table

# 2. Add missing info for some INVs 
## Extract missing INV info from unmerged VCFs
Rscript 01_scripts/utils/add_missing_INV_info.R $VCF_DIR/raw/"$(basename -s .sorted.vcf $RAW_SR_VCF)".table $VCF_DIR/raw/"$(basename -s .sorted.vcf $RAW_LR_VCF)".table $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".simpl.table
## 

## Add missing info to candidates VCF
bgzip $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".simpl.table.annot -f
tabix -s1 -b2 -e2 $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".simpl.table.annot.gz -f
## Important to set -e at 2 (same as POS, -s2), otherwise bcftools annotate will not be able to assign the right info to the right entry 

### For some reason, bcftools annotate does NOT carry over END from annotation file to VCF, so we bypass this issue by creating a whole new tag for the corrected END, which will differ from original END tag only for INVs (OK for SVLEN tag)
echo '##INFO=<ID=END_corr,Number=1,Type=Integer,Description="End position of structural variation">' > $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".END.hdr

bcftools annotate -a $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".simpl.table.annot.gz -c CHROM,POS,ID,SVTYPE,SVLEN,END_corr $CANDIDATES_VCF -h $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".END.hdr | bcftools sort -Oz > $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)"_correxINVs.vcf.gz 


# 3. Infer info about RAW genotyped SVs from candidates (merged, but unfiltered calls in 08_merged)
## Convert input candidates VCF to table
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SVTYPE\t%SVLEN\t%END_corr\t%SUPP\t%SUPP_VEC\t[%GT\t]\n' -H $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)"_correxINVs.vcf.gz > $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".table

## Convert input raw genotyped SVs VCF to table
#bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n" -H $GENO_VCF > "$MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.table"
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t%DP\t%GQ\t]\n" -H $GENO_VCF > "$MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.table"

## Match raw genotyped with candidates
Rscript 01_scripts/merge_compare_candidates_raw_geno.R $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".table $MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.table" 5 

# 4. Infer info about FILTERED genotyped SVs from candidates 
## Convert input filtered genotyped SVs VCF to table
#bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n" -H $FILT_GENO_VCF > "$FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".table"
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t%DP\t%GQ]\n" -H $FILT_GENO_VCF > "$FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".table"

## Match filtered genotyped with candidates
Rscript 01_scripts/merge_compare_candidates_filt_geno.R $VCF_DIR/candidates/"$(basename -s .vcf.gz $CANDIDATES_VCF)".table $MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.table" 5 

# Clean up