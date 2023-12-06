#!/bin/bash

# srun -p medium --time=7-00:00:00 -c 1 -J lastz_align_blocks -o log/lastz_align_blocks_%j.log /bin/sh 01_scripts/lastz_align_blocks.sh &

# VARIABLES
MERGED="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/merge_SVs_SRLR/06_filtered/merged_SUPP2.corrected.table"
#MERGED="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/merge_SVs_SRLR/06_filtered/merged_SUPP2.ready.vcf"

GENO_RAW="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/08_merged/merged_SUPP2_genotyped_matched_offset5bp.txt"
GENO_FILT="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5_matched_offset5bp.txt"

HOM_BLOCKS="XL_SVs/homolog_blocks_identity.bed"

WIN_SIZE=1000000
CHR_BED="02_infos/chrs.bed"


# LOAD REQUIRED MODULES
module load bedtools
module load bcftools/1.13
module load bedops


bedops --chop $WIN_SIZE -x $CHR_BED > 02_infos/chrs_win"$WIN_SIZE".bed


bedtools window -a  -b $HOM_BLOCKS -w 100 > XL_SVs/win100_filt_excl_SNPs_hom_regions.table

bedtools window -a XL_SVs/filtered_excluded_SNPs.bed -b "XL_SVs/genome.fasta.out.table" -w 100 > XL_SVs/win100_filt_excl_SNPs_RM.table