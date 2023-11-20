#!/bin/bash

# Compare genotypes called by LR callers with genotypes outputted by VG using short reads

# srun -p large -c 8 --mem=150G --time=21-00:00:00 -J addedN_mummer_self_synteny -o log/02_nucmer_self_synteny_addedN_%j.log /bin/sh 01_scripts/02_nucmer_symap_maxmatch_addedN.sh &


# VARIABLES
LR_SAMPLES="/project/lbernatchez/users/lalec31/RDC_Romaine/02_long_reads/SVs_long_reads/02_infos/ind_ONT.txt"
LR_VCF="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/merge_SVs_SRLR/04_vcf/LR/merged_sniffles_svim_nanovar_SUPP2.vcf"
SR_VCF="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/merge_SVs_SRLR/04_vcf/SR/merged_delly_manta_smoove_SUPP2.vcf"

RAW_VG_VCF="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/08_merged/merged_SUPP2_genotyped.tagged.vcf.gz"
FILT_VG_VCF="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"

OUT_DIR="GT_concordance"


# LOAD REQUIRED MODULES
module load bcftools/1.13


if [[ ! -d $OUT_DIR ]]
then
  mkdir $OUT_DIR
fi

# 1. From LR sample list, get corresponding IDs in each VCF
## For VG VCF
bcftools query -l $RAW_VG_VCF | grep -f $LR_SAMPLES > $OUT_DIR/raw_LR_IDs.txt 

bcftools query -l $FILT_VG_VCF | grep -f $LR_SAMPLES > $OUT_DIR/filt_LR_IDs.txt

## For the SR SV VCF
bcftools query -l $SR_VCF | grep -f $LR_SAMPLES > $OUT_DIR/filt_SR_IDs.txt

# 2. Extract required fields in each VCF
## LR SV calls
bcftools query -f '%CHROM\t%POS\t%END\t%ID\t[%GT\t]\n' -H $LR_VCF > $OUT_DIR/"$(basename -s .vcf $LR_VCF)"_GTs.table

## SR SV calls, only for LR samples
bcftools view -S $OUT_DIR/raw_LR_IDs.txt $SR_VCF | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t[%GT\t]\n' -H > $OUT_DIR/"$(basename -s .vcf.gz $SR_VCF)"_LR_samples_GTs.table

## Raw vg genotypes
bcftools view -S $OUT_DIR/raw_LR_IDs.txt $RAW_VG_VCF | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t[%GT\t]\n' -H > $OUT_DIR/"$(basename -s .vcf.gz $RAW_VG_VCF)"_GTs.table


## Filtered vg genotypes
bcftools view -S $OUT_DIR/filt_LR_IDs.txt $FILT_VG_VCF | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t[%GT\t]\n' -H > $OUT_DIR/"$(basename -s .vcf.gz $FILT_VG_VCF)"_GTs.table

