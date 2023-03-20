#!/bin/bash

# Re-map short reads to reference genome graph (giraffe), compute read support for variant sites (pack) and call genotypes (call). This is done independently for each sample, in parallel.

# manitou
# parallel -a 02_infos/ind_ALL.txt -j 4 srun -c 20 -p medium --time=7-00:00:00 --mem=100G -J 03_align_pack_call_{} -o log/03_align_pack_call_{}_%j.log /bin/sh ./01_scripts/03_align_pack_call.sh {} &

# valeria
# parallel -a 02_infos/ind_ALL.txt -j 3 srun -c 25 -p ibis_medium --time=7-00:00:00 --mem=100G -J 03_align_pack_call_{} -o log/03_align_pack_call_{}_%j.log /bin/sh ./01_scripts/03_align_pack_call.sh {} &


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

GENO_VCF="$MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.vcf"

MIN_MAF=0.05
MAX_MISS=0.5

FILT_GENO_VCF="$FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz"

OFFSET=5 # max distance allowed between candidate and SV for them to be merged, in bp

# LOAD REQUIRED MODULES
module load bcftools/1.13

# 1. Extract relevant fields from VCFs and export them as a table for easier handling
## candidates SVs VCF
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SVTYPE\t%END\t%SVLEN\t"%SUPP_VEC"\t[%GT\t]\n" -H $CANDIDATES_VCF >  "$VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.table"

## raw genotyped SVs VCF
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n" -H $GENO_VCF > "$MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.table"

## filtered genotyped SVs VCF
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n" -H $FILT_GENO_VCF > "$FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".table"


# 2. Merge candidates and genotyped SVs to determine genotyped SVs types
# Rscript script candidates_table geno_table filt_table allowed_pos_offset output_prefix
Rscript 01_scripts/utils/merge_candidates_genotyped.R "$VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.table" "$MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.table" "$FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".table" $OFFSET $FILT_DIR/geno_vs_cand_"$OFFSET"bp
