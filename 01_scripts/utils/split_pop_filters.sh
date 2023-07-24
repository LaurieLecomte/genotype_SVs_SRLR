#!/bin/bash

# Split pop filters for exploration purposes

# manitou
# srun -c 1 -p small --time=1-00:00:00 --mem=10G -J split_pop_filter -o log/split_pop_filter_%j.log /bin/sh ./01_scripts/utils/split_pop_filters.sh &

# valeria
# srun -c 1 -p ibis_small --time=1-00:00:00 --mem=10G -J split_pop_filter -o log/split_pop_filter_%j.log /bin/sh ./01_scripts/utils/split_pop_filters.sh &

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
MAX_MAF=0.95
MAX_MISS=0.5

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load htslib/1.13

# 1. Add tags : ALREADY DONE
#bcftools +fill-tags $MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.vcf.gz -Oz -- -t AC,AC_Hom,AC_Het,AC_Hemi,AF,AN,NS,ExcHet,HWE,MAF,F_MISSING,END > $MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.tagged.vcf.gz

# 2. Filter on proportion of missing genotypes only
#bcftools view --max-alleles 2 $MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.tagged.vcf.gz | bcftools filter -i "INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_FMISS"$MAX_MISS".vcf.gz
#tabix -p vcf $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_FMISS"$MAX_MISS".vcf.gz

#echo "number of SVs with F_MISSING <= $MAX_MISS : "$(zless $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_FMISS"$MAX_MISS".vcf.gz | grep -v ^# | wc -l)""

# 3. Filter on MAF only
#bcftools view --max-alleles 2 $MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF".vcf.gz
#tabix -p vcf $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF".vcf.gz

#echo "number of SVs with MAF >= $MIN_MAF and <= $MAX_MAF : "$(zless $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF".vcf.gz | grep -v ^# | wc -l)""


bcftools view --max-alleles 2 $MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.tagged.vcf.gz | bcftools filter -i "INFO/MAF > $MAX_MAF" | bcftools sort -Oz > $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_excessMAF"$MAX_MAF".vcf.gz

# Check overlap of both
#bcftools isec $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_FMISS"$MAX_MISS".vcf.gz $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_MAF"$MIN_MAF".vcf.gz -p $FILT_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_FMISS_vs_MAF
