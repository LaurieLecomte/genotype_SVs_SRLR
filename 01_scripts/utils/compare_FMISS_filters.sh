#!/bin/bash

# Filter merged and genotyped SVs set, final step ! 

# manitou
# srun -c 1 -p small --time=1-00:00:00 --mem=10G -J compare_FMISS_filters -o log/compare_FMISS_filters_%j.log /bin/sh ./01_scripts/utils/compare_FMISS_filters.sh &

# valeria
# srun -c 1 -p ibis_small --time=1-00:00:00 --mem=10G -J compare_FMISS_filters -o log/compare_FMISS_filters_%j.log /bin/sh ./01_scripts/utils/compare_FMISS_filters.sh &

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

#MAX_MISS=0.5

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load htslib/1.13


# 1. Filter on MAF and proportion of missing genotypes, with F_MISS = 0.125, 0.25, 0.375 

for MAX_MISS in 0.75 0.625 0.5 0.375 0.25 0.125 ; 
do 
  bcftools view --max-alleles 2 $MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.tagged.vcf.gz | bcftools filter -i "INFO/MAF >= $MIN_MAF & INFO/MAF <= $MAX_MAF & INFO/F_MISSING <= $MAX_MISS" | bcftools sort -Oz > $FILT_DIR/var_FMISS/"$(basename -s .ready.vcf $INPUT_VCF)"_"$MIN_DP"_GQ"$MIN_GQ"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz
  echo "with $MAX_MISS : $(zless $FILT_DIR/var_FMISS/"$(basename -s .ready.vcf $INPUT_VCF)"_"$MIN_DP"_GQ"$MIN_GQ"_MAF"$MIN_MAF"_FMISS"$MAX_MISS".vcf.gz | grep -v ^'#' | wc -l) variants" ;
done


