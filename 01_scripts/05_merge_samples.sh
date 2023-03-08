#!/bin/bash

# Merge genotyped SVs across samples to yield a single multisample VCF file

# manitou
# srun -c 1 -p small --time=1-00:00:00 -J --mem=10G 05_merge_samples -o log/05_merge_samples_%j.log /bin/sh ./01_scripts/05_merge_samples.sh &

# valeria
# srun -c 1 -p ibis_small --time=1-00:00:00 --mem=10G -J 05_merge_samples -o log/05_merge_samples_%j.log /bin/sh ./01_scripts/05_merge_samples.sh &

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
MERGED_DIR="08_MERGED"
FILT_DIR="09_filtered"

CANDIDATES_VCF="$VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.vcf.gz"

TMP_DIR="tmp"

REGIONS_EX="02_infos/excl_chrs.txt"

MIN_GQ=5
MAX_GQ=256 # deduct max allowed GQ from VCF, e.g. bcftools query -f "[%GQ\n]" $CALLS_DIR/raw/"$SAMPLE".vcf | sort -n | uniq | tail -n1

MIN_DP=4
MAX_DP=80 # arbitrary threshold of 5 * anticipated SR sequencing coverage, 5 * 16 


# LOAD REQUIRED MODULES
module load bcftools/1.13

# 1. Make a list of VCF files to merge
ls -1 $CALLS_DIR/filtered/*_DP"$MIN_DP"_GQ"$MIN_GQ".vcf.gz > 02_infos/samples_VCF.txt

# 2. Merge genotyped calls across samples
bcftools merge 02_infos/samples_VCF.txt | bcftools sort > $MERGED_DIR/"$(basename -s .ready.vcf $INPUT_VCF)"_genotyped.vcf