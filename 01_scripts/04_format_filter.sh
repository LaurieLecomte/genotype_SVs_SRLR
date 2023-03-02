#!/bin/bash

# Format each sample's VCF and filter genotype calls. This is done independently for each sample, in parallel.

# manitou
# parallel -a 02_infos/ind_ALL.txt -j 4 srun -c 1 -p small --time=1-00:00:00 -J --mem=10G 04_format_filter_{} -o log/04_format_filter_{}_%j.log /bin/sh ./01_scripts/04_format_filter.sh {} &

# valeria
# parallel -a 02_infos/ind_ALL.txt -j 4 srun -c 1 -p ibis_small --time=1-00:00:00 --mem=10G -J 04_format_filter_{} -o log/04_format_filter_{}_%j.log /bin/sh ./01_scripts/04_format_filter.sh {} &


# VARIABLES
GENOME="03_genome/genome.fasta"
FASTQ_DIR="04_reads"

VCF_DIR="05_candidates"
INPUT_VCF="$VCF_DIR/raw/merged_SUPP2.vcf"

GRAPH_DIR="06_graph"
INDEX_DIR="$GRAPH_DIR/index"
ALIGNED_DIR="$GRAPH_DIR/aligned"
SNARLS_DIR="$GRAPH_DIR/snarls"
PACKS_DIR="$GRAPH_DIR/packs"

CALLS_DIR="07_calls"
MERGED_DIR="08_MERGED"
FILT_DIR="09_filtered"

CPU=10
MEM="100G"

CANDIDATES_VCF="$VCF_DIR/candidates/"$(basename -s .vcf $INPUT_VCF)".candidates.vcf.gz"

TMP_DIR="tmp"

SAMPLE=$1

# LOAD REQUIRED MODULES
module load bcftools/1.13

# 1. Rename sample in output VCF, since vg call -s $ID does not output sample name in output VCF
echo -e "SAMPLE\t$SAMPLE" > 02_infos/"$SAMPLE".names
bcftools reheader -s 02_infos/"$SAMPLE".names $CALLS_DIR/raw/"$SAMPLE".vcf > $CALLS_DIR/raw/"$SAMPLE".vcf
rm 02_infos/"$SAMPLE".names

# 2. Compress 
#bgzip $CALLS_DIR/"$SAMPLE".vcf
#tabix $CALLS_DIR/"$SAMPLE".vcf.gz