#!/bin/bash

# Format each sample's VCF and filter genotype calls. This is done independently for each sample, in parallel.
# VERY IMPORTANT : deactivate conda env first because setGT will not work if older bcftools version is loaded by conda by default !

# manitou
# parallel -a 02_infos/ind_ALL.txt -j 4 srun -c 1 -p small --time=1-00:00:00 --mem=10G -J 04_format_filter_{} -o log/04_format_filter_{}_%j.log /bin/sh ./01_scripts/04_format_filter.sh {} &

# valeria
# parallel -a 02_infos/ind_ALL.txt -j 4 srun -c 1 -p ibis_small --time=1-00:00:00 --mem=10G -J 04_format_filter_{} -o log/04_format_filter_{}_%j.log /bin/sh ./01_scripts/04_format_filter.sh {} &


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

CPU=1
MEM="100G"

CANDIDATES_VCF="$VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.vcf.gz"

TMP_DIR="tmp"

SAMPLE=$1

REGIONS_EX="02_infos/excl_chrs.txt"

MIN_GQ=5
MAX_GQ=256 # deduct max allowed GQ from VCF, e.g. bcftools query -f "[%GQ\n]" $CALLS_DIR/raw/"$SAMPLE".vcf | sort -n | uniq | tail -n1

MIN_DP=4
MAX_DP=80 # arbitrary threshold of 5 * anticipated SR sequencing coverage, 5 * 16 

# LOAD REQUIRED MODULES
module load bcftools/1.15

# 1. Rename sample in output VCF, since vg call -s $ID does not output sample name in output VCF, and remove unwanted contigs from header
echo -e "SAMPLE\t$SAMPLE" > 02_infos/"$SAMPLE".names
bcftools reheader -s 02_infos/"$SAMPLE".names $CALLS_DIR/raw/"$SAMPLE".vcf | grep -vFf $REGIONS_EX > $CALLS_DIR/raw/"$SAMPLE".tmp

# 2. Remove non biallelic sites, and assign missing genotype where DP and GQ are too low or too extreme
#bcftools filter -i "FORMAT/GQ >= $MIN_GQ & FORMAT/DP >= $MIN_DP" $CALLS_DIR/raw/"$SAMPLE".tmp > $FILT_DIR/"$SAMPLE"_GQ"$MIN_GQ"_DP"$MIN_DP".vcf
bcftools view -m2 -M2 $CALLS_DIR/raw/"$SAMPLE".tmp | bcftools +setGT -- -t q -n . -e "FORMAT/DP >= $MIN_DP & FORMAT/DP < $MAX_DP & FORMAT/GQ >= $MIN_GQ & FORMAT/GQ < $MAX_GQ" | bcftools sort > $CALLS_DIR/filtered/"$SAMPLE"_DP"$MIN_DP"_GQ"$MIN_GQ".vcf

# 2. Compress 
bgzip $CALLS_DIR/filtered/"$SAMPLE"_DP"$MIN_DP"_GQ"$MIN_GQ".vcf -f
tabix -p vcf $CALLS_DIR/filtered/"$SAMPLE"_DP"$MIN_DP"_GQ"$MIN_GQ".vcf.gz -f

# Count SVs
## genotyped
zless $CALLS_DIR/filtered/"$SAMPLE"_DP"$MIN_DP"_GQ"$MIN_GQ".vcf.gz | grep -v ^'#' | wc -l 
## calls other than ./.
bcftools filter -i 'GT!="mis"' $CALLS_DIR/filtered/"$SAMPLE"_DP"$MIN_DP"_GQ"$MIN_GQ".vcf.gz | grep -v ^'#' | wc -l 

# Clean up 
rm 02_infos/"$SAMPLE".names
rm $CALLS_DIR/raw/"$SAMPLE".tmp