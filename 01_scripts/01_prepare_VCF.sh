#!/bin/bash

# Filter and format input VCF prior to building genome graph index

# manitou
# srun -c 1 -p small -J 01_prepare_VCF --mem=10G -o log/01_prepare_VCF_%j.log /bin/sh ./01_scripts/01_prepare_VCF.sh &

# valeria
# srun -c 1 -p ibis_small -J 01_prepare_VCF --mem=10G -o log/01_prepare_VCF_%j.log /bin/sh ./01_scripts/01_prepare_VCF.sh &

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

# LOAD REQUIRED MODULES
# module load bcftools/1.13 tabix 


# 1. Candidate SVs : Extract SVs with explicit ALT sequence : these CAN be genotyped and are candiates Remove SVs with no seq in ALT field
## Extract header
bcftools view -h $INPUT_VCF > $VCF_DIR/raw/"$(basename -s .ready.vcf $INPUT_VCF)".header

## Extract lines where 5th field does starts by a base or N
#grep -v ^'\#' $VCF_DIR/raw/"$(basename -s .vcf $INPUT_VCF)".simpl.vcf | awk 'BEGIN { OFS=FS="\t" } $5 ~ /^(A|a|T|t|G|g|C|c|N)/' > $VCF_DIR/raw/"$(basename -s .vcf $INPUT_VCF)".simpl.ALT.contents
grep -v ^'\#' $INPUT_VCF | awk 'BEGIN { OFS=FS="\t" } $5 ~ /^(A|a|T|t|G|g|C|c|N)/' > $VCF_DIR/raw/"$(basename -s .ready.vcf $INPUT_VCF)".ALT.contents

## Sort 'lexicographically'
less $VCF_DIR/raw/"$(basename -s .ready.vcf $INPUT_VCF)".ALT.contents | sort -k1,1 -k2,2n > $VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".ALT.contents.sorted

## Concatenate header and sorted contents
cat $VCF_DIR/raw/"$(basename -s .ready.vcf $INPUT_VCF)".header $VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".ALT.contents.sorted > $VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.vcf


# 2. Blacklisted SVs : Extract SVs with no explicit ALT seq, will not be genotyped
## Extract lines where 5th field does NOT start by a base or N
grep -v ^'\#' $INPUT_VCF | awk 'BEGIN { OFS=FS="\t" } $5 !~ /^(A|a|T|t|G|g|C|c|N)/' | awk 'BEGIN { OFS=FS="\t" } $5 ~ /^(<|[|]>|N|.)/' > $VCF_DIR/blacklisted/"$(basename -s .ready.vcf $INPUT_VCF)".noALT.contents

## Concatenate header and blacklisted contents
cat $VCF_DIR/raw/"$(basename -s .ready.vcf $INPUT_VCF)".header $VCF_DIR/blacklisted/"$(basename -s .ready.vcf $INPUT_VCF)".noALT.contents > $VCF_DIR/blacklisted/"$(basename -s .ready.vcf $INPUT_VCF)".blacklisted.vcf


# Remove STRANDS field : do before if possible
#sed -i -E "s/\;STRANDS\=\?\?//"  $VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.vcf


# 4. Compress and index 
bgzip $VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.vcf -f 
tabix $VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.vcf.gz -f 

# Clean up 
rm $VCF_DIR/raw/"$(basename -s .ready.vcf $INPUT_VCF)".header
rm $VCF_DIR/raw/"$(basename -s .ready.vcf $INPUT_VCF)".ALT.contents
rm $VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".ALT.contents.sorted
rm $VCF_DIR/blacklisted/"$(basename -s .ready.vcf $INPUT_VCF)".noALT.contents