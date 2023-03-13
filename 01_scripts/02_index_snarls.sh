#!/bin/bash

# Build variant-aware genome graph index and compute snarls, i.e. sites of variation. This is done ONCE for all samples.
# Correct CPU variable and srun command if required

# manitou
# srun -c 20 -p medium --time=7-00:00:00 -J 02_index_snarls --mem=200G -o log/02_index_snarls_%j.log /bin/sh ./01_scripts/02_index_snarls.sh &

# valeria
# srun -c 20 -p ibis_medium --time=7-00:00:00 -J 02_index_snarls --mem=200G -o log/02_index_snarls_%j.log /bin/sh ./01_scripts/02_index_snarls.sh &

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

CPU=20
MEM="200G"

CANDIDATES_VCF="$VCF_DIR/candidates/"$(basename -s .ready.vcf $INPUT_VCF)".candidates.vcf.gz"

TMP_DIR="tmp"

# LOAD REQUIRED MODULES

# 0. Create tmp directory if required
if [[ ! -d $TMP_DIR ]]
then
  mkdir $TMP_DIR
fi

# 1. Build the graph with genome + unpahsed vcf and index 
echo "starting vg index"
vg autoindex --workflow giraffe -R XG --prefix $INDEX_DIR/index --tmp-dir $TMP_DIR --target-mem $MEM --threads $CPU --ref-fasta $GENOME --vcf $CANDIDATES_VCF
echo "done building index"

# 2. Compute snarls
echo "starting vg snarls"
vg snarls -t $CPU $INDEX_DIR/index.xg -f $GENOME > $SNARLS_DIR/snarls.pb
echo "done computing snarls"