#!/bin/bash

# Re-map short reads to reference genome graph (giraffe), compute read support for variant sites (pack) and call genotypes (call). This is done independently for each sample, in parallel.

# manitou
# parallel -a 02_infos/ind_ALL.txt -j 4 srun -c 10 -p medium --time=7-00:00:00 --mem=100G -J 03_align_pack_call_{} -o log/03_align_pack_call_{}_%j.log /bin/sh ./01_scripts/03_align_pack_call.sh {} &

# valeria
# parallel -a 02_infos/ind_ALL.txt -j 4 srun -c 10 -p ibis_medium --time=7-00:00:00 --mem=100G -J 03_align_pack_call_{} -o log/03_align_pack_call_{}_%j.log /bin/sh ./01_scripts/003_align_pack_call.sh {} &


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
FASTQ1="$FASTQ_DIR/"$SAMPLE"_1.trimmed.fastq.gz"
FASTQ2="$FASTQ_DIR/"$SAMPLE"_2.trimmed.fastq.gz" 

# 1. Map paired short reads to the reference graph structure 
echo "Aligning $FASTQ1 and $FASTQ2 for $SAMPLE"
vg giraffe -t $CPU -x $INDEX_DIR/index.xg -Z $INDEX_DIR/index.giraffe.gbz -m $INDEX_DIR/index.min -d $INDEX_DIR/index.dist -f $FASTQ1 -f $FASTQ2 -N $SAMPLE -p > $ALIGNED_DIR/"$SAMPLE".gam
  
# 2. Pack the alignments
echo "Packing loops for $SAMPLE"
vg pack -t $CPU -Q 5 -x $INDEX_DIR/index.xg -g $ALIGNED_DIR/"$SAMPLE".gam -o $PACKS_DIR/"$SAMPLE".pack

# 3. Call
echo "Calling genotypes for $SAMPLE"
vg call -t $CPU -a $INDEX_DIR/index.xg -k $PACKS_DIR/"$SAMPLE".pack -r $SNARLS_DIR/snarls.pb -f $GENOME -s "$SAMPLE"  > $CALLS_DIR/raw/"$SAMPLE".vcf


