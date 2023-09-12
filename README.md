# Population-scale genotyping of SVs called from short AND/OR long reads using short reads

## Pipeline overview

1. `01_prepare_VCF.sh` : Filter input VCF to remove SVs that lack the proper info for genotyping (e.g., explicit ALT sequence) and format the VCF correctly.
2. `02_index_snarls.sh` : Build a variant-aware genome graph index and compute snarls, i.e. sites of variation. This is done ONCE for all samples.
3. `03_align_pack_call.sh` : Re-map short reads to reference genome graph (giraffe), compute read support for variant sites (pack) and call genotypes (call). This is done independently for each sample, in parallel.
4. `04_format_filter.sh` : Format each sample's genotyped VCF and filter genotype calls. This is done independently for each sample, in parallel.
5. `05_merge_samples.sh` : Merge genotyped SVs across samples to yield a single multisample VCF file.
6. `06_filter_pop.sh` : Filter the merged and genotyped SV set on missing data and minor allele frequency.

## Prerequisites

### Files
* A reference genome (`.fasta`) and its index (`.fai`) in `03_genome`

* A list of SVs to be genotyped, in the form of a VCF file. This pipeline is compatible with the output VCF from the [merge_SVs_SRLR pipeline](https://github.com/LaurieLecomte/merge_SVs_SRLR)


* A chromosomes list (or contigs, or sites) in `02_infos`. This list is used for parallelizing the SV calling step. It can be produced from the indexed genome file (`"$GENOME".fai`) : `less "$GENOME".fai | cut -f1 > 02_infos/chr.txt`. Warning : only add chromosomes or contigs for which variants are to be called (e.g. not unplaced contigs) 

### Software
