# Compute mean, min, max and sd of per-site total DP and F_MISSING, for raw and filtered data


# 1. Access file and import -----------------------------------------------
RAW_STATS <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/08_merged/merged_SUPP2_genotyped.tagged_DP_FMISS.table'
FILT_STATS <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5_DP_FMISS.table'

library(data.table)

raw_stats <- fread(RAW_STATS, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))
filt_stats <- fread(FILT_STATS, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))

# 2. Compute basic stats --------------------------------------------------
variants_stats <- function(x){
  # Compute mean
  print(paste('mean :', mean(x)))
  # Compute min and max
  print(paste('min :', min(x)))
  print(paste('max :', max(x)))
  # Compute sd
  print(paste('sd :', sd(x)))
}

# On raw, unfiltered data
variants_stats(raw_stats$DP)
variants_stats(raw_stats$F_MISSING)

# On filtered data
variants_stats(filt_stats$DP)
variants_stats(filt_stats$F_MISSING)
