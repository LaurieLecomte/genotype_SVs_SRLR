# Compute mean, min, max and sd of per-site total DP and F_MISSING, for raw and filtered data


# 1. Access file and import -----------------------------------------------
RAW_STATS <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/08_merged/merged_SUPP2_genotyped.tagged_DP_FMISS.table'
FILT_STATS <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5_DP_FMISS.table'

library(data.table)
library(dplyr)
library(ggplot2)

raw_stats <- fread(RAW_STATS, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))
filt_stats <- fread(FILT_STATS, col.names = c('CHROM', 'POS', 'DP', 'F_MISSING'))

# 2. Compute basic stats --------------------------------------------------
variants_stats <- function(x){
  # Total 
  print(paste('number of variants :', length(x)))
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


# 3. Explore F_MISSING ----------------------------------------------------
# Split F_MISSING into equal size bins
raw_stats$F_MISS_bins <- cut_interval(raw_stats$F_MISSING, 
                                    length = 0.1, right = FALSE)

# Assign 'kept' or 'filtered out' tag to each SV
raw_stats$F_MISS_groups <-
sapply(X = raw_stats$F_MISS_bins,
       FUN = function(x){
         ifelse(x %in% levels(raw_stats$F_MISS_bins)[1:5], 
                yes = 'passed',
                no = 'failed')
       })


# Plot
SVs_F_MISS_plot <- 
ggplot(data = raw_stats) +
  geom_bar(aes(F_MISS_bins, fill = F_MISS_groups)) + 
  theme(
    axis.text.x = element_text(angle = 45, size = 8, hjust = 1)
  ) + scale_fill_manual(values = c('red', 'grey60')) +
  labs(x = 'Proportion of missing genotypes',
       y = 'SV count',
       fill = 'F_MISSING filter')

saveRDS(SVs_F_MISS_plot,
        file = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5_DP_FMISS_FMISSplot.rds')


