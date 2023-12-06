# Check if difference in missing GT between SR and LR

library(ggplot2)

# Check missing data distribution in raw genotypes ------------------------
# Import
MATCHED_SV <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/08_merged/merged_SUPP2_genotyped_matched_offset5bp.txt"

MISSING_GT <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/08_merged/merged_SUPP2_genotyped.tagged.table"

matched_SVs <- read.table(MATCHED_SV, header = TRUE, colClasses = c(CHROM = 'character',
                                                                    POS = 'numeric',
                                                                    ID = 'character',
                                                                    CAND_ID = 'character',
                                                                    CAND_SVTYPE = 'character',
                                                                    CAND_SVLEN = 'numeric',
                                                                    CAND_SUPP_VEC = 'character'))

  
missing_GT <- read.table(MISSING_GT, col.names = c('CHROM', 'POS', 'END', 'ID', 'F_MISS'))

# Get platform support info by merging with genotyped SVs sucessfully matched to known candidates
matched_SVs_miss <- 
merge(x = matched_SVs, y = missing_GT,
      by = c('CHROM', 'POS', 'ID'))


# Add explicit platform info
matched_SVs_miss$platform <- sapply(X = as.character(matched_SVs_miss$CAND_SUPP_VEC), 
                        FUN = function(x){
                          switch(x,
                                 '10' = 'LR', 
                                 '01' = 'SR',
                                 '11' = 'SR + LR')
                        }
)

# Plot missing data in raw genotyped SVs by platform
ggplot(data = matched_SVs_miss) + 
  geom_boxplot(aes(x = platform, y = F_MISS))

ggplot(data = matched_SVs_miss) + 
  geom_violin(aes(x = platform, y = F_MISS, fill = platform)) + 
  theme_bw() +
  guides(fill = 'none') +
  labs(x = 'Platform support', y = 'Missing genotype proportion') +
  stat_summary(aes(x = platform, y = F_MISS), fun.y = "median", geom = "point") +
  geom_text(aes(y = 1, x = 'SR + LR'), label = 'A', hjust = -1.3, vjust = 1, size = 10)

# Save to external file
ggsave(filename = paste0(unlist(strsplit(MATCHED_SV, split = '.txt'))[1], 
                           '_FMISS_by_platform.png'),
         width = 2600,
         height = 2000,
         units = 'px',
         dpi = 600
  )
  

# Check missing data distribution in SR and LR seperatly ------------------
# split SR and LR
for (i in 1:nrow(matched_SVs_miss)){
  matched_SVs_miss$SR[i] <- 
  switch(matched_SVs_miss$CAND_SUPP_VEC[i],
         '01' = 1,
         '11' = 1,
         '10' = 0)
}
#for (i in 1:nrow(matched_SVs_miss)){
#  matched_SVs_miss$LR[i] <- 
#    switch(matched_SVs_miss$CAND_SUPP_VEC[i],
#           '01' = 0,
#           '11' = 1,
#           '10' = 1)
#}

  
ggplot(data = matched_SVs_miss) + 
  geom_boxplot(aes(x = factor(SR), y = F_MISS))

# By type
ggplot(data = matched_SVs_miss) + geom_boxplot(aes(x = factor(SR), y = F_MISS))



# Check missing data distribution in filtered SVs -------------------------
# Filtered vs non filtered
FILTERED_TABLE <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5.tagged.table"
FILTERED_MATCHED <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5_matched_offset5bp.txt"

filt_matched_SVs <- read.table(FILTERED_MATCHED, header = TRUE, colClasses = c(CHROM = 'character',
                                                                    POS = 'numeric',
                                                                    ID = 'character',
                                                                    CAND_ID = 'character',
                                                                    CAND_SVTYPE = 'character',
                                                                    CAND_SVLEN = 'numeric',
                                                                    CAND_SUPP_VEC = 'character'))


filt_missing_GT <- read.table(FILTERED_TABLE, col.names = c('CHROM', 'POS', 'END', 'ID', 'F_MISS'))

matched_filt_SVs_miss <- 
  merge(x = filt_matched_SVs, y = filt_missing_GT,
        by = c('CHROM', 'POS', 'ID'))

# Add explicit platform info
matched_filt_SVs_miss$platform <- sapply(X = as.character(matched_filt_SVs_miss$CAND_SUPP_VEC), 
                                    FUN = function(x){
                                      switch(x,
                                             '10' = 'LR', 
                                             '01' = 'SR',
                                             '11' = 'SR + LR')
                                    }
)

# Plot
ggplot(data = matched_filt_SVs_miss) + 
  geom_violin(aes(x = platform, y = F_MISS, fill = platform)) + ylim(0, 1) +
  theme_bw() +
  guides(fill = 'none') +
  labs(x = 'Platform support', y = 'Missing genotype proportion') +
  stat_summary(aes(x = platform, y = F_MISS), fun.y = "median", geom = "point") +
  geom_text(aes(y = 1, x = 'SR + LR'), label = 'B', hjust = -1.3, vjust = 1, size = 10)

# Save to external file
ggsave(filename = paste0(unlist(strsplit(FILTERED_MATCHED, split = '.txt'))[1], 
                         '_FMISS_by_platform.png'),
       width = 2600,
       height = 2000,
       units = 'px',
       dpi = 600
)

# By type
ggplot(data = matched_filt_SVs_miss) + 
  facet_wrap(~CAND_SVTYPE) + 
  geom_violin(aes(x = platform, y = F_MISS, fill = platform)) +
  theme_bw() +
  guides(fill = 'none') +
  labs(x = 'Platform support', y = 'Missing genotype proportion') +
  stat_summary(aes(x = platform, y = F_MISS), fun.y = "median", geom = "point")

ggsave(filename = paste0(unlist(strsplit(FILTERED_MATCHED, split = '.txt'))[1], 
                         '_FMISS_by_platform_type.png'),
       width = 2600,
       height = 2000,
       units = 'px',
       dpi = 600
)
