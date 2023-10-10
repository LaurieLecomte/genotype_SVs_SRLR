# Check if difference in missing GT between SR and LR

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


matched_SVs_miss <- 
merge(x = matched_SVs, y = missing_GT,
      by = c('CHROM', 'POS', 'ID'))


ggplot(data = matched_SVs_miss) + geom_boxplot(aes(x = CAND_SUPP_VEC, y = F_MISS))


# split SR and LR

sapply(X = matched_SVs_miss$CAND_SUPP_VEC, FUN = switch,
       '10' = 'LR', 
       '01' = 'SR', 
       )

for (i in 1:nrow(matched_SVs_miss)){
  matched_SVs_miss$SR[i] <- 
  switch(matched_SVs_miss$CAND_SUPP_VEC[i],
         '01' = 1,
         '11' = 1,
         '10' = 0)
}
for (i in 1:nrow(matched_SVs_miss)){
  matched_SVs_miss$LR[i] <- 
    switch(matched_SVs_miss$CAND_SUPP_VEC[i],
           '01' = 0,
           '11' = 1,
           '10' = 1)
}

ggplot(data = matched_SVs_miss) + geom_boxplot(aes(x = factor(SR), y = F_MISS))

# By type
ggplot(data = matched_SVs_miss) + geom_boxplot(aes(x = factor(SR), y = F_MISS))



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

ggplot(data = matched_filt_SVs_miss) + geom_boxplot(aes(x = CAND_SUPP_VEC, y = F_MISS)) + ylim(0, 1)

# By type
ggplot(data = matched_filt_SVs_miss) + facet_wrap(~CAND_SVTYPE) + geom_boxplot(aes(x = CAND_SUPP_VEC, y = F_MISS))