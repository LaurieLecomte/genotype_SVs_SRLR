setwd("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR")

library(data.table)

# import VCF
geno_VCF <- read.table('07_calls/raw/13070A.vcf', comment.char = "#", stringsAsFactors = FALSE,
                       )

geno_VCF <- as.data.table(geno_VCF[, 1:8])
colnames(geno_VCF) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')

cand_VCF <- read.table('05_candidates/candidates/merged_SUPP2.candidates.vcf.gz', 
                       comment.char = "#", stringsAsFactors = FALSE)


cand_VCF <- as.data.table(cand_VCF[, 1:8])
colnames(cand_VCF) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')


# merge
test_merge_V1_V2 <- merge(x = geno_VCF, cand_VCF, by = c('CHROM', 'POS'),
                          all = TRUE)

test_merge_V1_V2$SVTYPE <- sub(".*SVTYPE=([A-Z]+);.*", "\\1", test_merge_V1_V2$INFO.y)
test_merge_V1_V2$SVLEN  <- as.integer(sub(".*SVLEN=(-?[0-9]+).*", "\\1", test_merge_V1_V2$INFO.y))
test_merge_V1_V2$SUPP_VEC <- sub(".*;SUPP_VEC=([0-1]+);.*", "\\1", test_merge_V1_V2$INFO.y)
test_merge_V1_V2$SUPP <- as.integer(sub(".*;SUPP=([0-9]+).*", "\\1", test_merge_V1_V2$INFO.y))

not_genotyped <-test_merge_V1_V2[which(is.na(test_merge_V1_V2$ID.x)), ]
not_candidate <-test_merge_V1_V2[which(is.na(test_merge_V1_V2$ID.y)), ]
  
  
# merge with a window
WIN=10

cand_VCF$POS_max <- cand_VCF$POS + WIN
cand_VCF$POS_min <- cand_VCF$POS - WIN


test_merged_win10 <- 
  geno_VCF[cand_VCF, on = .(CHROM, POS <= POS_max, POS >= POS_min), # variables on which merging is done and conditions
    .(x.CHROM, x.POS, x.INFO, i.ID, i.POS, i.INFO)] # variables I need in merged output (i = cand_VCF) (x = geno_VCF)

test_merged_win10$SVTYPE <- sub(".*SVTYPE=([A-Z]+);.*", "\\1", test_merged_win10$i.INFO)
test_merged_win10$SVLEN  <- as.integer(sub(".*SVLEN=(-?[0-9]+).*", "\\1", test_merged_win10$i.INFO))
test_merged_win10$SUPP_VEC <- sub(".*;SUPP_VEC=([0-1]+);.*", "\\1", test_merged_win10$i.INFO)
test_merged_win10$SUPP <- as.integer(sub(".*;SUPP=([0-9]+).*", "\\1", test_merged_win10$i.INFO))

not_genotyped_win10 <- test_merged_win10[which(is.na(test_merged_win10$x.POS)),]

WIN=50

cand_VCF$POS_max <- cand_VCF$POS + WIN
cand_VCF$POS_min <- cand_VCF$POS - WIN


test_merged_win50 <- 
  geno_VCF[cand_VCF, on = .(CHROM, POS <= POS_max, POS >= POS_min), # variables on which merging is done and conditions
           .(x.CHROM, x.POS, x.INFO, i.ID, i.POS, i.INFO)] # variables I need in merged output (i = cand_VCF) (x = geno_VCF)

test_merged_win50$SVTYPE <- sub(".*SVTYPE=([A-Z]+);.*", "\\1", test_merged_win50$i.INFO)
test_merged_win50$SVLEN  <- as.integer(sub(".*SVLEN=(-?[0-9]+).*", "\\1", test_merged_win50$i.INFO))
test_merged_win50$SUPP_VEC <- sub(".*;SUPP_VEC=([0-1]+);.*", "\\1", test_merged_win50$i.INFO)
test_merged_win50$SUPP <- as.integer(sub(".*;SUPP=([0-9]+).*", "\\1", test_merged_win50$i.INFO))

not_genotyped_win50 <- test_merged_win50[which(is.na(test_merged_win50$x.POS)),]

# less duplicated i.IDs with WIN = 10 than 50
#

WIN=5

cand_VCF$POS_max <- cand_VCF$POS + WIN
cand_VCF$POS_min <- cand_VCF$POS - WIN


test_merged_win5 <- 
  geno_VCF[cand_VCF, on = .(CHROM, POS <= POS_max, POS >= POS_min), # variables on which merging is done and conditions
           .(x.CHROM, x.POS, x.INFO, i.ID, i.POS, i.INFO)] # variables I need in merged output (i = cand_VCF) (x = geno_VCF)

test_merged_win5$SVTYPE <- sub(".*SVTYPE=([A-Z]+);.*", "\\1", test_merged_win5$i.INFO)
test_merged_win5$SVLEN  <- as.integer(sub(".*SVLEN=(-?[0-9]+).*", "\\1", test_merged_win5$i.INFO))
test_merged_win5$SUPP_VEC <- sub(".*;SUPP_VEC=([0-1]+);.*", "\\1", test_merged_win5$i.INFO)
test_merged_win5$SUPP <- as.integer(sub(".*;SUPP=([0-9]+).*", "\\1", test_merged_win5$i.INFO))

not_genotyped_win5 <- test_merged_win5[which(is.na(test_merged_win5$x.POS)),]


WIN=2

cand_VCF$POS_max <- cand_VCF$POS + WIN
cand_VCF$POS_min <- cand_VCF$POS - WIN


test_merged_win2 <- 
  geno_VCF[cand_VCF, on = .(CHROM, POS <= POS_max, POS >= POS_min), # variables on which merging is done and conditions
           .(x.CHROM, x.POS, x.INFO, i.ID, i.POS, i.INFO)] # variables I need in merged output (i = cand_VCF) (x = geno_VCF)

test_merged_win2$SVTYPE <- sub(".*SVTYPE=([A-Z]+);.*", "\\1", test_merged_win2$i.INFO)
test_merged_win2$SVLEN  <- as.integer(sub(".*SVLEN=(-?[0-9]+).*", "\\1", test_merged_win2$i.INFO))
test_merged_win2$SUPP_VEC <- sub(".*;SUPP_VEC=([0-1]+);.*", "\\1", test_merged_win2$i.INFO)
test_merged_win2$SUPP <- as.integer(sub(".*;SUPP=([0-9]+).*", "\\1", test_merged_win2$i.INFO))

not_genotyped_win2 <- test_merged_win2[which(is.na(test_merged_win2$x.POS)),]


table(test_merge_V1_V2$SVTYPE)
table(test_merge_V1_V2$SVTYPE[test_merge_V1_V2$SUPP_VEC=='10'])
table(test_merge_V1_V2$SVTYPE[test_merge_V1_V2$SUPP_VEC=='01'])
table(test_merge_V1_V2$SVTYPE[test_merge_V1_V2$SUPP_VEC=='11'])
