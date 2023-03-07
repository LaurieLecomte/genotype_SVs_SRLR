# Explore relationships between allelic depth, genotype quality and genotype in genotyped SVs

argv <- commandArgs(T)
CAND_VCF <- argv[1]   # genotyped SVs VCF, for a single sample    
GENO_VCF <- argv[2]   # candidates SVs VCF, for a single sample


# Import ---------------------------------------------------------------
calls <- read.table(GENO_VCF,
                    header = TRUE, 
                    col.names = c('CHROM', 'POS', 'ID', 'DP', 'GT', 'GT_DP', 'GQ', 'GL{1}', 'GL{2}', 'GL{3}'))

all(calls$DP == calls$GT_DP)
## INFO/DP is equal to FORMAT/DP for a given sample

cand_VCF <- as.data.table(CAND_VCF[, c(1:8,25)])
colnames(cand_VCF) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'sample')


# merge
merged_geno_cand <- merge(x = calls, cand_VCF, by = c('CHROM', 'POS'),
                          all = TRUE)
merged_geno_cand$SVTYPE <- sub(".*SVTYPE=([A-Z]+);.*", "\\1", merged_geno_cand$INFO)
merged_geno_cand$SVLEN  <- as.integer(sub(".*SVLEN=(-?[0-9]+).*", "\\1", merged_geno_cand$INFO))
merged_geno_cand$SUPP_VEC <- sub(".*;SUPP_VEC=([0-1]+);.*", "\\1", merged_geno_cand$INFO)
merged_geno_cand$SUPP <- as.integer(sub(".*;SUPP=([0-9]+).*", "\\1", merged_geno_cand$INFO))


# import graph coverage
graph_depth <- read.table("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/07_calls/13066A.depth.bin200.txt", 
                          col.names = c('CHROM', 'start_bin', 'end_bin', 'mean_depth', 'stdev'))
graph_depth <- graph_depth[!grepl(pattern = 'CAK', graph_depth$CHROM), ]

# Plot DP distribution -------------------------------------------------
DP_distrib <- as.data.frame(table(calls$GT_DP))
colnames(DP_distrib) <- c('DP', 'freq')

# a LOT of SVs have no coverage
DP_distrib$freq[DP_distrib$DP == 0] / sum(DP_distrib$freq) # about 0.18

ggplot(data = calls) +
  geom_bar(aes(x = DP)) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      size = 9,
      hjust = 1
    ),
    plot.title = element_text(size = 11)
  ) 

# If we remove calls with no cov or DP > 100 : 
calls_DP <- subset(calls, DP > 0 & DP < 100)

ggplot(data = calls_DP) +
  geom_bar(aes(x = DP)) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      size = 9,
      hjust = 1
    ),
    plot.title = element_text(size = 11)
  ) 


# Explore GQ --------------------------------------------------------------

# How is GQ correlated with DP ?
ggplot(data = calls_DP) +
  geom_point(aes(x = DP, y = GQ)) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      size = 9,
      hjust = 1
    ),
    plot.title = element_text(size = 11)
  ) 


# What is going on with genotyped calls with DP = 0 ?
calls_DP0 <- subset(merged_geno_cand, DP == 0)

## are these SR or LR SVs ?
table(calls_DP0$SUPP_VEC) # mostly LR calls, so it's likely that no short reads could be mapped to these regions
table(merged_geno_cand$SUPP_VEC[merged_geno_cand$DP == 0])


## are they neighbors ?
### try on first chr 
calls_ssa01q <- subset(merged_geno_cand, CHROM == 'OV354429.1' & DP < 100)

ggplot(data = calls_ssa01q) + 
  facet_wrap(~SUPP_VEC) +
  geom_point(aes(x = POS, y = DP, color = SUPP_VEC), size = 0.08)

ggplot(data = subset(calls_ssa01q, DP < 1)) + 
  facet_wrap(~SUPP_VEC) +
  geom_point(aes(x = POS, y = DP, color = SUPP_VEC), size = 0.08)

## check graph depth
pl1 <- 
  ggplot(data = subset(graph_depth, CHROM == 'OV354429.1' & mean_depth < 500)) + 
  geom_point(aes(x = end_bin, y = mean_depth), size = 0.03, alpha = 0.5) 

pl1 + 
  geom_point(data = subset(graph_depth, CHROM == 'OV354429.1' & mean_depth < 3), 
             aes(x = end_bin, y = mean_depth), col = 'red', size = 0.05) +
  scale_x_continuous(label = function(x) x/1000) +
  labs(x = expression(paste("POS (x10"^3, ")"))) +
  
  geom_point(data = subset(calls_ssa01q, DP == 0), aes(x = POS, y = DP, col = SUPP_VEC), pch = 19, size = 0.05)

# how much of candidates LR SVs have SR support ? about 73 %
DP1_LR_SVs <- subset(merged_geno_cand, SUPP_VEC != '01' & DP > 0)

nrow(DP1_LR_SVs)/nrow(subset(merged_geno_cand, SUPP_VEC != '01'))


# Test effect of DP filter ------------------------------------------------

## DP > 0
geno_DP1 <- subset(merged_geno_cand, DP > 0)
nrow(geno_DP1)

geno_DP2 <- subset(merged_geno_cand, DP > 1)
nrow(geno_DP2)

geno_DP4 <- subset(merged_geno_cand, DP > 3)
nrow(geno_DP4)

# Check patterns asso with GQ = 0
geno_GQ0 <- subset(merged_geno_cand, GQ == 0)
table(geno_GQ0$SUPP_VEC)
table(geno_GQ0$GT) # most are missing GT ./.

## do GQ0 SVs include ALL DP0 SVs ? : YES
geno_DP0 <- subset(merged_geno_cand, DP == 0)
all(geno_DP0$ID.x %in% geno_GQ0$ID.x)
## so what about SVs with GQ0 that do NOT have DP = 0 ? 6838 SVs
geno_GQ0_wCOV <- subset(geno_GQ0, DP !=0)
## many were not among candidate SVs : 
length(is.na(geno_GQ0_wCOV$INFO))

# Test effect of GQ filter
## GQ > 0
geno_GQ1 <- subset(merged_geno_cand, GQ > 0)
nrow(geno_GQ1)
table(geno_GQ1$SUPP_VEC)
table(geno_GQ1$GT)

geno_GQ2 <- subset(merged_geno_cand, GQ > 1)
nrow(geno_GQ2)

geno_GQ4 <- subset(merged_geno_cand, GQ > 3)
nrow(geno_GQ4)
table(geno_GQ4$SUPP_VEC)
table(geno_GQ4$GT)



# Weird genotypes ---------------------------------------------------------
allowed_GTs <- c('0/0', '0/1', '1/0', '1/1', './.', NA)
weird_GTs <- subset(merged_geno_cand, ! GT %in% allowed_GTs)

nrow(weird_GTs)

table(weird_GTs$CHROM)
table(weird_GTs$SUPP_VEC) # nothing special
table(weird_GTs$DP) #nothing special

length(which(is.na(weird_GTs$INFO))) # many were not candidates


# Take home message -------------------------------------------------------
# GQ takes DP into account, i.e. all SVs with DP = 0 also have GQ = 0
# Many SVs with GQ = 0 but DP > 0 were not candidates (at least, not merged with known candidate)
# Weird genotypes have no special patterns, but many were not among candidates
# Filter on GQ, then on allowed genotypes


