# Compare genotyped SVs with original candidate SVs to determine type, since vg does not output SVTYPE field

library(data.table)
library(ggplot2)


# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
CANDIDATES <- argv[1] 
GENOTYPED <- argv[2]
GENOTYPED_FILT <- argv[3]
OFFSET <- argv[4]
OUT_PREFIX <- argv[5]


# 1.1 Candidate SVs -------------------------------------------------------
cand <- fread(CANDIDATES, header = TRUE, colClasses = 'character')
colnames(cand) <- c(sapply(X = strsplit(x = colnames(cand), split = ']'), FUN ="[", 2)[1:10],
                    sub(".*][0-9]_([0-9A-za-z]+_[a-z]+:GT).*", "\\1", colnames(cand)[11:length(colnames(cand))]))
## change END_corr to END
colnames(cand)[which(colnames(cand) == 'END_corr')] <- 'END'

cand <- as.data.frame(cand)
cand[, c('POS', 'SVLEN', 'END', 'SUPP')] <- sapply(cand[, c('POS', 'SVLEN', 'END', 'SUPP')],
                                                     as.numeric)

## remove weird extra column (likely a result from the tabs between samples)
cand <- cand[, 1:(ncol(cand)-1)]

# Convert candidate SVs length to num and bins
SVLEN_breaks <- c(-Inf, 50, 100, 250, 500, 1000, 2500, 5000, 10000, Inf)
SVLEN_names <- c('[0-50[',
                 '[50-100[',
                 '[100-250[',
                 '[250-500[',
                 '[500-1000[',
                 '[1000-2500[',
                 '[2500-5000[',
                 '[5000-10000[',
                 '[10000+')

cand$SVLEN_bin <-
  cut(abs(cand$SVLEN), breaks = SVLEN_breaks, labels = SVLEN_names, right = FALSE)

# Add explicit platform info for candidates
cand$platform <- sapply(X = as.character(cand$SUPP_VEC), 
                        FUN = function(x){
                          switch(x,
                                 '10' = 'LR', 
                                 '01' = 'SR',
                                 '11' = 'LR + SR')
                        }
)


#sapply(X = strsplit(x = colnames(cand), split = '1_'), FUN ="[", 2)[9:length(colnames(cand))]


# 1.2 Raw genotyped SVs ---------------------------------------------------
geno <- fread(GENOTYPED, header = TRUE)
colnames(geno) <- sapply(X = strsplit(x = colnames(geno), split = ']'), FUN ="[", 2)
#geno <- geno[, 1:(ncol(geno)-1)]

# Deduct genotyped SV type from length of REF and ALT alleles
geno$REF_len <- nchar(geno$REF)
geno$ALT_len <- nchar(geno$ALT)

geno$LEN <- 
  ifelse(test = geno$REF_len == 1,# means we have an INSs
         yes = geno$ALT_len,
         no = ifelse(test = geno$ALT_len == 1, # means we have a DEL
                     yes = (0 - geno$REF_len),
                     no = ifelse(test = geno$REF_len == geno$ALT_len,
                                 yes = geno$REF_len,
                                 no = NA)))
geno$absLEN <- 
  ifelse(test = geno$REF_len == 1,
         yes = geno$ALT_len,
         no = ifelse(test = geno$ALT_len == 1,
                     yes = geno$REF_len,
                     no = ifelse(test = geno$REF_len == geno$ALT_len,
                                 yes = geno$REF_len,
                                 no = NA)))

geno$VAR_TYPE <- 
  ifelse(test = abs(geno$LEN) >= 50,
         yes = 'SV',
         no = 'indel')

geno$TYPE <- 
  ifelse(test = geno$REF_len == 1 & geno$VAR_TYPE == 'SV',
         yes = 'INS or DUP',
         no = ifelse(test = geno$ALT_len == 1 & geno$VAR_TYPE == 'SV',
                     yes = 'DEL',
                     no = ifelse(test = geno$REF_len == geno$ALT_len & geno$VAR_TYPE == 'SV',
                                 yes = 'INV', 
                                 no = 'other')))
## others include INSs or DELs that have a REF or ALT allele larger than one, 
## or DUPs that have variable length relative to ref allele, 
## or INVs that REF and ALT seq are of sightly different length



# 1.3 Filtered genotyped SVs ----------------------------------------------
# Import filtered genotypes
geno_filt <- fread(GENOTYPED_FILT, header = TRUE)
colnames(geno_filt) <- sapply(X = strsplit(x = colnames(geno_filt), split = ']'), FUN ="[", 2)
#geno_filt <- geno_filt[, 1:(ncol(geno_filt)-1)]

# Deduct filtered genotyped SV type from length of REF and ALT alleles
geno_filt$REF_len <- nchar(geno_filt$REF)
geno_filt$ALT_len <- nchar(geno_filt$ALT)

geno_filt$LEN <- 
  ifelse(test = geno_filt$REF_len == 1,
         yes = geno_filt$ALT_len,
         no = ifelse(test = geno_filt$ALT_len == 1,
                     yes = (0 - geno_filt$REF_len),
                     no = ifelse(test = geno_filt$REF_len == geno_filt$ALT_len,
                                 yes = geno_filt$REF_len,
                                 no = NA)))
geno_filt$absLEN <- 
  ifelse(test = geno_filt$REF_len == 1,
         yes = geno_filt$ALT_len,
         no = ifelse(test = geno_filt$ALT_len == 1,
                     yes = (0 - geno_filt$REF_len),
                     no = ifelse(test = geno_filt$REF_len == geno_filt$ALT_len,
                                 yes = geno_filt$REF_len,
                                 no = NA)))

geno_filt$TYPE <- 
  ifelse(test = geno_filt$REF_len == 1,
         yes = 'INS or DUP',
         no = ifelse(test = geno_filt$ALT_len == 1,
                     yes = 'DEL',
                     no = ifelse(test = geno_filt$REF_len == geno_filt$ALT_len,
                                 yes = 'INV', 
                                 no = 'other')))
## others include INSs or DELs that have a REF or ALT allele larger than one, 
## or DUPs that have variable length relative to ref allele, 
## or INVs that REF and ALT seq are of sightly different length





# 2. Merge candidates with genotyped and filtered genotyped SVs -----------

# Merge with exact positions
## without genotypes, to make matters simpler
cand_cols <- which(! grepl(paste(c('GT', 'REF', 'ALT'), collapse = "|"), colnames(cand)))
      
#geno_cols <- which(!grepl('GT', colnames(geno)))

cand_geno <- merge(x = cand[, cand_cols], y = geno[geno$VAR_TYPE == 'SV', c('CHROM', 'POS', 'ID', 'LEN', 'VAR_TYPE', 'TYPE')], 
                   by = c('CHROM', 'POS'), all = TRUE)

# Merge with window around POS and SVLEN
## Add window around position and svlen of cand SVs
cand$POS_min <- cand$POS - OFFSET
cand$POS_max <- cand$POS + OFFSET

cand$absSVLEN_min <- abs(cand$SVLEN) - OFFSET
cand$absSVLEN_max <- abs(cand$SVLEN) + OFFSET

## Merge
cand_geno_offset5 <- 
  geno[cand, 
       on = .(CHROM, 
              POS <= POS_max, POS >= POS_min, # allow a window on POS (VARx <|>|>=|=< VARi)
              absLEN <= absSVLEN_max, absLEN >= absSVLEN_min), # allow window on SVLEN
       .(x.CHROM, x.POS, i.POS, # variables I need in merged output (i = cand) (x = geno)
         i.ID, x.ID, 
         i.SVLEN, x.LEN,
         i.SVTYPE, x.TYPE, x.VAR_TYPE,
         i.platform, i.SUPP, i.SUPP_VEC)] 

## How many candidates not genotypes ?
not_genotyped <- subset(cand_geno_offset5, is.na(x.ID))
nrow(not_genotyped) ## candidates that were unmatched to a genotyped SV

## How many genotyped SVs that were NOT candidates
### these SVs are in in geno, but not in cand_geno_offset5
not_candidates <- subset(geno, ID %in% setdiff(geno$ID, cand_geno_offset5$x.ID))


table(cand_geno_offset5$i.SVTYPE)



# Summarize SVTYPE by platform
table(cand_geno_offset5$x.TYPE)

# Plot by type, size and platform
ggplot(data = cand) +
  geom_bar(aes(x = SVLEN_bin, fill = platform)) +
  facet_wrap(~SVTYPE, scales = 'free_y') +
  theme(
    axis.text.x = element_text(
      angle = 45,
      size = 7,
      hjust = 1
    ),
    plot.title = element_text(size = 8)
  ) +
  labs(
    x = "SV size",
    y = "SV count",
    fill = "Data type",
    title = "SV count by type, size and sequencing data type"
  ) + 
  scale_fill_viridis_d(option = "D")


ggplot(data = cand) +
  geom_bar(aes(x = SVLEN_bin, fill = SVTYPE)) +
  facet_wrap(~platform, scales = 'free_y') +
  theme(
    axis.text.x = element_text(
      angle = 45,
      size = 7,
      hjust = 1
    ),
    plot.title = element_text(size = 8)
  ) +
  labs(
    x = "SV size",
    y = "SV count",
    fill = "SV type",
    title = "SV count by type, size and sequencing data type"
  ) + 
  scale_fill_viridis_d(option = "D")


# Add window around position

cand$POS_min <- cand$POS - OFFSET
cand$POS_max <- cand$POS + OFFSET

test_merged_win10 <- 
  geno_VCF[cand_VCF, on = .(CHROM, POS <= POS_max, POS >= POS_min), # variables on which merging is done and conditions
           .(x.CHROM, x.POS, x.INFO, i.ID, i.POS, i.INFO)] # variables I need in merged output (i = cand_VCF) (x = geno_VCF)


# Explore RAW genotyped SVs -----------------------------------------------


# merge with candidates, allowing position difference
merged_raw_offset <- 
  geno[cand, on = .(CHROM, POS <= POS_max, POS >= POS_min), # variables on which merging is done and conditions
           .(x.CHROM, x.POS, x.ID, i.ID, i.POS, i.SVTYPE, i.SVLEN, i.platform, i.END)] # variables I need in merged output (i = cand_VCF) (x = geno_VCF)

not_genotyped_raw_offset <- subset(merged_raw_offset, is.na(merged_raw_offset$x.ID))

cols_cand <- colnames(cand)[! grepl('GT', colnames(cand))]
cols_geno <- colnames(geno)[! grepl('GT', colnames(geno))]

merged_raw <- merge(x = cand[, ..cols_cand], y = geno[, ..cols_geno], 
                    by = c('CHROM', 'POS'), all = TRUE)
not_genotyped_raw <- subset(merged_raw, is.na(merged_raw$ID.y))
not_candidate_raw <- subset(merged_raw, is.na(merged_raw$ID.x))

# compare merged with offset with merge without offset
diff_offset_nooffset <- merged_raw_offset[which(merged_raw_offset$x.POS != merged_raw_offset$i.POS),]

## in the regular merged set, these IDs are not merged with anything, so appear as not genotyped
diff_offset_nooffset_IDs <- merged_raw[which(merged_raw$ID.y %in% diff_offset_nooffset$x.ID), c('CHROM', 'POS', 'ID.x', 'SVTYPE', 'SVLEN', 'platform')]

# Explore FILTERED genotyped SVs ------------------------------------------



