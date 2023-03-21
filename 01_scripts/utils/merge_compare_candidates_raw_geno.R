# Compare genotyped SVs with original candidate SVs to determine type, since vg does not output SVTYPE field

library(data.table)
library(ggplot2)
library(dplyr)


# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
CANDIDATES <- argv[1] 
GENOTYPED <- argv[2]
OFFSET <- argv[3]


# 1.1 Candidate SVs -------------------------------------------------------
# Import
cand <- fread(CANDIDATES, header = TRUE, 
              colClasses = 'character') # required to keep SUPP_VEC intact
colnames(cand) <- c(sapply(X = strsplit(x = colnames(cand), split = ']'), FUN ="[", 2)[1:10],
                    sub(".*][0-9]_([0-9A-za-z]+_[a-z]+:GT).*", "\\1", colnames(cand)[11:length(colnames(cand))]))

## Rename column END_corr to END
colnames(cand)[which(colnames(cand) == 'END_corr')] <- 'END'

## Convert back to data.frame for easier handling until merging
cand <- as.data.frame(cand)
cand[, c('POS', 'SVLEN', 'END', 'SUPP')] <- sapply(cand[, c('POS', 'SVLEN', 'END', 'SUPP')],
                                                     as.numeric)

## Remove weird extra column (likely a result from the tabs between samples)
cand <- cand[, 1:(ncol(cand)-1)]

# Compute info on each candidate
## Convert candidate SVs length to num and bins
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

## Add platform
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
# Import
geno <- fread(GENOTYPED, header = TRUE)
colnames(geno) <- sapply(X = strsplit(x = colnames(geno), split = ']'), FUN ="[", 2)
geno <- geno[, 1:(ncol(geno)-1)]

## Reconvert to data frame until merging
geno <- as.data.frame(geno)

# Compute info on each raw call
## Compute number of ALT alleles
nb_ALT <- unname(
  sapply(X = geno$ALT, 
         FUN = function(x){
           # multiple alleles are separated by commas
           (lengths(regmatches(x, gregexpr(",", x)))) + 1 # +1 because 0 commas = 1 allele only
           }
         )
  )

geno$num_alt_alleles <- nb_ALT

## Confirm that there is only 1 REF allele per site
nb_REF_alleles <- unname(
  sapply(X = geno$REF, 
         FUN = function(x){
           (lengths(regmatches(x, gregexpr(",", x)))) + 1
           }
         )
  )

table(nb_REF_alleles) ## should be 1 for all

# Deduct genotyped SV type from length of REF and ALT alleles
## Compute REF and ALT length from sequences
geno$REF_len <- nchar(geno$REF)
geno$ALT_len <- ifelse(geno$num_alt_alleles == 1,
                       yes = nchar(geno$ALT),
                       no = NA) ## we cannot compute ALT allele length if multiple ALT alleles

## Deduct SVTYPE from ALT vs REF alleles lengths
geno$LEN <- 
  ifelse(test = geno$REF_len == 1 & geno$num_alt_alleles == 1, 
         yes = geno$ALT_len, # means we have an INSs (or DUP), so SVLEN is ALT allele length
         no = ifelse(test = geno$ALT_len == 1 & geno$num_alt_alleles == 1, 
                     yes = (0 - geno$REF_len), # means we have a DEL, so SVLEN is REF allele length
                     no = ifelse(test = geno$REF_len == geno$ALT_len & geno$num_alt_alleles == 1, 
                                 yes = geno$REF_len, # possible INV
                                 no = ifelse(test = geno$num_alt_alleles == 1,
                                             yes = geno$ALT_len, # if REF and ALT are >1 bp and different length, we assign ALT_len
                                             no = NA)) # no length for non-biallelic sites
                                 
                     )
         )

geno$absLEN <- abs(geno$LEN)

## Assign variant type based on length (SV, indels or other SV, nothing if non-biallelic)
geno$VAR_TYPE <- 
  ifelse(test = geno$absLEN >= 50 & geno$num_alt_alleles == 1,
         yes = 'SV',
         no = ifelse(test = geno$absLEN < 50 & geno$num_alt_alleles == 1,
                     yes = 'indel',
                     no = ifelse(test = geno$num_alt_alleles == 1,
                                 yes = 'other',
                                 no = 'multiallelic')
                     ) # for multiallelic calls for which we cannot infer length
         )

## Assign a possible SV type for calls confidently labelled as SV
geno$TYPE <- 
  ifelse(test = geno$REF_len == 1 & geno$VAR_TYPE == 'SV', # SV and indel labels already take number of alt alleles into account
         yes = 'INS/DUP',
         no = ifelse(test = geno$ALT_len == 1 & geno$VAR_TYPE == 'SV',
                     yes = 'DEL',
                     no = ifelse(test = geno$REF_len == geno$ALT_len & geno$VAR_TYPE == 'SV',
                                 yes = 'INV', 
                                 no = ifelse(test = geno$VAR_TYPE == 'SV',
                                             yes = 'any SV type', 
                                             no = ifelse(test = geno$VAR_TYPE == 'indel',
                                                         yes = 'indel', 
                                                         no = 'multiallelic')
                                             )
                                 )
                     )
  )
         
         

## Compute mean DP and GQ per genotyped SV
geno_DP <- as.data.frame(select(geno, contains("DP")))
geno_DP <- replace(geno_DP, geno_DP == '.', NA)
geno_DP[, 1:ncol(geno_DP)] <- sapply(geno_DP[, 1:ncol(geno_DP)], as.numeric)

geno$mean_DP <- rowMeans(geno_DP, na.rm = TRUE)



# 2. Explore raw genotyped SVs --------------------------------------------
# How many raw calls indels, not SVs ?
table(geno$VAR_TYPE)
nrow(subset(geno, VAR_TYPE == 'indel'))

# How many are not biallelic calls ?
multiallelic <- (subset(geno, num_alt_alleles != 1))
nrow(multiallelic)


# How many are biallelic SVs ?
nrow(subset(geno, VAR_TYPE == 'SV' & num_alt_alleles == 1))
# How many are biallelic indels ?
indels <- subset(geno, VAR_TYPE == 'indel' & num_alt_alleles == 1)
nrow(indels)


# 2. Merge candidates with genotyped and filtered genotyped SVs -----------

# 2.1 Exact POS match -----------------------------------------------------
# Without genotypes, to make matters simpler
cand_cols <- which(! grepl(paste(c('GT', 'REF', 'ALT'), collapse = "|"), colnames(cand)))
      
geno_cols <- which(! grepl(paste(c(':GT', ':DP', ':GQ', 'ALT', 'REF'), collapse = "|"), colnames(geno)))

cand_geno <- merge(x = cand[, cand_cols], y = subset(geno[, geno_cols]), 
                   by = c('CHROM', 'POS'), all = TRUE)

# length(unique(cand_geno$ID.x[!is.na(cand_geno$ID.y)])) / length(unique(cand$ID))

# How many candidates not genotyped ?
nrow(subset(cand_geno, is.na(ID.y))) # candidate  length(unique(cand_geno$ID.y[is.na(cand_geno$ID.y)]))
length(unique(cand_geno$ID.x[is.na(cand_geno$ID.y)]))

# How many candidates genotyped ?
nrow(subset(cand_geno, !is.na(ID.y)))
length(unique(cand_geno$ID.x[!is.na(cand_geno$ID.y)]))
length(unique(cand_geno$ID.x[!is.na(cand_geno$ID.y)])) / length(unique(cand$ID))



# 2.2 With window around POS and SVLEN ------------------------------------
# Add window around position and svlen of cand SVs
cand$POS_min <- cand$POS - OFFSET
cand$POS_max <- cand$POS + OFFSET
cand$absSVLEN_min <- abs(cand$SVLEN) - OFFSET
cand$absSVLEN_max <- abs(cand$SVLEN) + OFFSET

# Merge
## Reconvert to data table
geno <- as.data.table(geno)
cand <- as.data.table(cand)

## Merge with offset
cand_geno_offset <- 
  geno[cand, 
       on = .(CHROM, 
              POS <= POS_max, POS >= POS_min, # allow a window on POS (VARx <|>|>=|=< VARi)
              absLEN <= absSVLEN_max, absLEN >= absSVLEN_min), # allow window on SVLEN
       .(x.CHROM, i.POS, x.POS, # variables I need in merged output (i = cand) (x = geno)
         i.ID, x.ID, 
         i.SVLEN_bin, i.SVLEN, x.LEN,
         i.SVTYPE, x.TYPE, x.VAR_TYPE,
         i.platform, i.SUPP, i.SUPP_VEC, x.num_alt_alleles)] 


# 3. Explore successfully MATCHED calls -----------------------------------
# How many assigned ?
cand_geno_offset_genotyped <- subset(cand_geno_offset, ! is.na(x.ID)) # vg ID is NOT NA if matched with a known candidate
length(unique(cand_geno_offset_genotyped$i.ID)) # unique candidates IDs assigned to a vg call
length(unique(cand_geno_offset_genotyped$x.ID)) # unique vg IDs assigned to a candidate SV

# How many assigned by candidate SVTYPE and platform?
table(cand_geno_offset_genotyped$i.SVTYPE)
## Reorder candidates' platform levels
reordered_platform <- c('SR', 'LR', 'LR + SR')
cand_geno_offset_genotyped$i.platform_reordered <- factor(cand_geno_offset_genotyped$i.platform, 
                                    levels = reordered_platform)
table(cand_geno_offset_genotyped$i.platform_reordered, cand_geno_offset_genotyped$i.SVTYPE)

# Plot assigned raw genotyped calls by platform and SVTYPE
ggplot(data = cand_geno_offset_genotyped) +
  facet_wrap(~i.SVTYPE, scales = 'free_y') +
  geom_bar(aes(x = i.SVLEN_bin, fill = i.platform_reordered)) + 
  theme(
    ## Plot title
    plot.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    ## Axis
    axis.text.x = element_text(angle = 45, size = 6, hjust = 1),
    axis.text.y = element_text(size = 6, hjust = 1),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    ## Legend
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 7),
    legend.key.size = unit(5, 'mm')
  ) +
  labs(
    x = "SV size (bp)",
    y = "SV count",
    fill = "Data type",
    title = "Raw genotyped SV count by size bins and platform"
  ) + 
  scale_fill_viridis_d(option = "B")


# 4. Explore UNMATCHED vg calls -------------------------------------------
# Which genotyped calls were not assigned ?
vg_not_assigned <- geno[geno$ID %in% setdiff(geno$ID, cand_geno_offset$x.ID), ]

# How many unique candidates not assigned ?
length(unique(cand_geno_offset$i.ID[is.na(cand_geno_offset$x.ID)])) # unique candidates not genotyped/assigned

# How many raw genotyped calls not assigned to a known candidate
length(unique(vg_not_assigned$ID))

# How many unmatched vg calls by putative variant type ?
table(vg_not_assigned$VAR_TYPE)
nrow(subset(vg_not_assigned, num_alt_alleles != 1)) ## most unmatched genotype calls are multiallelic

vg_not_assigned_biall <- (subset(vg_not_assigned, num_alt_alleles == 1))
table(vg_not_assigned_biall$VAR_TYPE)
table(vg_not_assigned_biall$TYPE)

# How many unassigned are duplicated CHROM-POS-ID entries ?
sum(duplicated(vg_not_assigned_biall[, c('CHROM', 'POS', 'ID')]))

# Convert estimated SVLEN to SVLEN_bins
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

vg_not_assigned_biall$SVLEN_bin_REF <-
  cut(abs(vg_not_assigned_biall$REF_len), breaks = SVLEN_breaks, labels = SVLEN_names, right = FALSE)
vg_not_assigned_biall$SVLEN_bin_ALT <-
  cut(abs(vg_not_assigned_biall$ALT_len), breaks = SVLEN_breaks, labels = SVLEN_names, right = FALSE)

# Plot unassigned raw genotyped calls by sizes
ggplot(data = vg_not_assigned_biall) +
  geom_bar(aes(x = SVLEN_bin_REF)) + 
  theme(
    ## Plot title
    plot.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    ## Axis
    axis.text.x = element_text(angle = 45, size = 6, hjust = 1),
    axis.text.y = element_text(size = 6, hjust = 1),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    ## Legend
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 7),
    legend.key.size = unit(5, 'mm')
  ) +
  labs(
    x = "SV size (bp)",
    y = "SV count",
    title = "Raw UNMATCHED genotyped SV count by REF allele size"
  ) 

ggplot(data = vg_not_assigned_biall) +
  geom_bar(aes(x = SVLEN_bin_ALT)) + 
  theme(
    ## Plot title
    plot.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    ## Axis
    axis.text.x = element_text(angle = 45, size = 6, hjust = 1),
    axis.text.y = element_text(size = 6, hjust = 1),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    ## Legend
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 7),
    legend.key.size = unit(5, 'mm')
  ) +
  labs(
    x = "SV size (bp)",
    y = "SV count",
    title = "Raw UNMATCHED genotyped SV count by ALT allele size"
  ) 


# 5. Export useful infos --------------------------------------------------
# Export matched calls to a list
write.table(x = cand_geno_offset_genotyped[, c('x.CHROM', 'x.POS', 'x.ID', 'i.ID', 'i.SVTYPE', 'i.SVLEN', 'i.SUPP_VEC')],
            file = paste0(strsplit(GENOTYPED, '.table')[[1]], 'matched_offset', OFFSET, 'bp.txt'),
            col.names = c('CHROM', 'POS', 'ID', 'CAND_ID', 'CAND_SVTYPE', 'CAND_SVLEN', 'CAND_SUPP_VEC'),
            quote = FALSE, row.names = FALSE, sep = "\t")

# Export UNmatched calls to a list
write.table(x = vg_not_assigned_biall[, c('CHROM', 'POS', 'ID')],
            file = paste0(strsplit(GENOTYPED, '.table')[[1]], 'unmatched_offset', OFFSET, 'bp.txt'),
            col.names = c('CHROM', 'POS', 'ID'),
            quote = FALSE, row.names = FALSE, sep = "\t")

