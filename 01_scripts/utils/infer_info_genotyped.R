# Compare filtered genotyped SVs with original candidate SVs to determine type, since vg does not output SVTYPE field

library(data.table)
library(ggplot2)
library(dplyr)
library(scales)


# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
CANDIDATES <- argv[1] 
#CANDIDATES <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/05_candidates/candidates/merged_SUPP2.candidates.table'
GENOTYPED_FILT <- argv[2]
#GENOTYPED_FILT <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5.table'
OFFSET <- as.numeric(argv[3])
#OFFSET <- 5

# 1.1 Candidate SVs -------------------------------------------------------
# Import
cand <- fread(CANDIDATES, header = TRUE, 
              colClasses = 'character') # required to keep SUPP_VEC intact

# Format column names
colnames(cand) <- c(sapply(X = strsplit(x = colnames(cand), split = ']'), FUN ="[", 2)[1:10],
                    sub(".*][0-9]_([0-9A-za-z]+_[a-z]+:GT).*", "\\1", colnames(cand)[11:length(colnames(cand))]))

# Rename column END_corr to END
colnames(cand)[which(colnames(cand) == 'END_corr')] <- 'END'


# Convert back to data.frame for easier handling until merging
cand <- as.data.frame(cand)

# Remove GQ and DP columns
cand <- select(cand, !grep(pattern = ':', colnames(cand)))

# Convert some cols to numeric
cand[, c('POS', 'SVLEN', 'END', 'SUPP')] <- sapply(cand[, c('POS', 'SVLEN', 'END', 'SUPP')],
                                                   as.numeric)

## Remove weird extra column (likely a result from the tabs between samples)
cand <- cand[, 1:(ncol(cand)-1)]



# 1.2 Filtered genotyped SVs ---------------------------------------------------
# Import
geno <- fread(GENOTYPED_FILT, header = TRUE)

# Format column names
colnames(geno) <- sapply(X = strsplit(x = colnames(geno), split = ']'), FUN ="[", 2)

# Reconvert to data frame until merging
geno <- as.data.frame(geno)

# Remove GQ and DP columns
geno <- select(geno, !grep(pattern = ':', colnames(geno)))

# Remove weird extra column (likely a result from the tabs between samples)
#geno <- geno[, 1:(ncol(geno)-1)]




# 2. Add and compute info on SV entries -----------------------------------
# 2.1 Putative SVs
# Compute REF and ALT length from sequences
cand$REF_len <- nchar(cand$REF)
cand$ALT_len <- nchar(cand$ALT)


# Convert candidate SVLEN to num and bins (for plotting)
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

# Add platform
cand$platform <- sapply(X = as.character(cand$SUPP_VEC), 
                        FUN = function(x){
                          switch(x,
                                 '10' = 'LR', 
                                 '01' = 'SR',
                                 '11' = 'LR + SR')
                        }
)

# 2.2 Filtered genotyped SVs
# Confirm that all genotyped SVs are biallelic
nb_ALT <- unname(
  sapply(X = geno$ALT, 
         FUN = function(x){
           # multiple alleles are separated by commas
           (lengths(regmatches(x, gregexpr(",", x)))) + 1 # +1 because 0 commas = 1 allele only
         }
  )
)

table(nb_ALT) # all are biallelic
geno$num_alt_alleles <- nb_ALT

# Confirm that there is only 1 REF allele per site
nb_REF_alleles <- unname(
  sapply(X = geno$REF, 
         FUN = function(x){
           (lengths(regmatches(x, gregexpr(",", x)))) + 1
         }
  )
)

table(nb_REF_alleles) ## should be 1 for all


# Deduct SV/indel type from alleles lengths
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
                     
         ))
  

#geno$LEN <- 
#  ifelse(test = geno$REF_len == 1 & geno$ALT_len > 1, 
#         yes = geno$ALT_len, # means we have an INSs (or DUP), so SVLEN is ALT allele length
#         no = ifelse(test = geno$ALT_len == 1 & geno$REF_len > 1, 
#                     yes = (0 - geno$REF_len), # means we have a DEL, so SVLEN is REF allele length
#                     no = ifelse(test = geno$REF_len == geno$ALT_len, 
#                                 yes = geno$REF_len, # possible INV
#                                 no = NA)) # no length for non-biallelic sites
#                     
#         )
#  

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

table(geno$VAR_TYPE) # we have 700 genotyped calls that are indels with <50 bp


# Assign a possible SV type for calls confidently labelled as SV
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
  )#





# 3. Match successfully genotyped SVs with a known putative SV ------------
# We match based on position, SV size and REF/ALT alleles lengths, 
# allowing a 5 bp window around each value

# Add window around POS, SVLEN, REF and ALT_len of each putative SV
cand$POS_min <- cand$POS - OFFSET
cand$POS_max <- cand$POS + OFFSET
cand$absSVLEN_min <- abs(cand$SVLEN) - OFFSET
cand$absSVLEN_max <- abs(cand$SVLEN) + OFFSET
cand$REF_len_min <- cand$REF_len - OFFSET
cand$REF_len_max <- cand$REF_len + OFFSET
cand$ALT_len_min <- cand$ALT_len - OFFSET
cand$ALT_len_max <- cand$ALT_len + OFFSET

# Reconvert to data table for merging
geno <- as.data.table(geno)
cand <- as.data.table(cand)

# Merge
cand_geno_offset <- 
  geno[cand, 
       on = .(CHROM, 
              POS <= POS_max, POS >= POS_min, # allow a window on POS (VARx <|>|>=|=< VARi)
              #absLEN <= absSVLEN_max, absLEN >= absSVLEN_min, # allow window on SVLEN
              REF_len <= REF_len_max, REF_len >= REF_len_min, # allow window on REF seq length
              ALT_len <= ALT_len_max, ALT_len >= ALT_len_min # allow window on ALT seq length
       ),
       .(x.CHROM, i.POS, x.POS, # variables I need in merged output (i = cand) (x = geno)
         i.ID, x.ID, 
         i.SVLEN_bin, i.SVLEN, x.LEN,
         i.SVTYPE, x.VAR_TYPE, x.TYPE,
         i.platform, i.SUPP, i.SUPP_VEC, x.num_alt_alleles,
         i.REF_len, x.REF_len,
         i.ALT_len, x.ALT_len),]



# 4. Explore successfully MATCHED calls -----------------------------------
# How many successful matches ?
cand_geno_offset_genotyped <- subset(cand_geno_offset, ! is.na(x.ID)) # vg ID is NOT NA if matched with a known candidate
length(unique(cand_geno_offset_genotyped$i.ID)) # unique candidates IDs assigned to a vg call
length(unique(cand_geno_offset_genotyped$x.ID)) # unique vg IDs assigned to a candidate SV


matched_geno_IDs <- geno$ID[which(geno$ID %in% cand_geno_offset$x.ID)]
paste(length(matched_geno_IDs), 'vg IDs matched to a known candidate')

unmatched_geno_IDs <- geno$ID[which(! geno$ID %in% cand_geno_offset$x.ID)]
paste(length(unmatched_geno_IDs), 'vg IDs NOT matched to a known candidate')

#cand_geno_offset_matched <- cand_geno_offset[cand_geno_offset$x.ID %in% matched_geno_IDs, ] # equivalent to cand_geno_offset_genotyped

## Remove vg IDs matched to more than 1 candidate
dup_vg_IDs <- cand_geno_offset_genotyped$x.ID[which(duplicated(cand_geno_offset_genotyped$x.ID))]
paste(length(dup_vg_IDs), 'vg IDs are duplicates and were matched to > 1 known candidate')

cand_geno_offset_genotyped_no_dups <- cand_geno_offset_genotyped[! which(cand_geno_offset_genotyped$x.ID %in% dup_vg_IDs), ]

## how many unmatched or duplicate (ambiguous) matches ? 
(length(dup_vg_IDs)) + length(unmatched_geno_IDs)
paste((length(dup_vg_IDs) + length(unmatched_geno_IDs)), 'vg IDs are duplicates OR unmatched')


# How many assigned by candidate SVTYPE and platform?
table(cand_geno_offset_genotyped_no_dups$i.SVTYPE)
table(cand_geno_offset_genotyped_no_dups$i.platform, cand_geno_offset_genotyped_no_dups$i.SVTYPE)


# Check that we have type correspondence
type_match <- vector(mode = 'character', length = nrow(cand_geno_offset_genotyped))
#type_match[i] <- 
for (i in 1:nrow(cand_geno_offset_genotyped)){
  cand_type <- cand_geno_offset_genotyped$i.SVTYPE[i]
  geno_type <- cand_geno_offset_genotyped$x.TYPE[i]
  #print(cand_type)
  
  if(cand_type == geno_type){
    type_match[i] <- 'yes'
    # } else if (cand_type != geno_type){
  } else if (cand_type == 'INS' & grepl(pattern = 'INS', geno_type)){
    type_match[i] <- 'yes'
  } else if (cand_type == 'DUP' & grepl(pattern = 'DUP', geno_type)){
    type_match[i] <-'yes'
  } else {
    type_match[i] <- 'no'
  }
}

cand_geno_offset_genotyped$type_match <- type_match
cand_geno_offset_genotyped$FINAL_TYPE <- ifelse(cand_geno_offset_genotyped$type_match == 'no' & cand_geno_offset_genotyped$x.VAR_TYPE == 'SV',
                                                yes = cand_geno_offset_genotyped$i.SVTYPE,
                                                no = cand_geno_offset_genotyped$x.TYPE)

cand_geno_offset_genotyped <- subset(cand_geno_offset_genotyped, FINAL_TYPE != 'indel')

# 5. Plot -----------------------------------------------------------------
## Reorder candidates' platform levels
reordered_platform <- c('SR', 'LR', 'LR + SR')
cand_geno_offset_genotyped_no_dups$i.platform_reordered <- factor(cand_geno_offset_genotyped_no_dups$i.platform, 
                                                                  levels = reordered_platform)
table(cand_geno_offset_genotyped_no_dups$i.platform_reordered, cand_geno_offset_genotyped_no_dups$i.SVTYPE)

# Plot assigned filtered genotyped calls by platform and SVTYPE
ggplot(data = cand_geno_offset_genotyped_no_dups) +
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
    title = "Filtered genotyped SV count matched with known candidates, by size bins and platform"
  ) + 
  scale_fill_viridis_d(option = "B")

# Plot assigned filtered genotyped calls by SVTYPE
# Get SVTYPE values
svtypes <- sort(unique(cand_geno_offset_genotyped_no_dups$i.SVTYPE)) 
#### we sort so that INV falls at the end of vector and 
#### is assigned the most divergent color from DELs, 
#### as INVs are rare and hard to distinguish bar plots

# Get hex code for as many colors as svtypes for a given viridis palette
hex_svtypes <- viridisLite::viridis(n = length(svtypes), option = 'D')
show_col(hex_svtypes)

# Assign a color to each svtype in a named vector
cols_svtypes <- vector(mode = 'character', length = length(svtypes))
for (i in 1:length(svtypes)) {
  names(cols_svtypes)[i] <- svtypes[i]
  cols_svtypes[i] <- hex_svtypes[i]
}


plot_matched_by_type <- 
  ggplot(data = cand_geno_offset_genotyped_no_dups) +
  #facet_wrap(~i.SVTYPE, scales = 'free_y') +
  geom_bar(aes(x = i.SVLEN_bin, fill = i.SVTYPE)) + 
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
    fill = "SV type",
    title = "Filtered genotyped SV count matched with known candidates"
  ) + 
  scale_fill_manual(values = cols_svtypes)

saveRDS(plot_matched_by_type, file = paste0(strsplit(GENOTYPED_FILT, '.table')[[1]], '_matched_offset', OFFSET, 'by_type.rds'))

# 4. Explore UNMATCHED vg calls -------------------------------------------
# Which genotyped calls were not assigned ?
vg_not_assigned <- geno[geno$ID %in% setdiff(geno$ID, cand_geno_offset$x.ID), ]

# How many unique candidates not assigned ?
length(unique(cand_geno_offset$i.ID[is.na(cand_geno_offset$x.ID)])) # unique candidates not genotyped/assigned

# How many filtered genotyped calls not assigned to a known candidate
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

# Plot unassigned filtered genotyped calls by sizes
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
    title = "filtered UNMATCHED genotyped SV count by REF allele size"
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
    title = "Filtered UNMATCHED genotyped SV count by ALT allele size"
  ) 


# 5. Export useful infos --------------------------------------------------
# Export matched calls to a list
write.table(x = cand_geno_offset_genotyped_no_dups[, c('x.CHROM', 'x.POS', 'x.ID', 'i.ID', 'i.SVTYPE', 'i.SVLEN', 'i.SUPP_VEC')],
            file = paste0(strsplit(GENOTYPED_FILT, '.table')[[1]], '_matched_offset', OFFSET, 'bp.txt'),
            col.names = c('CHROM', 'POS', 'ID', 'CAND_ID', 'CAND_SVTYPE', 'CAND_SVLEN', 'CAND_SUPP_VEC'),
            quote = FALSE, row.names = FALSE, sep = "\t")

# Export UNmatched calls to a list
write.table(x = vg_not_assigned_biall[, c('CHROM', 'POS', 'ID')],
            file = paste0(strsplit(GENOTYPED_FILT, '.table')[[1]], '_unmatched_offset', OFFSET, 'bp.txt'),
            col.names = c('CHROM', 'POS', 'ID'),
            quote = FALSE, row.names = FALSE, sep = "\t")

