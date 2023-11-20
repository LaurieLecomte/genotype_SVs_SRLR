# Check genotype concordance between SVs called by LR callers and these SVs when genotyped by vg using short reads

library(data.table)
library(dplyr)

# 1. Import and format ----------------------------------------------------
# Files
LR_SV <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/GT_concordance/merged_sniffles_svim_nanovar_SUPP2_GTs.table"
RAW_GT <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/GT_concordance/merged_SUPP2_genotyped.tagged_GTs.table"

MATCHED_RAW <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/08_merged/merged_SUPP2_genotyped_matched_offset5bp.txt"

FILT_GT <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/GT_concordance/merged_SUPP2_MAF0.05_FMISS0.5_GTs.table"
MATCHED_FILT <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5_matched_offset5bp.txt"

SR_SV <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/GT_concordance/merged_delly_manta_smoove_SUPP2.vcf_LR_samples_GTs.table"

samples <- c('13070', '14062', '14010', '14104')

# LR SVs
LR_SVs <- fread(LR_SV, colClasses = 'character')
LR_SVs <- LR_SVs[, 1:(ncol(LR_SVs)-1)]
colnames(LR_SVs) <- 
  c(
    sapply(X = strsplit(x = colnames(LR_SVs), split = ']'), FUN ="[", 2)[1:4],
    sub(".*]([0-9A-Z]+_[0-9A-Za-z]+):GT", "\\1", colnames(LR_SVs)[5:length(colnames(LR_SVs))])
  )

# SR SVs
SR_SVs <- fread(SR_SV, colClasses = 'character')
SR_SVs <- SR_SVs[, 1:(ncol(SR_SVs)-1)]
colnames(SR_SVs) <- 
  c(
    sapply(X = strsplit(x = colnames(SR_SVs), split = ']'), FUN ="[", 2)[1:4],
    sub(".*]([0-9A-Z]+_[0-9A-Za-z]+):GT", "\\1", colnames(SR_SVs)[5:length(colnames(SR_SVs))])
  )


# Raw vg genotypes
raw_GTs <- fread(RAW_GT, colClasses = 'character')
raw_GTs <- raw_GTs[, 1:(ncol(raw_GTs)-1)]

colnames(raw_GTs) <- 
  c(
    sapply(X = strsplit(x = colnames(raw_GTs), split = ']'), FUN ="[", 2)[1:4],
    sub(".*]([0-9A-Za-z]+):GT", "\\1", colnames(raw_GTs)[5:length(colnames(raw_GTs))])
  )

# Matched raw genotyped SVs
matched_raw <- fread(MATCHED_RAW, colClasses = 'character')


# 2. Check GT concordance for RAW vg genotypes -------------------------------------------------

# First extract genotyped SVs for which we KNOW platform support
raw_GTs_matched <- 
  merge(matched_raw,
        raw_GTs,
        by = c('CHROM', 'POS', 'ID'),
        sort = FALSE
  )

# Function to recode genotypes
recode_GTs <- function(x){
  switch(x,
         '1/1' = 2,
         '1/0' = 1,
         '0/1' = 1,
         '0/0' = 0,
         './.' = -1)
}

# Function to check GT concordance 
check_concord <- function(sample, matched_df, supp_type, type) {
  # First make sure that matched_df is a dataframe, not a datatable
  matched_table <- as.data.frame(matched_df)
  
  # Extract columns for given sample
  cols_to_keep <- c(1:8, grep(pattern = sample, x = colnames(matched_table)))
  GTs <- matched_table[, cols_to_keep]
  
  colnames(GTs) <- switch(supp_type,
                          'LR' = c(colnames(GTs)[1:8],
                                   'vg', 'sniffles', 'svim', 'nanovar'),
                          'SR' = c(colnames(GTs)[1:8],
                                   'vg', 'delly', 'manta', 'smoove'))
  
  # Check concordance within the 3 callers used for given platform
  for (i in 1:nrow(GTs)) {
    
    # Extract genotypes from the 3 callers
    GT_counts <- table(as.character(GTs[i, 10:12]))
    
    # Determine if these genotypes agree and get most frequent genotype per locus
    if(any(GT_counts == 3) ){ # All 3 callers agree
      GTs$caller_concord[i] <- names(which(GT_counts == 3))
      
    } else if (any(GT_counts == 2)) { # 2/3 callers agree
      GTs$caller_concord[i] <- names(which(GT_counts == 2))
      
    } else if (all(GT_counts == 1)){ # all callers differ 
      GTs$caller_concord[i] <- 'NC'
    }
  }
  
  
  # Check concordance between callers and vg
  for (i in 1:nrow(GTs)) {
    if (GTs$caller_concord[i] == 'NC'){ # if callers do not agree together, then they don't agree with vg
      GTs$callers_vg_concord[i] <- 'unknown'
    } else if (GTs$vg[i] == GTs$caller_concord[i] ) { # callers and vg both agree
      GTs$callers_vg_concord[i] <- 'concordant'
    } else if (GTs$vg[i] != GTs$caller_concord[i] && GTs$caller_concord[i] != 'NC') { # callers and vg do NOT agree
      GTs$callers_vg_concord[i] <- 'non-concordant'
    }
  }
  
  # Add sample info
  GTs$sample <- sample
  
  ## Check proportion of concordant genotypes
  #print('Proportion of concordant and non-concordant calls :')
  #return(table(GTs$callers_vg_concord) / nrow(GTs))
  
  assign(x = paste0('GTs_', supp_type, '_', type, '_', sample), value = GTs, envir = .GlobalEnv)
}




# 2.1 LR calls vs vg SR genotypes -----------------------------------------
# Extract KNOWN LR calls
raw_GTs_LR <- subset(raw_GTs_matched, CAND_SUPP_VEC == '10')

## Correct CAND_ID column to match
raw_GTs_LR$CAND_ID <- 
  sapply(X = raw_GTs_LR$CAND_ID,
         FUN = function(x) substr(x, start = 3, stop = nchar(x)))


# Merge with original LR info, to have indv LR genotypes
raw_GTs_LR_matched <- 
  merge(x = raw_GTs_LR, y = select(LR_SVs, -c('CHROM', 'POS', 'END')),
        by.x = 'CAND_ID', by.y = 'ID')

raw_GTs_LR_matched0 <- as.data.frame(raw_GTs_LR_matched)

# Recode genotypes as -1, 0, 1 and 2 for lighter processing
for(i in 9:24) {
  raw_GTs_LR_matched0[[i]] <- sapply(X = raw_GTs_LR_matched[[i]], FUN = recode_GTs)
}

# Check genotype concordance in each sample
for (i in samples){
  check_concord(i, raw_GTs_LR_matched0, supp_type = 'LR', type = 'raw')
}

#raw_GTs_LR_all_samples <- do.call('rbind', c())

all_raw_GTs_LR <- do.call(rbind, mget(paste0('GTs_LR_raw_', samples)))
all_raw_GTs_LR$platform <- 'LR'
all_raw_GTs_LR$type <- 'raw'

## Check proportion of concordant genotypes
table(all_raw_GTs_LR$callers_vg_concord, all_raw_GTs_LR$sample)/ nrow(raw_GTs_LR_matched0)

LR_raw <- as.data.frame(table(all_raw_GTs_LR$callers_vg_concord, all_raw_GTs_LR$sample)/ nrow(raw_GTs_LR_matched0))
colnames(LR_raw) <- c('GT', 'sample', 'prop')
LR_raw$platform <- 'LR'
LR_raw$type <- 'raw'


# 2.2 SR calls vs vg SR genotypes -----------------------------------------

# Extract known SR calls
raw_GTs_SR <- subset(raw_GTs_matched, CAND_SUPP_VEC == '01')

## Correct CAND_ID column to match
raw_GTs_SR$CAND_ID <- 
  sapply(X = raw_GTs_SR$CAND_ID,
         FUN = function(x) substr(x, start = 3, stop = nchar(x)))

# Merge with original SR info, to have indv LR genotypes
raw_GTs_SR_matched <- 
  merge(x = raw_GTs_SR, y = select(SR_SVs, -c('CHROM', 'POS', 'END')),
        by.x = 'CAND_ID', by.y = 'ID')

raw_GTs_SR_matched0 <- as.data.frame(raw_GTs_SR_matched)

# Recode genotypes as -1, 0, 1 and 2 for lighter processing
for(i in 9:24) {
  raw_GTs_SR_matched0[[i]] <- sapply(X = raw_GTs_SR_matched[[i]], FUN = recode_GTs)
}

# Check genotype concordance in each sample
for (i in samples){
  check_concord(i, raw_GTs_SR_matched0, supp_type = 'SR', type = 'raw')
}

all_raw_GTs_SR <- do.call(rbind, mget(paste0('GTs_SR_raw_', samples)))
all_raw_GTs_SR$platform <- 'SR'
all_raw_GTs_SR$type <- 'raw'

## Check proportion of concordant genotypes
table(all_raw_GTs_SR$callers_vg_concord, all_raw_GTs_SR$sample)/ nrow(raw_GTs_SR_matched0)

SR_raw <- as.data.frame(table(all_raw_GTs_SR$callers_vg_concord, all_raw_GTs_SR$sample)/ nrow(raw_GTs_SR_matched0))
colnames(SR_raw) <- c('GT', 'sample', 'prop')
SR_raw$platform <- 'SR'
SR_raw$type <- 'raw'

# 3. Check GT concordance for FILTERED vg genotypes --------------------
# First extract genotyped SVs for which we KNOW platform support
filt_GTs <- fread(FILT_GT, colClasses = 'character')
filt_GTs <- filt_GTs[, 1:(ncol(filt_GTs)-1)]

colnames(filt_GTs) <- 
  c(
    sapply(X = strsplit(x = colnames(filt_GTs), split = ']'), FUN ="[", 2)[1:4],
    sub(".*]([0-9A-Za-z]+):GT", "\\1", colnames(filt_GTs)[5:length(colnames(filt_GTs))])
  )


# Matched filt genotyped SVs
matched_filt <- fread(MATCHED_FILT, colClasses = 'character')

# First extract genotyped SVs for which we KNOW platform support
filt_GTs_matched <- 
  merge(matched_filt,
        filt_GTs,
        by = c('CHROM', 'POS', 'ID'),
        sort = FALSE
  )


# 3.1 LR calls vs FILTERED vg SR genotypes ----------------------------
# Extract KNOWN LR calls
filt_GTs_LR <- subset(filt_GTs_matched, CAND_SUPP_VEC == '10')

## Correct CAND_ID column to match
filt_GTs_LR$CAND_ID <- 
  sapply(X = filt_GTs_LR$CAND_ID,
         FUN = function(x) substr(x, start = 3, stop = nchar(x)))

# Merge with original LR info
filt_GTs_LR_matched <- 
  merge(x = filt_GTs_LR, y = select(LR_SVs, -c('CHROM', 'POS', 'END')),
        by.x = 'CAND_ID', by.y = 'ID')

filt_GTs_LR_matched0 <- as.data.frame(filt_GTs_LR_matched)


for(i in 9:24) {
  filt_GTs_LR_matched0[[i]] <- sapply(X = filt_GTs_LR_matched[[i]], FUN = recode_GTs)
}

# Check genotype concordance in each sample
for (i in samples){
  check_concord(i, filt_GTs_LR_matched0, supp_type = 'LR', type = 'filt')
}


all_filt_GTs_LR <- do.call(rbind, mget(paste0('GTs_LR_filt_', samples)))
all_filt_GTs_LR$platform <- 'LR'
all_filt_GTs_LR$type <- 'filt'


## Check proportion of concordant genotypes
table(all_filt_GTs_LR$callers_vg_concord, all_filt_GTs_LR$sample)/ nrow(filt_GTs_LR_matched0)

LR_filt <- as.data.frame(table(all_filt_GTs_LR$callers_vg_concord, all_filt_GTs_LR$sample)/ nrow(filt_GTs_LR_matched0))
colnames(LR_filt) <- c('GT', 'sample', 'prop')
LR_filt$platform <- 'LR'
LR_filt$type <- 'filt'

# 3.2 SR calls vs FILTERED vg SR genotypes ----------------------------
# Extract KNOWN LR calls
filt_GTs_SR <- subset(filt_GTs_matched, CAND_SUPP_VEC == '01')

## Correct CAND_ID column to match
filt_GTs_SR$CAND_ID <- 
  sapply(X = filt_GTs_SR$CAND_ID,
         FUN = function(x) substr(x, start = 3, stop = nchar(x)))

# Merge with original SR info
filt_GTs_SR_matched <- 
  merge(x = filt_GTs_SR, y = select(SR_SVs, -c('CHROM', 'POS', 'END')),
        by.x = 'CAND_ID', by.y = 'ID')

filt_GTs_SR_matched0 <- as.data.frame(filt_GTs_SR_matched)


for(i in 9:24) {
  filt_GTs_SR_matched0[[i]] <- sapply(X = filt_GTs_SR_matched[[i]], FUN = recode_GTs)
}

# Check genotype concordance in each sample
for (i in samples){
  check_concord(i, filt_GTs_SR_matched0, supp_type = 'SR', type = 'filt')
}


all_filt_GTs_SR <- do.call(rbind, mget(paste0('GTs_SR_filt_', samples)))
all_filt_GTs_SR$platform <- 'SR'
all_filt_GTs_SR$type <- 'filt'


## Check proportion of concordant genotypes
table(all_filt_GTs_SR$callers_vg_concord, all_filt_GTs_SR$sample)/ nrow(filt_GTs_SR_matched0)


SR_filt <- as.data.frame(table(all_filt_GTs_SR$callers_vg_concord, all_filt_GTs_SR$sample)/ nrow(filt_GTs_SR_matched0))
colnames(SR_filt) <- c('GT', 'sample', 'prop')
SR_filt$platform <- 'SR'
SR_filt$type <- 'filt'



# Get proportion of concordant, non-concordant and unknown across datasets

all_tables <- do.call(rbind, list(LR_raw, SR_raw, SR_filt, LR_filt))
aggregate(data = all_tables, prop ~ GT + platform + type, FUN = mean)



#all_tables_GTs <- do.call(rbind, list(all_raw_GTs_LR[, c(1:5,14:17)], all_raw_GTs_SR[, c(1:5,14:17)],
#                                      all_filt_GTs_LR[, c(1:5,14:17)], all_filt_GTs_SR[, c(1:5,14:17)]))






# Poubelle ----------------------------------------------------------------
# Recode genotypes
recode_GTs <- function(x){
  switch(x,
         '1/1' = 'hom-alt',
         '1/0' = 'het',
         '0/1' = 'het',
         '0/0' = 'hom-ref',
         './.' = 'mis')
}

