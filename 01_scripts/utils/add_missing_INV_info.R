# Add SVLEN and END info for INVs to candidates, since it was lost 
# when SR and LR calls were merged prior to genotyping

library(data.table)
library(ggplot2)

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
SR <- argv[1]
LR <- argv[2]
CANDIDATES <- argv[3]

# Short reads SVs
SR_SVs <- read.table(SR, header = FALSE, 
             col.names = c('CHROM', 'POS', 'ID', 'SVTYPE', 'SVLEN',
                           'END', 'SUPP', 'SUPP_VEC'),
             colClasses = c('character')
  )

SR_SVs[, c('POS', 'SVLEN', 'END', 'SUPP')] <- sapply(SR_SVs[, c('POS', 'SVLEN', 'END', 'SUPP')],
                                                     as.numeric)

# Long reads SVs
LR_SVs <- read.table(LR, header = FALSE, 
                     col.names = c('CHROM', 'POS', 'ID', 'SVTYPE', 'SVLEN',
                                   'END', 'SUPP', 'SUPP_VEC'),
                     colClasses = c('character')
)

LR_SVs[, c('POS', 'SVLEN', 'END', 'SUPP')] <- sapply(LR_SVs[, c('POS', 'SVLEN', 'END', 'SUPP')],
                                                     as.numeric)

# Candidate SVs
candidates_SVs <- read.table(CANDIDATES, header = FALSE, 
                     col.names = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'SVTYPE', 'SVLEN',
                                   'END', 'SUPP', 'SUPP_VEC'),
                     colClasses = c('character')
)
candidates_SVs[, c('POS', 'SVLEN', 'END', 'SUPP')] <- sapply(candidates_SVs[, c('POS', 'SVLEN', 'END', 'SUPP')],
                                                     as.numeric)

# 2. Get SVLEN and END for INVs with missing info -------------------------
# Which candidates INVs lost SVLEN and END info?
INVs_to_edit <- which(candidates_SVs$SVLEN == 0 & candidates_SVs$END == candidates_SVs$END)

# Copy original candidates
candidates_SVs_ed <- candidates_SVs

# Get info for each INVs with missing info
for (i in INVs_to_edit){
  # Get original ID prior to SR LR merge
  id <- sub("[0-9]_(.*)", "\\1", candidates_SVs$ID[i])
  # Search for SVLEN and END info in the right dataframe, based on SUPP_VEC
  candidates_SVs_ed$SVLEN[i] <- 
    ifelse(test = candidates_SVs[i, 'SUPP_VEC'] == '11' | candidates_SVs[i, 'SUPP_VEC'] == '10', # if LR or LR + SR
           yes = LR_SVs$SVLEN[LR_SVs$ID == id],
           no = SR_SVs$SVLEN[SR_SVs$ID == id]) # if SR
  candidates_SVs_ed$END[i] <- 
    ifelse(test = candidates_SVs[i, 'SUPP_VEC'] == '11' | candidates_SVs[i, 'SUPP_VEC'] == '10',
           yes = LR_SVs$END[LR_SVs$ID == id],
           no = SR_SVs$END[SR_SVs$ID == id])
}

# Convert SVLEN and END to integer
candidates_SVs_ed$END <- as.integer(candidates_SVs_ed$END)
candidates_SVs_ed$SVLEN <- as.integer(candidates_SVs_ed$SVLEN)



# 3. Write to output annotation file --------------------------------------
write.table(candidates_SVs_ed[, c('CHROM', 'POS', 'ID','SVTYPE', 'SVLEN', 'END')], 
            file = paste0(CANDIDATES, '.annot'),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


