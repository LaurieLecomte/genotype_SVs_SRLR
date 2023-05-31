# Compute mean number of genotyped sites per sample

# 1. Access files from command line ---------------------------------------
argv <- commandArgs(T)
COUNTS <- argv[1]


# 2. Import and compute mean ----------------------------------------------
samples_calls <- read.table(COUNTS, col.names = c('SAMPLE', 'CALLS'))

print(mean(samples_calls$CALLS))

