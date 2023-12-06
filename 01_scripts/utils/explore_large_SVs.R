XL_CUTOFF <- 30000



options(scipen=999)
# Very large SVs in raw, merged dataset -----------------------------------
# Import
merged_SRLR <- read.delim("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/merge_SVs_SRLR/06_filtered/merged_SUPP2.corrected.table", header=FALSE,
                          col.names = c('CHROM', 'POS', 'ID', 'SVTYPE', 'SVLEN', 'END', 'SUPP', 'SUPP_VEC'))

# Extract very large SVs
very_large_merged_SRLR <- subset(merged_SRLR, abs(SVLEN) > XL_CUTOFF)
nrow(very_large_merged_SRLR)

# Add explicit platform info
very_large_merged_SRLR$platform <- sapply(X = as.character(very_large_merged_SRLR$SUPP_VEC), 
                                FUN = function(x){
                                  switch(x,
                                         '10' = 'LR', 
                                         '1' = 'SR',
                                         '11' = 'SR + LR')
                                }
)

# Get number of large SVs by platform
table(very_large_merged_SRLR$platform)
table(very_large_merged_SRLR$platform)/nrow(very_large_merged_SRLR)
## 54.4 % of theses are unique to SR
table(very_large_merged_SRLR$platform, very_large_merged_SRLR$SVTYPE)/nrow(very_large_merged_SRLR)
## 29.5% are very large INVs unique to SR

# Export, for subsetting candidate merged VCF
write.table(very_large_merged_SRLR,
            "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/XL_SVs/merged_SUPP2.corrected_XL_SVs_30kb.table",
            quote = FALSE, sep = "\t", row.names = FALSE)


# Find the 2.5-Mb DEL on ssa10 supported by SR only
subset(very_large_merged_SRLR, CHROM == 'OV354439.1' & SVTYPE == 'DEL' & abs(SVLEN) > 2500000)


# Very large SVs in raw, genotyped dataset --------------------------------
## we do not know genotyped SV length, so we need to use genotyped SVs that were
## matched to a known candidate

# Import
geno_raw <- read.delim("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/08_merged/merged_SUPP2_genotyped_matched_offset5bp.txt")
very_large_geno_raw_SVs <- subset(geno_raw, abs(CAND_SVLEN) > XL_CUTOFF)


# Add explicit platform info
very_large_geno_raw_SVs$platform <- sapply(X = as.character(very_large_geno_raw_SVs$CAND_SUPP_VEC), 
                                          FUN = function(x){
                                            switch(x,
                                                   '10' = 'LR', 
                                                   '1' = 'SR',
                                                   '11' = 'SR + LR')
                                          }
)

# Get number of large SVs by platform
table(very_large_geno_raw_SVs$platform)
table(very_large_geno_raw_SVs$platform)/nrow(very_large_geno_raw_SVs)
## only 25 % of these are unique to SR
table(very_large_geno_raw_SVs$platform, very_large_geno_raw_SVs$CAND_SVTYPE)/nrow(very_large_geno_raw_SVs)
## 63.2 % are very large DUPs unique to LR, whereas only 4.9% are SR INVs, in contrast with merged SV dataset



# Very large SVs in filtered, genotyped dataset --------------------------------
# Import
geno_filtered <- read.delim("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5_matched_offset5bp.txt")
very_large_geno_filt_SVs <- subset(geno_filtered, abs(CAND_SVLEN) > XL_CUTOFF)

# Add explicit platform info
very_large_geno_filt_SVs$platform <- sapply(X = as.character(very_large_geno_filt_SVs$CAND_SUPP_VEC), 
                                           FUN = function(x){
                                             switch(x,
                                                    '10' = 'LR', 
                                                    '1' = 'SR',
                                                    '11' = 'SR + LR')
                                           }
)

# Get number of large SVs by platform
table(very_large_geno_filt_SVs$platform)
table(very_large_geno_filt_SVs$platform)/nrow(very_large_geno_filt_SVs)
## only 33.3 % of these are unique to SR
table(very_large_geno_filt_SVs$platform, very_large_geno_filt_SVs$CAND_SVTYPE)/nrow(very_large_geno_filt_SVs)
## all remaining very large SVs are DELs, mostly from LR


very_large_geno_filt_SVs$END <- very_large_geno_filt_SVs$POS + abs(very_large_geno_filt_SVs$CAND_SVLEN)
# Export, for subsetting candidate merged VCF
write.table(very_large_geno_filt_SVs,
            paste0("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/XL_SVs/merged_SUPP2_MAF0.05_FMISS0.5_XL_SVs_",
            XL_CUTOFF/1000, "kb.table"),
            quote = FALSE, sep = "\t", row.names = FALSE)





# Plot density of candidate XL SVs ----------------------------------------


WIN_CHUNKS <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/02_infos/chrs_win1000000.bed"


win_chunks <- as.data.table(read.table(WIN_CHUNKS,
                                       col.names = c('CHROM', 'START', 'STOP')))


# Get corresponding chunk for each SV
very_large_merged_SRLR <- as.data.table(very_large_merged_SRLR)

very_large_merged_SRLR_win <- 
  very_large_merged_SRLR[win_chunks, 
                    on = .(CHROM, 
                           POS >= START, END <= STOP),
                    .(x.CHROM, x.POS, x.END, x.ID, 
                      i.START, i.STOP # variables I need in merged output (i = overlap_all_bed) (x = win_chunks)
                    )] 




# remove windows where there is nothing
very_large_merged_SRLR_win <- subset(very_large_merged_SRLR_win, !is.na(x.CHROM))


# Count by group  
very_large_merged_SRLR_win_dens <-
  very_large_merged_SRLR_win %>% count(x.CHROM, i.START, i.STOP, 
                               sort = TRUE)
  


# Convert OV to SSA chrom
OV_2_SSA <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/02_infos/OV_to_ssa.txt"
OV_2_ssa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))
## simplify for plotting
OV_2_ssa$CHROM_NUM <- sapply(X = OV_2_ssa$CHROM_SSA, FUN = function(x){
  unlist(strsplit(x, split = 'ssa'))[2]}
)

## merge with full dataset
very_large_merged_SRLR_win_dens <- merge(x = very_large_merged_SRLR_win_dens, y = OV_2_ssa, 
                      by.x = 'x.CHROM', by.y = 'CHROM_OV',
                      sort = FALSE)

## add midPOS
very_large_merged_SRLR_win_dens$midPOS <- (very_large_merged_SRLR_win_dens$i.START + 500000)

ggplot(data = very_large_merged_SRLR_win_dens) +
  facet_grid(. ~ CHROM_NUM, 
             scales = 'free', space = 'free_x') +
  geom_point(aes(x = midPOS, y = n, color = CHROM_NUM), size =0.5
             ) +
  theme_bw() +
  theme(
    # Panels and background
    panel.spacing.x = unit(0.6, 'points'),
    panel.spacing.y = unit(3, 'points'),
    #panel.background = element_rect(color = 'black', linewidth = 0.1),
    #panel.border = element_rect(color = 'black', linewidth = 0.1, fill = NA),
    panel.grid = element_blank(),
    
    # Strips
    strip.text.x.top = element_text(size = 3.5, 
                                    margin = margin(3,0,3,0, 'pt')),
    strip.text.y.right = element_text(size = 4,
                                      margin = margin(0,1,0,1, 'pt')),
    strip.background.y = element_rect(color = 'black', linewidth = 0.1),
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    
    # Axis
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08)
    
  ) + 
  
  #scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  
  guides(color = 'none') +
  
  scale_color_manual(values = rep(c('black', 'grey60'), 
                                  length(unique(very_large_merged_SRLR_win_dens$CHROM_NUM))/2)) +
  labs(x = 'Position along each chromosome',
       y = 'Variant density')






