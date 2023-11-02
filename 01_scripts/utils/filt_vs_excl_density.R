# Compute mean variant density by chromosome and plot

options(scipen=999)

library(ggplot2)

# 1. Access files in command line, import and format ----------------------
#argv <- commandArgs(T)
EXCL_DENSITY_TABLE <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_excluded_density_win100000bp.txt" 
FILT_DENSITY_TABLE <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_filtered_density_win100000bp.txt"
OV_2_SSA <- "02_infos/OV_to_ssa.txt"

WIN_SIZE <- 100000

excl_var_density <-  read.delim(EXCL_DENSITY_TABLE, header = FALSE,
                           col.names = c('CHROM', 'BIN_START', 'BIN_STOP', 'DENSITY'))

filt_var_density <-  read.delim(FILT_DENSITY_TABLE, header = FALSE,
                                col.names = c('CHROM', 'BIN_START', 'BIN_STOP', 'DENSITY'))
## Convert OV chromosomes to ssa notation
OV_2_sasa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))
## simplify for plotting
OV_2_sasa$CHROM_NUM <- sapply(X = OV_2_sasa$CHROM_SSA, FUN = function(x){
  unlist(strsplit(x, split = 'ssa'))[2]}
)


excl_var_density_OV <- merge(x = excl_var_density, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV', sort = FALSE)
filt_var_density_OV <- merge(x = filt_var_density, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV', sort = FALSE)

excl_var_density_OV$type<- 'excl'
filt_var_density_OV$type <- 'kept'


var_density_OV <- rbind(excl_var_density_OV, filt_var_density_OV)

var_density_OV$BIN_MIDDLE <- var_density_OV$BIN_START + WIN_SIZE/2

#density_plot <-
  ggplot(data = var_density_OV) + 
  #facet_grid(.~CHROM_NUM, scales = 'free_x', space = 'free_x', nrow = 2) +
  facet_wrap(~CHROM_NUM, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = DENSITY, color = type), size = 0.8, alpha = 0.5) +
  theme(# Panels and background
    panel.spacing.x = unit(0.6, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    
        strip.text.x = element_text(size = 4),
        axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
        strip.placement = "inside",
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08),
    
  ) + 

  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = paste('Density by', WIN_SIZE/1000000, 'Mb window')
       )
  



# Plot proportions of excluded/filtered variants --------------------------
library(tidyr)
  

  
excl_var_density$type<- 'excl'
filt_var_density$type <- 'kept'
  
var_density <- rbind(filt_var_density, excl_var_density)  

var_density <- 
merge(x = var_density, y = OV_2_sasa, 
      by.x = 'CHROM', by.y = 'CHROM_OV',
      all = TRUE, sort = FALSE)

var_density$BIN_MIDDLE <- var_density$BIN_START + WIN_SIZE/2



var_density_splitted <- 
  pivot_wider(data = var_density,
            values_from = DENSITY,
            names_from = type)

var_density_splitted$prop_kept <- var_density_splitted$kept / (var_density_splitted$kept + var_density_splitted$excl)
var_density_splitted$prop_excl <- var_density_splitted$excl / (var_density_splitted$kept + var_density_splitted$excl)

# Prop excl plot
ggplot(data = var_density_splitted) + 
  facet_wrap(~CHROM_NUM, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = prop_excl, color = CHROM_NUM), 
             size = 0.6, alpha = 0.8) +
  theme(
    panel.spacing.x = unit(1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    
    strip.text.x = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "inside",
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08)
  ) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  
  guides(color = 'none') +
  
  scale_color_manual(values = rep(c('red4', 'snow4'), 
                                  length(unique(var_density_splitted$CHROM_NUM))/2)) +
  labs(x = 'Position along chromosomes',
       y = 'Variant density'
  )

ggsave(filename = paste0(unlist(strsplit(EXCL_DENSITY_TABLE, split = '.txt')), '_propEXCL.png'),
       width = 3100,
       height = 2800,
       units = 'px',
       dpi = 700,
)


# Abs excl plot
ggplot(data = var_density_splitted) + 
  facet_wrap(~CHROM_NUM, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = excl, color = CHROM_NUM), 
             size = 0.6, alpha = 0.8) +
  theme(
    panel.spacing.x = unit(1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    
    strip.text.x = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "inside",
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08)
  ) + 
  guides(color = 'none') +
  
  scale_color_manual(values = rep(c('red4', 'snow4'), 
                                  length(unique(var_density_splitted$CHROM_NUM))/2)) +
  labs(x = 'Position along chromosomes',
       y = 'Variant density'
  )

ggsave(filename = paste0(unlist(strsplit(EXCL_DENSITY_TABLE, split = '.txt')), '_absEXCL.png'),
       width = 3100,
       height = 2800,
       units = 'px',
       dpi = 700,
)

# Prop kept plot
ggplot(data = var_density_splitted) + 
  facet_wrap(~CHROM_NUM, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = prop_kept, color = CHROM_NUM), 
             size = 0.6, alpha = 0.8) +
  theme(
    panel.spacing.x = unit(1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    
    strip.text.x = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "inside",
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08)
  ) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  
  guides(color = 'none') +
  
  scale_color_manual(values = rep(c('darkgreen', 'snow4'), 
                                  length(unique(var_density_splitted$CHROM_NUM))/2)) +
  labs(x = 'Position along chromosomes',
       y = 'Variant density'
  )

ggsave(filename = paste0(unlist(strsplit(FILT_DENSITY_TABLE, split = '.txt')), '_propKEPT.png'),
       width = 3100,
       height = 2800,
       units = 'px',
       dpi = 700,
)

# Abs kept plot
ggplot(data = var_density_splitted) + 
  facet_wrap(~CHROM_NUM, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = kept, color = CHROM_NUM), 
             size = 0.6, alpha = 0.8) +
  theme(
    panel.spacing.x = unit(1, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    
    strip.text.x = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "inside",
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08)
  ) + 
  guides(color = 'none') +
  
  scale_color_manual(values = rep(c('darkgreen', 'snow4'), 
                                  length(unique(var_density_splitted$CHROM_NUM))/2)) +
  labs(x = 'Position along chromosomes',
       y = 'Variant density'
  )

ggsave(filename = paste0(unlist(strsplit(FILT_DENSITY_TABLE, split = '.txt')), '_absKEPT.png'),
       width = 3100,
       height = 2800,
       units = 'px',
       dpi = 700,
)


# Plot both kept and excl proportions together in different panels
var_density_splitted_long <- 
  pivot_longer(data = var_density_splitted,
               cols = c('prop_kept', 'prop_excl'),
               names_to = 'type',
               values_to = 'prop')


ggplot(data = var_density_splitted_long) + 
  facet_grid(type ~ CHROM_NUM, scales = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = prop, color = interaction(CHROM_NUM, type)), 
             size = 0.1, alpha = 0.5) +
  theme(
      # Panels and background
      panel.spacing.x = unit(0.6, 'points'),
      panel.spacing.y = unit(3, 'points'),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      
      # Strips
      strip.text.x.top = element_text(size = 3, 
                                      margin = margin(3,0,3,0, 'pt')),
      strip.text.y.right = element_blank(),
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
  
  guides(color = 'none') +
  
  scale_color_manual(values = c(
    rep(c('red4', 'snow4'), length(unique(var_density_splitted$CHROM_NUM))/2),
    rep(c('darkgreen', 'snow4'), length(unique(var_density_splitted$CHROM_NUM))/2))
    )+
  labs(x = 'Position along chromosomes',
       y = 'Proportion of variants'
  )

ggsave(filename = paste0(unlist(strsplit(FILT_DENSITY_TABLE, split = '.txt')), '_kept_excl_props.png'),
       width = 2800,
       height = 3000,
       units = 'px',
       dpi = 700,
)

# Plot both kept and excl absolute numbers together in different panels
var_density_splitted_long_abs <- 
  pivot_longer(data = var_density_splitted,
               cols = c('kept', 'excl'),
               names_to = 'type',
               values_to = 'number')

ggplot(data = var_density_splitted_long_abs) + 
  facet_grid(type ~ CHROM_NUM, scales = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = number, color = interaction(CHROM_NUM, type)), 
             size = 0.1, alpha = 0.5) +
  theme(
    # Panels and background
    panel.spacing.x = unit(0.6, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    # Strips
    strip.text.x.top = element_text(size = 3, 
                                    margin = margin(3,0,3,0, 'pt')),
    strip.text.y.right = element_blank(),
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
  
  guides(color = 'none') +
  
  scale_color_manual(values = c(
    rep(c('red4', 'snow4'), length(unique(var_density_splitted$CHROM_NUM))/2),
    rep(c('darkgreen', 'snow4'), length(unique(var_density_splitted$CHROM_NUM))/2))
  )+
  labs(x = 'Position along chromosomes',
       y = 'Number of variants'
  )

ggsave(filename = paste0(unlist(strsplit(FILT_DENSITY_TABLE, split = '.txt')), '_kept_excl_abs.png'),
       width = 2800,
       height = 3000,
       units = 'px',
       dpi = 700,
)
