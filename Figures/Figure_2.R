options(stringsAsFactors = FALSE)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)

source('Figures/functions_plots.R')

## Data ------------------------------------------------------------------------
message('Loading data')

all <- readRDS('/home/acorraliza/TOFA_data/20220222_TOFAS_23/00_annotation_process/00_anotadas/todas_2023.RDS')

## Figure 2 A: UMAP subset -----------------------------------------------------

fig <- DimPlot(all, group.by = 'subset') +
  scale_color_manual(values = colors_subset) +
  theme_umap()

save_sizes(plot = fig, filename = 'Figure_2A', device = 'tiff')
save_sizes(plot = fig, filename = 'Figure_2A', device = 'svg')
save_sizes(plot = fig, filename = 'Figure_2A', device = 'jpeg')
save_sizes(plot = fig, filename = 'Figure_2A', device = 'pdf')


