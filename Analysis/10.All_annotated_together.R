#
# Libraries --------------------------------------------------------------------
#
message('Loading libraries')

library(Seurat)
library(plyr)
library(ggplot2)
library(viridis)
#
# Extra functions
#
message('Loading functions')
source('Analysis/extra_functions/functions_rnaseq.R')

colors_subset <- c('epi' = '#332f6f',
                   'stroma' =  '#e28c40',
                   'tcells' = '#eec35f',
                   'plasmas' = '#985f9d',
                   'myeloids' = '#75bccb',
                   'cycling' = '#7d9f56')

#
# Data loading -----------------------------------------------------------------
#
plasmas <- readRDS(file = 'Analysis/00_annotation_process/00_anotadas/plasmas.RDS')
tcells <- readRDS(file = 'Analysis/00_annotation_process/00_anotadas/tcells.RDS')
epi <- readRDS(file = 'Analysis/00_annotation_process/00_anotadas/epi.RDS')
myeloids <- readRDS(file = 'Analysis/00_annotation_process/00_anotadas/myeloids.RDS')
stroma <- readRDS(file = 'Analysis/00_annotation_process/00_anotadas/stroma.RDS')
cycling <- readRDS(file = 'Analysis/00_annotation_process/00_anotadas/cycling.RDS')

#
# Data cleaning ----------------------------------------------------------------
#
plasmas$response <- factor(plasmas$response, levels = c('R', 'NR'))
tcells$response <- factor(tcells$response, levels = c('R', 'NR'))
epi$response <- factor(epi$response, levels = c('R', 'NR'))
myeloids$response <- factor(myeloids$response, levels = c('R', 'NR'))
stroma$response <- factor(stroma$response, levels = c('R', 'NR'))
cycling$response <- factor(cycling$response, levels = c('R', 'NR'))

plasmas$week_3 <- factor(mapvalues(x = as.character(plasmas$week_2), from = c('W0', 'W8','W16'), to = c('W0', 'POST', 'POST')), levels = c('W0', 'POST'))
tcells$week_3 <- factor(mapvalues(x = as.character(tcells$week_2), from = c('W0', 'W8','W16'), to = c('W0', 'POST', 'POST')), levels = c('W0', 'POST'))
epi$week_3 <- factor(mapvalues(x = as.character(epi$week_2), from = c('W0', 'W8','W16'), to = c('W0', 'W8', 'POST')), levels = c('W0', 'POST'))
myeloids$week_3 <- factor(mapvalues(x = as.character(myeloids$week_2), from = c('W0', 'W8','W16'), to = c('W0', 'POST', 'POST')), levels = c('W0', 'POST'))
stroma$week_3 <- factor(mapvalues(x = as.character(stroma$week_2), from = c('W0', 'W8','W16'), to = c('W0', 'POST', 'POST')), levels = c('W0', 'POST'))
cycling$week_3 <- factor(mapvalues(x = as.character(cycling$week_2), from = c('W0', 'W8','W16'), to = c('W0', 'POST', 'POST')), levels = c('W0', 'POST'))


plasmas$subset <- 'plasmas'
tcells$subset <- 'tcells'
epi$subset <- 'epi'
myeloids$subset <- 'myeloids'
stroma$subset <- 'stroma'
cycling$subset <- 'cycling'

#
# Subset data saving -----------------------------------------------------------
#
saveRDS(plasmas, file = 'Analysis/00_annotation_process/00_anotadas/plasmas.RDS')
saveRDS(tcells, file = 'Analysis/00_annotation_process/00_anotadas/tcells.RDS')
saveRDS(epi, file = 'Analysis/00_annotation_process/00_anotadas/epi.RDS')
saveRDS(myeloids, file = 'Analysis/00_annotation_process/00_anotadas/myeloids.RDS')
saveRDS(stroma, file = 'Analysis/00_annotation_process/00_anotadas/stroma.RDS')
saveRDS(cycling, file = 'Analysis/00_annotation_process/00_anotadas/cycling.RDS')


#
# Data joining -----------------------------------------------------------------
#
todas <- merge(plasmas, tcells)
todas <- merge(todas, epi)
todas <- merge(todas, myeloids)
todas <- merge(todas, stroma)
todas <- merge(todas, cycling)

#
# UMAP Re-generation
#
todas <- seurat_to_pca(todas)
PCS <- select_pcs(todas, 2)
PCS2 <- select_pcs(todas, 1.6)

ElbowPlot(todas, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - merge all'))

todas <- FindNeighbors(todas, dims = 1:36)
todas <- RunUMAP(todas, dims = 1:36)

DimPlot(todas, group.by = 'subset', label=T,shuffle = T, cols = colors_subset)+
  labs(title='All 69813 cells')

#
# Harmony
#
todas <- RunHarmony(todas, group.by = 'sample', dims.use = 1:34)
ElbowPlot(todas, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 35, linetype = 2) +
  geom_vline(xintercept = 25, linetype = 2) +
  labs(title = paste0('Harmony - 34PCS'))
todas <- FindNeighbors(todas, reduction = "harmony", dims = 1:35)
todas <- RunUMAP(todas, dims=1:35, reduction= "harmony")

DimPlot(todas, group.by = 'subset', label=T,shuffle = T, cols = colors_subset)+
  labs(title='All 69813 cells - Harmony')

saveRDS(todas, file = 'Analysis/00_annotation_process/00_anotadas/todas.RDS')
