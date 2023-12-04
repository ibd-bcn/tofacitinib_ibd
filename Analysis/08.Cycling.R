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

#
# Data filtering  --------------------------------------------------------------
#
cycling <- readRDS(file = 'Analysis/data/00_annotation_process/cycling.RDS')

#
# %MT filtering
#
VlnPlot(cycling, features = 'percent.mt')
cycling <- cycling[,cycling$percent.mt < 25]

#
# Genes with no counts filtering
#
counts <- cycling@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 9072
xx <-setdiff(rownames(cycling), names(pp))
cycling <- subset(cycling, features = xx)
cycling
# 19391 features across 2385 samples within 1 assay


#
# Normalization, scaling and UMAP generation -----------------------------------
#
cycling <- NormalizeData(cycling)
cycling <- FindVariableFeatures(cycling,selection.method = "vst", nfeatures = 2000)
cycling <- ScaleData(cycling)
cycling <- RunPCA(cycling, npcs = 100)

PCS <- select_pcs(cycling, 2)
PCS2 <- select_pcs(cycling, 1.6)
ElbowPlot(cycling, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - cycling'))

cycling<- FindNeighbors(cycling,  dims = 1:20, reduction = 'pca')
cycling<-RunUMAP(cycling, dims=1:20, reduction = 'pca')

DimPlot(cycling, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - 20PCs')


#
# Louvain clustering without batch correction ----------------------------------
#
dir.create('Analysis/data/00_annotation_process/cycling')
dir.create('Analysis/data/00_annotation_process/cycling_no_harmony')
cycling <- resolutions(cycling,
                       workingdir = 'Analysis/data/00_annotation_process/cycling_no_harmony',
                       title = 'cycling_no_harmony')
saveRDS(cycling, file = 'Analysis/data/00_annotation_process/cycling/cycling_no_harmony/cycling_filtered_20.RDS')

#
# Louvain clustering with batch correction -------------------------------------
#
cycling <- RunHarmony(cycling, group.by = 'sample', dims.use = 1:20)
ElbowPlot(cycling, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 23, linetype = 2) +
  labs(title = paste0('Harmony - 20PCS'))
cycling<- FindNeighbors(cycling, reduction = "harmony", dims = 1:23)
cycling<-RunUMAP(cycling, dims=1:23, reduction= "harmony")

DimPlot(cycling, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - cycling - Harmony 20-23')


dir.create('Analysis/data/00_annotation_process/cycling_harmony_20_23/')
cycling <- resolutions(cycling,
                       workingdir = 'Analysis/data/00_annotation_process/cycling_harmony_20_23',
                       title = 'Togehter_cycling_filt25_20_23')
saveRDS(cycling, file = 'Analysis/data/00_annotation_process/cycling/cycling_harmony_20_23/cycling_filtered_20_23.RDS')

