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
plasmas <- readRDS(file = 'Analysis/data/00_annotation_process/plasmas.RDS')

FeaturePlot(plasmas, features= c('MS4A1', 'DERL3', 'EPCAM'), order = T)
FeaturePlot(plasmas, features= c('CD3E', 'CD3D', 'CD3G', 'EPCAM'), order = T)

#
# cell doublets from other subsets
#
counts <- plasmas@assays$RNA@counts
p <- grep("CD3E$|CD3D$|CD3G$|EPCAM$",rownames(plasmas))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 1008
xx <-setdiff(colnames(plasmas), names(pp))

plasmas <- subset(plasmas,cells = xx)
plasmas
# 28451 features across 28474 samples within 1 assay

#
# %MT filtering
#
plasmas <- plasmas[,plasmas$percent.mt < 25]
plasmas
# 28451 features across 27882 samples within 1 assay

#
# Genes with no counts filtering
#
counts <- plasmas@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp) #5046
xx <-setdiff(rownames(plasmas), names(pp))
plasmas <- subset(plasmas, features = xx)
plasmas
# 23405 features across 27882 samples within 1 assay

#
# Normalization, scaling and UMAP generation -----------------------------------
#
plasmas <- NormalizeData(plasmas)
plasmas <- FindVariableFeatures(plasmas,selection.method = "vst", nfeatures = 2000)
plasmas <- ScaleData(plasmas)
plasmas <- RunPCA(plasmas, npcs = 100)

PCS <- select_pcs(plasmas, 2)
PCS2 <- select_pcs(plasmas, 1.6)
ElbowPlot(plasmas, ndims = 100) +
  geom_vline(xintercept = 45, linetype = 2)+
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - plasmas '))


plasmas<- FindNeighbors(plasmas,  dims = 1:45, reduction = 'pca')
plasmas<-RunUMAP(plasmas, dims=1:45, reduction = 'pca')

DimPlot(plasmas, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - 45PCs')


#
# Louvain clustering without batch correction ----------------------------------
#
setwd('Analysis/00_annotation_process/')
dir.create('plasmas')
dir.create('plasmas/plasmas_no_harmony')

plasmas <- resolutions(plasmas,
                       workingdir = 'Analysis/00_annotation_process/plasmas/plasmas_no_harmony',
                       title = 'plasmas_no_harmony')
saveRDS(plasmas, file = 'Analysis/00_annotation_process/plasmas/plasmas_no_harmony/plasmas_filtered_45.RDS')

#
# Louvain clustering with batch correction -------------------------------------
#
plasmas <- RunHarmony(plasmas, group.by = 'sample', dims.use = 1:25)
ElbowPlot(plasmas, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept =27, linetype = 2) +
  labs(title = paste0('Harmony - 45PCS'))
plasmas<- FindNeighbors(plasmas, reduction = "harmony", dims = 1:27)
plasmas<-RunUMAP(plasmas, dims=1:27, reduction= "harmony")


dir.create('Analysis/00_annotation_process/plasmas/plasmas_harmony_45_27/')
plasmas <- resolutions(plasmas,
                       workingdir = 'Analysis/00_annotation_process/plasmas/plasmas_harmony_45_27',
                       title = 'plasmas_harmony_45_27')
saveRDS(plasmas, file = 'Analysis/00_annotation_process/plasmas/plasmas_harmony_45_27/plasmas_filtered_45_27.RDS')

#
# Extra cleaning after analysis ------------------------------------------------
#
plasmas <- readRDS( 'Analysis/00_annotation_process/plasmas/plasmas_harmony_45_27/plasmas_filtered_45_27.RDS')
DimPlot(plasmas, label=T, group.by = 'RNA_snn_res.0.9')
plasmas <- plasmas[,!(plasmas$RNA_snn_res.0.9 %in% c(14))]

#
# Remove genes with no counts
#
counts <- plasmas@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 4
xx <-setdiff(rownames(plasmas), names(pp))
plasmas <- subset(plasmas, features = xx)
plasmas
# 23401 features across 27808 samples within 1 assay

#
# Reanalysis with Harmony ------------------------------------------------------
#
plasmas <- NormalizeData(plasmas)
plasmas <- FindVariableFeatures(plasmas,selection.method = "vst", nfeatures = 2000)
plasmas <- ScaleData(plasmas)
plasmas <- RunPCA(plasmas, npcs = 100)

PCS <- select_pcs(plasmas, 2)
PCS2 <- select_pcs(plasmas, 1.6)
ElbowPlot(plasmas, ndims = 100) +
  geom_vline(xintercept = 25, linetype = 'dotted') +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS, y=6, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - plasmas'))

plasmas <- FindNeighbors(plasmas,  dims = 1:25, reduction = 'pca')
plasmas <- RunUMAP(plasmas, dims=1:25, reduction = 'pca')

plasmas <- RunHarmony(plasmas, group.by = 'sample', dims.use = 1:25)
ElbowPlot(plasmas, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 17, linetype = 2) +
  geom_vline(xintercept = 22, linetype = 2) +
  labs(title = paste0('Harmony - 25PCS'))
plasmas<- FindNeighbors(plasmas, reduction = "harmony", dims = 1:22)
plasmas<-RunUMAP(plasmas, dims=1:22, reduction= "harmony")

DimPlot(plasmas, group.by = 'sample') + labs(title = 'plasmas - Harmony 25_22')

dir.create('Analysis/00_annotation_process/plasmas/plasmas_reanalysis_harmony_25_22')
plasmas <- resolutions(plasmas,
                       workingdir = 'Analysis/00_annotation_process/plasmas/plasmas_reanalysis_harmony_25_22',
                       title = 'plasmas_reanalysis_harmony_25_22')
saveRDS(plasmas, file = 'Analysis/00_annotation_process/plasmas/plasmas_reanalysis_harmony_25_22/plasmas_reanalysis_harmony_25_22.RDS')


