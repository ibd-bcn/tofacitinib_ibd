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
myeloids <- readRDS(file = 'Analysis/data/00_annotation_process/myeloids.RDS')
# 28451 features across 11271 samples within 1 assay

#
# remove cells with genes from other subsets
#
FeaturePlot(myeloids, features= c('CD3E', 'CD3G', 'MS4A1', 'DERL3'), order = T)
FeaturePlot(myeloids, features= c('IGHA1', 'IGHG1', 'EPCAM', 'COL3A1'), order = T)
counts <- myeloids@assays$RNA@counts
p <- grep("CD3E$|CD3D$|CD3G$|MS4A1$|^DERL3$|^EPCAM$|^COL3A1$",rownames(myeloids))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 836
xx <-setdiff(colnames(myeloids), names(pp))
myeloids <- subset(myeloids,cells = xx)
myeloids
#28451 features across 10435 samples within 1 assay

#
# remove genes from IGs
#
gg <- rownames(myeloids)[c(grep("^IGH",rownames(myeloids)),
                           grep("^IGK", rownames(myeloids)),
                           grep("^IGL", rownames(myeloids)))]
genes <- setdiff(rownames(myeloids),gg)
myeloids <- subset(myeloids,features = genes)
myeloids
# 28185 features across 10435 samples within 1 assay

#
# %MT filtering
#
VlnPlot(myeloids, feature = 'percent.mt')
myeloids <- myeloids[,myeloids$percent.mt < 25]
# 28185 features across 9902 samples within 1 assay

counts <- myeloids@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 7547
xx <-setdiff(rownames(myeloids), names(pp))
myeloids <- subset(myeloids, features = xx)
myeloids
#20638 features across 9902 samples within 1 assay


#
# Normalization, scaling and UMAP generation -----------------------------------
#
myeloids <- NormalizeData(myeloids)
myeloids <- FindVariableFeatures(myeloids,selection.method = "vst", nfeatures = 2000)
myeloids <- ScaleData(myeloids)
myeloids <- RunPCA(myeloids, npcs = 100)

PCS <- select_pcs(myeloids, 2)
PCS2 <- select_pcs(myeloids, 1.6)
ElbowPlot(myeloids, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 30, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - myeloids'))

myeloids<- FindNeighbors(myeloids,  dims = 1:30, reduction = 'pca')
myeloids<-RunUMAP(myeloids, dims=1:30, reduction = 'pca')

DimPlot(myeloids, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - 30PCs')


#
# Louvain clustering without batch correction ----------------------------------
#
dir.create('Analysis/00_annotation_process/myeloids')
dir.create('Analysis/00_annotation_process/myeloids/myeloids_no_harmony')
myeloids <- resolutions(myeloids,
                        workingdir = 'Analysis/00_annotation_process/myeloids/myeloids_no_harmony',
                        title = 'myeloids_no_harmony')
saveRDS(myeloids, file = 'Analysis/00_annotation_process/myeloids/myeloids_no_harmony/myeloids_filtered_30.RDS')


#
# Louvain clustering with batch correction -------------------------------------
#
myeloids <- RunHarmony(myeloids, group.by = 'sample', dims.use = 1:PCS2)
ElbowPlot(myeloids, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 25, linetype = 2) +
  labs(title = paste0('Harmony - 25PCS'))
myeloids<- FindNeighbors(myeloids, reduction = "harmony", dims = 1:25)
myeloids<-RunUMAP(myeloids, dims=1:25, reduction= "harmony")

DimPlot(myeloids, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - myeloids - Harmony 25-25')

dir.create('Analysis/00_annotation_process/myeloids/myeloids_harmony_25_25')
myeloids <- resolutions(myeloids,
                        workingdir = 'Analysis/00_annotation_process/myeloids/myeloids_harmony_25_25',
                        title = 'Togehter_myeloids_filt25_25_25')
saveRDS(myeloids, file = 'Analysis/00_annotation_process/myeloids/myeloids_harmony_25_25/myeloids_filtered_25_25.RDS')
