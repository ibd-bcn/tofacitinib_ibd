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
stroma <- readRDS(file = 'Analysis/data/00_annotation_process/stroma.RDS')

#
# cell doublets from other subsets
#
FeaturePlot(stroma, features = grep("CD3E$|CD3D$|CD3G$|MS4A1$|^DERL3$|^EPCAM$",
                                    rownames(stroma), value=T))
counts <- stroma@assays$RNA@counts
p <- grep("CD3E$|CD3D$|CD3G$|MS4A1$|^DERL3$|^EPCAM$",rownames(stroma))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 467
xx <-setdiff(colnames(stroma), names(pp))

stroma <- subset(stroma,cells = xx)
stroma
# An object of class Seurat
# 28451 features across 6522 samples within 1 assay
# Active assay: RNA (28451 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

#
# remove genes from IGs
#
gg <- rownames(stroma)[c(grep("^IGH",rownames(stroma)),
                         grep("^IGK", rownames(stroma)),
                         grep("^IGL", rownames(stroma)))]
genes <- setdiff(rownames(stroma),gg)
stroma <- subset(stroma,features = genes)
stroma
# 228185 features across 6522 samples within 1 assay


#
# Filter by %MT genes
#
stroma <- stroma[,stroma$percent.mt < 25]


#
# genes with no expression out
#
counts <- stroma@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 4446
xx <-setdiff(rownames(stroma), names(pp))
stroma <- subset(stroma, features = xx)
stroma
#23739 features across 5726 samples within 1 assay


#
# Normalization, scaling and UMAP generation -----------------------------------
#
stroma <- seurat_to_pca(stroma)

PCS <- select_pcs(stroma, 2)
PCS2 <- select_pcs(stroma, 1.6)
ElbowPlot(stroma, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS+13, y=10, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2+15, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - stroma'))

stroma<- FindNeighbors(stroma,  dims = 1:42, reduction = 'pca')
stroma<-RunUMAP(stroma, dims=1:42, reduction = 'pca')
DimPlot(stroma, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - 42PCs')

#
# Louvain clustering without batch correction ----------------------------------
#
dir.create('Analysis/00_annotation_process/stroma/no_harmony')
stroma <- resolutions(stroma,
                      workingdir = 'Analysis/00_annotation_process/stroma/no_harmony',
                      title = 'stroma_no_h')
saveRDS(stroma, file = 'Analysis/00_annotation_process/stroma/no_harmony/stroma.RDS')

#
# Louvain clustering with batch correction -------------------------------------
#
stroma <- RunHarmony(stroma, group.by = 'sample', dims.use = 1:42)
ElbowPlot(stroma, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 27, linetype = 2) +
  labs(title = paste0('Harmony - 42PCS'))
stroma<- FindNeighbors(stroma, reduction = "harmony", dims = 1:27)
stroma<-RunUMAP(stroma, dims=1:27, reduction= "harmony")

DimPlot(stroma, group.by = 'sample') + labs(title = 'Stroma - Harmony 42_27')


dir.create('Analysis/00_annotation_process/stroma/harmony_42_27')
stroma <- resolutions(stroma,
                      workingdir = 'Analysis/00_annotation_process/stroma/harmony_42_27',
                      title = 'Togehter_stroma_filt27_42_27')
saveRDS(stroma, file = 'Analysis/00_annotation_process/stroma/harmony_42_27/stroma_filtered_42_27.RDS')


#
# Louvain clustering with batch correction (2) ---------------------------------
#
stroma <- RunHarmony(stroma, group.by = 'sample', dims.use = 1:25)
ElbowPlot(stroma, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 27, linetype = 2) +
  labs(title = paste0('Harmony - 25PCS'))
stroma<- FindNeighbors(stroma, reduction = "harmony", dims = 1:27)
stroma<-RunUMAP(stroma, dims=1:27, reduction= "harmony")

DimPlot(stroma, group.by = 'sample') + labs(title = 'Stroma - Harmony 25_27')


dir.create('Analysis/00_annotation_process/stroma/harmony_25_27')
stroma <- resolutions(stroma,
                      workingdir = 'Analysis/00_annotation_process/stroma/harmony_25_27',
                      title = 'stroma_harmony_25_27')
saveRDS(stroma, file = 'Analysis/00_annotation_process/stroma/harmony_25_27/stroma_harmony_25_27.RDS')

#
# Extra cleaning after analysis ------------------------------------------------
#
stroma <- readRDS( 'Analysis/00_annotation_process/stroma/harmony_25_27/stroma_harmony_25_27.RDS')

stroma <- stroma[,!(stroma$RNA_snn_res.1.5 %in% c(16))]

#
# Remove genes with no counts
#
counts <- stroma@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 7
xx <-setdiff(rownames(stroma), names(pp))
stroma <- subset(stroma, features = xx)
stroma
#23732 features across 5674 samples within 1 assay


#
# Reanalysis harmony -----------------------------------------------------------
#
stroma <- NormalizeData(stroma)
stroma <- FindVariableFeatures(stroma,selection.method = "vst", nfeatures = 2000)
stroma <- ScaleData(stroma)
stroma <- RunPCA(stroma, npcs = 100)

PCS <- select_pcs(stroma, 2)
PCS2 <- select_pcs(stroma, 1.6)
ElbowPlot(stroma, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS, y=10, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - stroma'))

stroma <- FindNeighbors(stroma,  dims = 1:41, reduction = 'pca')
stroma <- RunUMAP(stroma, dims=1:41, reduction = 'pca')

DimPlot(stroma, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - 41PCs')

stroma <- RunHarmony(stroma, group.by = 'sample', dims.use = 1:41)

ElbowPlot(stroma, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 25, linetype = 2) +
  geom_vline(xintercept = 15, linetype = 2) +
  labs(title = paste0('Harmony - 41PCS'))

stroma <- FindNeighbors(stroma, reduction = "harmony", dims = 1:27)
stroma <- RunUMAP(stroma, dims=1:27, reduction= "harmony")

DimPlot(stroma, group.by = 'sample') + labs(title = 'Stroma - Harmony 41_27')

dir.create('Analysis/00_annotation_process/stroma/clean_harmony_41_27')
stroma <- resolutions(stroma,
                      workingdir = 'Analysis/00_annotation_process/stroma/clean_harmony_41_27',
                      title = 'clean_harmony_41_27')
saveRDS(stroma, file = 'Analysis/00_annotation_process/stroma/clean_harmony_41_27/clean_harmony_41_27.RDS')

