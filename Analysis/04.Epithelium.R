#
# Libraries --------------------------------------------------------------------
#
message('Loading libraries')

library(Seurat)
library(plyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(scater)
#
# Extra functions
#
message('Loading functions')
source('Analysis/extra_functions/functions_rnaseq.R')

#
# Data filtering  --------------------------------------------------------------
#
epi <- readRDS(file = 'Analysis/data/00_annotation_process/epi.RDS')

sce <- as.SingleCellExperiment(epi)
# Mitochondrial genes
mito_genes <- rownames(sce)[grep("^MT-", rownames(sce))]

# Ribosomal genes
ribo_genes <- rownames(sce)[grep("^RP[SL]", rownames(sce))]

# Hemoglobin genes - includes all genes starting with HB except HBP.
hb_genes <- rownames(sce)[grep("^HB[^(P)]", rownames(sce))]

sce <- scuttle::addPerCellQC(sce, flatten = T, subsets = list(mt = mito_genes, hb = hb_genes,
                                                     ribo = ribo_genes))
sce@colData$sample <- 'nuclei'

plot_grid(plotColData(sce, y = "detected", x = "total"),
          plotColData(sce, y = "detected", x = "RNA_snn_res.0.5"),
          plotColData(sce, y = "total", x = "RNA_snn_res.0.5"),
          plotColData(sce, y = "subsets_mt_percent",  x = "RNA_snn_res.0.5"),
          plotColData(sce, y = "subsets_ribo_percent", x = 'RNA_snn_res.0.5'),
          plotColData(sce, y = "subsets_hb_percent",  x = "RNA_snn_res.0.5"),
          ncol = 3)


#
# cell doublets from other subsets
#
FeaturePlot(epi, features = c('CD79A','MS4A1', 'DERL3', 'FCGR3B'), order =T, cols = c('lightgray', 'red'), pt.size = 3 )
FeaturePlot(epi, features = c('CD3E','CD3D', 'CD3G', 'CD8A'), order =T, cols = c('lightgray', 'red'), pt.size = 3 )
VlnPlot(epi, features = c('CD3E','DERL3', 'MS4A1', 'CD79A','EPCAM', 'MKI67'))

counts <- epi@assays$RNA@counts
p <- grep("CD3E$|CD3D$|CD3G$|MS4A1$|^DERL3$|^CD79A$|^FCGR3B$",rownames(epi))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 1394
xx <-setdiff(colnames(epi), names(pp))

epi <- subset(epi, cells = xx)


#
# IGs genes out
#
gg <- rownames(epi)[c(grep("^IGH",rownames(epi)),
                      grep("^IGK", rownames(epi)),
                      grep("^IGL", rownames(epi)))]
genes <- setdiff(rownames(epi),gg)
epi <- subset(epi,features = genes)
epi
# 28185 features across 9680 samples within 1 assay

#
# genes with no expression out
#
counts <- epi@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 3752
xx <-setdiff(rownames(epi), names(pp))
epi <- subset(epi, features = xx)
epi
# 24433 features across 9680 samples within 1 assay


#
# Normalization, scaling and UMAP generation -----------------------------------
#
epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi,selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi)
epi <- RunPCA(epi, npcs = 100)

PCS <- select_pcs(epi, 2)
PCS2 <- select_pcs(epi, 1.6)
ElbowPlot(epi, ndims = 100) +
  geom_vline(xintercept = PCS) +
  # geom_vline(xintercept = 22, linetype = 2) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS+29, y=4, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2+19, y=6, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - Epi'))


epi<- FindNeighbors(epi,  dims = 1:PCS2, reduction = 'pca')
epi<-RunUMAP(epi, dims=1:PCS2, reduction = 'pca')

DimPlot(epi, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - 33PCs')

#
# Louvain clustering without batch correction ----------------------------------
#
setwd('Analysis/00_annotation_process/')
dir.create('epi')
dir.create('epi/epi_no_harmony')

epi <- resolutions(epi,
                   workingdir = 'Analysis/00_annotation_process/epi/epi_no_harmony',
                   title = 'epi_no_harmony')
saveRDS(epi, file = 'Analysis/00_annotation_process/epi/epi_no_harmony/epi_filtered_33.RDS')


#
# Louvain clustering with batch correction -------------------------------------
#
epi <- RunHarmony(epi, group.by = 'sample', dims.use = 1:PCS2)
ElbowPlot(epi, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 26, linetype = 2) +
  labs(title = paste0('Harmony - 33PCS'))
epi<- FindNeighbors(epi, reduction = "harmony", dims = 1:26)
epi<-RunUMAP(epi, dims=1:26, reduction= "harmony")

DimPlot(epi, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - epi - Harmony 33-26')

dir.create('Analysis/00_annotation_process/epi/epi_harmony_33_26')
epi <- resolutions(epi,
                   workingdir = 'Analysis/00_annotation_process/epi/epi_harmony_33_26',
                   title = 'epi_harmony_33_26')

saveRDS(epi, file = 'Analysis/00_annotation_process/epi/epi_harmony_33_26/epi_filtered_33_26.RDS')

#
# Extra cleaning after analysis ------------------------------------------------
#
# Since epithelium has a less strict filter for MT genes, we have to make 5 rounds
# of correction until we get the cells that have key information.
#
# 1rst
#
epi <- readRDS('Analysis/00_annotation_process/epi/epi_harmony_33_26/epi_filtered_33_26.RDS')
DimPlot(epi, label=T)
epi <- epi[,!(epi$RNA_snn_res.1.5 %in% c(27, 26, 28, 31))]

counts <- epi@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 90
xx <-setdiff(rownames(epi), names(pp))
epi <- subset(epi, features = xx)
epi
# 24343 features across 9456 samples within 1 assay

epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi,selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi)
epi <- RunPCA(epi, npcs = 100)

PCS <- select_pcs(epi, 2)
PCS2 <- select_pcs(epi, 1.6)
ElbowPlot(epi, ndims = 100) +
  geom_vline(xintercept = 50, linetype = 'dotted') +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS, y=6, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - epi'))

epi<- FindNeighbors(epi,  dims = 1:33, reduction = 'pca')
epi<-RunUMAP(epi, dims=1:33, reduction = 'pca')

DimPlot(epi, group.by = 'sample') +
  labs(title = 'Tofacitinib 23 samples - 33PCs')

epi <- RunHarmony(epi, group.by = 'sample', dims.use = 1:33)
ElbowPlot(epi, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 19, linetype = 2) +
  geom_vline(xintercept = 36, linetype = 2) +
  labs(title = paste0('Harmony - 33PCS'))
epi<- FindNeighbors(epi, reduction = "harmony", dims = 1:36)
epi<-RunUMAP(epi, dims=1:36, reduction= "harmony")

DimPlot(epi, group.by = 'sample') + labs(title = 'epi - Harmony 33_36')

setwd('Analysis/00_annotation_process/epi/')
dir.create('epi_reanalysis_harmony_33_36')
epi <- resolutions(epi,
                   workingdir = 'epi_reanalysis_harmony_33_36',
                   title = 'epi_reanalysis_harmony_33_36')
saveRDS(epi, file = 'epi_reanalysis_harmony_33_36/epi_reanalysis_33_36.RDS')
epi <- readRDS('Analysis/00_annotation_process/epi/epi_reanalysis_harmony_33_36/epi_reanalysis_33_36.RDS')


#
# 2nd time
#
epi <- readRDS('Analysis/00_annotation_process/epi/epi_reanalysis_harmony_33_36/epi_reanalysis_33_36.RDS')
DimPlot(epi, label=T)
epi <- epi[,!(epi$RNA_snn_res.1.5 %in% c(21, 23, 24, 25, 27))] #

counts <- epi@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 226
xx <-setdiff(rownames(epi), names(pp))
epi <- subset(epi, features = xx)
epi
# 24117 features across 8867 samples within 1 assay

epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi,selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi)
epi <- RunPCA(epi, npcs = 100)

PCS <- select_pcs(epi, 2)
PCS2 <- select_pcs(epi, 1.6)
ElbowPlot(epi, ndims = 100) +
  geom_vline(xintercept = 50, linetype = 'dotted') +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS, y=6, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - epi'))

epi<- FindNeighbors(epi,  dims = 1:50, reduction = 'pca')
epi<-RunUMAP(epi, dims=1:50, reduction = 'pca')

DimPlot(epi, group.by = 'sample') +
  labs(title = 'Tofacitinib 23 samples - 50PCs')

epi <- RunHarmony(epi, group.by = 'sample', dims.use = 1:50)
ElbowPlot(epi, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 19, linetype = 2) +
  geom_vline(xintercept = 35, linetype = 2) +
  labs(title = paste0('Harmony - 50PCS'))
epi<- FindNeighbors(epi, reduction = "harmony", dims = 1:35)
epi<-RunUMAP(epi, dims=1:35, reduction= "harmony")

DimPlot(epi, group.by = 'sample') + labs(title = 'epi - Harmony 50_35')

setwd('Analysis/00_annotation_process/epi/')
dir.create('epi_reanalysis_2_harmony_50_35')
epi <- resolutions(epi,
                   workingdir = 'epi_reanalysis_2_harmony_50_35',
                   title = 'epi_reanalysis_2_harmony_50_35')
saveRDS(epi, file = 'epi_reanalysis_2_harmony_50_35/epi_reanalysis_2_50_35.RDS')

#
# 3rd time
#
epi <- readRDS('Analysis/00_annotation_process/epi/epi_reanalysis_2_harmony_50_35/epi_reanalysis_2_50_35.RDS')
DimPlot(epi, label=T)
epi <- epi[,!(epi$RNA_snn_res.1.5 %in% c(27))] #

counts <- epi@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 23
xx <-setdiff(rownames(epi), names(pp))
epi <- subset(epi, features = xx)
epi
# 24094 features across 8844 samples within 1 assay

epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi,selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi)
epi <- RunPCA(epi, npcs = 100)

PCS <- select_pcs(epi, 2)
PCS2 <- select_pcs(epi, 1.6)
ElbowPlot(epi, ndims = 100) +
  geom_vline(xintercept = 50, linetype = 'dotted') +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS, y=6, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - epi'))

epi<- FindNeighbors(epi,  dims = 1:50, reduction = 'pca')
epi<-RunUMAP(epi, dims=1:50, reduction = 'pca')

DimPlot(epi, group.by = 'sample') +
  labs(title = 'Tofacitinib 23 samples - 50PCs')

epi <- RunHarmony(epi, group.by = 'sample', dims.use = 1:50)
ElbowPlot(epi, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 19, linetype = 2) +
  geom_vline(xintercept = 34, linetype = 2) +
  labs(title = paste0('Harmony - 50PCS'))
epi<- FindNeighbors(epi, reduction = "harmony", dims = 1:34)
epi<-RunUMAP(epi, dims=1:34, reduction= "harmony")

DimPlot(epi, group.by = 'sample') + labs(title = 'epi - Harmony 50_34')

setwd('Analysis/00_annotation_process/epi/')
dir.create('epi_reanalysis_3_harmony_50_34')
epi <- resolutions(epi,
                   workingdir = 'epi_reanalysis_3_harmony_50_34',
                   title = 'epi_reanalysis_3_harmony_50_34')
saveRDS(epi, file = 'epi_reanalysis_3_harmony_50_34/epi_reanalysis_3_50_34.RDS')
epi <- readRDS('Analysis/00_annotation_process/epi/epi_reanalysis_3_harmony_50_34/epi_reanalysis_3_50_34.RDS')
FeaturePlot(epi, features = c('KRT34', 'COL3A1', 'VWF', 'CHI3L1'), order =T, cols = c('lightgray', 'red'), )
FeaturePlot(epi, features = c('COL3A1'), order =T, cols = c('lightgray', 'red'), pt.size = 3 )

#
# 4th time
#
epi <- readRDS('Analysis/00_annotation_process/epi/epi_reanalysis_3_harmony_50_34/epi_reanalysis_3_50_34.RDS')
epi <- epi[,!(epi$RNA_snn_res.1.5 %in% c(25))]

counts <- epi@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 19
xx <-setdiff(rownames(epi), names(pp))
epi <- subset(epi, features = xx)
epi
# 24075 features across 8827 samples within 1 assay

epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi,selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi)
epi <- RunPCA(epi, npcs = 100)

PCS <- select_pcs(epi, 2)
PCS2 <- select_pcs(epi, 1.6)
ElbowPlot(epi, ndims = 100) +
  geom_vline(xintercept = 50, linetype = 'dotted') +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS, y=6, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - epi'))

epi<- FindNeighbors(epi,  dims = 1:50, reduction = 'pca')
epi<-RunUMAP(epi, dims=1:50, reduction = 'pca')

DimPlot(epi, group.by = 'sample') +
  labs(title = 'Tofacitinib 23 samples - 50PCs')

epi <- RunHarmony(epi, group.by = 'sample', dims.use = 1:50)
ElbowPlot(epi, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 19, linetype = 2) +
  geom_vline(xintercept = 25, linetype = 2) +
  labs(title = paste0('Harmony - 50PCS'))
epi<- FindNeighbors(epi, reduction = "harmony", dims = 1:25)
epi<-RunUMAP(epi, dims=1:25, reduction= "harmony")

DimPlot(epi, group.by = 'sample') + labs(title = 'epi - Harmony 50_25')

setwd('Analysis/00_annotation_process/epi/')
dir.create('epi_reanalysis_4_harmony_50_25')
epi <- resolutions(epi,
                   workingdir = 'epi_reanalysis_4_harmony_50_25',
                   title = 'epi_reanalysis_4_harmony_50_25')
saveRDS(epi, file = 'epi_reanalysis_4_harmony_50_25/epi_reanalysis_4_50_25.RDS')
epi <- readRDS('Analysis/00_annotation_process/epi/epi_reanalysis_4_harmony_50_25/epi_reanalysis_4_50_25.RDS')
FeaturePlot(epi, features = c('KRT25', 'COL3A1', 'VWF', 'CHI3L1'), order =T, cols = c('lightgray', 'red'), )
FeaturePlot(epi, features = c('COL3A1'), order =T, cols = c('lightgray', 'red'), pt.size = 3 )


#
# 5th time
#

epi <- readRDS('Analysis/00_annotation_process/epi/epi_reanalysis_4_harmony_50_25/epi_reanalysis_4_50_25.RDS')
DimPlot(epi, label=T)
epi <- epi[,!(epi$RNA_snn_res.1.5 %in% c(25))] #

counts <- epi@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 8
xx <-setdiff(rownames(epi), names(pp))
epi <- subset(epi, features = xx)
epi
# 24067 features across 8788 samples within 1 assay

epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi,selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi)
epi <- RunPCA(epi, npcs = 100)

PCS <- select_pcs(epi, 2)
PCS2 <- select_pcs(epi, 1.6)
ElbowPlot(epi, ndims = 100) +
  geom_vline(xintercept = 50, linetype = 'dotted') +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS, y=6, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - epi'))

epi<- FindNeighbors(epi,  dims = 1:28, reduction = 'pca')
epi<-RunUMAP(epi, dims=1:28, reduction = 'pca')

DimPlot(epi, group.by = 'sample') +
  labs(title = 'Tofacitinib 23 samples - 28PCs')

epi <- RunHarmony(epi, group.by = 'sample', dims.use = 1:28)
ElbowPlot(epi, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 19, linetype = 2) +
  geom_vline(xintercept = 25, linetype = 2) +
  labs(title = paste0('Harmony - 28PCS'))
epi<- FindNeighbors(epi, reduction = "harmony", dims = 1:25)
epi<-RunUMAP(epi, dims=1:25, reduction= "harmony")

DimPlot(epi, group.by = 'sample') + labs(title = 'epi - Harmony 28_25')

setwd('Analysis/00_annotation_process/epi/')
dir.create('epi_reanalysis_5_harmony_28_25')
epi <- resolutions(epi,
                   workingdir = 'epi_reanalysis_5_harmony_28_25',
                   title = 'epi_reanalysis_5_harmony_28_25')
saveRDS(epi, file = 'epi_reanalysis_5_harmony_28_25/epi_reanalysis_5_28_25.RDS')
