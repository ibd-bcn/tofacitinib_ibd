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
tcells <- readRDS(file = 'Analysis/data/00_annotation_process/tcells.RDS')

FeaturePlot(tcells, features= c('MS4A1', 'DERL3', 'EPCAM', 'FCGR3B', 'HBB'), order = T)

#
# remove cells with genes from other subsets
#
counts <- tcells@assays$RNA@counts
p <- grep("MS4A1$|DERL3$|EPCAM$|^HBB$|^FCGR3B$",rownames(tcells))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 786
xx <-setdiff(colnames(tcells), names(pp))
tcells <- subset(tcells,cells = xx)
tcells
#28451 features across 15833 samples within 1 assay


#
# remove genes from IGs
#
gg <- rownames(tcells)[c(grep("^IGH",rownames(tcells)),
                         grep("^IGK", rownames(tcells)),
                         grep("^IGL", rownames(tcells)))]
genes <- setdiff(rownames(tcells),gg)
tcells <- subset(tcells,features = genes)
tcells
# 28185 features across 15833 samples within 1 assay


#
# Filter by %MT genes
#
VlnPlot(tcells, features = 'percent.mt')
tcells <- tcells[,tcells$percent.mt < 25]
# 28185 features across 15256 samples within 1 assay

#
# genes with no expression out
#
counts <- tcells@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp) #6406
xx <-setdiff(rownames(tcells), names(pp))

tcells <- subset(tcells, features = xx)
tcells
# 21779 features across 15256 samples within 1 assay


#
# Normalization, scaling and UMAP generation -----------------------------------
#
tcells <- NormalizeData(tcells)
tcells <- FindVariableFeatures(tcells,selection.method = "vst", nfeatures = 2000)
tcells <- ScaleData(tcells)
tcells <- RunPCA(tcells, npcs = 100)

PCS <- select_pcs(tcells, 2)
PCS2 <- select_pcs(tcells, 1.6)
ElbowPlot(tcells, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = 40, linetype = 2) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - tcells'))

tcells<- FindNeighbors(tcells,  dims = 1:40, reduction = 'pca')
tcells<-RunUMAP(tcells, dims=1:40, reduction = 'pca')

DimPlot(tcells, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - 40PCs')

#
# Louvain clustering without batch correction ----------------------------------
#
dir.create('Analysis/00_annotation_process/tcells')
setwd('Analysis/00_annotation_process/tcells')
dir.create('tcells_no_harmony')
tcells <- resolutions(tcells,
                      workingdir = 'tcells_no_harmony',
                      title = 'tcells_no_harmony')
saveRDS(tcells, file = 'tcells_no_harmony/tcells_filtered_40.RDS')
tcells <- readRDS( 'tcells_no_harmony/tcells_filtered_40.RDS')
FeaturePlot(tcells, features = c('EPCAM', 'percent.mt', 'AQP8', 'OLFM4'), order =T, cols = c('lightgray', 'red'), )


#
# Louvain clustering with batch correction -------------------------------------
#
tcells <- RunHarmony(tcells, group.by = 'sample', dims.use = 1:40)
ElbowPlot(tcells, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 32, linetype = 2) +
  labs(title = paste0('Harmony - 40PCS'))
tcells<- FindNeighbors(tcells, reduction = "harmony", dims = 1:32)
tcells<-RunUMAP(tcells, dims=1:32, reduction= "harmony")

DimPlot(tcells, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - tcells - Harmony 40-32')


FeaturePlot(tcells, features= c('CD3E','CD3D', 'MS4A1', 'DERL3'), order = T)

dir.create('tcells_harmony_40_32/')
tcells <- resolutions(tcells,
                      workingdir = 'tcells_harmony_40_32',
                      title = 'tcells_harmony_40_32')
saveRDS(tcells, file = 'tcells_harmony_40_32/tcells_filtered_40_32.RDS')
tcells <- readRDS( 'tcells_harmony_40_32/tcells_filtered_40_32.RDS')
FeaturePlot(tcells, features = c('percent.mt', 'C1QA', 'NRG1', 'FCGR3B'), order =T, cols = c('lightgray', 'red'), )
DimPlot(tcells, group.by='sample' )

a <- DimPlot(tcells, group.by = 'RNA_snn_res.0.1', label=T) + NoLegend()
b <- DimPlot(tcells, group.by = 'RNA_snn_res.0.3', label=T) + NoLegend()
c <- DimPlot(tcells, group.by = 'RNA_snn_res.0.5', label=T) + NoLegend()
d <- DimPlot(tcells, group.by = 'RNA_snn_res.0.7', label=T) + NoLegend()
e <- DimPlot(tcells, group.by = 'RNA_snn_res.0.9', label=T) + NoLegend()
f <- DimPlot(tcells, group.by = 'RNA_snn_res.1.1', label=T) + NoLegend()
g <- DimPlot(tcells, group.by = 'RNA_snn_res.1.3', label=T) + NoLegend()
h <- DimPlot(tcells, group.by = 'RNA_snn_res.1.5', label=T) + NoLegend()

a+b+c+d
e+f+g+h
