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
# Data processing --------------------------------------------------------------
#
message('Data processing')
seudata_f <- readRDS('~/TOFA_data/TOFAS_21_SAMPLES/seudata_f.RDS')
seudata_f <- seurat_to_pca(seudata_f)

PCS <- select_pcs(seudata_f, 2)
PCS2 <- select_pcs(seudata_f, 1.6) # 35
p <- ElbowPlot(seudata_f, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('Tofacitinib - merge 23 samples'))

seudata_f <- FindNeighbors(seudata_f,  dims = 1:PCS2, reduction = 'pca')
seudata_f <- RunUMAP(seudata_f, dims=1:PCS2, reduction = 'pca')

# plot the output UMAP
DimPlot(seudata_f, group.by = 'sample') + labs(title = 'Tofacitinib 23 samples - 35PCs')

# explore gene expression
FeaturePlot(seudata_f, features = c('EPCAM', 'CD3E', 'C1QB', 'DERL3'), order =T, cols = c('lightgray', 'red'))
FeaturePlot(seudata_f, features = c('LYZ', 'C1QA', 'C1QB', 'AIF1'), order =T, cols = c('lightgray', 'red'))
FeaturePlot(seudata_f, features = c('FCGR3B', 'CMTM2', 'MS4A1', 'MKI67'), order =T, cols = c('lightgray', 'red'))
FeaturePlot(seudata_f, features = c('KRT17', 'COL3A1', 'VWF', 'CHGA', 'TPSAB1', 'CHI3L1'), order =T, cols = c('lightgray', 'red'))
FeaturePlot(seudata_f, features = c('percent.mt'), order =T, cols = c('lightgray', 'red'))

message('Saving the analysis')
saveRDS(seudata_f,'Analysis/data/seudata_f.RDS')

#
# Louvain clustering -----------------------------------------------------------
#
message('Louvain clustering')
if(!exists('seudata_f')){seudata_f <- readRDS('data/seudata_f.RDS')}

dir.create('Analysis/data/00_annotation_process')
dir.create('Analysis/data/00_annotation_process/together_filter50')
seudata_f <- resolutions(seudata_f, resolutions = c(0.1,0.3,0.5),
                         workingdir = 'data/00_annotation_process/together_filter50/',
                         title = '23 samples filt50')

DimPlot(seudata_f, label=T, group.by = 'RNA_snn_res.0.5')
saveRDS(seudata_f,'Analysis/data/00_annotation_process/together_filter50/seudata_f.RDS')

#
# Splitting the cells into subsets ---------------------------------------------
#
if(!exists('seudata_f')){seudata_f <- readRDS('Analysis/data/00_annotation_process/together_filter50/seudata_f.RDS')}

# annotation based on cluster markers:
message('Cluster annotation and splitting')
subsets <- tibble::tribble(
  ~cluster,    ~subset,
  0L,   "tcells",
  1L,  "plasmas",
  2L,      "epi",
  3L,  "plasmas",
  4L,  "plasmas",
  5L,   "stroma",
  6L,   "tcells",
  7L,  "plasmas",
  8L,  "plasmas",
  9L, "myeloids",
  10L, "myeloids",
  11L, "myeloids",
  12L,  "cycling",
  13L,  "plasmas",
  14L,   "tcells",
  15L,      "epi",
  16L,   "stroma",
  17L, "myeloids",
  18L,      "epi",
  19L,   "stroma",
  20L,   "stroma",
  21L,      "epi",
  22L,  "cycling"
)

tcells <- seudata_f[,seudata_f$RNA_snn_res.0.5 %in% subsets$cluster[subsets$subset == 'tcells']]
plasmas <- seudata_f[,seudata_f$RNA_snn_res.0.5 %in% subsets$cluster[subsets$subset == 'plasmas']]
myeloids <- seudata_f[,seudata_f$RNA_snn_res.0.5 %in% subsets$cluster[subsets$subset == 'myeloids']]
epi <- seudata_f[,seudata_f$RNA_snn_res.0.5 %in% subsets$cluster[subsets$subset == 'epi']]
stroma <- seudata_f[,seudata_f$RNA_snn_res.0.5 %in% subsets$cluster[subsets$subset == 'stroma']]
cycling <- seudata_f[,seudata_f$RNA_snn_res.0.5 %in% subsets$cluster[subsets$subset == 'cycling']]

saveRDS(tcells, file = 'Analysis/data/00_annotation_process/tcells.RDS')
saveRDS(plasmas, file = 'Analysis/data/00_annotation_process/plasmas.RDS')
saveRDS(tcells, file = 'Analysis/data/00_annotation_process/tcells.RDS')
saveRDS(tcells, file = 'Analysis/data/00_annotation_process/tcells.RDS')
saveRDS(tcells, file = 'Analysis/data/00_annotation_process/tcells.RDS')
saveRDS(tcells, file = 'Analysis/data/00_annotation_process/tcells.RDS')

