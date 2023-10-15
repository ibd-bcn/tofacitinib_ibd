#
# Libraries
#
message('Loading libraries')

library(Seurat)
library(plyr)
library(ggplot2)
library(DropletUtils)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scDblFinder)
library(viridis)

#
# Extra functions
#
message('Loading functions')
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
source('~/TOFA_data/functions_scrseq.R')
setwd('~/TOFA_data/output_cellranger/')

## 01. Loading raw data --------------------------------------------------------
messaage('01. Loading raw data')

list.files <-c("TOF-005-W0",
             "TOF-005-W8","TOF-005-W24","TOF-009-W0",
             "TOF-009-W48","TOF-010-W0","TOF-011-W0",
             "TOF-011-W8","TOF-012-W0","TOF-012-W8","TOF-012-W16",
             "TOF-013-W0","TOF-013-W8","TOF-015-W0",
             "TOF-015-W8","TOF-016-W0","TOF-016-W8","TOF-019-W0",
             "TOF-019-W8","TOF-019-W16","TOF-022-W0",
             "TOF-023-W0","TOF-023-W8")
meta <- data.frame(
  stringsAsFactors = FALSE,
  PATIENT = c("TOF_005",
              "TOF_005","TOF_005","TOF_009","TOF_009",
              "TOF_010","TOF_011","TOF_011","TOF_012","TOF_012",
              "TOF_012","TOF_013","TOF_013","TOF_015",
              "TOF_015","TOF_016","TOF_016","TOF_019",
              "TOF_019","TOF_019","TOF_022","TOF_023","TOF_023"),
  TIMEPOINT = c("W0","W8",
                "W24","W0","w48","W0","W0","W8","W0","W8",
                "W16","W0","W8","W0","W8","W0","W8","W0",
                "W8","W16","W0","W0","W8"),
  week_2 = c("W0","W8",
             "W16","W0","W16","W0","W0","W8","W0","W8",
             "W16","W0","W8","W0","W8","W0","W8","W0",
             "W8","W16","W0","W0","W8"),
  SAMPLE = c("TOF-005-W0",
             "TOF-005-W8","TOF-005-W24","TOF-009-W0",
             "TOF-009-W48","TOF-010-W0","TOF-011-W0",
             "TOF-011-W8","TOF-012-W0","TOF-012-W8","TOF-012-W16",
             "TOF-013-W0","TOF-013-W8","TOF-015-W0",
             "TOF-015-W8","TOF-016-W0","TOF-016-W8","TOF-019-W0",
             "TOF-019-W8","TOF-019-W16","TOF-022-W0",
             "TOF-023-W0","TOF-023-W8"),
  Response = c("NR","NR",
               "NR","R","R","NR","R","R","R","R","R",
               "NR","NR","R","R","NR","NR","NR","NR","NR",
               "NR","R","R"),
  BIOPSIES.LOCATION = c("Recto",
                        "Recto","Recto-Sigma","Recto","Recto",
                        "Recto-Sigma","Sigma","Sigma","Sigma","Sigma","Sigma",
                        "Sigma","Sigma","Sigma","Sigma","Recto",
                        "Recto","Sigma","Sigma","Sigma",
                        "Descendent colon","Sigma","Sigma")
)

list_data <- list()

for(i in list_files){
  print(i)
  sce2 <- Read10X(i)

  keep_feature <- rowSums(sce2 > 0) > 0
  sce2 <- sce2[keep_feature, ]

  sce2 <- scDblFinder(sce2, verbose=FALSE)

  m4 <- data.frame('sample' = rep(i, ncol(sce2)),
                   'doublet' = sce2$scDblFinder.class,
                   'week' = rep(meta[meta$SAMPLE == i, 'TIMEPOINT'][1], ncol(sce2)),
                   'week_2' = rep(meta[meta$SAMPLE == i, 'week_2'][1], ncol(sce2)),
                   'patient' = rep(meta[meta$SAMPLE == i, 'PATIENT'][1], ncol(sce2)),
                   'BIOPSIES_LOCATION' = rep(meta[meta$SAMPLE == i, 'BIOPSIES.LOCATION'][1], ncol(sce2)),
                   'response' = rep(meta[meta$SAMPLE == i, 'Response'][1], ncol(sce2))
  )

  colnames(sce2) <- paste(i, colnames(sce2), sep='_')
  rownames(m4) <- colnames(sce2)


  data <- CreateSeuratObject(
    counts(sce2),
    min.features = 100,
    project = i,
    assay = "RNA",
    meta.data = m4
  )

  data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = "^MT-")

  list_data[[i]] <- data
}

seudata <- list_data[[1]]
for(i in 2:length(list_files)){
  print(i)
  seudata <- merge(seudata, list_data[[i]])
}

dir.create('~/TOFA_data/20220222_TOFAS_23')
setwd('~/TOFA_data/20220222_TOFAS_23')
VlnPlot(seudata, features = c('percent.mt', 'nFeature_RNA'), group.by = 'sample')
VlnPlot(seudata, features = c( 'nCount_RNA'), group.by = 'sample')

# doublet singlet
#    6679   82700

seudata <- seudata[,seudata$doublet == 'singlet']
counts <- seudata@assays$RNA@counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 220
xx <-setdiff(rownames(seudata), names(pp))
seudata <- subset(seudata, features = xx)
seudata


meta <- seudata@meta.data
meta$density <- get_density(meta$percent.mt, meta$nFeature_RNA, n = 100)
jpeg(filename = '~/TOFA_data/20220222_TOFAS_23/density.jpeg', width = 1500, height = 1500, res = 150)
ggplot(meta) +
  geom_point(aes(percent.mt, nFeature_RNA, color = density)) +
  geom_hline(yintercept = 100, color = 'gray', linetype="dashed") +
  geom_vline(xintercept = 50, color = 'gray', linetype="dashed") +
  geom_vline(xintercept = 25, color = 'gray', linetype="dashed") +
  scale_color_viridis() + theme_classic() + theme(plot.title = element_text(size = 25))+
  labs(title = 'Tofacitinib - 23 samples together')
dev.off()

meta$density <- get_density(meta$nCount_RNA, meta$nFeature_RNA, n = 1000)
jpeg(filename = '~/TOFA_data/20220222_TOFAS_23/feature_counts.jpeg', width = 1500, height = 1500, res = 150)
ggplot(meta) +
  geom_point(aes(nCount_RNA, nFeature_RNA, color = density)) +
  scale_color_viridis() +
  theme_classic() +
  theme(plot.title = element_text(size = 25))+
  labs(title = 'Tofacitinib - 23 samples together')
dev.off()

seudata_f <- seudata[, seudata$percent.mt < 50 & seudata$nFeature_RNA > 100]
saveRDS(seudata_f,'~/TOFA_data/20220222_TOFAS_23/seudata_f.RDS')
saveRDS(seudata,'~/TOFA_data/20220222_TOFAS_23/seudata.RDS')
