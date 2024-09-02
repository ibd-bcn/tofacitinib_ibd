print(Sys.time())

library(Seurat)
library(plyr)
library(ggplot2)
library(DropletUtils)
library(beepr)
library(celda)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scDblFinder)
library(viridis)
library(MASS)
library(harmony)


#
# working in paral·lel
#
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize= 9891289600)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
source('~/TOFA_data/functions_scrseq.R')


plasmas <- readRDS(file = 'Analysis/data/00_annotation_process/00_anotadas/plasmas.RDS')
tcells <- readRDS(file = 'Analysis/data/00_annotation_process/00_anotadas/tcells.RDS')
epi <- readRDS(file = 'Analysis/data/00_annotation_process/00_anotadas/epi.RDS')
myeloids <- readRDS(file = 'Analysis/data/00_annotation_process/00_anotadas/myeloids.RDS')
myeloids$annotation_2 <- myeloids$annotation_refined
myeloids$annotation_2[myeloids$annotation_2 %in% clusters_Azu] <- 'Macrophages'

stroma <- readRDS(file = 'Analysis/data/00_annotation_process/00_anotadas/stroma.RDS')
cycling <- readRDS(file = 'Analysis/data/00_annotation_process/00_anotadas/cycling.RDS')



listdata <- list(epi,
                 stroma,
                 plasmas,
                 tcells,
                 cycling,
                 myeloids)

names(listdata) <- c('epi',
                     'stroma',
                     'plasmas',
                     'tcells',
                     'cycling',
                     'myeloids')


#### loop DE


dir.create('/Analysis/temp/')
for(name in names(listdata)){
  data <- listdata[[name]]
  DES <- NULL
  rm(dataset_de, dataset_de0, dataset_de16, dataset_de8, dataset_deNR, dataset_deR, DE)
  for( ANOT in c('annotation_intermediate', 'annotation_refined', 'annotation_2')){
    if(ANOT %in% colnames(data@meta.data)){
      print(ANOT)
      for(cluster in levels(as.factor(data@meta.data[,ANOT]))){
        print(cluster)
        dataset_de <- SetIdent(data, value = ANOT)
        dataset_de <- dataset_de[,dataset_de@active.ident == cluster]
        dataset_de <- SetIdent(dataset_de, value = 'week_2')
        if('W0' %in% dataset_de@active.ident  & 'W8' %in% dataset_de@active.ident){
          if(sum(dataset_de@active.ident == 'W0') > 3 &
             sum(dataset_de@active.ident == 'W8') > 3){
            DE <- FindMarkers(dataset_de, ident.1 = 'W8',
                              ident.2 ='W0', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            head(DE)
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w0_vs_w8'
            DE$cluster <- cluster
            DE$annotation <- ANOT

            DES <- rbind(DES, DE)

            print('w0_vs_w8')
          }
        }

        if('W0' %in% dataset_de@active.ident  & 'W16' %in% dataset_de@active.ident){
          if(sum(dataset_de@active.ident == 'W0') > 3 & sum(dataset_de@active.ident == 'W16') > 3){
            DE <- FindMarkers(dataset_de, ident.1 = 'W16',
                              ident.2 ='W0', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w0_vs_w16'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w0_vs_w16')
          }
        }


        if('W8' %in% dataset_de@active.ident  & 'W16' %in% dataset_de@active.ident){
          if(sum(dataset_de@active.ident == 'W8') > 3 & sum(dataset_de@active.ident == 'W16') > 3){
            DE <- FindMarkers(dataset_de, ident.1 = 'W16',
                              ident.2 ='W8', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w8_vs_w16'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w8_vs_w16')
          }
        }

        #
        # dataset_deR
        #

        dataset_deR <- dataset_de[,dataset_de$response == 'R']
        if('W0' %in% dataset_deR@active.ident  & 'W8' %in% dataset_deR@active.ident){
          if(sum(dataset_deR@active.ident == 'W0') > 3 & sum(dataset_deR@active.ident == 'W8') > 3){
            DE <- FindMarkers(dataset_deR, ident.1 = 'W8',
                              ident.2 ='W0', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w0R_vs_w8R'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w0R_vs_w8R')
          }
        }

        if('W0' %in% dataset_deR@active.ident  & 'W16' %in% dataset_deR@active.ident){
          if(sum(dataset_deR@active.ident == 'W0') > 3 & sum(dataset_deR@active.ident == 'W16') > 3){
            DE <- FindMarkers(dataset_deR, ident.1 = 'W16',
                              ident.2 ='W0', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w0R_vs_w16R'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w0R_vs_w16R')
          }
        }

        if('W8' %in% dataset_deR@active.ident  & 'W16' %in% dataset_deR@active.ident){
          if(sum(dataset_deR@active.ident == 'W8') > 3 & sum(dataset_deR@active.ident == 'W16') > 3){
            DE <- FindMarkers(dataset_deR, ident.1 = 'W16',
                              ident.2 ='W8', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w8R_vs_w16R'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w8R_vs_w16R')
          }
        }

        dataset_deNR <- dataset_de[,dataset_de$response == 'NR']

        if('W0' %in% dataset_deNR@active.ident  & 'W8' %in% dataset_deNR@active.ident){
          if(sum(dataset_deNR@active.ident == 'W0') > 3 & sum(dataset_deNR@active.ident == 'W8') > 3){
            DE <- FindMarkers(dataset_deNR, ident.1 = 'W8',
                              ident.2 ='W0', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w0NR_vs_w8NR'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w0NR_vs_w8NR')
          }
        }

        if('W0' %in% dataset_deNR@active.ident  & 'W16' %in% dataset_deNR@active.ident){
          if(sum(dataset_deNR@active.ident == 'W0') > 3 & sum(dataset_deNR@active.ident == 'W16') > 3){
            DE <- FindMarkers(dataset_deNR, ident.1 = 'W16',
                              ident.2 ='W0', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w0NR_vs_w16NR'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w0NR_vs_w16NR')
          }
        }


        if('W8' %in% dataset_deNR@active.ident  & 'W16' %in% dataset_deNR@active.ident){
          if(sum(dataset_deNR@active.ident == 'W8') > 3 & sum(dataset_deNR@active.ident == 'W16') > 3){
            DE <- FindMarkers(dataset_deNR, ident.1 = 'W16',
                              ident.2 ='W8', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w8NR_vs_w16NR'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w8NR_vs_w16NR')
          }
        }

        dataset_de <- SetIdent(data, value = ANOT)
        dataset_de <- dataset_de[,dataset_de@active.ident == cluster]
        dataset_de <- SetIdent(dataset_de, value = 'week_2')
        dataset_de0 <- dataset_de[,dataset_de$week_2 == 'W0']
        dataset_de0 <- SetIdent(dataset_de0, value = 'response')

        if('NR' %in% dataset_de0@active.ident  & 'R' %in% dataset_de0@active.ident){
          if(sum(dataset_de0@active.ident == 'NR') > 3 & sum(dataset_de0@active.ident == 'R') > 3){
            DE <- FindMarkers(dataset_de0, ident.1 = 'NR',
                              ident.2 ='R', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'Rw0_NRw0'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('Rw0_NRw0')
          }}

        dataset_de8 <- dataset_de[,dataset_de$week_2 == 'W8']
        dataset_de8 <- SetIdent(dataset_de8, value = 'response')

        if('NR' %in% dataset_de8@active.ident  & 'R' %in% dataset_de8@active.ident){
          if(sum(dataset_de8@active.ident == 'NR') > 3 & sum(dataset_de8@active.ident == 'R') > 3){
            DE <- FindMarkers(dataset_de8, ident.1 = 'NR',
                              ident.2 ='R', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'Rw8_NRw8'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('Rw8_NRw8')
          }}

        if(sum(dataset_de$week_2 == 'W16')>3){
          dataset_de16 <- dataset_de[,dataset_de$week_2 == 'W16']
          dataset_de16 <- SetIdent(dataset_de16, value = 'response')

          if('NR' %in% dataset_de16@active.ident  &
             'R' %in% dataset_de16@active.ident){
            if(sum(dataset_de16@active.ident == 'NR') > 3 &
               sum(dataset_de16@active.ident == 'R') > 3){
              DE <- FindMarkers(dataset_de16, ident.1 = 'NR',
                                ident.2 ='R', logfc.threshold = 0.01,
                                verbose = FALSE,
                                test.use = 'wilcox')
              DE$gene <- rownames(DE)
              DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
              DE$comp <- 'Rw16_NRw16'
              DE$cluster <- cluster
              DE$annotation <- ANOT
              DES <- rbind(DES, DE)
              print('Rw16_NRw16')
            }

          }

        }

        #
        # week_3
        #
        dataset_de <- SetIdent(dataset_de, value = 'week_3')
        if('W0' %in% dataset_de@active.ident  & 'POST' %in% dataset_de@active.ident){
          if(sum(dataset_de@active.ident == 'W0') > 3 &
             sum(dataset_de@active.ident == 'POST') > 3){
            DE <- FindMarkers(dataset_de, ident.1 = 'POST',
                              ident.2 ='W0', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            head(DE)
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w0_vs_POST'
            DE$cluster <- cluster
            DE$annotation <- ANOT

            DES <- rbind(DES, DE)

            print('w0_vs_POST')
          }
        }
        dataset_deR <- dataset_de[,dataset_de$response == 'R']
        if('W0' %in% dataset_deR@active.ident  & 'POST' %in% dataset_deR@active.ident){
          if(sum(dataset_deR@active.ident == 'W0') > 3 & sum(dataset_deR@active.ident == 'POST') > 3){
            DE <- FindMarkers(dataset_deR, ident.1 = 'POST',
                              ident.2 ='W0', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w0R_vs_POSTR'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w0R_vs_POSTR')
          }
        }

        dataset_deNR <- dataset_de[,dataset_de$response == 'NR']

        if('W0' %in% dataset_deNR@active.ident  & 'POST' %in% dataset_deNR@active.ident){
          if(sum(dataset_deNR@active.ident == 'W0') > 3 & sum(dataset_deNR@active.ident == 'POST') > 3){
            DE <- FindMarkers(dataset_deNR, ident.1 = 'POST',
                              ident.2 ='W0', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'w0NR_vs_POSTNR'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('w0NR_vs_POSTNR')
          }
        }

        dataset_depost <- dataset_de[,dataset_de$week_3 == 'POST']
        dataset_depost <- SetIdent(dataset_depost, value = 'response')

        if('NR' %in% dataset_depost@active.ident  & 'R' %in% dataset_depost@active.ident){
          if(sum(dataset_depost@active.ident == 'NR') > 3 & sum(dataset_depost@active.ident == 'R') > 3){
            DE <- FindMarkers(dataset_depost, ident.1 = 'NR',
                              ident.2 ='R', logfc.threshold = 0.01,
                              verbose = FALSE,
                              test.use = 'wilcox')
            DE$gene <- rownames(DE)
            DE$sign <- funfun(DE$p_val, DE$avg_log2FC, DE$p_val_adj)
            DE$comp <- 'RPOST_NRPOST'
            DE$cluster <- cluster
            DE$annotation <- ANOT
            DES <- rbind(DES, DE)
            print('RPOST_NRPOST')
          }}

        print('next cluster')
      }
      print('ALL CLUSTERS DONE! WRITE TABLE!')
      write.table(DES, sep = ';', dec = ',', row.names = F, col.names = T,
                  file = paste('Analysis/temp/',name, ANOT, '.csv', sep='_'))
      DES <- NULL
    }
  }
}

# DE together unique file -----------------------------------------------

setwd("Analysis/temp/")

files <-list.files(path = "Analysis/temp/",
                   full.names = T)
complete <- data.frame()
for(file in files){
  print(file)
  subset <- stringr::str_split_fixed(file, '_',9)[,6]
  annotation <- stringr::str_split_fixed(file, '_',9)[,8]
  a <- read.delim(file, header=T, sep = ';', dec = ',')
  a$subset <- subset
  a$annotation <- paste0('annotation_',annotation)
  complete <- rbind(complete, a)
}

write.table(complete, file = 'Figures/extra_data/new_complete_DE.csv',
            sep = ';', dec = ',', row.names = F, col.names = T)
saveRDS(complete, file = 'Figures/extra_data/new_complete.RDS')

print(Sys.time())
