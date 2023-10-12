options(stringsAsFactors = FALSE)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readr)
library(dplyr)

source('Figures/functions_plots.R')

## Data ------------------------------------------------------------------------
message('Loading data')

myeloids <- readRDS('/home/acorraliza/TOFA_data/20220222_TOFAS_23/00_annotation_process/00_anotadas/myeloids.RDS')
myeloids$pre_post <- plyr::mapvalues(myeloids$week_3, from = c('W0', 'POST'), to = c('PRE', 'POST'))
myeloids$annotation_intermediate <- gsub('Macrophage NRG1', 'IDA macrophages', myeloids$annotation_intermediate)
myeloids$pre_post <- factor(myeloids$pre_post, levels = c('PRE', 'POST'))
myeloids$response <- factor(myeloids$response, levels = c('R', 'NR'))
myeloids$subset <- factor(myeloids$subset, levels = c('epi', 'stroma', 'plasmas', 'myeloids', 'cycling', 'tcells'))

## Figure 3A  ------------------------------------------------------------------
# UMAP split by response and treatment, with some celltypes in colour
# HS, Inflammatory monocytes, M0, M1, M2, IDA macrophages, Neutrophils
