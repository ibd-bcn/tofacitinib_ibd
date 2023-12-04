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
# Data loading -----------------------------------------------------------------
#
annotation <- data.frame(
  stringsAsFactors = FALSE,
  check.names = FALSE,
  Cluster = c(0L,1L,2L,
              3L,4L,5L,6L,7L,
              8L,9L,10L,11L,
              12L,13L,14L,15L,
              16L,17L,18L,19L,
              20L,21L),
  `Together_myeloids_30_25_1,3_intermediate` = c("Neutrophils","Mast",
                                                 "Neutrophils","M1","M2",
                                                 "M0",
                                                 "Macrophage NRG1",
                                                 "Inflammatory monocytes","M0",
                                                 "Rib hi myeloids",
                                                 "Neutrophils","HS","DCs",
                                                 "DCs","DCs",
                                                 "Cycling myeloids",
                                                 "Eosinophils","pDC",NA,
                                                 NA,NA,NA),
  `Plasmas_25_22_1,3_intermediate` = c("Plasmablast_IgG","PC_IgG",
                                       "B_cell",
                                       "Plasmablast_IgG","PC_IgA",
                                       "B_cell",
                                       "Plasmablast_IgA",
                                       "Plasmablast_IgG","PC_IGLL5",
                                       "Cycling_B_cell",
                                       "PC_IGLV6-57",
                                       "Plasmablast_IgA",
                                       "PC_heat_shock","PC_IER",
                                       "PC_IgA","PC_PSAT1",
                                       "B_cell",
                                       "Plasmablast_IGKC","B_cell",'PC_IFIT1',
                                       NA,NA),
  `Stroma_harmony_clean_41_27_1,3_intermediate` = c("Endothelium","IER_fibroblasts",
                                                    "S3","S1","S2",
                                                    "Inflammatory_fibroblasts",
                                                    "Endothelium","MT fibroblasts",
                                                    "S1","Perycites",
                                                    "Inflammatory_fibroblasts",
                                                    "Myofibroblasts","S2","Glia",
                                                    "Perycites",
                                                    "Cycling_fibroblasts",
                                                    "Endothelium",NA,NA,
                                                    NA,NA,NA),
  `epi_clean_5_harmony_28_25_1,1_intermediate` = c("Colonocytes","Colonocytes",
                                                   "Cycling TA",
                                                   "Undifferentiated epithelium","IER Epithelium",
                                                   "Colonocytes",
                                                   "Secretory progenitor",
                                                   "Goblet",
                                                   "Epithlium Rib hi",
                                                   "Cycling TA","Colonocytes",
                                                   "Secretory progenitor",
                                                   "Epithlium Rib hi","Stem",
                                                   "Epithlium Rib hi",
                                                   "Mature goblet","Tuft",
                                                   "Colonocytes",
                                                   "Colonocytes","Colonocytes",
                                                   "Enteroendocrines",
                                                   "APOA4"),
  `tcells_harmony_40_32_1,3_intermediate` = c("CD8","CD4",
                                              "Ribhi T cells",
                                              "Tregs","Thf",
                                              "Naïve T cells","Tregs",
                                              "Ribhi T cells",
                                              "Th17","Thf","CD8",
                                              "CD4","MT hi IER",
                                              "HS","NK","CD8",
                                              "CD4 CD8 IFIT3",
                                              "CD8","CD4","ILC3",
                                              "DN EOMES",NA),
  Together_cycling_filt25_20_23_intermediate = c("Cycling_B_cell","Cycling_PC",
                                                 "Cycling_PC",
                                                 "Cycling_T_cell",
                                                 "Cycling_B_cell",
                                                 "Cycling_T_cell","Cycling_PC",
                                                 "Cycling_Myeloid",
                                                 NA,NA,NA,NA,NA,
                                                 NA,NA,NA,NA,NA,
                                                 NA,NA,NA,NA)
)

myeloids <- readRDS( 'Analysis/data/00_annotation_process/myeloids/myeloids_harmony_30_25/myeloids_filtered_30_25.RDS')
plasmas <- readRDS( 'Analysis/data/00_annotation_process/plasmas/plasmas_reanalysis_harmony_25_22/plasmas_reanalysis_harmony_25_22.RDS')
stroma <- readRDS( 'Analysis/data/00_annotation_process/stroma/clean_harmony_41_27/clean_harmony_41_27.RDS')
epi <- readRDS( 'Analysis/data/00_annotation_process/epi/epi_reanalysis_5_harmony_28_25/epi_reanalysis_5_28_25.RDS')
tcells <- readRDS( 'Analysis/data/00_annotation_process/tcells/tcells_harmony_40_32/tcells_filtered_40_32.RDS')
cycling <- readRDS( 'Analysis/data/00_annotation_process/cycling/cycling_harmony_20_23/cycling_filtered_20_23.RDS')

#
# Data annotation --------------------------------------------------------------
#
plasmas$annotation_intermediate <- mapvalues(x = plasmas$RNA_snn_res.1.3,
                                             from = annotation$Cluster,
                                             to = annotation$`Plasmas_25_22_1,3_intermediate`)
tcells$annotation_intermediate <- mapvalues(x = tcells$RNA_snn_res.1.3,
                                            from = annotation$Cluster,
                                            to = annotation$`tcells_harmony_40_32_1,3_intermediate`)
epi$annotation_intermediate <- mapvalues(x = epi$RNA_snn_res.1.1,
                                         from = annotation$Cluster,
                                         to = annotation$`epi_clean_5_harmony_28_25_1,1_intermediate`)
myeloids$annotation_intermediate <- mapvalues(x = myeloids$RNA_snn_res.1.3,
                                              from = annotation$Cluster,
                                              to = annotation$`Together_myeloids_30_25_1,3_intermediate`)
stroma$annotation_intermediate <- mapvalues(x = stroma$RNA_snn_res.1.3,
                                            from = annotation$Cluster,
                                            to = annotation$`Stroma_harmony_clean_41_27_1,3_intermediate`)
cycling$annotation_intermediate <- mapvalues(x = cycling$RNA_snn_res.0.3,
                                             from = annotation$Cluster,
                                             to = annotation$Together_cycling_filt25_20_23_intermediate)

#
# Data saving ------------------------------------------------------------------
#

dir.create('Analysis/data/00_annotation_process/00_anotadas')
saveRDS(plasmas, file = 'Analysis/data/00_annotation_process/00_anotadas/plasmas.RDS')
saveRDS(tcells, file = 'Analysis/data/00_annotation_process/00_anotadas/tcells.RDS')
saveRDS(epi, file = 'Analysis/data/00_annotation_process/00_anotadas/epi.RDS')
saveRDS(myeloids, file = 'Analysis/data/00_annotation_process/00_anotadas/myeloids.RDS')
saveRDS(stroma, file = 'Analysis/data/00_annotation_process/00_anotadas/stroma.RDS')
saveRDS(cycling, file = 'Analysis/data/00_annotation_process/00_anotadas/cycling.RDS')
