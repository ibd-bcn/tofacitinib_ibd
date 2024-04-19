options(stringsAsFactors = FALSE,bitmapType = "cairo")
library(Seurat)
library(ggplot2)
library(patchwork)
library(grid)


## Supplementary Figure 2A: Dotplot -------------------------------------------
# All annotated together

todas <- readRDS('Analysis/data//00_annotation_process/00_anotadas/todas.RDS')
todas$subset  <- gsub('myeloids', 'Myeloids', todas$subset)
todas$subset  <- gsub('stroma', 'Stroma', todas$subset)
todas$subset  <- gsub('epi', 'Epithelium', todas$subset)
todas$subset  <- gsub('tcells', 'T cells', todas$subset)
todas$subset  <- gsub('plasmas', 'Plasma and B cells', todas$subset)
todas$subset  <- gsub('cycling', 'Cycling', todas$subset)
todas$subset <- factor(todas$subset, levels = c("Myeloids", "Stroma", "Epithelium", "T cells", "Plasma and B cells", "Cycling"))


dotplot_together <- DotPlot(todas, features = c("TYROBP", "FCER1G", "S100A8", "PLAUR", "LYZ","NNMT", "COL3A1", "COL1A1", "DCN", "LUM","KRT8", "PIGR",
                            "EPCAM", "MUC12", "TFF3","CD2", "CD3E", "TRBC2", "CD7", "KLRB1",  "IGHA1", "CD79A", "MZB1", "DERL3",
                            "IGHG1", "PCLAF", "PTTG1", "TOP2A", "CDKN3", "MKI67"),
        group.by = 'subset', cols = 'RdYlBu',
        cluster.idents = F, dot.scale = 2.5) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 9)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


save_sizes(plot = dotplot_together, filename = 'sup2_dotplot_together', device = 'jpeg')
save_sizes(plot = dotplot_together, filename = 'sup2_dotplot_together', device = 'tiff')
save_sizes(plot = dotplot_together, filename = 'sup2_dotplot_together', device = 'svg')
save_sizes(plot = dotplot_together, filename = 'sup2_dotplot_together', device = 'pdf')



## Supplementary Figure 2B: Dotplot per subset----------------------------------
# Myeloid

myeloids <- readRDS('Analysis/data/00_annotation_process/00_anotadas/myeloids.RDS')
myeloids$annotation_intermediate  <- gsub('Macrophage NRG1', 'IDA macrophages', myeloids$annotation_intermediate)
myeloids$annotation_intermediate <- factor(myeloids$annotation_intermediate,
                                           levels = c("Cycling myeloids", "DCs", "Eosinophils", "HS",
                                                      "Inflammatory monocytes", "M0", "M1", "M2",
                                                      "IDA macrophages", "Mast", "Neutrophils", "pDC",
                                                      "Rib hi myeloids"))




dotplot_myeloids <- DotPlot(myeloids, features = c("TOP2A", "PCLAF", "CLEC9A", "CD1C", "CLC", "IL4", "DNAJB1", "HSPA1B", "VCAN",
                               "FCN1", "C1QB", "SELENOP", "SPP1", "INHBA", "CD209", "ATF3", "NRG1", "TPSB2", "TPSAB1",
                               "CMTM2", "PROK2", "GZMB", "LILRA4", "RPS19", "RPS18"),
        group.by = 'annotation_intermediate', cols = 'RdYlBu',
        cluster.idents = F, dot.scale = 2.5) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 9)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))



umap_myeloids <- DimPlot(myeloids, group.by = "annotation_intermediate", pt.size = 0.1) +
  theme_void() +
  theme(legend.position = "left",
        text = element_text(family = "Helvetica", size = 12),
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.size = unit(0.25, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  ggtitle(NULL)




save_sizes(plot = dotplot_myeloids, filename = 'sup2_dotplot_myeloids', device = 'jpeg')
save_sizes(plot = dotplot_myeloids, filename = 'sup2_dotplot_myeloids', device = 'tiff')
save_sizes(plot = dotplot_myeloids, filename = 'sup2_dotplot_myeloids', device = 'svg')
save_sizes(plot = dotplot_myeloids, filename = 'sup2_dotplot_myeloids', device = 'pdf')

save_sizes(plot = umap_myeloids, filename = 'sup_2umap_myeloids', device = 'pdf')
save_sizes(plot = umap_myeloids, filename = 'sup_2umap_myeloids', device = 'jpeg')
save_sizes(plot = umap_myeloids, filename = 'sup_2umap_myeloids', device = 'tiff')
save_sizes(plot = umap_myeloids, filename = 'sup_2umap_myeloids', device = 'svg')


# Stroma

stroma <- readRDS('Analysis/data/00_annotation_process/00_anotadas/stroma.RDS')
stroma$annotation_intermediate  <- gsub('_', ' ', stroma$annotation_intermediate)
stroma$annotation_intermediate  <- gsub('MT fibroblasts', 'MT hi fibroblasts', stroma$annotation_intermediate)
stroma$annotation_intermediate <- factor(stroma$annotation_intermediate,
                                           levels = rev(c("S3", "S2", "S1", "Perycites", "Myofibroblasts",
                                                      "MT hi fibroblasts", "Inflammatory fibroblasts", "IER fibroblasts",
                                                      "Glia", "Endothelium", "Cycling fibroblasts")))




dotplot_stroma <- DotPlot(stroma, features = c("TOP2A", "PCLAF", "PLVAP", "VWF", "NRXN1", "S100B", "EGR1", "FOSB", "IL13RA2",
                                          "IL11", "MT-ND3", "MT-CO2", "ACTG2", "MYH11", "NOTCH3", "COX4I2", "ADAMDEC1", "DCN", "SOX6",
                                          "F3", "OGN", "GREM1"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2.5) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 9)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap_stroma <- DimPlot(stroma, group.by = "annotation_intermediate", pt.size = 0.1) +
  theme_void() +
  theme(legend.position = "left",
        text = element_text(family = "Helvetica", size = 12),
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.size = unit(0.25, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  ggtitle(NULL)





save_sizes(plot = dotplot_stroma, filename = 'sup2_dotplot_stroma', device = 'jpeg')
save_sizes(plot = dotplot_stroma, filename = 'sup2_dotplot_stroma', device = 'tiff')
save_sizes(plot = dotplot_stroma, filename = 'sup2_dotplot_stroma', device = 'svg')
save_sizes(plot = dotplot_stroma, filename = 'sup2_dotplot_stroma', device = 'pdf')

save_sizes(plot = umap_stroma, filename = 'sup_2umap_stroma', device = 'pdf')
save_sizes(plot = umap_stroma, filename = 'sup_2umap_stroma', device = 'jpeg')
save_sizes(plot = umap_stroma, filename = 'sup_2umap_stroma', device = 'tiff')
save_sizes(plot = umap_stroma, filename = 'sup_2umap_stroma', device = 'svg')


# Tcells

tcells <- readRDS('Analysis/data/00_annotation_process/00_anotadas/tcells.RDS')
tcells$annotation_intermediate  <- gsub('MT hi IER', 'MT hi T cells', tcells$annotation_intermediate)
tcells$annotation_intermediate  <- gsub('CD4 CD8 IFIT3', 'IFN-activated', tcells$annotation_intermediate)
tcells$annotation_intermediate  <- gsub('CD8', 'CD8 CTL', tcells$annotation_intermediate)

tcells$annotation_intermediate <- factor(tcells$annotation_intermediate,
                                         levels = rev(c("Tregs", "Thf", "Th17", "Ribhi T cells", "NK", "Naïve T cells","MT hi T cells", "ILC3",
                                                        "HS", "DN EOMES", "CD8 CTL", "IFN-activated", "CD4")))




dotplot_tcells <- DotPlot(tcells, features = c("ANXA1", "IL7R", "IFIT3", "IFIT1", "NKG7", "GZMK", "MYC", "PTK7", "HSPA1B",
                                        "PCDH9", "KIT", "MALAT1", "MT-ATP6", "SELL", "CCR7", "GNLY", "FCER1G",
                                        "RPS4X", "RPS18", "CCL20", "IL17A", "CXCL13", "MAGEH1", "TNFRSF4",
                                        "FOXP3"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2.5) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 9)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap_tcells <- DimPlot(tcells, group.by = "annotation_intermediate", pt.size = 0.1) +
  theme_void() +
  theme(legend.position = "left",
        text = element_text(family = "Helvetica", size = 12),
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.size = unit(0.25, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  ggtitle(NULL)





save_sizes(plot = dotplot_tcells, filename = 'sup2_dotplot_tcells', device = 'jpeg')
save_sizes(plot = dotplot_tcells, filename = 'sup2_dotplot_tcells', device = 'tiff')
save_sizes(plot = dotplot_tcells, filename = 'sup2_dotplot_tcells', device = 'svg')
save_sizes(plot = dotplot_tcells, filename = 'sup2_dotplot_tcells', device = 'pdf')

save_sizes(plot = umap_tcells, filename = 'sup_2umap_tcells', device = 'pdf')
save_sizes(plot = umap_tcells, filename = 'sup_2umap_tcells', device = 'jpeg')
save_sizes(plot = umap_tcells, filename = 'sup_2umap_tcells', device = 'tiff')
save_sizes(plot = umap_tcells, filename = 'sup_2umap_tcells', device = 'svg')

# Epi


epi <- readRDS('Analysis/data/00_annotation_process/00_anotadas/epi.RDS')
epi$annotation_intermediate  <- gsub('Epithlium Rib hi', 'Epithelium Rib hi', epi$annotation_intermediate)
epi$annotation_intermediate <- factor(epi$annotation_intermediate,
                                         levels = rev(c("Undifferentiated epithelium", "Tuft", "Stem", "Secretory progenitor",
                                                        "Mature goblet", "IER Epithelium", "Goblet", "Epithelium Rib hi",
                                                        "Enteroendocrines", "Cycling TA", "Colonocytes", "APOA4")))

dotplot_epi <- DotPlot(epi, features = c("APOA4", "APOA1", "GUCA2A", "CLCA4", "MKI67", "CENPF", "PYY", "CHGA", "RPLP1",
                                     "RPL41", "MUC2", "REP15", "PLCG2", "CSKMT", "FER1L6", "ZG16", "SPINK4", "TFF3",
                                     "LGR5", "SMOC2", "LRMP", "TRPM5", "MRPL11", "TOMM22"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2.5) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 9)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap_epi <- DimPlot(epi, group.by = "annotation_intermediate", pt.size = 0.1) +
  theme_void() +
  theme(legend.position = "left",
        text = element_text(family = "Helvetica", size = 12),
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.size = unit(0.25, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  ggtitle(NULL)





save_sizes(plot = dotplot_epi, filename = 'sup2_dotplot_epi', device = 'jpeg')
save_sizes(plot = dotplot_epi, filename = 'sup2_dotplot_epi', device = 'tiff')
save_sizes(plot = dotplot_epi, filename = 'sup2_dotplot_epi', device = 'svg')
save_sizes(plot = dotplot_epi, filename = 'sup2_dotplot_epi', device = 'pdf')

save_sizes(plot = umap_epi, filename = 'sup_2umap_epi', device = 'pdf')
save_sizes(plot = umap_epi, filename = 'sup_2umap_epi', device = 'jpeg')
save_sizes(plot = umap_epi, filename = 'sup_2umap_epi', device = 'tiff')
save_sizes(plot = umap_epi, filename = 'sup_2umap_epi', device = 'svg')

# Plasmas


plasmas <- readRDS('Analysis/data/00_annotation_process/00_anotadas/plasmas.RDS')
plasmas$annotation_intermediate  <- gsub('_', ' ', plasmas$annotation_intermediate)
plasmas$annotation_intermediate <- factor(plasmas$annotation_intermediate,
                                      levels = rev(c("Plasmablast IGKC", "Plasmablast IgG","Plasmablast IgA","PC PSAT1","PC IGLV6-57",
                                                     "PC IGLL5", "PC IgG", "PC IgA", "PC IFIT1", "PC IER", "PC heat shock", "Cycling B cell",
                                                      "B cell")))

dotplot_plasmas <- DotPlot(plasmas, features = c("MS4A1", "CXCR4", "TMSB4X", "ACTB", "HSPA1A", "HIST1H2BG", "HIST1H2AC", "IFIT1", "RSAD2", "IGHA2",
                                         "IGHA1", "IGHGP", "IGHG4", "IGLL5", "NEAT1", "IGLC3", "IGLC2", "PSAT1", "CHAC1", "JCHAIN", "IGHG1",
                                         "IGHG3", "IGKC", "MZB1"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2.5) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 9)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap_plasmas <- DimPlot(plasmas, group.by = "annotation_intermediate", pt.size = 0.1) +
  theme_void() +
  theme(legend.position = "left",
        text = element_text(family = "Helvetica", size = 12),
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.size = unit(0.25, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  ggtitle(NULL)





save_sizes(plot = dotplot_plasmas, filename = 'sup2_dotplot_plasmas', device = 'jpeg')
save_sizes(plot = dotplot_plasmas, filename = 'sup2_dotplot_plasmas', device = 'tiff')
save_sizes(plot = dotplot_plasmas, filename = 'sup2_dotplot_plasmas', device = 'svg')
save_sizes(plot = dotplot_plasmas, filename = 'sup2_dotplot_plasmas', device = 'pdf')

save_sizes(plot = umap_plasmas, filename = 'sup_2umap_plasmas', device = 'pdf')
save_sizes(plot = umap_plasmas, filename = 'sup_2umap_plasmas', device = 'jpeg')
save_sizes(plot = umap_plasmas, filename = 'sup_2umap_plasmas', device = 'tiff')
save_sizes(plot = umap_plasmas, filename = 'sup_2umap_plasmas', device = 'svg')

# Cycling
cycling <- readRDS('Analysis/data/00_annotation_process/00_anotadas/cycling.RDS')

cycling$annotation_intermediate  <- gsub('_', ' ', cycling$annotation_intermediate)
cycling$annotation_intermediate <- factor(cycling$annotation_intermediate,
                                          levels = rev(c("Cycling T cell", "Cycling PC",  "Cycling Myeloid","Cycling B cell")))

dotplot_cycling <- DotPlot(cycling, features = c("RPL28", "MS4A1", "AIF1", "LYZ", "MZB1", "DERL3", "IL32", "CD3D"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2.5) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 9)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))

umap_cycling <- DimPlot(cycling, group.by = "annotation_intermediate", pt.size = 0.1) +
  theme_void() +
  theme(legend.position = "left",
        text = element_text(family = "Helvetica", size = 12),
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.size = unit(0.25, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  ggtitle(NULL)





save_sizes(plot = dotplot_cycling, filename = 'sup2_dotplot_cycling', device = 'jpeg')
save_sizes(plot = dotplot_cycling, filename = 'sup2_dotplot_cycling', device = 'tiff')
save_sizes(plot = dotplot_cycling, filename = 'sup2_dotplot_cycling', device = 'svg')
save_sizes(plot = dotplot_cycling, filename = 'sup2_dotplot_cycling', device = 'pdf')

save_sizes(plot = umap_cycling, filename = 'sup_2umap_cycling', device = 'pdf')
save_sizes(plot = umap_cycling, filename = 'sup_2umap_cycling', device = 'jpeg')
save_sizes(plot = umap_cycling, filename = 'sup_2umap_cycling', device = 'tiff')
save_sizes(plot = umap_cycling, filename = 'sup_2umap_cycling', device = 'svg')

