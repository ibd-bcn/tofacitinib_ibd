

# Dotplot figure by subset

library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
# Myeloid

myeloids <- readRDS('Analysis/data/00_annotation_process/00_anotadas/myeloids.RDS')
myeloids$annotation_intermediate  <- gsub('Macrophage NRG1', 'IDA macrophages', myeloids$annotation_intermediate)
myeloids$annotation_intermediate <- factor(myeloids$annotation_intermediate,
                                           levels = c("Cycling myeloids", "DCs", "Eosinophils", "HS",
                                                      "Inflammatory monocytes", "M0", "M1", "M2",
                                                      "IDA macrophages", "Mast", "Neutrophils", "pDC",
                                                      "Rib hi myeloids"))




dotplot <- DotPlot(myeloids, features = c("TOP2A", "PCLAF", "CLEC9A", "CD1C", "CLC", "IL4", "DNAJB1", "HSPA1B", "VCAN",
                               "FCN1", "C1QB", "SELENOP", "SPP1", "INHBA", "CD209", "ATF3", "NRG1", "TPSB2", "TPSAB1",
                               "CMTM2", "PROK2", "GZMB", "LILRA4", "RPS19", "RPS18"),
        group.by = 'annotation_intermediate', cols = 'RdYlBu',
        cluster.idents = F, dot.scale = 2) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 6)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap <- DimPlot(myeloids, group.by = "annotation_intermediate", pt.size = 0.1) + theme_void() + theme(legend.position = "none")+ theme(text = element_text(family = "Helvetica")) + ggtitle(NULL)

sup2_myeloids <- umap + dotplot

save_sizes(plot = sup2_myeloids, filename = 'sup2_myeloids', device = 'jpeg')
save_sizes(plot = sup2_myeloids, filename = 'sup2_myeloids', device = 'tiff')
save_sizes(plot = sup2_myeloids, filename = 'sup2_myeloids', device = 'svg')
save_sizes(plot = sup2_myeloids, filename = 'sup2_myeloids', device = 'pdf')

legend <- umap & theme(legend.position = 'left', legend.text = element_text(size = 6), text = element_text(family = "Helvetica"))
leg <- get_legend(legend)

plot_legend_sup2myeloids <- patchwork::wrap_plots(leg)
save_sizes(plot = plot_legend_sup2myeloids, filename = 'legend_sup2myeloids', device = 'pdf')
save_sizes(plot = plot_legend_sup2myeloids, filename = 'legend_sup2myeloids', device = 'jpeg')
save_sizes(plot = plot_legend_sup2myeloids, filename = 'legend_sup2myeloids', device = 'tiff')
save_sizes(plot = plot_legend_sup2myeloids, filename = 'legend_sup2myeloids', device = 'svg')


# Stroma

stroma <- readRDS('Analysis/data/00_annotation_process/00_anotadas/stroma.RDS')
stroma$annotation_intermediate  <- gsub('_', ' ', stroma$annotation_intermediate)
stroma$annotation_intermediate  <- gsub('MT fibroblasts', 'MT hi fibroblasts', stroma$annotation_intermediate)
stroma$annotation_intermediate <- factor(stroma$annotation_intermediate,
                                           levels = rev(c("S3", "S2", "S1", "Perycites", "Myofibroblasts",
                                                      "MT hi fibroblasts", "Inflammatory fibroblasts", "IER fibroblasts",
                                                      "Glia", "Endothelium", "Cycling fibroblasts")))




dotplot <- DotPlot(stroma, features = c("TOP2A", "PCLAF", "PLVAP", "VWF", "NRXN1", "S100B", "EGR1", "FOSB", "IL13RA2",
                                          "IL11", "MT-ND3", "MT-CO2", "ACTG2", "MYH11", "NOTCH3", "COX4I2", "ADAMDEC1", "DCN", "SOX6",
                                          "F3", "OGN", "GREM1"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 6)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap <- DimPlot(stroma, group.by = "annotation_intermediate", pt.size = 0.1) + theme_void() + theme(legend.position = "none")+ theme(text = element_text(family = "Helvetica")) + ggtitle(NULL)

sup2_stroma <- umap + dotplot

save_sizes(plot = sup2_stroma, filename = 'sup2_stroma', device = 'jpeg')
save_sizes(plot = sup2_stroma, filename = 'sup2_stroma', device = 'tiff')
save_sizes(plot = sup2_stroma, filename = 'sup2_stroma', device = 'svg')
save_sizes(plot = sup2_stroma, filename = 'sup2_stroma', device = 'pdf')

legend <- umap & theme(legend.position = 'left', legend.text = element_text(size = 6), text = element_text(family = "Helvetica"))
leg <- get_legend(legend)

plot_legend_sup2stroma <- patchwork::wrap_plots(leg)
save_sizes(plot = plot_legend_sup2stroma, filename = 'legend_sup2stroma', device = 'pdf')
save_sizes(plot = plot_legend_sup2stroma, filename = 'legend_sup2stroma', device = 'jpeg')
save_sizes(plot = plot_legend_sup2stroma, filename = 'legend_sup2stroma', device = 'tiff')
save_sizes(plot = plot_legend_sup2stroma, filename = 'legend_sup2stroma', device = 'svg')


# Tcells

tcells <- readRDS('Analysis/data/00_annotation_process/00_anotadas/tcells.RDS')
tcells$annotation_intermediate  <- gsub('MT hi IER', 'MT hi T cells', tcells$annotation_intermediate)
tcells$annotation_intermediate  <- gsub('CD4 CD8 IFIT3', 'IFN-activated', tcells$annotation_intermediate)
tcells$annotation_intermediate  <- gsub('CD8', 'CD8 CTL', tcells$annotation_intermediate)

tcells$annotation_intermediate <- factor(tcells$annotation_intermediate,
                                         levels = rev(c("Tregs", "Thf", "Th17", "Ribhi T cells", "NK", "Naïve T cells","MT hi T cells", "ILC3",
                                                        "HS", "DN EOMES", "CD8 CTL", "IFN-activated", "CD4")))




dotplot <- DotPlot(tcells, features = c("ANXA1", "IL7R", "IFIT3", "IFIT1", "NKG7", "GZMK", "MYC", "PTK7", "HSPA1B",
                                        "PCDH9", "KIT", "MALAT1", "MT-ATP6", "SELL", "CCR7", "GNLY", "FCER1G",
                                        "RPS4X", "RPS18", "CCL20", "IL17A", "CXCL13", "MAGEH1", "TNFRSF4",
                                        "FOXP3"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 6)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap <- DimPlot(tcells, group.by = "annotation_intermediate", pt.size = 0.1) + theme_void() + theme(legend.position = "none")+ theme(text = element_text(family = "Helvetica")) + ggtitle(NULL)

sup2_tcells <- umap + dotplot

save_sizes(plot = sup2_tcells, filename = 'sup2_tcells', device = 'jpeg')
save_sizes(plot = sup2_tcells, filename = 'sup2_tcells', device = 'tiff')
save_sizes(plot = sup2_tcells, filename = 'sup2_tcells', device = 'svg')
save_sizes(plot = sup2_tcells, filename = 'sup2_tcells', device = 'pdf')

legend <- umap & theme(legend.position = 'left', legend.text = element_text(size = 6), text = element_text(family = "Helvetica"))
leg <- get_legend(legend)

plot_legend_sup2tcells <- patchwork::wrap_plots(leg)
save_sizes(plot = plot_legend_sup2tcells, filename = 'legend_sup2tcells', device = 'pdf')
save_sizes(plot = plot_legend_sup2tcells, filename = 'legend_sup2tcells', device = 'jpeg')
save_sizes(plot = plot_legend_sup2tcells, filename = 'legend_sup2tcells', device = 'tiff')
save_sizes(plot = plot_legend_sup2tcells, filename = 'legend_sup2tcells', device = 'svg')

# Epi


Epi <- readRDS('Analysis/data/00_annotation_process/00_anotadas/Epi.RDS')
Epi$annotation_intermediate  <- gsub('Epithlium Rib hi', 'Epithelium Rib hi', Epi$annotation_intermediate)
Epi$annotation_intermediate <- factor(Epi$annotation_intermediate,
                                         levels = rev(c("Undifferentiated epithelium", "Tuft", "Stem", "Secretory progenitor",
                                                        "Mature goblet", "IER Epithelium", "Goblet", "Epithelium Rib hi",
                                                        "Enteroendocrines", "Cycling TA", "Colonocytes", "APOA4")))

dotplot <- DotPlot(Epi, features = c("APOA4", "APOA1", "GUCA2A", "CLCA4", "MKI67", "CENPF", "PYY", "CHGA", "RPLP1",
                                     "RPL41", "MUC2", "REP15", "PLCG2", "CSKMT", "FER1L6", "ZG16", "SPINK4", "TFF3",
                                     "LGR5", "SMOC2", "LRMP", "TRPM5", "MRPL11", "TOMM22"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 6)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap <- DimPlot(Epi, group.by = "annotation_intermediate", pt.size = 0.1) + theme_void() + theme(legend.position = "none")+ theme(text = element_text(family = "Helvetica")) + ggtitle(NULL)

sup2_Epi <- umap + dotplot

save_sizes(plot = sup2_Epi, filename = 'sup2_Epi', device = 'jpeg')
save_sizes(plot = sup2_Epi, filename = 'sup2_Epi', device = 'tiff')
save_sizes(plot = sup2_Epi, filename = 'sup2_Epi', device = 'svg')
save_sizes(plot = sup2_Epi, filename = 'sup2_Epi', device = 'pdf')

legend <- umap & theme(legend.position = 'left', legend.text = element_text(size = 6), text = element_text(family = "Helvetica"))
leg <- get_legend(legend)

plot_legend_sup2Epi <- patchwork::wrap_plots(leg)
save_sizes(plot = plot_legend_sup2Epi, filename = 'legend_sup2Epi', device = 'pdf')
save_sizes(plot = plot_legend_sup2Epi, filename = 'legend_sup2Epi', device = 'jpeg')
save_sizes(plot = plot_legend_sup2Epi, filename = 'legend_sup2Epi', device = 'tiff')
save_sizes(plot = plot_legend_sup2Epi, filename = 'legend_sup2Epi', device = 'svg')



# Plasmas


plasmas <- readRDS('Analysis/data/00_annotation_process/00_anotadas/plasmas.RDS')
plasmas$annotation_intermediate  <- gsub('_', ' ', plasmas$annotation_intermediate)
plasmas$annotation_intermediate <- factor(plasmas$annotation_intermediate,
                                      levels = rev(c("Plasmablast IGKC", "Plasmablast IgG","Plasmablast IgA","PC PSAT1","PC IGLV6-57",
                                                     "PC IGLL5", "PC IgG", "PC IgA", "PC IFIT1", "PC IER", "PC heat shock", "Cycling B cell",
                                                      "B cell")))

dotplot <- DotPlot(plasmas, features = c("MS4A1", "CXCR4", "TMSB4X", "ACTB", "HSPA1A", "HIST1H2BG", "HIST1H2AC", "IFIT1", "RSAD2", "IGHA2",
                                         "IGHA1", "IGHGP", "IGHG4", "IGLL5", "NEAT1", "IGLC3", "IGLC2", "PSAT1", "CHAC1", "JCHAIN", "IGHG1",
                                         "IGHG3", "IGKC", "MZB1"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 6)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap <- DimPlot(plasmas, group.by = "annotation_intermediate", pt.size = 0.1) + theme_void() + theme(legend.position = "none")+ theme(text = element_text(family = "Helvetica")) + ggtitle(NULL)

sup2_plasmas <- umap + dotplot

save_sizes(plot = sup2_plasmas, filename = 'sup2_plasmas', device = 'jpeg')
save_sizes(plot = sup2_plasmas, filename = 'sup2_plasmas', device = 'tiff')
save_sizes(plot = sup2_plasmas, filename = 'sup2_plasmas', device = 'svg')
save_sizes(plot = sup2_plasmas, filename = 'sup2_plasmas', device = 'pdf')

legend <- umap & theme(legend.position = 'left', legend.text = element_text(size = 6), text = element_text(family = "Helvetica"))
leg <- get_legend(legend)

plot_legend_sup2plasmas <- patchwork::wrap_plots(leg)
save_sizes(plot = plot_legend_sup2plasmas, filename = 'legend_sup2plasmas', device = 'pdf')
save_sizes(plot = plot_legend_sup2plasmas, filename = 'legend_sup2plasmas', device = 'jpeg')
save_sizes(plot = plot_legend_sup2plasmas, filename = 'legend_sup2plasmas', device = 'tiff')
save_sizes(plot = plot_legend_sup2plasmas, filename = 'legend_sup2plasmas', device = 'svg')

# Cycling
cycling <- readRDS('Analysis/data/00_annotation_process/00_anotadas/cycling.RDS')

cycling$annotation_intermediate  <- gsub('_', ' ', cycling$annotation_intermediate)
cycling$annotation_intermediate <- factor(cycling$annotation_intermediate,
                                          levels = rev(c("Cycling T cell", "Cycling PC",  "Cycling Myeloid","Cycling B cell")))

dotplot <- DotPlot(cycling, features = c("RPL28", "MS4A1", "AIF1", "LYZ", "MZB1", "DERL3", "IL32", "CD3D"),
                   group.by = 'annotation_intermediate', cols = 'RdYlBu',
                   cluster.idents = F, dot.scale = 2) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 6)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))


umap <- DimPlot(cycling, group.by = "annotation_intermediate", pt.size = 0.1) + theme_void() + theme(legend.position = "none")+ theme(text = element_text(family = "Helvetica")) + ggtitle(NULL)

sup2_cycling <- umap + dotplot

save_sizes(plot = sup2_cycling, filename = 'sup2_cycling', device = 'jpeg')
save_sizes(plot = sup2_cycling, filename = 'sup2_cycling', device = 'tiff')
save_sizes(plot = sup2_cycling, filename = 'sup2_cycling', device = 'svg')
save_sizes(plot = sup2_cycling, filename = 'sup2_cycling', device = 'pdf')

legend <- umap & theme(legend.position = 'left', legend.text = element_text(size = 6), text = element_text(family = "Helvetica"))
leg <- get_legend(legend)

plot_legend_sup2cycling <- patchwork::wrap_plots(leg)
save_sizes(plot = plot_legend_sup2cycling, filename = 'legend_sup2cycling', device = 'pdf')
save_sizes(plot = plot_legend_sup2cycling, filename = 'legend_sup2cycling', device = 'jpeg')
save_sizes(plot = plot_legend_sup2cycling, filename = 'legend_sup2cycling', device = 'tiff')
save_sizes(plot = plot_legend_sup2cycling, filename = 'legend_sup2cycling', device = 'svg')
