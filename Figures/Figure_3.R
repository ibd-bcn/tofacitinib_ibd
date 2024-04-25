options(stringsAsFactors = FALSE,)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readr)
library(dplyr)
library(ggrepel)
library(cowplot)
library(ggpubr)
source('Figures/functions_plots.R')

## Data ------------------------------------------------------------------------
message('Loading data')

myeloids <- readRDS('Analysis/data/00_annotation_process/00_anotadas/myeloids.RDS')
myeloids$pre_post <- plyr::mapvalues(myeloids$week_3, from = c('W0', 'POST'), to = c('PRE', 'POST'))
myeloids$annotation_intermediate <- gsub('Macrophage NRG1', 'IDA macrophages', myeloids$annotation_intermediate)
myeloids$pre_post <- factor(myeloids$pre_post, levels = c('PRE', 'POST'))
myeloids$response <- factor(myeloids$response, levels = c('R', 'NR'))
myeloids$subset <- factor(myeloids$subset, levels = c('epi', 'stroma', 'plasmas', 'myeloids', 'cycling', 'tcells'))


stroma <- readRDS('Analysis/data/00_annotation_process/00_anotadas/stroma.RDS')
stroma$annotation_intermediate <- gsub('Inflammatory_fibroblasts', 'Inflammatory fibroblasts', stroma$annotation_intermediate)
stroma$annotation_intermediate <- gsub('IER_fibroblasts', 'IER fibroblasts', stroma$annotation_intermediate)
stroma$pre_post <- plyr::mapvalues(stroma$week_3, from = c('W0', 'POST'), to = c('PRE', 'POST'))
stroma$pre_post <- factor(stroma$pre_post, levels = c('PRE', 'POST'))
stroma$response <- factor(stroma$response, levels = c('R', 'NR'))
stroma$subset <- factor(stroma$subset, levels = c('epi', 'stroma', 'plasmas', 'myeloids', 'cycling', 'tcells'))
## Figure 3A  ------------------------------------------------------------------
# UMAP split by response and treatment, with some celltypes in colour
# HS, Inflammatory monocytes, M0, M1, M2, IDA macrophages, Neutrophils

colors <- c('HS' = '#bed0e8',
            'Inflammatory monocytes' = '#eb6c92',
            'M0' = '#c1d665',
            'M1' = '#053370',
            'M2' = '#018d7e',
            'IDA macrophages' = '#f3bf5d',
            'Neutrophils' = '#ba4042')

resp_pre <- DimPlot(myeloids[,myeloids@meta.data$response == 'R' & myeloids@meta.data$pre_post == 'PRE'],
                    group.by = 'annotation_intermediate', pt.size = 0.1) +
  scale_color_manual(values = colors, na.value = '#e6e5e5') +
  theme_umap()

nresp_pre <- DimPlot(myeloids[,myeloids@meta.data$response == 'NR' & myeloids@meta.data$pre_post == 'PRE'],
                    group.by = 'annotation_intermediate', pt.size = 0.1) +
  scale_color_manual(values = colors, na.value = '#e6e5e5') +
  theme_umap()

resp_post <- DimPlot(myeloids[,myeloids@meta.data$response == 'R' & myeloids@meta.data$pre_post == 'POST'],
                    group.by = 'annotation_intermediate', pt.size = 0.1) +
  scale_color_manual(values = colors, na.value = '#e6e5e5') +
  theme_umap()

nresp_post <- DimPlot(myeloids[,myeloids@meta.data$response == 'NR' & myeloids@meta.data$pre_post == 'POST'],
                    group.by = 'annotation_intermediate', pt.size = 0.1) +
  scale_color_manual(values = colors, na.value = '#e6e5e5') +
  theme_umap()+
  theme(legend.position = 'left')

fig3a <- patchwork::wrap_plots(resp_pre, nresp_pre, resp_post, nresp_post, guides = 'collect') &
  theme(text = element_text(family = 'Helvetica', size = 8))

fig3a_nl <- fig3a & theme(legend.position = 'none')

save_sizes(plot = fig3a, filename = 'Figure_3A', device = 'jpeg')
save_sizes(plot = fig3a, filename = 'Figure_3A', device = 'tiff')
save_sizes(plot = fig3a, filename = 'Figure_3A', device = 'svg')
save_sizes(plot = fig3a, filename = 'Figure_3A', device = 'pdf')

save_sizes(plot = fig3a_nl, filename = 'Figure_3A_no_legend', device = 'jpeg')
save_sizes(plot = fig3a_nl, filename = 'Figure_3A_no_legend', device = 'tiff')
save_sizes(plot = fig3a_nl, filename = 'Figure_3A_no_legend', device = 'svg')
save_sizes(plot = fig3a_nl, filename = 'Figure_3A_no_legend', device = 'pdf')

## Figure 3B -------------------------------------------------------------------
de_data <- readRDS('Analysis/data/01_DE/REPASO/new_complete.RDS')

# Volcano plots M2

# M2 responders
cluster <- "M2"
comp <- "w0R_vs_POSTR"

genes_up <-  c("IGF1", "CLEC10A", "CD163L1", "IL10RA", "AHR", "MAF")
genes_down <- c("S100A9", "GBP1", "MMP12", "HLA-A","IFITM3", "STAT1", "HLA-B", "FCGR3A", "FCGR2A", "CCL13")


volcano_plot <- function(cluster, comp, genes_up, genes_down) {

  labels <- c(
    "UPP" = "UPP",
    "UP" = "UP",
    "0" = "0",
    "DW" = "DW",
    "DWW" = "DWW"
  )

  de_data2 <- de_data[de_data$cluster == cluster &
                        de_data$annotation == 'annotation_intermediate' &
                        de_data$comp == comp &
                        de_data$sign %in% names(labels), c("p_val", "avg_log2FC", "sign", "comp", "gene")]

  data <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))
  data$sign <- factor(data$sign, levels = c("UPP", "UP", "0", "DW", "DWW"))
  data_dw <- data[data$gene %in% genes_down, ]
  data_up <- data[data$gene %in% genes_up, ]
  fig <- ggplot(data = data, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
    geom_point(size = 1) +
    scale_color_manual(values = colors_volcano, labels = labels) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica")) +
    guides(color = guide_legend(override.aes = list(shape = 1))) +
    theme(legend.position = "none") +
    scale_y_continuous(breaks = c(seq(0, 25, 5)), limits = c(0, 25))

  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", fill = "#911704" ) +
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38") +
    geom_label_repel(data = data_dw, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold') +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold')+
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )



  return(fig)
}


fig3b_M2R <- volcano_plot(cluster, comp, genes_up, genes_down)

print(fig3b_M2R)

save_sizes(plot =fig3b_M2R , filename = 'fig3b_M2R', device = 'jpeg')
save_sizes(plot = fig3b_M2R, filename = 'fig3b_M2R', device = 'tiff')
save_sizes(plot = fig3b_M2R, filename = 'fig3b_M2R', device = 'svg')
save_sizes(plot = fig3b_M2R, filename = 'fig3b_M2R', device = 'pdf')

# M2 non-responders

cluster <- "M2"
comp <- "w0NR_vs_POSTNR"
genes_up <- c("MMP9", "INHBA", "SPP1", "IDO1", "PLAUR", "IL1RN", "IL7R", "MT2A",
                         "SOD2", "IL1B", "CCL5", "CXCL8", "IL6", "CXCL1", "S100A9")
genes_down <- c("TMSB4X", "TPT1", "FUCA1")
clec <- "CLEC5A"

volcano_plot <- function(cluster, comp, genes_up, genes_down, clec) {

  labels <- c(
    "UPP" = "UPP",
    "UP" = "UP",
    "0" = "0",
    "DW" = "DW",
    "DWW" = "DWW"
  )

  de_data2 <- de_data[de_data$cluster == cluster &
                        de_data$annotation == 'annotation_intermediate' &
                        de_data$comp == comp &
                        de_data$sign %in% names(labels), c("p_val", "avg_log2FC", "sign", "comp", "gene")]


  data <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))
  data$sign <- factor(data$sign, levels = c("UPP", "UP", "0", "DW", "DWW"))
  data_down <- data[data$gene %in% genes_down, ]
  data_up <- data[data$gene %in% genes_up, ]
  data_clec <- data[data$gene %in% clec,]
  fig <- ggplot(data = data, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
    geom_point(size = 1) +
    scale_fill_manual(values = colors_volcano, labels = labels) +
    scale_color_manual(values = colors_volcano, labels = labels) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica")) +
    guides(color = guide_legend(override.aes = list(shape = 1))) +
    theme(legend.position = "none") +
    scale_y_continuous(breaks = c(seq(0, 25, 5)), limits = c(0, 25)) +
    scale_x_continuous(breaks = c(seq(-2, 7, 1)), limits = c(-2, 7))


  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", aes(fill = sign )) +
    geom_point(data = data_down,shape = 21, color = "black", aes(fill = sign )) +
    geom_point(data = data_clec,shape = 21,color = "black", aes(fill = sign )) +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold', max.overlaps = Inf,
                     position = position_nudge_repel(x = 3, y = 1)) +
    geom_label_repel(data = data_down, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold', max.overlaps = Inf)+
    geom_label_repel(data = data_clec, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold', max.overlaps = Inf,
                     position = position_nudge_repel(x = 0.5, y = 0.5)) +

    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"))
  return(fig)
}
fig3b_M2NR <- volcano_plot(cluster, comp, genes_up, genes_down, clec)
print(fig3b_M2NR)

save_sizes(plot =fig3b_M2NR , filename = 'fig3b_M2NR', device = 'jpeg')
save_sizes(plot = fig3b_M2NR, filename = 'fig3b_M2NR', device = 'tiff')
save_sizes(plot = fig3b_M2NR, filename = 'fig3b_M2NR', device = 'svg')
save_sizes(plot = fig3b_M2NR, filename = 'fig3b_M2NR', device = 'pdf')


# 3B legend

fill_scale <- scale_fill_manual(values = colors_volcano, guide = guide_legend(override.aes = list(shape = 21, color = NULL, size = 3)))
color_scale <- scale_color_manual(values = colors_volcano, guide = guide_legend(override.aes = list(shape = 21, color = NULL, size = 3)))

fig <- ggplot(data = data, aes(x = avg_log2FC, y = -log10(p_val), fill = sign, color = sign), show.legend = F) +
  geom_point(size = 1) +
  fill_scale +
  color_scale +
  theme_classic() +
  theme(text = element_text(family = "Helvetica")) +

fig <- fig + guides(color = "none")
leg <- get_legend(fig)
as_ggplot(leg)
save_sizes(plot =leg , filename = 'legvolcano', device = 'jpeg')
save_sizes(plot = leg, filename = 'legvolcano', device = 'tiff')
save_sizes(plot = leg, filename = 'legvolcano', device = 'svg')
save_sizes(plot = leg, filename = 'legvolcano', device = 'pdf')


## Figure 3C -------------------------------------------------------------------
# healthy Macrophages from Garrido-Trigo et al.
myeloids_HC <- readRDS("~/data_Albas/HC_13122021/01_piezas_anotadas/myeloids.RDS")

#
# get markers from this project and Garrido-Trigo et al.
#
myeloids@active.ident <- myeloids$annotation_refined
myeloids_mk <- FindAllMarkers(myeloids,  only.pos = TRUE,
                              min.pct = 0.25, thresh.use = 0.25)

myeloids_HC@active.ident <- myeloids_HC$annotation_refined
myeloids_HC_mk <- FindAllMarkers(myeloids_HC,  only.pos = TRUE,
                                 min.pct = 0.25, thresh.use = 0.25)

#
# Get M1 and M2 markers
#
m2 <- myeloids_HC_mk[myeloids_HC_mk$cluster == 'M2' & myeloids_HC_mk$p_val_adj < 0.001,'gene']
m1 <- myeloids_mk[myeloids_mk$cluster == 'M1' & myeloids_mk$p_val_adj < 0.05,'gene']

m2_ok <- m2[!(m2 %in% m1)]
m1_ok <- m1[!(m1 %in% m2)]

gene_lists <- read_delim("~/TOFA_data/20220222_TOFAS_23/03_M1_M2_W0_POST/gene_lists_M2.csv",
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)

# W0 vs POST

rlist <- gene_lists$ALL_w0R_vs_POSTR[gene_lists$ALL_w0R_vs_POSTR_UPP_DWW == 'UPP']
nrlist <- gene_lists$ALL_w0NR_vs_POSTNR[gene_lists$ALL_w0NR_vs_POSTNR_UPP_DWW == 'UPP']

rlist_ok <- rlist[!(rlist %in% nrlist)]
nrlist_ok <- nrlist[!(nrlist %in% rlist)]

the_list <- vector(mode = 'list', length = 4)
the_list[[1]] <- rlist_ok
the_list[[2]] <- m2_ok
the_list[[3]] <- m1_ok
the_list[[4]] <- nrlist_ok[!is.na(nrlist_ok)]
names(the_list) <- c('UPP Post-tx\nRESPONDERS',
                     'M2 markers',
                     'M1 markers',
                     'UPP Post-tx\nNO RESPONDERS')

#
# jaccard
#

library(matchSCore2)
ms <- matchSCore2(gene_cl.ref = the_list[2:3],gene_cl.obs = the_list[c(1,4)],
                  ylab = "Markers", xlab = "DE genes")
ms$ggplot +
  theme_figure()

df <- as.data.frame(ms$JI.mat)
df$rownames <- rownames(df)
matrix <- tidyr::pivot_longer(df, cols = 1:2)

fig3d <- ggplot(matrix,
       aes(x = name, y = rownames, fill= value)) +
  geom_raster() +
  geom_text(mapping = aes(label = round(value, digits = 2)),  size = 10/.pt) +
  scale_fill_gradient(low  = 'white',
                      high = 'red',
                      limits = c(0,0.13),
                      name="Jaccard\nIndex",
                      breaks=c(0, 0.13),
                      labels=c("Min", 'Max')) +
  theme_figure() +
  theme(text = element_text(family = 'Helvetica', size = 10),
        legend.position = 'right',
        legend.title = element_text(family = 'Helvetica', size = 9))

save_sizes(plot = fig3d, filename = 'Figure_3D', device = 'jpeg')
save_sizes(plot = fig3d, filename = 'Figure_3D', device = 'tiff')
save_sizes(plot = fig3d, filename = 'Figure_3D', device = 'svg')
save_sizes(plot = fig3d, filename = 'Figure_3D', device = 'pdf')

fig3d_nt <- fig3d +
  theme(
    text = element_text(color = 'white'),
    axis.text = element_text(color = 'white'))


save_sizes(plot = fig3d_nt, filename = 'Figure_3D_no_text', device = 'jpeg')
save_sizes(plot = fig3d_nt, filename = 'Figure_3D_no_text', device = 'tiff')
save_sizes(plot = fig3d_nt, filename = 'Figure_3D_no_text', device = 'svg')
save_sizes(plot = fig3d_nt, filename = 'Figure_3D_no_text', device = 'pdf')


# Figure 3D ----------------------------------------------------------------------------------
# UMAP split by response and treatment, with some celltypes in colour
color_list <- c( "S3" = "#C5D5EA", "Inflammatory fibroblasts" = "#EC769A", "S2" = "#C5D86D", "S1" = "#1F487E", "Myofibroblasts" = "#1B998B", "IER fibroblasts" = "#BC4749" )


resp_pre <- DimPlot(stroma[,stroma@meta.data$response == 'R' & stroma@meta.data$pre_post == 'PRE'],
                    group.by = 'annotation_intermediate', pt.size = 0.1) +
  scale_color_manual(values = color_list, na.value = '#e6e5e5') +
  theme_umap()

nresp_pre <- DimPlot(stroma[,stroma@meta.data$response == 'NR' & stroma@meta.data$pre_post == 'PRE'],
                     group.by = 'annotation_intermediate', pt.size = 0.1) +
  scale_color_manual(values = color_list, na.value = '#e6e5e5') +
  theme_umap()

resp_post <- DimPlot(stroma[,stroma@meta.data$response == 'R' & stroma@meta.data$pre_post == 'POST'],
                     group.by = 'annotation_intermediate', pt.size = 0.1) +
  scale_color_manual(values = color_list, na.value = '#e6e5e5') +
  theme_umap()

nresp_post <- DimPlot(stroma[,stroma@meta.data$response == 'NR' & stroma@meta.data$pre_post == 'POST'],
                      group.by = 'annotation_intermediate', pt.size = 0.1) +
  scale_color_manual(values = color_list, na.value = '#e6e5e5') +
  theme_umap()+
  theme(legend.position = 'left')

fig4a <- patchwork::wrap_plots(resp_pre, nresp_pre, resp_post, nresp_post, guides = 'collect') &
  theme(text = element_text(family = 'Helvetica', size = 8))

fig4a_nl <- fig4a & theme(legend.position = 'none')

save_sizes(plot = fig4a, filename = 'Figure_4A', device = 'jpeg')
save_sizes(plot = fig4a, filename = 'Figure_4A', device = 'tiff')
save_sizes(plot = fig4a, filename = 'Figure_4A', device = 'svg')
save_sizes(plot = fig4a, filename = 'Figure_4A', device = 'pdf')

save_sizes(plot = fig4a_nl, filename = 'Figure_4A_no_legend', device = 'jpeg')
save_sizes(plot = fig4a_nl, filename = 'Figure_4A_no_legend', device = 'tiff')
save_sizes(plot = fig4a_nl, filename = 'Figure_4A_no_legend', device = 'svg')
save_sizes(plot = fig4a_nl, filename = 'Figure_4A_no_legend', device = 'pdf')

# Figure 3E  -------------------------------------------------------------------

de_data <- readRDS('Analysis/data/01_DE/REPASO/new_complete.RDS')

# Volcano plots S1

# S1 Responders

cluster <- "S1"
comp <- "w0R_vs_POSTR"

genes_up <- c("ADIRF", "ADAMDEC1", "FOS")
genes_down <- c("CHI3L1", "WNT5A", "IFITM3", "GBP1", "CCL19", "JAK1",  "SOCS1",
                "OSMR", "IL13RA2")
volcano_plot <- function(cluster, comp, genes_up, genes_down) {

  labels <- c(
    "UPP" = "UPP",
    "UP" = "UP",
    "0" = "0",
    "DW" = "DW",
    "DWW" = "DWW"
  )

  de_data2 <- de_data[de_data$cluster == cluster &
                        de_data$annotation == 'annotation_intermediate' &
                        de_data$comp == comp &
                        de_data$sign %in% names(labels), c("p_val", "avg_log2FC", "sign", "comp", "gene")]

  data <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))
  data$sign <- factor(data$sign, levels = c("UPP", "UP", "0", "DW", "DWW"))
  data_dw <- data[data$gene %in% genes_down, ]
  data_up <- data[data$gene %in% genes_up, ]


  fig <- ggplot(data = data, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
    geom_point(size = 1) +
    scale_color_manual(values = colors_volcano, labels = labels) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica")) +
    guides(color = guide_legend(override.aes = list(shape = 1))) +
    theme(legend.position = "none") +
    scale_y_continuous(breaks = c(seq(0, 55, 5)), limits = c(0, 55))+
    scale_x_continuous(breaks = c(seq(-8, 8, 1)), limits = c(-8, 8))


  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", fill = "#911704" , size = 1) +
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38", size = 1) +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold',
                     position = position_nudge_repel(x = 3, y = 4)) +
    geom_label_repel(data = data_dw, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold',
                     position = position_nudge_repel(x = -2, y = 5))+
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}


fig4c_S1R <- volcano_plot(cluster, comp, genes_up, genes_down)

print(fig4c_S1R)


save_sizes(plot =fig4c_S1R , filename = 'fig4c_S1R', device = 'jpeg')
save_sizes(plot = fig4c_S1R, filename = 'fig4c_S1R', device = 'tiff')
save_sizes(plot = fig4c_S1R, filename = 'fig4c_S1R', device = 'svg')
save_sizes(plot = fig4c_S1R, filename = 'fig4c_S1R', device = 'pdf')

# S1 non-responders
cluster <- "S1"
comp <- "w0NR_vs_POSTNR"
genes_up <- c("POSTN", "FTH1", "RPS4Y1", "TIMP1", "MT2A","CHI3L1", "DDT", "RARRES2", "LGALS3")
genes_down <- c("IFI27", "MT-CO3")


volcano_plot <- function(cluster, comp, genes_up, genes_down) {

  labels <- c(
    "UPP" = "UPP",
    "UP" = "UP",
    "0" = "0",
    "DW" = "DW",
    "DWW" = "DWW"
  )

  de_data2 <- de_data[de_data$cluster == cluster &
                        de_data$annotation == 'annotation_intermediate' &
                        de_data$comp == comp &
                        de_data$sign %in% names(labels), c("p_val", "avg_log2FC", "sign", "comp", "gene")]
  data <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))
  data$sign <- factor(data$sign, levels = c("UPP", "UP", "0", "DW", "DWW"))
  data_dw <- data[data$gene %in% genes_down, ]
  data_up <- data[data$gene %in% genes_up, ]

  fig <- ggplot(data = data, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
    geom_point(size = 1) +
    scale_color_manual(values = colors_volcano, labels = labels) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica")) +
    guides(color = guide_legend(override.aes = list(shape = 1))) +
    theme(legend.position = "none") +
    scale_y_continuous(breaks = c(seq(0, 20, 5)), limits = c(0, 20))+
    scale_x_continuous(breaks = c(seq(-4, 7, 1)), limits = c(-4, 7))




  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", fill = "#911704" , size = 1) +
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38", size = 1) +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold',
                     position = position_nudge_repel(x = 1.5, y = 0.5))+

    geom_label_repel(data = data_dw, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold')+

    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)

}


fig4c_S1NR <- volcano_plot(cluster, comp, genes_up, genes_down)

print(fig4c_S1NR)

save_sizes(plot =fig4c_S1NR , filename = 'fig4c_S1NR', device = 'jpeg')
save_sizes(plot = fig4c_S1NR, filename = 'fig4c_S1NR', device = 'tiff')
save_sizes(plot = fig4c_S1NR, filename = 'fig4c_S1NR', device = 'svg')
save_sizes(plot = fig4c_S1NR, filename = 'fig4c_S1NR', device = 'pdf')

# Figure 3F  -------------------------------------------------------------------

# healthy stroma from Garrido-Trigo et al.
stroma_HC <- readRDS("~/data_Albas/HC_13122021/01_piezas_anotadas/stroma.RDS")

#
# get markers from this project and Garrido-Trigo et al.
#
stroma@active.ident <- stroma$annotation_refined
stroma_mk <- FindAllMarkers(stroma,  only.pos = TRUE,
                              min.pct = 0.25, thresh.use = 0.25)

stroma_HC@active.ident <- stroma_HC$annotation_refined
stroma_HC_mk <- FindAllMarkers(stroma_HC,  only.pos = TRUE,
                                 min.pct = 0.25, thresh.use = 0.25)

#
# Get S1 and IF markers
#
S1_UC <- stroma_HC_mk[stroma_HC_mk$cluster == 'S1' & stroma_HC_mk$p_val_adj < 0.001,'gene']
IF_UC <- stroma_mk[stroma_mk$cluster == 'Inflammatory fibroblasts' & stroma_mk$p_val_adj < 0.05,'gene']





# W0 vs POST
de_data <- readRDS('/home/acorraliza/TOFA_data/20220222_TOFAS_23/01_DE/REPASO/new_complete.RDS')
de <- de_data[de_data$cluster == 'S1' &
                de_data$annotation == 'annotation_refined' &
                de_data$comp %in% c('w0R_vs_POSTR', 'w0NR_vs_POSTNR', 'RPOST_NRPOST') &
                de_data$sign %in% c('UPP', 'DWW', 'UP', 'DW'),]
x <- list(w0R_vs_POSTR = de$gene[de$comp == 'w0R_vs_POSTR'],
          w0NR_vs_POSTNR = de$gene[de$comp == 'w0NR_vs_POSTNR'],
          RPOST_NRPOST = de$gene[de$comp == 'RPOST_NRPOST']
)

gene_list <- list()
gene_list$ALL_w0R_vs_POSTR <- de$gene[de$comp == 'w0R_vs_POSTR']
gene_list$ALL_w0R_vs_POSTR_UPP_DWW <- de$sign[de$comp == 'w0R_vs_POSTR']
gene_list$ALL_w0NR_vs_POSTNR <- de$gene[de$comp == 'w0NR_vs_POSTNR']
gene_list$ALL_w0NR_vs_POSTNR_UPP_DWW <- de$sign[de$comp == 'w0NR_vs_POSTNR']
gene_list$ALL_RPOST_NRPOST <- de$gene[de$comp == 'RPOST_NRPOST']
gene_list$ALL_RPOST_NRPOST_UPP_DWW <- de$sign[de$comp == 'RPOST_NRPOST']
gene_list$INTERSECT_w0R_vs_POSTR__RPOST_NRPOST <- intersect(de$gene[de$comp == 'RPOST_NRPOST'],
                                                            de$gene[de$comp == 'w0R_vs_POSTR'])
gene_list$INTERSECT_w0NR_vs_POSTNR__RPOST_NRPOST <- intersect(de$gene[de$comp == 'RPOST_NRPOST'],
                                                              de$gene[de$comp == 'w0NR_vs_POSTNR'])
gene_list$INTERSECT_w0R_vs_POSTR__w0NR_vs_POSTNR <- intersect(de$gene[de$comp == 'w0NR_vs_POSTNR'],
                                                              de$gene[de$comp == 'w0R_vs_POSTR'])
gene_list$INTERSECT_ALL <- intersect(intersect(de$gene[de$comp == 'w0NR_vs_POSTNR'],
                                               de$gene[de$comp == 'w0R_vs_POSTR']),
                                     de$gene[de$comp == 'RPOST_NRPOST'])
df <- data.frame(Reduce(cbind, gene_list))
colnames(df) <- names(gene_list)
for(i in c(1,3,5,7:10)){
  df[,i][duplicated(df[,i])] <- NA
}
head(df)

S1_1list <- gene_list$ALL_w0R_vs_POSTR[gene_list$ALL_w0R_vs_POSTR_UPP_DWW == 'UPP' | gene_list$ALL_w0R_vs_POSTR_UPP_DWW == 'UP' ]
S1_2list <- gene_list$ALL_w0NR_vs_POSTNR[gene_list$ALL_w0NR_vs_POSTNR_UPP_DWW == 'UPP' | gene_list$ALL_w0NR_vs_POSTNR_UPP_DWW == 'UP' ]

the_list <- vector(mode = 'list', length = 4)
the_list[[1]] <- S1_1list
the_list[[2]] <- S1_2list
the_list[[3]] <- S1_UC
the_list[[4]] <- IF_UC
names(the_list) <- c('S1_W0R_POSTR_UP',
                     'S1_W0NR_POSTNR_UP',
                     'S1_markers',
                     'IF_markers')
#
# jaccard
#

library(matchSCore2)
ms <- matchSCore2(gene_cl.ref = the_list[2:3],gene_cl.obs = the_list[c(1,4)],
                  ylab = "Markers", xlab = "DE genes")
ms$ggplot +
  theme_figure()

df <- as.data.frame(ms$JI.mat)
df$rownames <- rownames(df)
matrix <- tidyr::pivot_longer(df, cols = 1:2)

fig3d <- ggplot(matrix,
                aes(x = name, y = rownames, fill= value)) +
  geom_raster() +
  geom_text(mapping = aes(label = round(value, digits = 2)),  size = 10/.pt) +
  scale_fill_gradient(low  = 'white',
                      high = 'red',
                      limits = c(0,0.13),
                      name="Jaccard\nIndex",
                      breaks=c(0, 0.13),
                      labels=c("Min", 'Max')) +
  theme_figure() +
  theme(text = element_text(family = 'Helvetica', size = 10),
        legend.position = 'right',
        legend.title = element_text(family = 'Helvetica', size = 9))

save_sizes(plot = fig3d, filename = 'Figure_3D', device = 'jpeg')
save_sizes(plot = fig3d, filename = 'Figure_3D', device = 'tiff')
save_sizes(plot = fig3d, filename = 'Figure_3D', device = 'svg')
save_sizes(plot = fig3d, filename = 'Figure_3D', device = 'pdf')

fig3d_nt <- fig3d +
  theme(
    text = element_text(color = 'white'),
    axis.text = element_text(color = 'white'))


save_sizes(plot = fig3d_nt, filename = 'Figure_3D_no_text', device = 'jpeg')
save_sizes(plot = fig3d_nt, filename = 'Figure_3D_no_text', device = 'tiff')
save_sizes(plot = fig3d_nt, filename = 'Figure_3D_no_text', device = 'svg')
save_sizes(plot = fig3d_nt, filename = 'Figure_3D_no_text', device = 'pdf')
