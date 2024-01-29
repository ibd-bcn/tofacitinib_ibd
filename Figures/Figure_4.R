options(stringsAsFactors = FALSE)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readr)
library(dplyr)
library(ggrepel)
source('Figures/functions_plots.R')


## Data ------------------------------------------------------------------------
message('Loading data')

stroma <- readRDS('Analysis/data/00_annotation_process/00_anotadas/stroma.RDS')
stroma$annotation_intermediate <- gsub('Inflammatory_fibroblasts', 'Inflammatory fibroblasts', stroma$annotation_intermediate)
stroma$annotation_intermediate <- gsub('IER_fibroblasts', 'IER fibroblasts', stroma$annotation_intermediate)
stroma$pre_post <- plyr::mapvalues(stroma$week_3, from = c('W0', 'POST'), to = c('PRE', 'POST'))
stroma$pre_post <- factor(stroma$pre_post, levels = c('PRE', 'POST'))
stroma$response <- factor(stroma$response, levels = c('R', 'NR'))
stroma$subset <- factor(stroma$subset, levels = c('epi', 'stroma', 'plasmas', 'myeloids', 'cycling', 'tcells'))
# Figure 4A ----------------------------------------------------------------------------------
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

# Figure 4 B

# TO DO!

# Figure 4C  -------------------------------------------------------------------------

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

