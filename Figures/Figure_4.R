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

# Figure 4C Volcano plots stromal cells -------------------------------------------------------------------------

de_data <- readRDS('Analysis/data/01_DE/REPASO/new_complete.RDS')

# Volcano plots  myeloid cells

# S1 responders
cluster <- "S1"
comp <- "w0R_vs_POSTR"
filtered_genes_down <- c("ADIRF", "ADAMDEC1", "FOS")
filtered_genes_up <- c("CHI3L1", "HIF1A", "WNT5A", "IFI27", "IFITM3", "ISG15", "GBP1", "CCL19", "JAK1", "STAT1", "SOCS1",
                       "IL7R", "IL4R", "OSMR", "IL15RA", "IL13RA2")

volcano_plot <- function(cluster, comp, filtered_genes_down, filtered_genes_up) {

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

  response <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))

  if (cluster != "IDA macrophages") {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = c(seq(0, 50, 5)), limits = c(0, 50)) +
      scale_x_continuous(breaks = c(seq(-2, 2, 1)), limits = c(-2, 2))
  } else {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))
  }



  filtered_data <- response[response$gene %in% filtered_genes_down, ]
  filtered_data2 <- response[response$gene %in% filtered_genes_up, ]

  fig <- fig+ geom_label_repel(data = filtered_data, aes(label = gene, group = gene, col = sign), size = 6/.pt,
                               segment.color = "black",
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines"),
                               position = position_nudge_repel(x = 0, y = 0)) +
    geom_label_repel(data = filtered_data2, aes(label = gene, group = gene, col = sign), size = 6/.pt,
                     segment.color = "black",
                     # fill = colors_volcano[filtered_data$sign],
                     # color = "white",
                     # segment.color = "black",
                     fontface = 'bold',
                     # box.padding = unit(0.2, "lines"),
                     # point.padding = unit(0.5, "lines"),
                     position = position_nudge_repel(x = 0, y = 0)) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}


fig4c_S1 <- volcano_plot(cluster, comp, filtered_genes_down, filtered_genes_up)
print(fig4c_S1)

save_sizes(plot =fig4c_S1 , filename = 'fig4c_S1R', device = 'jpeg')
save_sizes(plot = fig4c_S1, filename = 'fig4c_S1R', device = 'tiff')
save_sizes(plot = fig4c_S1, filename = 'fig4c_S1R', device = 'svg')
save_sizes(plot = fig4c_S1, filename = 'fig4c_S1R', device = 'pdf')

# S1 non-responders
cluster <- "S1"
comp <- "w0NR_vs_POSTNR"
filtered_genes_down <- c("CHI3L1", "POSTN", "FTH1", "RPS4Y1", "TIMP1", "MT-CO3")
filtered_genes_up <- c("MT2A", "IFI27")


volcano_plot <- function(cluster, comp, filtered_genes_down, filtered_genes_up) {

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

  response <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))

  if (cluster != "IDA macrophages") {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = c(seq(0, 30, 5)), limits = c(0, 30))
  } else {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))
  }



  filtered_data <- response[response$gene %in% filtered_genes_down, ]
  filtered_data2 <- response[response$gene %in% filtered_genes_up, ]

  fig <- fig+ geom_label_repel(data = filtered_data, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                               segment.color = "black",
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines"),
                               position = position_nudge_repel(x = 0, y = 2)) +
    geom_label_repel(data = filtered_data2, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                     segment.color = "black",
                     # fill = colors_volcano[filtered_data$sign],
                     # color = "white",
                     # segment.color = "black",
                     fontface = 'bold',
                     # box.padding = unit(0.2, "lines"),
                     # point.padding = unit(0.5, "lines"),
                     position = position_nudge_repel(x = 0, y = 1)) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}

fig4c_S1NR <- volcano_plot(cluster, comp, filtered_genes_down, filtered_genes_up)
print(fig4c_S1NR)

save_sizes(plot =fig4c_S1NR , filename = 'fig4c_S1NR', device = 'jpeg')
save_sizes(plot = fig4c_S1NR, filename = 'fig4c_S1NR', device = 'tiff')
save_sizes(plot = fig4c_S1NR, filename = 'fig4c_S1NR', device = 'svg')
save_sizes(plot = fig4c_S1NR, filename = 'fig4c_S1NR', device = 'pdf')





# S2 Responders

cluster <- "S2"
comp <- "w0R_vs_POSTR"
filtered_genes_down <- c("HMGB1", "EGR1", "IGFBP3", "HIF1A", "ISG15", "SOCS1", "OSMR")
filtered_genes_up <- c("CHI3LI", "STAT1", "IL7R", "IL15RA")


volcano_plot <- function(cluster, comp, filtered_genes_down, filtered_genes_up) {

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

  response <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))

  if (cluster != "IDA macrophages") {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = c(seq(0, 30, 5)), limits = c(0, 30))
  } else {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))
  }



  filtered_data <- response[response$gene %in% filtered_genes_down, ]
  filtered_data2 <- response[response$gene %in% filtered_genes_up, ]

  fig <- fig+ geom_label_repel(data = filtered_data, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                               segment.color = "black",
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines"),
                               position = position_nudge_repel(x = 0, y = 2)) +
    geom_label_repel(data = filtered_data2, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                     segment.color = "black",
                     # fill = colors_volcano[filtered_data$sign],
                     # color = "white",
                     # segment.color = "black",
                     fontface = 'bold',
                     # box.padding = unit(0.2, "lines"),
                     # point.padding = unit(0.5, "lines"),
                     position = position_nudge_repel(x = 0, y = 1)) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}

fig4c_S2R <- volcano_plot(cluster, comp, filtered_genes_down, filtered_genes_up)
print(fig4c_S2R)

save_sizes(plot =fig4c_S2R , filename = 'fig4c_S2R', device = 'jpeg')
save_sizes(plot = fig4c_S2R, filename = 'fig4c_S2R', device = 'tiff')
save_sizes(plot = fig4c_S2R, filename = 'fig4c_S2R', device = 'svg')
save_sizes(plot = fig4c_S2R, filename = 'fig4c_S2R', device = 'pdf')

#S2 non-responders

cluster <- "S2"
comp <- "w0NR_vs_POSTNR"
filtered_genes_down <- c("PLCG2", "FTH1", "TMEM176A", "CXCL14")
filtered_genes_up <- c("GBP3", "GBP1", "MT-ND3", "MT-CO3")

volcano_plot <- function(cluster, comp, filtered_genes_down, filtered_genes_up) {

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

  response <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))

  if (cluster != "IDA macrophages") {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = c(seq(0, 30, 5)), limits = c(0, 30))
  } else {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))
  }



  filtered_data <- response[response$gene %in% filtered_genes_down, ]
  filtered_data2 <- response[response$gene %in% filtered_genes_up, ]

  fig <- fig+ geom_label_repel(data = filtered_data, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                               segment.color = "black",
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines"),
                               position = position_nudge_repel(x = 0, y = 2)) +
    geom_label_repel(data = filtered_data2, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                     segment.color = "black",
                     # fill = colors_volcano[filtered_data$sign],
                     # color = "white",
                     # segment.color = "black",
                     fontface = 'bold',
                     # box.padding = unit(0.2, "lines"),
                     # point.padding = unit(0.5, "lines"),
                     position = position_nudge_repel(x = 0, y = 1)) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}

fig4c_S2NR <- volcano_plot(cluster, comp, filtered_genes_down, filtered_genes_up )
print(fig4c_S2NR)

save_sizes(plot =fig4c_S2NR , filename = 'fig4c_S2NR', device = 'jpeg')
save_sizes(plot = fig4c_S2NR, filename = 'fig4c_S2NR', device = 'tiff')
save_sizes(plot = fig4c_S2NR, filename = 'fig4c_S2NR', device = 'svg')
save_sizes(plot = fig4c_S2NR, filename = 'fig4c_S2NR', device = 'pdf')

# S3 Responders

cluster <- "S3"
comp <- "w0R_vs_POSTR"
filtered_genes_down <- c("ADAMDEC1", "ABCA8", "CXCL12", "SELENOP",
                    "HIF1A", "ISG15", "STAT1", "IL15RA")
filtereD_genes_up <- c("CHI3L1", "CCL19")


volcano_plot <- function(cluster, comp, filtered_genes_down, filtered_genes_up) {

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

  response <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))

  if (cluster != "IDA macrophages") {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = c(seq(0, 30, 5)), limits = c(0, 30))
  } else {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))
  }



  filtered_data <- response[response$gene %in% filtered_genes_down, ]
  filtered_data2 <- response[response$gene %in% filtered_genes_up, ]

  fig <- fig+ geom_label_repel(data = filtered_data, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                               segment.color = "black",
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines"),
                               position = position_nudge_repel(x = 0, y = 2)) +
    geom_label_repel(data = filtered_data2, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                     segment.color = "black",
                     # fill = colors_volcano[filtered_data$sign],
                     # color = "white",
                     # segment.color = "black",
                     fontface = 'bold',
                     # box.padding = unit(0.2, "lines"),
                     # point.padding = unit(0.5, "lines"),
                     position = position_nudge_repel(x = 0, y = 1)) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}
fig4c_S3R <- volcano_plot(cluster, comp, filtered_genes_down, filtered_genes_up)
print(fig4c_S3R)

save_sizes(plot =fig4c_S3R , filename = 'fig4c_S3R', device = 'jpeg')
save_sizes(plot = fig4c_S3R, filename = 'fig4c_S3R', device = 'tiff')
save_sizes(plot = fig4c_S3R, filename = 'fig4c_S3R', device = 'svg')
save_sizes(plot = fig4c_S3R, filename = 'fig4c_S3R', device = 'pdf')

#S3 non-responders

cluster <- "S3"
comp <- "w0NR_vs_POSTNR"
filtered_genes_down <- c("CXCL14", "PLCG2", "RPS4Y1")
filtered_genes_up <- c("SOCS3", "MT-CO3", "MT-ATP6", "RPL9")

volcano_plot <- function(cluster, comp, filtered_genes_down, filtered_genes_up) {

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

  response <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))

  if (cluster != "IDA macrophages") {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = c(seq(0, 30, 5)), limits = c(0, 30))
  } else {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))
  }



  filtered_data <- response[response$gene %in% filtered_genes_down, ]
  filtered_data2 <- response[response$gene %in% filtered_genes_up, ]

  fig <- fig+ geom_label_repel(data = filtered_data, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                               segment.color = "black",
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines"),
                               position = position_nudge_repel(x = 0, y = 2)) +
    geom_label_repel(data = filtered_data2, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                     segment.color = "black",
                     # fill = colors_volcano[filtered_data$sign],
                     # color = "white",
                     # segment.color = "black",
                     fontface = 'bold',
                     # box.padding = unit(0.2, "lines"),
                     # point.padding = unit(0.5, "lines"),
                     position = position_nudge_repel(x = 0, y = 1)) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}

fig4c_S3NR <- volcano_plot(cluster, comp, filtered_genes_down, filtered_genes_up)
print(fig4c_S3NR)

save_sizes(plot =fig4c_S3NR , filename = 'fig4c_S3NR', device = 'jpeg')
save_sizes(plot = fig4c_S3NR, filename = 'fig4c_S3NR', device = 'tiff')
save_sizes(plot = fig4c_S3NR, filename = 'fig4c_S3NR', device = 'svg')
save_sizes(plot = fig4c_S3NR, filename = 'fig4c_S3NR', device = 'pdf')


# Inflammatory fibroblasts responders

cluster <- "Inflammatory_fibroblasts"
comp <- "w0R_vs_POSTR"
filtered_genes_up <- c("CHI3L1", "PDPN", "PLAU")

volcano_plot <- function(cluster, comp, filtered_genes_down, filtered_genes_up) {

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

  response <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))

  if (cluster != "IDA macrophages") {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = c(seq(0, 30, 5)), limits = c(0, 30))
  } else {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))
  }



  filtered_data <- response[response$gene %in% filtered_genes_down, ]
  filtered_data2 <- response[response$gene %in% filtered_genes_up, ]

  fig <- fig+ geom_label_repel(data = filtered_data, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                               segment.color = "black",
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines"),
                               position = position_nudge_repel(x = 0, y = 2)) +
    geom_label_repel(data = filtered_data2, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                     segment.color = "black",
                     # fill = colors_volcano[filtered_data$sign],
                     # color = "white",
                     # segment.color = "black",
                     fontface = 'bold',
                     # box.padding = unit(0.2, "lines"),
                     # point.padding = unit(0.5, "lines"),
                     position = position_nudge_repel(x = 0, y = 1)) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}

fig4c_IFR <- volcano_plot(cluster, comp, filtered_genes_up)
print(fig4c_IFR)

save_sizes(plot =fig4c_IFR , filename = 'fig4c_IFR', device = 'jpeg')
save_sizes(plot = fig4c_IFR, filename = 'fig4c_IFR', device = 'tiff')
save_sizes(plot = fig4c_IFR, filename = 'fig4c_IFR', device = 'svg')
save_sizes(plot = fig4c_IFR, filename = 'fig4c_IFR', device = 'pdf')

#Inflammatory_fibroblasts non-responders

cluster <- "Inflammatory_fibroblasts"
comp <- "w0NR_vs_POSTNR"
filtered_genes_down <- c("PLCG2", "IL32", "MT2A", "TNFRSF12A")
filtered_genes_up <- c("MT-ND1", "MT-CO3", "SOCS3", "CXCL13")

volcano_plot <- function(cluster, comp, filtered_genes_down, filtered_genes_up) {

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

  response <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))

  if (cluster != "IDA macrophages") {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = c(seq(0, 30, 5)), limits = c(0, 30))
  } else {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))
  }



  filtered_data <- response[response$gene %in% filtered_genes_down, ]
  filtered_data2 <- response[response$gene %in% filtered_genes_up, ]

  fig <- fig+ geom_label_repel(data = filtered_data, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                               segment.color = "black",
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines"),
                               position = position_nudge_repel(x = 0, y = 2)) +
    geom_label_repel(data = filtered_data2, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                     segment.color = "black",
                     # fill = colors_volcano[filtered_data$sign],
                     # color = "white",
                     # segment.color = "black",
                     fontface = 'bold',
                     # box.padding = unit(0.2, "lines"),
                     # point.padding = unit(0.5, "lines"),
                     position = position_nudge_repel(x = 0, y = 1)) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}

fig4c_IFNR <- volcano_plot(cluster, comp, filtered_genes_down, filtered_genes_up)
print(fig4c_IFNR)

save_sizes(plot =fig4c_IFNR , filename = 'fig4c_IFNR', device = 'jpeg')
save_sizes(plot = fig4c_IFNR, filename = 'fig4c_IFNR', device = 'tiff')
save_sizes(plot = fig4c_IFNR, filename = 'fig4c_IFNR', device = 'svg')
save_sizes(plot = fig4c_IFNR, filename = 'fig4c_IFNR', device = 'pdf')
