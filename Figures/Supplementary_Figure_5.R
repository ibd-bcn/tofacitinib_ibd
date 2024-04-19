options(stringsAsFactors = FALSE,bitmapType = "cairo")

library(ggplot2)
library(readxl)
source('Figures/functions_plots.R')

## Data ------------------------------------------------------------------------
# Volcano plots stroma subset

de_data <- readRDS('Analysis/data/01_DE/REPASO/new_complete.RDS')

## Supplementary Figure 5A-------------------------------------------------------

# S2 Responders

cluster <- "S2"
comp <- "w0R_vs_POSTR"
genes_up <- c("HMGB1", "EGR1", "IGFBP3")
genes_down <- c("CHI3LI", "STAT1", "IL7R", "IL15RA", "OSMR","ISG15", "SOCS1", "HIF1A")


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
    scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))

  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", fill = "#911704" ) +
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38") +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold') +
    geom_label_repel(data = data_dw, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold',
                     position = position_nudge_repel(x = -0.5, y =0))+
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )


  return(fig)
}


sup5_S2R <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup5_S2R)

save_sizes(plot =sup5_S2R , filename = 'sup5_S2R', device = 'jpeg')
save_sizes(plot = sup5_S2R, filename = 'sup5_S2R', device = 'tiff')
save_sizes(plot = sup5_S2R, filename = 'sup5_S2R', device = 'svg')
save_sizes(plot = sup5_S2R, filename = 'sup5_S2R', device = 'pdf')

# S3 Responders

cluster <- "S3"
comp <- "w0R_vs_POSTR"
genes_up <- c("ADAMDEC1", "ABCA8", "CXCL12", "SELENOP")
genes_down <- c("CHI3L1", "CCL19", "STAT1", "HIF1A", "ISG15", "IL15RA")


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
    scale_y_continuous(breaks = c(seq(0, 40, 5)), limits = c(0, 40)) +
    scale_x_continuous(breaks = c(seq(-5, 5, 1)), limits = c(-5, 5))


  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", fill = "#911704" ) +
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38") +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold',
                     position = position_nudge_repel(x = 1, y = 0)) +

    geom_label_repel(data = data_dw, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold',
                     position = position_nudge_repel(x = -1, y = 0))+
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}


sup5_S3R <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup5_S3R)

save_sizes(plot =sup5_S3R , filename = 'sup5_S3R', device = 'jpeg')
save_sizes(plot = sup5_S3R, filename = 'sup5_S3R', device = 'tiff')
save_sizes(plot = sup5_S3R, filename = 'sup5_S3R', device = 'svg')
save_sizes(plot = sup5_S3R, filename = 'sup5_S3R', device = 'pdf')


# Inflammatory fibroblasts responders

cluster <- "Inflammatory_fibroblasts"
comp <- "w0R_vs_POSTR"
genes_down <- c("CHI3L1", "PDPN", "PLAU")

volcano_plot <- function(cluster, comp, genes_down) {

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
    scale_y_continuous(breaks = c(seq(0, 10, 5)), limits = c(0, 10)) +
    scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))




  fig <- fig+
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38") +
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


sup5_IFR <- volcano_plot(cluster, comp, genes_down)

print(sup5_IFR)

save_sizes(plot =sup5_IFR , filename = 'sup5_IFR', device = 'jpeg')
save_sizes(plot = sup5_IFR, filename = 'sup5_IFR', device = 'tiff')
save_sizes(plot = sup5_IFR, filename = 'sup5_IFR', device = 'svg')
save_sizes(plot = sup5_IFR, filename = 'sup5_IFR', device = 'pdf')

## Supplementary Figure 5B------------------------------------------------------

#S2 non-responders

cluster <- "S2"
comp <- "w0NR_vs_POSTNR"
genes_up <- c("PLCG2", "FTH1", "TMEM176A", "CXCL14")
genes_down <- c("GBP3", "GBP1", "MT-ND3", "MT-CO3")

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
    scale_y_continuous(breaks = c(seq(0, 30, 5)), limits = c(0, 30))


  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", fill = "#911704" ) +
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38") +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold') +
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


sup5_S2NR <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup5_S2NR)

save_sizes(plot =sup5_S2NR , filename = 'sup5_S2NR', device = 'jpeg')
save_sizes(plot = sup5_S2NR, filename = 'sup5_S2NR', device = 'tiff')
save_sizes(plot = sup5_S2NR, filename = 'sup5_S2NR', device = 'svg')
save_sizes(plot = sup5_S2NR, filename = 'sup5_S2NR', device = 'pdf')


#S3 non-responders

cluster <- "S3"
comp <- "w0NR_vs_POSTNR"
genes_up <- c("CXCL14", "PLCG2", "RPS4Y1")
genes_down <- c("SOCS3", "MT-CO3", "MT-ATP6", "RPL9")

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
    scale_y_continuous(breaks = c(seq(0, 15, 5)), limits = c(0, 15)) +
    scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))


  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", fill = "#911704" ) +
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38") +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold') +
    geom_label_repel(data = data_dw, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold',
                     position = position_nudge_repel(x = -1, y = 0.2))+

    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}


sup5_S3NR <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup5_S3NR)

save_sizes(plot =sup5_S3NR , filename = 'sup5_S3NR', device = 'jpeg')
save_sizes(plot = sup5_S3NR, filename = 'sup5_S3NR', device = 'tiff')
save_sizes(plot = sup5_S3NR, filename = 'sup5_S3NR', device = 'svg')
save_sizes(plot = sup5_S3NR, filename = 'sup5_S3NR', device = 'pdf')


#Inflammatory_fibroblasts non-responders

cluster <- "Inflammatory_fibroblasts"
comp <- "w0NR_vs_POSTNR"
genes_up <- c("PLCG2", "IL32", "MT2A", "TNFRSF12A")
genes_down <- c("MT-ND1", "MT-CO3", "SOCS3", "CXCL13")

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
    scale_y_continuous(breaks = c(seq(0, 40, 5)), limits = c(0, 40)) +
    scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))


  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", fill = "#911704" ) +
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38") +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                     segment.color = "black",
                     fontface = 'bold') +

    geom_label_repel(data = data_dw, aes(label = gene, group = gene, col = sign), size = 8/.pt,
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


sup5_IFNR <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup5_IFNR)

save_sizes(plot =sup5_IFNR , filename = 'sup5_IFNR', device = 'jpeg')
save_sizes(plot = sup5_IFNR, filename = 'sup5_IFNR', device = 'tiff')
save_sizes(plot = sup5_IFNR, filename = 'sup5_IFNR', device = 'svg')
save_sizes(plot = sup5_IFNR, filename = 'sup5_IFNR', device = 'pdf')

