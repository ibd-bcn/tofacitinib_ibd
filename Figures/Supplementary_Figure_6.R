options(stringsAsFactors = FALSE,bitmapType = "cairo")

library(ggplot2)
library(readxl)
source('Figures/functions_plots.R')

## Data ------------------------------------------------------------------------
# Volcano plots myeloid subset

de_data <- readRDS('Figures/extra_data/new_complete.RDS')
de_data$cluster <- gsub('Macrophage NRG1', 'IDA macrophages', de_data$cluster)

# Supplementary Figure 5A-------------------------------------------------------

#Responders
# M0 Responders

cluster <- "M0"
comp <- "w0R_vs_POSTR"
genes_down <- c("S100A9","HLA-A", "HLA-C", "STAT1", "ISG15",
                "IFI27" )
genes_up <- c("RPL39", "RPS24", "FABP1")

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


sup6_M0R <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup6_M0R)

save_sizes(plot =sup6_M0R , filename = 'sup6_M0R', device = 'jpeg')
save_sizes(plot = sup6_M0R, filename = 'sup6_M0R', device = 'tiff')
save_sizes(plot = sup6_M0R, filename = 'sup6_M0R', device = 'svg')
save_sizes(plot = sup6_M0R, filename = 'sup6_M0R', device = 'pdf')


# IDA macrophages

cluster <- "IDA macrophages"
comp <- "w0R_vs_POSTR"
genes_down <- c("S100A9", "HLA-B", "GBP1", "IFITM3", "ISG15")
genes_up <- c( "SELENOP", "C1QA", "C1QC", "RPS28")

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
    scale_x_continuous(breaks = c(seq(-4, 4, 1)), limits = c(-4, 4)) +
    scale_y_continuous(breaks = c(seq(0, 15, 5)), limits = c(0, 15))


  fig <- fig+ geom_point(data = data_up,shape = 21,color = "black", fill = "#911704" ) +
    geom_point(data = data_dw,shape = 21, color = "black", fill = "#376D38") +
    geom_label_repel(data = data_up, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold') +

    geom_label_repel(data = data_dw, aes(label = gene, group = gene, col = sign), size = 9/.pt,
                     segment.color = "black",
                     fontface = 'bold',
                     position = position_nudge_repel(x = -0.6, y=0))+

    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}


sup6_IDAR <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup6_IDAR)

save_sizes(plot =sup6_IDAR , filename = 'sup6_IDAR', device = 'jpeg')
save_sizes(plot = sup6_IDAR, filename = 'sup6_IDAR', device = 'tiff')
save_sizes(plot = sup6_IDAR, filename = 'sup6_IDAR', device = 'svg')
save_sizes(plot = sup6_IDAR, filename = 'sup6_IDAR', device = 'pdf')


# Supplementary Figure 5B-------------------------------------------------------
#Non responders

#M0 non-responders

cluster <- "M0"
comp <- "w0NR_vs_POSTNR"
genes_up <- c("MMP9", "MALAT1", "NEAT1", "C15orf48", "FABP5")
genes_down <- c("CYBA", "MTRNR2L8", "MTRNR2L12", "RPS2")

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
                     position = position_nudge_repel(x = -0.5, y = 0.2))+

    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}


sup6_M0NR <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup6_M0NR)

save_sizes(plot =sup6_M0NR , filename = 'sup6_M0NR', device = 'jpeg')
save_sizes(plot = sup6_M0NR, filename = 'sup6_M0NR', device = 'tiff')
save_sizes(plot = sup6_M0NR, filename = 'sup6_M0NR', device = 'svg')
save_sizes(plot = sup6_M0NR, filename = 'sup6_M0NR', device = 'pdf')

# IDA macrophages non-responders

cluster <- "IDA macrophages"
comp <- "w0NR_vs_POSTNR"
genes_up <- c("PLGC2", "C15orf48", "MT-ATP8", "MYL6", "S100A11")
genes_down <- c( "TPT1","RPL9", "HERPUD1", "RETN")

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
    scale_x_continuous(breaks = c(seq(-2, 2, 1)), limits = c(-2, 2))

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


sup6_IDANR <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup6_IDANR)

save_sizes(plot =sup6_IDANR , filename = 'sup6_IDANR', device = 'jpeg')
save_sizes(plot = sup6_IDANR, filename = 'sup6_IDANR', device = 'tiff')
save_sizes(plot = sup6_IDANR, filename = 'sup6_IDANR', device = 'svg')
save_sizes(plot = sup6_IDANR, filename = 'sup6_IDANR', device = 'pdf')



# INHBA + Macrophages

# non-responders
cluster <- "M1"
comp <- "w0NR_vs_POSTNR"
genes_up <- c("PLCG2", "MIF", "FTH1", "GAPDH")
genes_down <- c("MT-CO3", "MT-CO2", "IFI6", "OAS1", "S100A9", "MX1")


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


sup6_M1NR <- volcano_plot(cluster, comp, genes_up, genes_down)

print(sup6_M1NR)

save_sizes(plot =sup6_M1NR , filename = 'sup6_M1NR', device = 'jpeg')
save_sizes(plot = sup6_M1NR, filename = 'sup6_M1NR', device = 'tiff')
save_sizes(plot = sup6_M1NR, filename = 'sup6_M1NR', device = 'svg')
save_sizes(plot = sup6_M1NR, filename = 'sup6_M1NR', device = 'pdf')










