options(stringsAsFactors = FALSE,bitmapType = "cairo")
library(openxlsx)
library(ggplot2)
library(readr)
library(nlcor)
library(ggpubr)
library(patchwork)
library(readr)
library(ComplexHeatmap)
library(dplyr)
library(plyr)
library(circlize)
library(gridExtra)
library(grid)
library(cowplot)
source('Figures/functions_plots.R')


## Figure 1 C: Single-Cell UMAP-------------------------------------------------

replacements <- list(
  "M0" = "Macrophages",
  "M1" = "Macrophages",
  "M2" = "Macrophages",
  "Macrophage NRG1" = "Macrophages",
  "Rib hi myeloids" = "Macrophages",
  "pDC" = "DCs",
  "Plasmablast_IgG" = "Plasma cell",
  "PC_IgG" = "Plasma cell",
  "B_cell" = "Plasma cell",
  "PC_IgA" = "Plasma cell",
  "Plasmablast_IgA" = "Plasma cell",
  "PC_IGLL5" = "Plasma cell",
  "Cycling_B_cell" = "Plasma cell",
  "PC_IGLV6-57" = "Plasma cell",
  "PC_heat_shock" = "Plasma cell",
  "PC_IER" = "Plasma cell",
  "PC_PSAT1" = "Plasma cell",
  "Plasmablast_IGKC" = "Plasma cell",
  "PC_IFIT1" = "Plasma cell",
  "IER_fibroblasts" = "Fibroblasts",
  "S1" = "Fibroblasts",
  "S2" = "Fibroblasts",
  "S3" = "Fibroblasts",
  "Inflammatory_fibroblasts" = "Fibroblasts",
  "MT fibroblasts" = "Fibroblasts",
  "Perycites" = "Fibroblasts",
  "Myofibroblasts" = "Fibroblasts",
  "Cycling_fibroblasts" = "Fibroblasts",
  "Colonocytes" = "Epithelium",
  "Cycling TA" = "Epithelium",
  "Undifferentiated epithelium" = "Epithelium",
  "IER Epithelium" = "Epithelium",
  "Secretory progenitor" = "Epithelium",
  "Goblet" = "Epithelium",
  "Epithlium Rib hi" = "Epithelium",
  "Stem" = "Epithelium",
  "Mature goblet" = "Epithelium",
  "Tuft" = "Epithelium",
  "Enteroendocrines" = "Epithelium",
  "APOA4" = "Epithelium",
  "CD8" = "T cell",
  "CD4" = "T cell",
  "Ribhi T cells" = "T cell",
  "Tregs" = "T cell",
  "Thf" = "T cell",
  "Naïve T cells" = "T cell",
  "Th17" = "T cell",
  "MT hi IER" = "T cell",
  "HS" = "T cell",
  "NK" = "T cell",
  "CD4 CD8 IFIT3" = "T cell",
  "ILC3" = "T cell",
  "DN EOMES" = "T cell",
  "Cycling_B_cell" = "Cycling",
  "Cycling_PC" = "Cycling",
  "Cycling_T_cell" = "Cycling",
  "Cycling myeloids" = "Cycling",
  "Cycling_Plasma cell" = "Cycling",
  "T cell T cell IFIT3" = "T cell",
  "Inflammatory monocytes" = "Inf mono"
)

for (pattern in names(replacements)) {
  todas$annotation_intermediate <- gsub(pattern, replacements[[pattern]], todas$annotation_intermediate)
}

color_vector <- c(
  "Neutrophils" = "#BC4749",
  "Plasma cell" = "#965e97",
  "Fibroblasts" = "#93b7db",
  "Epithelium" = "#2f2a57",
  "T cell" = "#ecc45d",
  "Cycling" = "#7d9f56",
  "Macrophages" = "#C5D86D",
  "Endothelium" = "#D99697",
  "Glia" = "#F4C366",
  "Inf mono" = "#1B998B",
  "DCs" = "#EC769A",
  "Eosinophils" = "#6B5E62",
  "Mast" = "#D98560"
)


fig1c <- DimPlot(todas, group.by = 'annotation_intermediate') +
  scale_color_manual(values = color_vector) +
  theme_umap()

save_sizes(plot = fig1c, filename = 'fig1b', device = 'tiff')
save_sizes(plot = fig1c, filename = 'fig1b', device = 'svg')
save_sizes(plot = fig1c, filename = 'fig1b', device = 'jpeg')
save_sizes(plot = fig1c, filename = 'fig1b', device = 'pdf')

fig1c_legend <- fig1c +
  theme(legend.position = 'right')

save_sizes(plot = fig1c_legend, filename = 'fig1c_legend', device = 'tiff')
save_sizes(plot = fig1c_legend, filename = 'fig1c_legend', device = 'svg')
save_sizes(plot = fig1c_legend, filename = 'fig1c_legend', device = 'jpeg')
save_sizes(plot = fig1c_legend, filename = 'fig1c_legend', device = 'pdf')



## Figure 1 D: PROGENy UMAP-----------------------------------------------------
prog_sub <- readRDS("Analysis/PROGENy/data/progeny.RDS")
png(
  filename = "output/NFkB.png",
  width = 6,
  height = 5,
  units = "in",
  res = 800
)
feature_plot_progeny(todas = prog_sub, feat = "NFkB")
dev.off()
png(
  filename = "output//JAK.STAT.png",
  width = 6,
  height = 5,
  units = "in",
  res = 800
)
feature_plot_progeny(todas = prog_sub, feat = "JAK-STAT")
dev.off()


## Figure 1 E: PROGENy Heatmap--------------------------------------------------
png(
  filename = "output/heatmap_NFkB.png",
  width = 8,
  height = 4,
  units = "in",
  res = 800
)
heatmap_progeny(stimuli = "NFkB")
dev.off()
png(
  filename = "output/heatmap_JAK.STAT.png",
  width = 8,
  height = 4,
  units = "in",
  res = 800
)
heatmap_progeny(stimuli = "JAK-STAT")
dev.off()










# ## Figure 1 E: Correlation plot-------------------------------------------------
# corr_tofa <- read_delim("Figures/extra_data/corr_tofa_sangre.csv",
#                         delim = ";", escape_double = FALSE,
#                         locale = locale(decimal_mark = ",",
#                                         grouping_mark = "."),
#                         na = c("ND","n.d"), trim_ws = TRUE)
#
# SOCS1 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, SOCS1)) +
#   geom_point(aes(color = Response)) +
#   geom_smooth(formula = y ~ log(x), method = 'lm') +
#   scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
#   labs(x = 'Tofa conc (ng/ml)') +
#   theme_classic()+
#   theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))
#
# SOCS3 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, SOCS3)) +
#   geom_point(aes(color = Response)) +
#   geom_smooth(formula = y ~ log(x), method = 'lm') +
#   scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
#   labs(x = 'Tofa conc (ng/ml)') +
#   theme_classic()+
#   theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))
#
# IRF1 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, IRF1)) +
#   geom_point(aes(color = Response)) +
#   geom_smooth(formula = y ~ log(x), method = 'lm') +
#   scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
#   labs(x = 'Tofa conc (ng/ml)') +
#   theme_classic()+
#   theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))
#
#
# patchwork::wrap_plots(SOCS1, SOCS3, IRF1, guides = 'collect') &
#   theme(legend.position = 'top', axis.text = element_text(colour = 'black', size = 8))
