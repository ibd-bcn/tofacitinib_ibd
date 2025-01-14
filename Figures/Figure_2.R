options(stringsAsFactors = FALSE,bitmapType = "cairo")
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readr)
library(dplyr)
library(readxl)


source('Figures/functions_plots.R')

## Data ------------------------------------------------------------------------
message('Loading data')

all <- readRDS('Analysis/data/00_annotation_process/todas_2023.RDS')
all$pre_post <- plyr::mapvalues(all$week_3, from = c('W0', 'POST'), to = c('PRE', 'POST'))
all$pre_post <- factor(all$pre_post, levels = c('PRE', 'POST'))
all$response <- factor(all$response, levels = c('R', 'NR'))
all$subset <- factor(all$subset, levels = c('epi', 'stroma', 'plasmas', 'myeloids', 'cycling', 'tcells'))

#
# Abundances data
#
cluster_per_heatmap <- read_delim("Figures/extra_data/clusters.csv",
                                  delim = ";", escape_double = FALSE, trim_ws = TRUE)
cluster_per_heatmap$Cluster_ok <- make.unique(cluster_per_heatmap$Cluster)
cluster_per_heatmap$Cluster_ok <- gsub(pattern = "Macrophage NRG1", replacement = 'IDA macrophage', x = cluster_per_heatmap$Cluster_ok)

abundances_data_immune <- readr::read_csv2(file = 'Figures/extra_data/Figure_2C_IMMUNE.csv')
abundances_data_immune$Cluster <- gsub(pattern = "Macrophage NRG1", replacement = 'IDA macrophage', x = abundances_data_immune$Cluster)
abundances_data_immune$Cluster <- factor(abundances_data_immune$Cluster,
                                  levels = rev(cluster_per_heatmap$Cluster_ok[cluster_per_heatmap$Plot == "IMMUNE" ]))
abundances_data_immune$Response <- factor(abundances_data_immune$Response,
                                          levels = c("PRE R VS NR", "R POST",  "NR POST"))

abundances_data_noimmune <- readr::read_csv2(file = 'Figures/extra_data/Figure_2C_NO_IMMUNE.csv')
abundances_data_noimmune$Cluster <- factor(abundances_data_noimmune$Cluster,
                                  levels = rev(cluster_per_heatmap$Cluster_ok[cluster_per_heatmap$Plot == "NO_IMMUNE" ]))
abundances_data_noimmune$Response <- factor(abundances_data_noimmune$Response,
                                          levels = c("PRE R VS NR", "R POST",  "NR POST"))


## Figure 2 A_1: UMAP AND BARPLOTS pre-post-response -----------------------------------------------------
all_pre <- all[,all$pre_post == "PRE"]
umap2a1 <- DimPlot(all_pre, group.by = 'response', split.by = 'pre_post') +
  scale_color_manual(values = colors_response) +
  theme_umap()

save_sizes(plot = umap2a1, filename = 'Figure_2A1', device = 'tiff')
save_sizes(plot = umap2a1, filename = 'Figure_2A1', device = 'svg')
save_sizes(plot = umap2a1, filename = 'Figure_2A1', device = 'jpeg')
save_sizes(plot = umap2a1, filename = 'Figure_2A1', device = 'pdf')

fig_legend <- umap2a1 +
  theme(legend.position = 'right')

## Figure 2 A_2: UMAP AND BARPLOTS pre-post-response -----------------------------

#
# UMAP
#
all_post <- all[,all$pre_post == "POST"]
umap2a2 <- DimPlot(all_post, group.by = 'response', split.by = 'pre_post') +
  scale_color_manual(values = colors_response) +
  theme_umap()

save_sizes(plot = umap2a2, filename = 'Figure_2B_UMAP', device = 'tiff')
save_sizes(plot = umap2a2, filename = 'Figure_2B_UMAP', device = 'svg')
save_sizes(plot = umap2a2, filename = 'Figure_2B_UMAP', device = 'jpeg')
save_sizes(plot = umap2a2, filename = 'Figure_2B_UMAP', device = 'pdf')

umap2a2_legend <- umap2a2 +
  theme(legend.position = 'right')
#
# UMAPs 2a1 & 2a2 together
#
umaps_tog <- umap2a1 + umap2a2 + patchwork::plot_layout(nrow = 2, widths = c(0.66,0.66))
umaps_tog_legend <- umap2a1 + umap2a2 +
  theme(legend.position = 'right')
save_sizes(plot = umaps_tog, filename = 'Figure_2AB_UMAPS_together', device = 'tiff')
save_sizes(plot = umaps_tog, filename = 'Figure_2AB_UMAPS_together', device = 'svg')
save_sizes(plot = umaps_tog, filename = 'Figure_2AB_UMAPS_together', device = 'jpeg')
save_sizes(plot = umaps_tog, filename = 'Figure_2AB_UMAPS_together', device = 'pdf')

#
# Barplot 2a3
#
list_plots <- vector(mode = 'list')

for(subset in levels(all$subset)){

  ppl <- ggplot(all@meta.data[all@meta.data$subset==subset,]) +
    geom_bar(aes(fill = response, x = pre_post),
             position = position_dodge()) +
    facet_grid(cols = vars(response)) +
    scale_fill_manual(values = colors_response)

  list_plots[[subset]]<-ppl

}

pat <- patchwork::wrap_plots(list_plots, ncol=2, guides='collect') &
  theme_figure()

pat2 <- patchwork::wrap_plots(list_plots, ncol=2, guides='collect') &
  theme_figure_wo_text()

pat3 <- umaps_tog / pat
pat4 <- umaps_tog / pat2

save_sizes(plot = pat, filename = 'Figure_2a3_BARPLOTS_text', device = 'tiff')
save_sizes(plot = pat, filename = 'Figure_2a3_BARPLOTS_text', device = 'svg')
save_sizes(plot = pat, filename = 'Figure_2a3_BARPLOTS_text', device = 'jpeg')
save_sizes(plot = pat, filename = 'Figure_2a3_BARPLOTS_text', device = 'pdf')

save_sizes(plot = pat2, filename = 'Figure_2a3_BARPLOTS', device = 'tiff')
save_sizes(plot = pat2, filename = 'Figure_2a3_BARPLOTS', device = 'svg')
save_sizes(plot = pat2, filename = 'Figure_2a3_BARPLOTS', device = 'jpeg')
save_sizes(plot = pat2, filename = 'Figure_2a3_BARPLOTS', device = 'pdf')

save_sizes(plot = pat3, filename = 'Figure_2a3_BARPLOTS_AND_UMAP_text', device = 'tiff')
save_sizes(plot = pat3, filename = 'Figure_2a3_BARPLOTS_AND_UMAP_text', device = 'svg')
save_sizes(plot = pat3, filename = 'Figure_2a3_BARPLOTS_AND_UMAP_text', device = 'jpeg')
save_sizes(plot = pat3, filename = 'Figure_2a3_BARPLOTS_AND_UMAP_text', device = 'pdf')

save_sizes(plot = pat4, filename = 'Figure_2a3_BARPLOTS_AND_UMAP', device = 'tiff')
save_sizes(plot = pat4, filename = 'Figure_2a3_BARPLOTS_AND_UMAP', device = 'svg')
save_sizes(plot = pat4, filename = 'Figure_2a3_BARPLOTS_AND_UMAP', device = 'jpeg')
save_sizes(plot = pat4, filename = 'Figure_2a3_BARPLOTS_AND_UMAP', device = 'pdf')

# Figure 2a1, 2a2 and 2a3 together ---------------------------------------------------

design <- "#AA#
           #AA#
           BBBB
           BBBB
           CCCC
           CCCC
           CCCC"

fin <-  umaps_tog  + pat + patchwork::plot_layout(design = design) & theme(text = element_text(family = 'Helvetica', size = 6))

ggsave(filename = 'Figure_2a123_text.jpeg', plot = fin, dpi = 300, device = 'jpeg', path = 'Figures/output/')

# Figure 2 B: Abundances -------------------------------------------------------


immune <- ggplot(abundances_data_immune[abundances_data_immune$Response != "PRE R VS NR", ],
            aes(Response, Cluster, fill= e.score)) +
  geom_raster() +
  geom_text(mapping = aes(label = ast),  size = 10/.pt, nudge_y = -0.25) +
  scale_fill_gradient2(low = '#023fa5', mid = '#FFFCFC',
                       high = '#8e063b',limits = c(-30,30),
                       midpoint = 0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(family = 'Helvetica',
                            size = 6),
        axis.title = element_blank(),
        legend.position = 'none') +
  labs(fill = 'Enrichment\nscore\nvs\nPre-tx')

plot(immune)

no_immune <- ggplot(abundances_data_noimmune[abundances_data_noimmune$Response != "PRE R VS NR", ],
                 aes(Response, Cluster, fill= e.score)) +
  geom_raster() +
  geom_text(mapping = aes(label = ast), size = 10/.pt, nudge_y = -0.15) +
  scale_fill_gradient2(low = '#023fa5', mid = '#FFFCFC',
                       high = '#8e063b',limits = c(-30,30),
                       midpoint = 0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(family = 'Helvetica',
                            size = 6),
        axis.title = element_blank(),
        legend.position = 'none') +
  labs(fill = 'Enrichment\nscore\nvs\nPre-tx')

plot(no_immune)

fig2b <- immune + no_immune &
  theme_figure() &
  theme(panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none')

save_sizes(plot = fig2b, filename = 'Figure_2b', device = 'tiff', max = 7)
save_sizes(plot = fig2b, filename = 'Figure_2b', device = 'svg', max = 7)
save_sizes(plot = fig2b, filename = 'Figure_2b', device = 'jpeg', max = 7)
save_sizes(plot = fig2b, filename = 'Figure_2b', device = 'pdf', max = 7)

a <- get_legend(fig2b + theme(legend.position = 'right'))

ggsave(a, filename = 'Figures/output/Figure_2b_legend.jpeg', dpi = 300)
ggsave(a, filename = 'Figures/output/Figure_2b_legend.pdf', dpi = 300)
ggsave(a, filename = 'Figures/output/Figure_2b_legend.svg', dpi = 300)
ggsave(a, filename = 'Figures/output/Figure_2b_legend.tiff', dpi = 300)

# Vertical figure

fig2b_legend <- fig2b +
  theme(panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        legend.position = 'right')

fig2b_legend & theme(text = element_text(color = 'white'), axis.text = element_text(color = 'white'))

# Horizontal figure

t_immune <- ggplot(abundances_data_immune[abundances_data_immune$Response != "PRE R VS NR", ],
                 aes(Cluster,Response, fill= e.score)) +
  geom_raster() +
  geom_text(mapping = aes(label = ast),  size = 10/.pt,angle = 90, vjust = 0.85) +
  scale_fill_gradient2(low = '#023fa5', mid = '#FFFCFC',
                       high = '#8e063b',limits = c(-30,30),
                       midpoint = 0) +
  theme_figure() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust= 1),
        legend.position = 'none') +
  labs(fill = 'Enrichment\nscore\nvs\nPre-tx')

t_no_immune <- ggplot(abundances_data_noimmune[abundances_data_noimmune$Response != "PRE R VS NR", ],
                    aes(Cluster, Response, fill= e.score)) +
  geom_raster() +
  geom_text(mapping = aes(label = ast), size = 10/.pt, angle = 90, vjust = 0.85) +
  scale_fill_gradient2(low = '#023fa5', mid = '#FFFCFC',
                       high = '#8e063b',limits = c(-30,30),
                       midpoint = 0) +
  theme_figure() +
  theme(panel.border = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
        axis.title = element_blank(),
        legend.position = 'right',
        legend.title = element_text(family = 'Helvetica', size = 6)) +
  labs(fill = 'Enrichment\nscore\nvs\nPre-tx')

fig2bt <- t_immune / t_no_immune  &
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank()) &
  patchwork::plot_layout(guides = 'collect')

# Figure 2 A y B -----------------------------------------------------------------


design <- "#AAAA#
           #AAAA#
           #AAAA#
           BBBBBB
           BBBBBB
           BBBBBB
           CCCCCC
           CCCCCC
           CCCCCC
           CCCCCC
           DDDDDD
           DDDDDD
           DDDDDD
           DDDDDD"

fin <- fig_legend + umaps_tog_legend + pat + fig2ct +
  patchwork::plot_layout(design = design, guides = 'collect') &
  theme(text = element_text(family = 'Helvetica', size = 7), legend.text = element_text(family = 'Helvetica', size = 6))

ggsave(filename = 'Figure_2AB_text_larga.pdf', plot = fin, dpi = 300, device = 'pdf',
       path = 'Figures/output/', width = 5, height = 9, units = 'in')


# Figure 2C --------------------------------------------------------------------


#
# qPCR Boxplots data
#

qpcr <- read_excel("Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
                   col_types = c("text", "text", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric"))

colnames(qpcr) <- gsub(" \\(AU\\)","",colnames(qpcr))
qpcr <- head(qpcr, -2)

#
# qPCR Boxplots
#

qpcr$Response <- factor(qpcr$Response, levels = c("R", "NR"))
qpcr$Time <- factor(qpcr$Time, levels = c("Pre-tx", "Post-tx"))
qpcr$facet_group <- factor(qpcr$Response)


qpcr_r <- qpcr[qpcr$Response=="R",]
qpcr_nr <- qpcr[qpcr$Response == "NR",]
list_genes <- c("CDH1","DERL3","PDGFRA","CHI3L1","PROK2","CD3E")


# Iterate over genes to obtain qPCR boxplots
for(gene in list_genes){

  qpcr_plot <- boxplot_plot(qpcr_r,qpcr_nr,gene)
  print(qpcr_plot)

  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'jpeg')
  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'tiff')
  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'svg')
  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'pdf')

}

