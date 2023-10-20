options(stringsAsFactors = FALSE)
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

all <- readRDS('/home/acorraliza/TOFA_data/20220222_TOFAS_23/00_annotation_process/00_anotadas/todas_2023.RDS')
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


## Figure 2 A: UMAP subset -----------------------------------------------------

fig <- DimPlot(all, group.by = 'subset') +
  scale_color_manual(values = colors_subset) +
  theme_umap()

save_sizes(plot = fig, filename = 'Figure_2A', device = 'tiff')
save_sizes(plot = fig, filename = 'Figure_2A', device = 'svg')
save_sizes(plot = fig, filename = 'Figure_2A', device = 'jpeg')
save_sizes(plot = fig, filename = 'Figure_2A', device = 'pdf')

fig_legend <- fig +
  theme(legend.position = 'right')

## Figure 2 B: UMAP AND BARPLOTS pre-post-response -----------------------------

#
# UMAP
#
umap2b <- DimPlot(all, group.by = 'response', split.by = 'pre_post') +
  scale_color_manual(values = colors_response) +
  theme_umap()

save_sizes(plot = umap2b, filename = 'Figure_2B_UMAP', device = 'tiff')
save_sizes(plot = umap2b, filename = 'Figure_2B_UMAP', device = 'svg')
save_sizes(plot = umap2b, filename = 'Figure_2B_UMAP', device = 'jpeg')
save_sizes(plot = umap2b, filename = 'Figure_2B_UMAP', device = 'pdf')

umap2b_legend <- umap2b +
  theme(legend.position = 'right')
#
# UMAPs 2a & 2b together
#
umaps_tog <- fig + umap2b + patchwork::plot_layout(ncol = 2, widths = c(0.33,0.66))

save_sizes(plot = umaps_tog, filename = 'Figure_2AB_UMAPS_together', device = 'tiff')
save_sizes(plot = umaps_tog, filename = 'Figure_2AB_UMAPS_together', device = 'svg')
save_sizes(plot = umaps_tog, filename = 'Figure_2AB_UMAPS_together', device = 'jpeg')
save_sizes(plot = umaps_tog, filename = 'Figure_2AB_UMAPS_together', device = 'pdf')

#
# Barplot
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

pat3 <- umap2b / pat
pat4 <- umap2b / pat2

save_sizes(plot = pat, filename = 'Figure_2B_BARPLOTS_text', device = 'tiff')
save_sizes(plot = pat, filename = 'Figure_2B_BARPLOTS_text', device = 'svg')
save_sizes(plot = pat, filename = 'Figure_2B_BARPLOTS_text', device = 'jpeg')
save_sizes(plot = pat, filename = 'Figure_2B_BARPLOTS_text', device = 'pdf')

save_sizes(plot = pat2, filename = 'Figure_2B_BARPLOTS', device = 'tiff')
save_sizes(plot = pat2, filename = 'Figure_2B_BARPLOTS', device = 'svg')
save_sizes(plot = pat2, filename = 'Figure_2B_BARPLOTS', device = 'jpeg')
save_sizes(plot = pat2, filename = 'Figure_2B_BARPLOTS', device = 'pdf')

save_sizes(plot = pat3, filename = 'Figure_2B_BARPLOTS_AND_UMAP_text', device = 'tiff')
save_sizes(plot = pat3, filename = 'Figure_2B_BARPLOTS_AND_UMAP_text', device = 'svg')
save_sizes(plot = pat3, filename = 'Figure_2B_BARPLOTS_AND_UMAP_text', device = 'jpeg')
save_sizes(plot = pat3, filename = 'Figure_2B_BARPLOTS_AND_UMAP_text', device = 'pdf')

save_sizes(plot = pat4, filename = 'Figure_2B_BARPLOTS_AND_UMAP', device = 'tiff')
save_sizes(plot = pat4, filename = 'Figure_2B_BARPLOTS_AND_UMAP', device = 'svg')
save_sizes(plot = pat4, filename = 'Figure_2B_BARPLOTS_AND_UMAP', device = 'jpeg')
save_sizes(plot = pat4, filename = 'Figure_2B_BARPLOTS_AND_UMAP', device = 'pdf')

# Figure 2 A and 2B together ---------------------------------------------------

design <- "#AA#
           #AA#
           BBBB
           BBBB
           CCCC
           CCCC
           CCCC"

fin <- fig + umap2b + pat + patchwork::plot_layout(design = design) & theme(text = element_text(family = 'Helvetica', size = 6))

ggsave(filename = 'Figure_2AB_text.jpeg', plot = fin, dpi = 300, device = 'jpeg', path = 'Figures/output/')

# Figure 2 C: Abundances -------------------------------------------------------


immune <- ggplot(abundances_data_immune[abundances_data_immune$Response != "PRE R VS NR", ],
            aes(Response, Cluster, fill= e.score)) +
  geom_raster() +
  geom_text(mapping = aes(label = ast_Chisq),  size = 10/.pt, nudge_y = -0.25) +
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
  labs(fill = 'Enrichment\nscore\nvs\nPRE')

plot(immune)

no_immune <- ggplot(abundances_data_noimmune[abundances_data_noimmune$Response != "PRE R VS NR", ],
                 aes(Response, Cluster, fill= e.score)) +
  geom_raster() +
  geom_text(mapping = aes(label = ast_Chisq), size = 10/.pt, nudge_y = -0.15) +
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
  labs(fill = 'Enrichment\nscore\nvs\nPRE')

plot(no_immune)

fig2c <- immune + no_immune &
  theme_figure() &
  theme(panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none')

save_sizes(plot = fig2c, filename = 'Figure_2C', device = 'tiff', max = 7)
save_sizes(plot = fig2c, filename = 'Figure_2C', device = 'svg', max = 7)
save_sizes(plot = fig2c, filename = 'Figure_2C', device = 'jpeg', max = 7)
save_sizes(plot = fig2c, filename = 'Figure_2C', device = 'pdf', max = 7)

a <- get_legend(fig2c + theme(legend.position = 'right'))

ggsave(a, filename = 'Figures/output/Figure_2C_legend.jpeg', dpi = 300)
ggsave(a, filename = 'Figures/output/Figure_2C_legend.pdf', dpi = 300)
ggsave(a, filename = 'Figures/output/Figure_2C_legend.svg', dpi = 300)
ggsave(a, filename = 'Figures/output/Figure_2C_legend.tiff', dpi = 300)


fig2c_legend <- fig2c +
  theme(panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        legend.position = 'right')

fig2c_legend & theme(text = element_text(color = 'white'), axis.text = element_text(color = 'white'))


t_immune <- ggplot(abundances_data_immune[abundances_data_immune$Response != "PRE R VS NR", ],
                 aes(Cluster,Response, fill= e.score)) +
  geom_raster() +
  geom_text(mapping = aes(label = ast_Chisq),  size = 10/.pt,angle = 90, vjust = 0.85) +
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
  labs(fill = 'Enrichment\nscore\nvs\nPRE')

t_no_immune <- ggplot(abundances_data_noimmune[abundances_data_noimmune$Response != "PRE R VS NR", ],
                    aes(Cluster, Response, fill= e.score)) +
  geom_raster() +
  geom_text(mapping = aes(label = ast_Chisq), size = 10/.pt, angle = 90, vjust = 0.85) +
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
  labs(fill = 'Enrichment\nscore\nvs\nPRE')

fig2ct <- t_immune / t_no_immune  &
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank()) &
  patchwork::plot_layout(guides = 'collect')

# Figure 2 ABC -----------------------------------------------------------------


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

fin <- fig_legend + umap2b_legend + pat + fig2ct +
  patchwork::plot_layout(design = design, guides = 'collect') &
  theme(text = element_text(family = 'Helvetica', size = 7), legend.text = element_text(family = 'Helvetica', size = 6))

ggsave(filename = 'Figure_2ABC_text_larga.pdf', plot = fin, dpi = 300, device = 'pdf',
       path = 'Figures/output/', width = 5, height = 9, units = 'in')

# Figure 3C --------------------------------------------------------------------


#
# qPCR Boxplots data
#

qpcr <- read_excel("Figures/extra_data/Base_de_datos_bx_reals_bcn_bram_boxplots.xlsx",
                   col_types = c("text", "text", "text",
                                 "text", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric","numeric","numeric","numeric"))

#
# qPCR Boxplots
#

qpcr$Response <- factor(qpcr$Response, levels = c("R", "NR"))
qpcr$Cohort <- factor(qpcr$Cohort, levels = c("BCN", "LEU"))
qpcr$Time <- factor(qpcr$Time, levels = c("Pre-tx", "Post-tx"))
qpcr$facet_group <- factor(qpcr$Response)


qpcr_r <- qpcr[qpcr$Response=="R",]
qpcr_nr <- qpcr[qpcr$Response == "NR",]
list_genes <- colnames(qpcr)[6:length(colnames(qpcr))-1]


# Iterate over genes to obtain qPCR boxplots
for(gene in list_genes){

  qpcr_plot <- boxplot_plot(qpcr_r,qpcr_nr,gene)
  print(qpcr_plot)

  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'jpeg')
  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'tiff')
  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'svg')
  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'pdf')

}
