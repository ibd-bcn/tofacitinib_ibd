options(stringsAsFactors = FALSE)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)

source('Figures/functions_plots.R')

## Data ------------------------------------------------------------------------
message('Loading data')

all <- readRDS('/home/acorraliza/TOFA_data/20220222_TOFAS_23/00_annotation_process/00_anotadas/todas_2023.RDS')
all$pre_post <- plyr::mapvalues(all$week_3, from = c('W0', 'POST'), to = c('PRE', 'POST'))
all$pre_post <- factor(all$pre_post, levels = c('PRE', 'POST'))
all$response <- factor(all$response, levels = c('R', 'NR'))
all$subset <- factor(all$subset, levels = c('epi', 'stroma', 'plasmas', 'myeloids', 'cycling', 'tcells'))


## Figure 2 A: UMAP subset -----------------------------------------------------

fig <- DimPlot(all, group.by = 'subset') +
  scale_color_manual(values = colors_subset) +
  theme_umap()

save_sizes(plot = fig, filename = 'Figure_2A', device = 'tiff')
save_sizes(plot = fig, filename = 'Figure_2A', device = 'svg')
save_sizes(plot = fig, filename = 'Figure_2A', device = 'jpeg')
save_sizes(plot = fig, filename = 'Figure_2A', device = 'pdf')

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

