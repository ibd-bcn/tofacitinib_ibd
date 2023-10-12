options(stringsAsFactors = FALSE)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readr)
library(dplyr)

source('Figures/functions_plots.R')

## Data ------------------------------------------------------------------------
message('Loading data')

myeloids <- readRDS('/home/acorraliza/TOFA_data/20220222_TOFAS_23/00_annotation_process/00_anotadas/myeloids.RDS')
myeloids$pre_post <- plyr::mapvalues(myeloids$week_3, from = c('W0', 'POST'), to = c('PRE', 'POST'))
myeloids$annotation_intermediate <- gsub('Macrophage NRG1', 'IDA macrophages', myeloids$annotation_intermediate)
myeloids$pre_post <- factor(myeloids$pre_post, levels = c('PRE', 'POST'))
myeloids$response <- factor(myeloids$response, levels = c('R', 'NR'))
myeloids$subset <- factor(myeloids$subset, levels = c('epi', 'stroma', 'plasmas', 'myeloids', 'cycling', 'tcells'))

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
                    group.by = 'annotation_intermediate', pt.size = 0.5) +
  scale_color_manual(values = colors, na.value = '#e6e5e5') +
  theme_umap()

nresp_pre <- DimPlot(myeloids[,myeloids@meta.data$response == 'NR' & myeloids@meta.data$pre_post == 'PRE'],
                    group.by = 'annotation_intermediate', pt.size = 0.5) +
  scale_color_manual(values = colors, na.value = '#e6e5e5') +
  theme_umap()

resp_post <- DimPlot(myeloids[,myeloids@meta.data$response == 'R' & myeloids@meta.data$pre_post == 'POST'],
                    group.by = 'annotation_intermediate', pt.size = 0.5) +
  scale_color_manual(values = colors, na.value = '#e6e5e5') +
  theme_umap()

nresp_post <- DimPlot(myeloids[,myeloids@meta.data$response == 'NR' & myeloids@meta.data$pre_post == 'POST'],
                    group.by = 'annotation_intermediate', pt.size = 0.5) +
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

#
# Pendiente!
#

## Figure 3C -------------------------------------------------------------------

fig3c <- feature_plot(myeloids,
             feature = c('MT2A','CLEC5A','IDO1', 'MMP9'),
             split.by = c('pre_post', 'response'),
             size = 0.1) +
  patchwork::plot_layout(ncol = 2) &
  theme_figure() &
  theme(
    title = element_text(family = 'Helvetica', size = 10, colour = 'black', hjust = 0),
    plot.title = element_text(family = 'Helvetica', size = 10, colour = 'black'),
    axis.text = element_text(family = 'Helvetica', size = 6, colour = 'black'),
    strip.text.x = element_text(family = 'Helvetica', size = 8, colour = "black"),
    strip.text.y = element_text(family = 'Helvetica', size = 8, colour = "black"),
    legend.position = 'none',
    strip.background = element_blank()
  )

fig3c_nt <- feature_plot(myeloids,
                      feature = c('MT2A','CLEC5A','IDO1', 'MMP9'),
                      split.by = c('pre_post', 'response'),
                      size = 0.1) +
  patchwork::plot_layout(ncol = 2) &
  theme_figure() &
  theme(
    title = element_text(family = 'Helvetica', size = 10, colour = 'white', hjust = 0),
    plot.title = element_text(family = 'Helvetica', size = 10, colour = 'white'),
    axis.text = element_text(family = 'Helvetica', size = 6, colour = 'white'),
    strip.text.x = element_text(family = 'Helvetica', size = 8, colour = "white"),
    strip.text.y = element_text(family = 'Helvetica', size = 8, colour = "white"),
    legend.position = 'none',
    strip.background = element_blank()
  )

save_sizes(plot = fig3c, filename = 'Figure_3C', device = 'jpeg')
save_sizes(plot = fig3c, filename = 'Figure_3C', device = 'tiff')
save_sizes(plot = fig3c, filename = 'Figure_3C', device = 'svg')
save_sizes(plot = fig3c, filename = 'Figure_3C', device = 'pdf')

save_sizes(plot = fig3c_nt, filename = 'Figure_3C_no_text', device = 'jpeg')
save_sizes(plot = fig3c_nt, filename = 'Figure_3C_no_text', device = 'tiff')
save_sizes(plot = fig3c_nt, filename = 'Figure_3C_no_text', device = 'svg')
save_sizes(plot = fig3c_nt, filename = 'Figure_3C_no_text', device = 'pdf')

legend <- fig3c & theme(legend.position = 'right')
leg <- get_legend(legend)

plot_legend <- patchwork::wrap_plots(leg)
save_sizes(plot = plot_legend, filename = 'Figure_3C_legend', device = 'pdf')
save_sizes(plot = plot_legend, filename = 'Figure_3C_legend', device = 'jpeg')
save_sizes(plot = plot_legend, filename = 'Figure_3C_legend', device = 'tiff')
save_sizes(plot = plot_legend, filename = 'Figure_3C_legend', device = 'svg')

