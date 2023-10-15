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

## Figure 3D -------------------------------------------------------------------
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
# DE entre los clusters de Elisa. Genes UP en NRR que no son UPP en R i viceversa
# los queremos comparar con los marcadores M1 de tofa i los M2 de HC
# i Jaccard entre todo.


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

