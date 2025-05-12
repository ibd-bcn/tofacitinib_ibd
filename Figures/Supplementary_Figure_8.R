options(stringsAsFactors = FALSE,bitmapType = "cairo")

library(ggplot2)
library(readxl)
source('Figures/functions_plots.R')

## Data ------------------------------------------------------------------------

pathway_RESP_DWW_enrichGO_BP_fig4 <- read_excel("Figures/extra_data/pathway_RESP_DWW_enrichGO_BP_fig4.xlsx")

## Supplementary Figure 8-------------------------------------------------------

filas_seleccionadas <- c(1,4,7,8,9,11,13,16,21,27,31,47)

pathway <- pathway_RESP_DWW_enrichGO_BP_fig4[filas_seleccionadas,]

pathway$S1 <- rep("S1", nrow(pathway))



pathway <- pathway %>%
  arrange(pvalue)

pathway <- pathway %>%
  dplyr::mutate(Description = factor(Description,
                                     levels = unique(Description))) %>%
  dplyr::arrange(pvalue)

pathway$gene_ratio_ok <- as.numeric(sub("\\/.*", "", pathway$GeneRatio)) / as.numeric(sub(".*\\/", "", pathway$GeneRatio))
pathway$log10pval <- -log10(pathway$pvalue)

#
# Polar plot ----
#
library(plotrix)
colorsElena<-colorRampPalette(colors=c('#C1D4E8', '#83A9D3', '#0652A9'))
fun <- function(str, n) {gsub(paste0("([^ ]+( +[^ ]+){",n-1,"}) +"),
                              "\\1\n", str)}

taula <- arrange(pathway, gene_ratio_ok) %>% mutate(Color = colorsElena(n = nrow(taula)))
taula <- arrange(taula, desc(log10pval))
taula$nom <- fun(taula$Description, n = 3)

svg(filename = 'polar_plot.svg', onefile = T, family = 'helvetica')
par(cex.lab=0.6, mar=c(5.1,5,5,2.1), cex.axis = 0.8, cex = 1)
polar.plot(taula$log10pval, start = 90, rp.type = 'p',show.centroid = F,
           boxed.radial=F, labels = taula$nom, radial.lim=c(0,8),
           poly.col = NULL)

polar.plot(taula$log10pval, start = 90,
           rp.type = 's',show.centroid = F,
           boxed.radial=F, labels = NULL , radial.lim=c(0,360),
           add=T, point.symbols=19,
           point.col = as.character(taula$Color), cex=2.5)

dev.off()

taula <- taula %>% arrange(desc(gene_ratio_ok))
legend("topleft", legend=unique(taula$GeneRatio), pch=16, col=unique(taula$Color))

# no labels ------
plot.new()
taula <- taula %>% arrange(desc(log10pval))
polar.plot(taula$log10pval, start = 90, rp.type = 'p',show.centroid = F,
           boxed.radial=F,
           radial.lim=c(0,8), labels = rep('', 12),
           poly.col = NULL)
polar.plot(taula$log10pval, start = 90,
           rp.type = 's',show.centroid = F,
           boxed.radial=F, labels = NULL , radial.lim=c(0,8),
           add=T, point.symbols=19,
           point.col = as.character(taula$Color), cex=2.5)

## ggplot legend -----

a <- ggplot(data = taula, aes(x = S1 , y = Description,  color = gene_ratio_ok)) +
  geom_point() +
  scale_color_stepsn(colours = c('#C1D4E8', '#83A9D3', '#0652A9'), guide = 'colourbar') +
  theme(legend.title = element_blank(), legend.key.size = unit(1, units = 'cm'))

legend <- cowplot::get_legend(a)

plot(legend)
