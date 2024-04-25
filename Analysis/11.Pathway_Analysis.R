# Pathway analysis for S1 cells

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

de_data <- readRDS('Analysis/data/01_DE/REPASO/new_complete.RDS')
de_data_resp <- de_data[de_data$cluster == 'S1' & de_data$annotation == 'annotation_intermediate' & de_data$comp == 'w0R_vs_POSTR',] %>%
  arrange(desc(avg_log2FC))
de_data_nresp <- de_data[de_data$cluster == 'S1' & de_data$annotation == 'annotation_intermediate' & de_data$comp == 'w0NR_vs_POSTNR',] %>%
  arrange(desc(avg_log2FC))

# upp --------------------------------------------------------------------------
# Extraer resultados significativos
signif_res <- de_data_resp[de_data_resp$sign == 'UPP', ]
# Vector con genes
signif_genes <- as.character(signif_res$gene)
ego <- enrichGO(gene = signif_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db,
                ont = "BP", # CC y MF disponibles
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = F)

write.csv(ego@result, file = 'Analysis/pathway_RESP_UPP_enrichGO_BP.csv')

# DWW --------------------------------------------------------------------------
# Extraer resultados significativos
signif_res <- de_data_resp[de_data_resp$sign == 'DWW', ]
# Vector con genes
signif_genes <- as.character(signif_res$gene)
ego <- enrichGO(gene = signif_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db,
                ont = "BP", # CC y MF disponibles
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = F)

write.csv(ego@result, file = 'Analysis/pathway_RESP_DWW_enrichGO_BP.csv')

# gse todo --------------------------------------------------------------------
geneList <- de_data_resp$avg_log2FC
names(geneList) <- de_data_resp$gene
# Ejecutar análisis de enriquecimiento GO
ego3 <- gseGO(geneList     = geneList,
              keyType      = 'SYMBOL',
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
write.csv(ego3, file = 'Analysis/pathway_GSE_BP.csv')



# upp NR --------------------------------------------------------------------------
# Extraer resultados significativos
signif_res <- de_data_nresp[de_data_nresp$sign == 'UPP', ]
# Vector con genes
signif_genes <- as.character(signif_res$gene)
ego <- enrichGO(gene = signif_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db,
                ont = "BP", # CC y MF disponibles
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = F)

write.csv(ego@result, file = 'Analysis/pathway_NRESP_UPP_enrichGO_BP.csv')

# DWW NR--------------------------------------------------------------------------
# Extraer resultados significativos
signif_res <- de_data_nresp[de_data_nresp$sign == 'DWW', ]
# Vector con genes
signif_genes <- as.character(signif_res$gene)
ego <- enrichGO(gene = signif_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db,
                ont = "BP", # CC y MF disponibles
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = F)

write.csv(ego@result, file = 'Analysis/pathway_NRESP_DWW_enrichGO_BP.csv')

# gse todo --------------------------------------------------------------------
geneList <- de_data_nresp$avg_log2FC
names(geneList) <- de_data_nresp$gene
# Ejecutar análisis de enriquecimiento GO
ego3 <- gseGO(geneList     = geneList,
              keyType      = 'SYMBOL',
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
write.csv(ego3, file = 'Analysis/pathway_NRESP_GSE_BP.csv')



## Figure ----------------------------------------------------------------------



library(readxl)
pathway_RESP_DWW_enrichGO_BP_fig4 <- read_excel("Figures/extra_data/pathway_RESP_DWW_enrichGO_BP_fig4.xlsx")


## 2,5,8,9,10, 12. 14, 17 , 22, 28, 32, 48 hay que quitarle 1 porque es excel
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
colorsElena<-colorRampPalette(colors=c('Forestgreen', 'white', 'Firebrick4'))
fun <- function(str, n) {gsub(paste0("([^ ]+( +[^ ]+){",n-1,"}) +"),
                              "\\1\n", str)}

taula <- arrange(pathway, gene_ratio_ok) %>% mutate(Color = colorsElena(n = nrow(taula)))
taula <- arrange(taula, desc(log10pval))
taula$nom <- fun(taula$Description, n = 3)
theta=seq(0, 360, length = nrow(taula)+1)[-1] -15
# theta10=seq(0, 360, length = 10 + 1)[-(10+1)]
# theta15=seq(0, 360, length = 15 + 1)[-(15+1)]

svg(filename = 'polar_plot.svg', onefile = T, family = 'helvetica')
par(cex.lab=0.6, mar=c(5.1,5,5,2.1), cex.axis = 0.8, cex = 1)
polar.plot(taula$log10pval, start = 90, rp.type = 'p',show.centroid = F,
           boxed.radial=F, labels = taula$nom, radial.lim=c(0,8),
           # label.pos = 1:12,
           poly.col = NULL)

polar.plot(taula$log10pval, start = 90,
           rp.type = 's',show.centroid = F,
           boxed.radial=F, labels = NULL , radial.lim=c(0,360),
           add=T, point.symbols=19,
           point.col = as.character(taula$Color), cex=2.5)

dev.off()

plot.new()
polar.plot(taula$log10pval, start = 90, rp.type = 'p',show.centroid = F,
           boxed.radial=F,
           # labels = taula$nom,
           radial.lim=c(0,8), labels = rep('', 12),
           # label.pos = 1:12,
           poly.col = NULL)
polar.plot(taula$log10pval, start = 90,
           rp.type = 's',show.centroid = F,
           boxed.radial=F, labels = NULL , radial.lim=c(0,8),
           add=T, point.symbols=19,
           point.col = as.character(taula$Color), cex=2.5)

## dotplot ---


ggplot(pathway, aes(x = S1 , y = Description, size = gene_ratio_ok, color = pvalue)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  labs(title = "", x = "", y = "Pathway description", size = "Gene ratio", color = "p-value")+
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica"),  # Set font to Helvetica
    plot.title = element_text(face = "bold", size = 14),  # Customize title font
    axis.title.x = element_text(face = "bold", size = 12),  # Customize x-axis label font
    axis.title.y = element_text(face = "bold", size = 12),  # Customize y-axis label font
    axis.text = element_text(size = 10),  # Customize axis text font
    legend.text = element_text(size = 10),  # Customize legend text font
    legend.title = element_text(face = "bold", size = 10)  # Customize legend title font
  )


