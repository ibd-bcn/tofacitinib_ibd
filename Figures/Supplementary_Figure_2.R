

# Dotplot figure by subset

library(Seurat)
# Myeloid

myeloids <- readRDS('Analysis/data/00_annotation_process/00_anotadas/myeloids.RDS')
myeloids$annotation_intermediate  <- gsub('Macrophage NRG1', 'IDA macrophages', myeloids$annotation_intermediate)
myeloids$annotation_intermediate <- factor(myeloids$annotation_intermediate,
                                           levels = c("Cycling myeloids", "DCs", "Eosinophils", "HS",
                                                      "Inflammatory monocytes", "M0", "M1", "M2",
                                                      "IDA macrophages", "Mast", "Neutrophils", "pDC",
                                                      "Rib hi myeloids"))




DotPlot(myeloids, features = c(),
        group.by = 'annotation_intermediate', cols = 'RdYlBu',
        cluster.idents = F) +
  theme(axis.text.x = element_text(angle=90,  hjust = 1, vjust = 0.2, size = 9),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5), legend.key.size = unit(0.25, 'cm')) +
  theme(axis.text.y = element_text(angle=0, size = 6)) +
  theme(axis.title = element_blank()) +
  guides(color=guide_legend(title="Exp"),
         size = guide_legend(title = '%'))

