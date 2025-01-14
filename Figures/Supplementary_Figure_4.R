options(stringsAsFactors = FALSE,bitmapType = "cairo")

library(ggplot2)
library(readxl)
source('Figures/functions_plots.R')

## Supplementary Figure 4-------------------------------------------------------

#
# qPCR Boxplots data
#

qpcr <- read_excel("Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
                   sheet = "Sup Figure 4", col_types = c("text",
                                                         "text", "numeric", "numeric", "numeric",
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
list_genes <- c("TFF3","FCGR3A","CLEC5A","CXCL10","S100A8")


# Iterate over genes to obtain qPCR boxplots
for(gene in list_genes){

  qpcr_plot <- boxplot_plot(qpcr_r,qpcr_nr,gene)
  print(qpcr_plot)

  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'jpeg')
  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'tiff')
  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'svg')
  save_sizes(plot = qpcr_plot, filename = paste0(gene,"_qpcr",sep = ""), device = 'pdf')

}

