options(stringsAsFactors = FALSE)

library(ggplot2)
library(readxl)
source('Figures/functions_plots.R')


## Supplementary Figure 3-------------------------------------------------------

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

