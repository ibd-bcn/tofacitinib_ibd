options(stringsAsFactors = FALSE,bitmapType = "cairo")
library(ggplot2)
library(readxl)
source('Figures/functions_plots.R')

# Correlations found in Analysis folder


## Data ------------------------------------------------------------------------
fig1<- read_excel("Figures/extra_data/## Supplementary Figure 2A: Dotplot -------------------------------------------Analiticas_sup_fig1.xlsx",
                   col_types = c("text", "text", "text",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric"))



## Supplementary Figure 1A: Boxplots -------------------------------------------

cells <- c("neutros","linfos","mono","leuco")

for(c in cells){

  anal_plot <- boxplot_analitica(fig1,c)
  print(anal_plot)

  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'jpeg')
  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'tiff')
  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'svg')
  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'pdf')

}

## Supplementary Figure 1B: Correlations ---------------------------------------
