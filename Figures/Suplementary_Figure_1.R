options(stringsAsFactors = FALSE)

library(ggplot2)
library(readxl)
source('Figures/functions_plots.R')


## Data ------------------------------------------------------------------------
fig1<- read_excel("Figures/extra_data/Analiticas_sup_fig1.xlsx",
                   col_types = c("text", "text", "text",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric"))



## Supplementary Figure 1: Boxplots --------------------------------------------

cells <- c("neutros","linfos","mono","leuco")

for(c in cells){

  anal_plot <- boxplot_analitica(fig1,c)
  print(anal_plot)

  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'jpeg')
  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'tiff')
  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'svg')
  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'pdf')

}

