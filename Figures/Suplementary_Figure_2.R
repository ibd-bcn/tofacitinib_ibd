options(stringsAsFactors = FALSE,bitmapType = "cairo")
library(ggplot2)
library(readxl)
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(nlcor))
suppressPackageStartupMessages(library(ggpubr))
library(patchwork)

source('Figures/functions_plots.R')

# Correlations found in Analysis folder


## Data ------------------------------------------------------------------------
fig2<- read_excel("Figures/extra_data/Analiticas_sup_fig1.xlsx",
                   col_types = c("text", "text", "text",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric"))



## Supplementary Figure 2A: Boxplots -------------------------------------------

cells <- c("neutros","linfos","mono","leuco")

for(c in cells){

  anal_plot <- boxplot_analitica(fig2,c)
  print(anal_plot)

  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'jpeg')
  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'tiff')
  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'svg')
  save_sizes(plot = anal_plot, filename = paste0(c,"_anal",sep = ""), device = 'pdf')

}

## Supplementary Figure 2B: Correlations ---------------------------------------


corr_tofa <- read_delim("/Figures/extra_data/corr_tofa_sangre.csv",
                        delim = ";", escape_double = FALSE,
                        locale = locale(decimal_mark = ",",
                                        grouping_mark = "."),
                        na = c("ND","n.d"), trim_ws = TRUE)

corr_tofa <- read_delim("/Figures/extra_data/corr_tofa.csv",
                        delim = ";", escape_double = FALSE,
                        locale = locale(decimal_mark = ",",
                                        grouping_mark = "."),
                        na = c("ND","n.d"), trim_ws = TRUE)


SOCS1 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, SOCS1)) +
  geom_point(aes(color = Response), size = 1) +
  geom_smooth(formula = y ~ log(x), method = 'lm') +
  scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
  labs(x = 'Tofa conc (ng/ml)') +
  theme_classic()+
  theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))

SOCS3 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, SOCS3)) +
  geom_point(aes(color = Response), size = 1) +
  geom_smooth(formula = y ~ log(x), method = 'lm') +
  scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
  labs(x = 'Tofa conc (ng/ml)') +
  theme_classic()+
  theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))

IRF1 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, IRF1)) +
  geom_point(aes(color = Response), size = 1) +
  geom_smooth(formula = y ~ log(x), method = 'lm') +
  scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
  labs(x = 'Tofa conc (ng/ml)') +
  theme_classic()+
  theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))



patchwork::wrap_plots(SOCS1, SOCS3, IRF1, guides = 'collect') &
  theme(legend.position = 'top', axis.text = element_text(colour = 'black', size = 8))

combined_plot <- patchwork::wrap_plots(SOCS1, SOCS3, IRF1, guides = 'collect') &
  theme(legend.position = 'top', axis.text = element_text(colour = 'black', size = 8))

ggsave("combined_plot_biopsia.png", plot = combined_plot, width = 8, height = 2.5, units = "in")
