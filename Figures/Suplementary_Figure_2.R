options(stringsAsFactors = FALSE,bitmapType = "cairo")
library(ggplot2)
library(readxl)
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(nlcor))
suppressPackageStartupMessages(library(ggpubr))
library(patchwork)

source('Figures/functions_plots.R')

# Correlations found in Analysis folder

## Functions -------------------------------------------------------------------
boxplot_analitica <- function(df, cp) {
  #Dictionary to grab column
  figures <- c(
    "neutros" = 4,
    "linfos" = 5,
    "mono" = 6,
    "leuco" = 3
  )
  #Grab columns
  fig_cell <- df[, c(1, 2, figures[[cp]])]
  #Convert week_2 column to factor
  fig_cell$Time <-
    factor(fig_cell$Time,
           levels = c("w0", "w2", "≥w8"))
  #Change column name of the celltype
  colnames(fig_cell)[3] <- cp
  #Convert column response to factor
  fig_cell$Response <-
    factor(fig_cell$Response, levels = c("R", "NR"))
  fig_cell$black <- "black"

  #Grab Max value of datapoints
  max <- max(fig_cell[, 3], na.rm = TRUE)

  #Jitter width selection
  if (cp == "mono") {
    jw <- 0.5
  } else{
    jw <- 0.2
  }

  #Plot
  p <- ggplot(na.omit(fig_cell), aes(x = Time, y = .data[[cp]])) +
    geom_boxplot(
      aes(fill = Response),
      color = "black",
      alpha = 0.9,
      size = 0.5,
      outlier.shape = NA
    ) +
    geom_jitter(
      aes(color = black),
      position = position_jitterdodge(dodge.width = 0.4, jitter.width = jw),
      alpha = 1,
      size = 1
    ) +
    facet_grid(~ Response, scales = "free", space = "free") +
    theme_bw(base_rect_size = 1, base_size = 20) +
    theme(
      legend.position = 'none',
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      strip.text.y = element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(0.1, "cm")
    ) + scale_fill_manual(values = c("#70ADE6", "#FF8E47")) + scale_color_manual(values = c("#000000")) + ylim(c(0, round(max)))

}


## Data ------------------------------------------------------------------------
fig2 <- read_excel("Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
                   sheet = "Sup figure 2A", col_types = c("text",
                                                          "text", "numeric", "numeric", "numeric",
                                                          "numeric"))
colnames(fig2) <- gsub(" \\(week\\)", "", colnames(fig2))


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
