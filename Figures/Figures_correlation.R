suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(nlcor))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
library(patchwork)

corr_tofa <- read_delim("Figures/extra_data/corr_tofa_sangre.csv",
                        delim = ";", escape_double = FALSE,
                        locale = locale(decimal_mark = ",",
                                        grouping_mark = "."),
                        na = c("ND","n.d"), trim_ws = TRUE)

SOCS1 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, SOCS1)) +
  geom_point(aes(color = Response)) +
  geom_smooth(formula = y ~ log(x), method = 'lm') +
  scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
  labs(x = 'Tofa conc (ng/ml)') +
  theme_classic()+
  theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))

SOCS3 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, SOCS3)) +
  geom_point(aes(color = Response)) +
  geom_smooth(formula = y ~ log(x), method = 'lm') +
  scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
  labs(x = 'Tofa conc (ng/ml)') +
  theme_classic()+
  theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))

IRF1 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, IRF1)) +
  geom_point(aes(color = Response)) +
  geom_smooth(formula = y ~ log(x), method = 'lm') +
  scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
  labs(x = 'Tofa conc (ng/ml)') +
  theme_classic()+
  theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))



patchwork::wrap_plots(SOCS1, SOCS3, IRF1, guides = 'collect') &
  theme(legend.position = 'top', axis.text = element_text(colour = 'black', size = 8))
