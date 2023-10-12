library(ggplot2)
library(ggpubr)
#
# colors -----------------------------------------------------------------------
#
message('Loading colors')

colors_subset <- c('epi' = '#332f6f',
                   'stroma' =  '#e28c40',
                   'tcells' = '#eec35f',
                   'plasmas' = '#985f9d',
                   'myeloids' = '#75bccb',
                   'cycling' = '#7d9f56')
colors_response <- c('R' = '#70ADE6',
                     'NR' = '#FF8E47')

#
# functions --------------------------------------------------------------------
#
message('Loading functions')

save_sizes <- function(plot,
                       filename,
                       path = 'Figures/output',
                       device = 'svg',
                       max = 5){
  ggsave(filename = paste0(filename, '_larga.', device),
         path = path,
         plot = plot,
         dpi = 'print',
         device = device,
         height = max,
         width = max - 2,
         units = 'in')
  ggsave(filename = paste0(filename, '_ancha.', device),
         path = path,
         plot = plot,
         dpi = 'print',
         device = device,
         units = 'in',
         height = max - 2,
         width = max)
  ggsave(filename = paste0(filename, '_cuadrada.', device),
         path = path,
         plot = plot,
         dpi = 'print',
         device = device,
         units = 'in',
         height = max -2 ,
         width = max -2)
}


theme_figure <- function(){
  theme(
    plot.title = element_blank(),
    legend.position = 'none',
    axis.line = element_line(colour = 'black', size = 0.5),
    axis.title = element_blank(),
    axis.text = element_text(colour = 'black'),
    title =  element_text(family = 'Helvetica', size = 12),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    text = element_text(family = 'Helvetica', size = 10)
  )
}

theme_figure_wo_text <- function(){
  theme(
    plot.title = element_blank(),
    legend.position = 'none',
    axis.line = element_line(colour = 'black', size = 0.5),
    axis.title = element_blank(),
    axis.text = element_text(colour = 'white'),
    title = element_text(family = 'Helvetica', size = 12, color = 'white'),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    text = element_text(family = 'Helvetica', size = 10)
  )
}

theme_umap <- function(){
  theme(
    line = element_blank(),
    rect = element_blank(),
    plot.title = element_blank(),
    legend.position = 'none',
    strip.text.x = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    title = element_blank()
    )
}
