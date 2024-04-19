## Figure 3C -------------------------------------------------------------------

fig3c <- feature_plot(myeloids,
                      feature = c('SPP1'), #,'CLEC5A', 'MMP9', 'SPP1', 'IDO1'
                      split.by = c('response','pre_post'),
                      size = 0.1) +
  # patchwork::plot_layout(ncol = 2) &
  theme_figure() &
  theme(
    title = element_text(family = 'Helvetica', size = 12, colour = 'black', hjust = 0),
    plot.title = element_text(family = 'Helvetica', size = 12, colour = 'black'),
    axis.text = element_text(family = 'Helvetica', size = 9, colour = 'black'),
    strip.text.x = element_text(family = 'Helvetica', size = 10, colour = "black"),
    strip.text.y = element_text(family = 'Helvetica', size = 10, colour = "black"),
    legend.position = 'none',
    strip.background = element_blank()
  )

fig3c_nt <- feature_plot(myeloids,
                         feature = c('CLEC5A','IDO1', 'MMP9', 'SPP1'),
                         split.by = c('response','pre_post'),
                         size = 0.1) +
  patchwork::plot_layout(ncol = 2) &
  theme_figure() &
  theme(
    title = element_text(family = 'Helvetica', size = 10, colour = 'white', hjust = 0),
    plot.title = element_text(family = 'Helvetica', size = 10, colour = 'white'),
    axis.text = element_text(family = 'Helvetica', size = 6, colour = 'white'),
    strip.text.x = element_text(family = 'Helvetica', size = 8, colour = "white"),
    strip.text.y = element_text(family = 'Helvetica', size = 8, colour = "white"),
    legend.position = 'none',
    strip.background = element_blank()
  )

save_sizes(plot = fig3c, filename = 'Figure_3C', device = 'jpeg')
save_sizes(plot = fig3c, filename = 'Figure_3C', device = 'tiff')
save_sizes(plot = fig3c, filename = 'Figure_3C', device = 'svg')
save_sizes(plot = fig3c, filename = 'Figure_3C', device = 'pdf')

save_sizes(plot = fig3c_nt, filename = 'Figure_3C_no_text', device = 'jpeg')
save_sizes(plot = fig3c_nt, filename = 'Figure_3C_no_text', device = 'tiff')
save_sizes(plot = fig3c_nt, filename = 'Figure_3C_no_text', device = 'svg')
save_sizes(plot = fig3c_nt, filename = 'Figure_3C_no_text', device = 'pdf')

legend <- fig3c & theme(legend.position = 'right')
leg <- get_legend(legend)

plot_legend <- patchwork::wrap_plots(leg)
save_sizes(plot = plot_legend, filename = 'Figure_3C_legend', device = 'pdf')
save_sizes(plot = plot_legend, filename = 'Figure_3C_legend', device = 'jpeg')
save_sizes(plot = plot_legend, filename = 'Figure_3C_legend', device = 'tiff')
save_sizes(plot = plot_legend, filename = 'Figure_3C_legend', device = 'svg')
