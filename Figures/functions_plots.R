library(ggplot2)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
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


colors_volcano <- c("UPP" = "#911704",
                    "UP" = "#C38683",
                    "0" = "#D1D1D1",
                    "DW" = "#98AE99",
                    'DWW' = "#376D38")
#
# functions --------------------------------------------------------------------
#
message('Loading functions')

save_sizes <- function(plot,
                       filename,
                       path = 'Figures/output',
                       device = 'svg',
                       max = 5){
  #
  # save plot in different sizes to try the best position
  #
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
  #
  # theme to use in plots by defect.
  #
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
  #
  # theme of figures wo text
  #
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
  #
  # Theme to use in umap plots
  #
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

feature_plot <- function(object, features, size = 3, cols = c('lightgray', 'red'), split.by = NULL, ncol = NULL){
  #
  # feature plot that allows 2 splits.
  #
  require(ggplot2)
  require(Seurat)
  require(patchwork)
  require(formulaic)
  plots <- NULL
  if (sum(features %in% rownames(object)) != length(features)){
    print(paste('Feature(s)', features[!(features %in% rownames(object))], 'not present in the object'))
  }
  if(sum(features %in% rownames(object)) == length(features)){
    data <- as.data.frame(object@reductions$umap@cell.embeddings)
    feat <- 2+length(features)
    for(j in 1:length(features)){
      data <- cbind(data, object@assays[['RNA']]@data[features[j],])}
    colnames(data)[3:feat] <- features
    if(is.null(split.by)){
      for(i in 1:length(features)){
        data <- data[order(data[,features[i]]),]
        if(sum(data[,features[i]]) != 0){
          plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color = add.backtick(features[i], include.backtick = 'all'))) +
            geom_point(size = size) +
            scale_color_gradient(low=cols[1], high=cols[2]) +
            theme_classic() +
            labs(title = features[i], x = element_blank(), y = element_blank()) +
            theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))
          plots[[i]] <- plot
        }else{
          plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color =add.backtick(features[i], include.backtick = 'all'))) +
            geom_point(size = size)   +
            scale_color_gradient(low=cols[1], high=cols[1]) +
            theme_classic() +
            labs(title = features[i], x = element_blank(), y = element_blank()) +
            theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))
          plots[[i]] <- plot
        }
      }
      if (is.null(x = ncol)) {
        ncol <- 2
        if(length(x = features) == 1) {
          ncol <- 1
        }
        if(length(x = features) > 6) {
          ncol <- 3
        }
        if(length(x = features) > 9) {
          ncol <- 4
        }
      }

      plots <- wrap_plots(plots, ncol = ncol)
      return(plots)
    }


    if(!is.null(split.by)){
      splitings <- split.by
      if(length(splitings) == 1){
        plots <- NULL
        metadata <- FetchData(object, split.by)
        data <- cbind(data, metadata)
        colnames(data)[ncol(data)]<- 'split.by'
        ncol <- 1
        for(i in 1:length(features)){
          if(i == 1){
            data <- data[order(data[,features[i]]),]
            if(sum(data[,features[i]]) != 0){
              plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color = add.backtick(features[i], include.backtick = 'all'))) +
                geom_point(size = size) +
                scale_color_gradient(low=cols[1], high=cols[2]) +
                theme_classic() + facet_grid(. ~ split.by) +
                labs(title = features[i], x = element_blank(), y = '') +
                theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))+
                theme(strip.background = element_rect(colour="white", fill="white",
                                                      size=1.5, linetype="solid"))+
                theme(strip.text.x = element_text(size=12, color="black",
                                                  face="bold"))+
                theme(panel.spacing = unit(1, "lines"))
              plots[[i]] <- plot
            }else{
              plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color = add.backtick(features[i], include.backtick = 'all'))) +
                geom_point(size = size)  +
                scale_color_gradient(low=cols[1], high=cols[1]) +
                theme_classic() + facet_grid(. ~ split.by) +
                labs(title = features[i], x = element_blank(), y = '') +
                theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))+
                theme(strip.background = element_rect(colour="white", fill="white",
                                                      size=1.5, linetype="solid"))+
                theme(strip.text.x = element_text(size=12, color="black",
                                                  face="bold"))+
                theme(panel.spacing = unit(1, "lines"))
              plots[[i]] <- plot
            }
          }
          if(i != 1){
            data <- data[order(data[,features[i]]),]
            if(sum(data[,features[i]]) != 0){
              plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color = add.backtick(features[i], include.backtick = 'all'))) +
                geom_point(size = size) +
                scale_color_gradient(low=cols[1], high=cols[2]) +
                theme_classic() + facet_grid(. ~ split.by) +
                labs(title = features[i], x = element_blank(), y ='') +
                theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))+
                theme(strip.background = element_rect(colour="white", fill="white",
                                                      size=1.5, linetype="solid"))+
                theme(strip.text.x = element_text(size=1, color="white",
                                                  face="bold"))+
                theme(panel.spacing = unit(1, "lines"))
              plots[[i]] <- plot
            }else{
              plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color = add.backtick(features[i], include.backtick = 'all'))) +
                geom_point(size = size)  +
                scale_color_gradient(low=cols[1], high=cols[1]) +
                theme_classic() + facet_grid(. ~ split.by) +
                labs(title = features[i], x = element_blank(), y ='') +
                theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))+
                theme(strip.background = element_rect(colour="white", fill="white",
                                                      size=1.5, linetype="solid"))+
                theme(strip.text.x = element_text(size=1, color="white",
                                                  face="bold"))+
                theme(panel.spacing = unit(1, "lines"))
              plots[[i]] <- plot
            }
          }
        }
        plots <-  wrap_plots(plots, ncol = ncol)
        return(plots)
      }
      if(length(splitings) == 2){
        plots <- NULL
        metadata <- FetchData(object, split.by)
        data <- cbind(data, metadata)
        colnames(data)[ncol(data) -1]<- 'split.by.1'
        colnames(data)[ncol(data)]<- 'split.by.2'
        ncol <- 1
        for(i in 1:length(features)){
          if(i == 1){
            data <- data[order(data[,features[i]]),]
            if(sum(data[,features[i]]) != 0){
              plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color = add.backtick(features[i], include.backtick = 'all'))) +
                geom_point(size = size) +
                scale_color_gradient(low=cols[1], high=cols[2]) +
                theme_classic() + facet_grid(cols = vars(split.by.1), rows = vars(split.by.2)) +
                labs(title = features[i], x = element_blank(), y = '') +
                theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))+
                theme(strip.background = element_rect(colour="white", fill="white",
                                                      size=1.5, linetype="solid"))+
                theme(strip.text.x = element_text(size=12, color="black",
                                                  face="bold"))+
                theme(panel.spacing = unit(1, "lines"))
              plots[[i]] <- plot
            }else{
              plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color = add.backtick(features[i], include.backtick = 'all'))) +
                geom_point(size = size)  +
                scale_color_gradient(low=cols[1], high=cols[1]) +
                theme_classic() + facet_grid(cols = vars(split.by.1), rows = vars(split.by.2)) +
                labs(title = features[i], x = element_blank(), y = '') +
                theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))+
                theme(strip.background = element_rect(colour="white", fill="white",
                                                      size=1.5, linetype="solid"))+
                theme(strip.text.x = element_text(size=12, color="black",
                                                  face="bold"))+
                theme(panel.spacing = unit(1, "lines"))
              plots[[i]] <- plot
            }
          }
          if(i != 1){
            data <- data[order(data[,features[i]]),]
            if(sum(data[,features[i]]) != 0){
              plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color = add.backtick(features[i], include.backtick = 'all'))) +
                geom_point(size = size) +
                scale_color_gradient(low=cols[1], high=cols[2]) +
                theme_classic() + facet_grid(cols = vars(split.by.1), rows = vars(split.by.2)) +
                labs(title = features[i], x = element_blank(), y ='') +
                theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))+
                theme(strip.background = element_rect(colour="white", fill="white",
                                                      size=1.5, linetype="solid"))+
                theme(strip.text.x = element_text(size=1, color="white",
                                                  face="bold"))+
                theme(panel.spacing = unit(1, "lines"))
              plots[[i]] <- plot
            }else{
              plot <- ggplot(data, aes_string(x='UMAP_1', y='UMAP_2', color = add.backtick(features[i], include.backtick = 'all'))) +
                geom_point(size = size)  +
                scale_color_gradient(low=cols[1], high=cols[1]) +
                theme_classic() + facet_grid(cols = vars(split.by.1), rows = vars(split.by.2)) +
                labs(title = features[i], x = element_blank(), y ='') +
                theme(title = element_text(size = 15), axis.text = element_text(size = 12, colour = 'black'))+
                theme(strip.background = element_rect(colour="white", fill="white",
                                                      size=1.5, linetype="solid"))+
                theme(strip.text.x = element_text(size=1, color="white",
                                                  face="bold"))+
                theme(panel.spacing = unit(1, "lines"))
              plots[[i]] <- plot
            }
          }
        }
        plots <-  wrap_plots(plots, ncol = ncol)
        return(plots)
      }
    }
    rm(plot, plots, data, metadata, i, j, ncol)
  }

}


# Volcano_plot function
volcano_plot <- function(cluster, comp, filtered_genes_down, filtered_genes_up) {

  labels <- c(
    "UPP" = "UPP",
    "UP" = "UP",
    "0" = "0",
    "DW" = "DW",
    "DWW" = "DWW"
  )

  de_data2 <- de_data[de_data$cluster == cluster &
                        de_data$annotation == 'annotation_intermediate' &
                        de_data$comp == comp &
                        de_data$sign %in% names(labels), c("p_val", "avg_log2FC", "sign", "comp", "gene")]

  response <- subset(de_data2, comp == comp, select = c("avg_log2FC", "p_val", "gene", "sign"))

  if (cluster != "IDA macrophages") {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = c(seq(0, 30, 5)), limits = c(0, 30))
  } else {
    fig <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_volcano, labels = labels) +
      theme_classic() +
      theme(text = element_text(family = "Helvetica")) +
      guides(color = guide_legend(override.aes = list(shape = 1))) +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3))
  }



  filtered_data <- response[response$gene %in% filtered_genes_down, ]
  filtered_data2 <- response[response$gene %in% filtered_genes_up, ]

  fig <- fig+ geom_label_repel(data = filtered_data, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                               segment.color = "black",
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines"),
                               position = position_nudge_repel(x = 0, y = 2)) +
    geom_label_repel(data = filtered_data2, aes(label = gene, group = gene, col = sign), size = 8/.pt,
                     segment.color = "black",
                     # fill = colors_volcano[filtered_data$sign],
                     # color = "white",
                     # segment.color = "black",
                     fontface = 'bold',
                     # box.padding = unit(0.2, "lines"),
                     # point.padding = unit(0.5, "lines"),
                     position = position_nudge_repel(x = 0, y = 1)) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")

    )

  return(fig)
}

# Boxplot function

boxplot_plot <- function(qpcr_r,qpcr_nr,gene) {

  #Compute max Y label
  value1 <- max(na.omit(qpcr_r[[gene]]))
  value2 <- max(na.omit(qpcr_nr[[gene]]))
  max_value <- max(c(value1,value2))
  label_y_max <- max_value + (max_value*0.2)
  label_y_stat <- max_value + (max_value*0.11)

  #Responder plot
  p <- ggplot(qpcr_r, aes(x = Time, y = .data[[gene]])) +
    geom_boxplot(aes(fill = facet_group), color = "black", alpha = 0.9,size = 0.5,outlier.shape = NA) +
    geom_jitter(aes(color = facet_group),
                position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
                alpha = 1, size = 1) +
    theme_bw(base_rect_size = 1,base_size = 20) +  facet_grid(. ~ "") +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          plot.title = element_blank(),
          axis.ticks=element_line(size=0.5),
          axis.ticks.length = unit(0.1, "cm")
    ) +
    ylim(c(0,round(label_y_max)))


  p <- p + scale_fill_manual(values = c("#70ADE6")) + scale_color_manual(values = c("#000000"))


  #Non-responder plot
  comparison_result <- pairwise.wilcox.test(qpcr_nr[[gene]], qpcr_nr$Time, p.adjust.method = "none")
  p_value <- comparison_result$p.value[1]

  q <- ggplot(qpcr_nr, aes(x = Time, y = .data[[gene]])) +
    geom_boxplot(aes(fill = facet_group), color = "black", alpha = 0.9,size = 0.5, outlier.shape = NA) +
    geom_jitter(aes(color = facet_group),
                position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
                alpha = 1,size = 1) +
    theme_bw( base_rect_size = 1,base_size = 20)  + facet_grid(. ~ "") +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          plot.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks=element_line(size=0.5),
          axis.ticks.length = unit(0.1, "cm"))  +
    ylim(c(0,round(label_y_max)))


  q <- q + scale_fill_manual(values = c("#FF8E47")) + scale_color_manual(values = c("#000000"))

  #Join two plots
  combined_plot <- p + q

  return(combined_plot)

}

# Heatmap function


heatmap_plot <- function(input_hm, df) {
  df <- df %>%
    mutate(across(where(is.numeric), log2))


  colnames(df)[2] <- "Sample"
  colnames(df)[3] <- "Condition"

  if (input_hm == "DMSO") {
    #Color
    # col_fun = colorRamp2(c(-5, 0, 5, 15), c("green", "black", "darkred", "red"))
    col_fun = colorRamp2(c(-5, 0, 15), c("green", "black", "red"))
    #Get rid of IL4
    df <- df[df$Condition != "M-DMSO-IL4",]
  } else{
    #Color
    # col_fun = colorRamp2(c(-4, 0, 1, 3.5), c("green", "black", "darkred", "red"))

    col_fun = colorRamp2(c(-4,0,3.5), c("green","black","red"))
    #Get rid of IL4
    df <- df[df$Condition != "M-TOFA-IL4",]
  }

  ht_opt$TITLE_PADDING = unit(c(10, 10), "points")

  #Gene vector
  gene_vector <-
    c(
      "TNF",
      "IDO1",
      "SOCS3",
      "IRF1",
      "OAS1",
      "MX1",
      "CXCL1",
      "CXCL8",
      "CXCL10",
      "CCL5",
      "IL1B",
      "IL6",
      "IL23A",
      "IL10",
      "CD209",
      "MMP9",
      "CLEC5A",
      "SPP1",
      "ACOD1",
      "INHBA"
    )

  #Group rows
  row_groups <- factor(c(
    rep("IFNg", 6),
    rep("Inflammatory cytokines", 5),
    rep("JAK dependent", 3),
    rep("M2", 2),
    rep("M1", 4)
  ))

  #Colnames
  col_names <-
    c(
      "LPS1",
      "LPS2",
      "LPS3",
      "LPS4",
      "LPS5",
      "TNF1",
      "TNF2",
      "TNF3",
      "TNF4",
      "TNF5" ,
      "IFNG1",
      "IFNG2",
      "IFNG3",
      "IFNG4",
      "IFNG5"
    )

  #Create col groups
  col_groups <-
    factor(
      substring(col_names, 1, 3),
      levels = c("LPS", "TNF", "IFN"),
      labels = c("LPS Group", "TNF Group", "IFNG Group")
    )

  # note how we set the width of this empty annotation
  ha = rowAnnotation(foo = anno_empty(border = FALSE,
                                      width = unit(0.5, "mm")))

  #Column title
  column_title <-
    gpar(
      fill = "lightgrey",
      col = "black",
      border = "white",
      lwd = 3,
      fontsize = 25
    )

  #Row title
  row_title <- gpar(fontsize = 14)

  # Identify numeric columns in the elisa dataframe
  numeric_cols <- sapply(df, is.numeric)

  if (input_hm == "DMSO") {
    #LPS
    LPS <- subset(df, Condition == "M-DMSO-LPS")
    LPS <- data.frame(LPS, row.names = NULL)
    rownames(LPS) <- c("LPS1", "LPS2", "LPS3", "LPS4", "LPS5")
    LPS <- LPS[, 4:length(colnames(LPS))]
    #IFNg
    TNFa <- subset(df, Condition == "M-DMSO-TNF?")
    TNFa <- data.frame(TNFa, row.names = NULL)
    rownames(TNFa) <- c("TNF1", "TNF2", "TNF3", "TNF4", "TNF5")
    TNFa <- TNFa[, 4:length(colnames(TNFa))]
    #TNFa
    IFNG <- subset(df, Condition == "M-DMSO-IFN?")
    IFNG <- data.frame(IFNG, row.names = NULL)
    rownames(IFNG) <- c("IFNG1", "IFNG2", "IFNG3", "IFNG4", "IFNG5")
    IFNG <- IFNG[, 4:length(colnames(IFNG))]
  } else{
    #LPS
    LPS <- subset(df, Condition == "M-TOFA-LPS")
    LPS <- data.frame(LPS, row.names = NULL)
    rownames(LPS) <- c("LPS1", "LPS2", "LPS3", "LPS4", "LPS5")
    LPS <- LPS[, 4:length(colnames(LPS))]
    #IFNg
    TNFa <- subset(df, Condition == "M-TOFA-TNF?")
    TNFa <- data.frame(TNFa, row.names = NULL)
    rownames(TNFa) <- c("TNF1", "TNF2", "TNF3", "TNF4", "TNF5")
    TNFa <- TNFa[, 4:length(colnames(TNFa))]
    #TNFa
    IFNG <- subset(df, Condition == "M-TOFA-IFN?")
    IFNG <- data.frame(IFNG, row.names = NULL)
    rownames(IFNG) <- c("IFNG1", "IFNG2", "IFNG3", "IFNG4", "IFNG5")
    IFNG <- IFNG[, 4:length(colnames(IFNG))]
  }

  t_df <- rbind(LPS, TNFa, IFNG)
  t_df <- t_df[gene_vector]
  t_df <- t(t_df)

  desired_order <-
    c(
      "LPS1",
      "LPS2",
      "LPS3",
      "LPS4",
      "LPS5",
      "TNF1",
      "TNF2",
      "TNF3",
      "TNF5",
      "TNF4",
      "IFNG1",
      "IFNG2",
      "IFNG3",
      "IFNG4",
      "IFNG5"
    )
  ## Obtain database

  heatmap <- Heatmap(
    t_df,
    na_col = "white",
    name = "Legend",
    col = col_fun,
    column_title_gp = column_title,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    show_heatmap_legend = FALSE,
    column_split = col_groups,
    row_split = row_groups,
    row_title = NULL,
    column_title = c("LPS", "TNF", "IFN"),
    show_column_names = FALSE,
    right_annotation = ha,
    row_names_gp = gpar(fontsize = 20),

  )

  draw(heatmap)

  #Row groups colors
  heatmap_colors <-
    c("#B4DC7F", "#FF6392", "#F9A03F", "#8BAEC7", "#516C7B")

  for (i in 1:5) {
    decorate_annotation("foo", slice = i, {
      grid.rect(
        x = 0,
        width = unit(2, "mm"),
        gp = gpar(fill = heatmap_colors[i], col = NA),
        just = "left"
      )
    })
  }

  dev.off()
}


#Legend function

heatmap_lgd <- function(input_hm) {

  if(input_hm == "DMSO"){
    #Color
    col_fun = colorRamp2(c(-5,0,15), c("green", "black","red"))
    #At
    at_vec <- c(-5,0,15)
    #Label
    lab_vec <- c("-5","0","15")
  }else{
    #Color
    col_fun = colorRamp2(c(-4,0,3.5), c("green","black","red"))
    #At
    at_vec <- c(-4,0,3.5)
    #Label
    lab_vec <- c("-4","0","3.5")
  }

  lgd = Legend(col_fun = col_fun, direction = "horizontal", legend_width = unit(7, "cm"),at = at_vec,labels = lab_vec)
  draw(lgd)

  dev.off()
}

#Boxplot analitics function

boxplot_analitica <- function(df, cp) {
  #Dictionary to grab column
  figures <- c(
    "neutros" = 4,
    "linfos" = 7,
    "mono" = 8,
    "leuco" = 9
  )
  #Grab columns
  fig_cell <- df[, c(1, 2, 3, figures[[cp]])]
  #Change colum name
  colnames(fig_cell)[3] <- "WEEK_2"
  #Convert week_2 column to factor
  fig_cell$`WEEK_2` <-
    factor(fig_cell$`WEEK_2`,
           levels = c("Pre-tx", "Week 2", "Week 8 and later"))
  #Change column name of the celltype
  colnames(fig_cell)[4] <- cp
  #Convert column response to factor
  fig_cell$RESPONSE <-
    factor(fig_cell$RESPONSE, levels = c("R", "NR"))
  fig_cell$black <- "black"

  #Grab Max value of datapoints
  max <- max(fig_cell[, 4], na.rm = TRUE)

  #Jitter width selection
  if (cp == "mono") {
    jw <- 0.5
  } else{
    jw <- 0.2
  }

  #Plot
  p <- ggplot(na.omit(fig_cell), aes(x = WEEK_2, y = .data[[cp]])) +
    geom_boxplot(
      aes(fill = RESPONSE),
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
    facet_grid( ~ RESPONSE, scales = "free", space = "free") +
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

