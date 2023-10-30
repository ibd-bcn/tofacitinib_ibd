library(ggplot2)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
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

volcano_plot <- function(cluster, comp, filtered_genes) {

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
      scale_y_continuous(breaks = c(seq(0, 20, 5)), limits = c(0, 25))
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



  filtered_data <- response[response$gene %in% filtered_genes, ]

  fig <- fig+ geom_text_repel(data = filtered_data, aes(label = gene, color = "gray4"), size = 14/.pt
                               # fill = colors_volcano[filtered_data$sign],
                               # color = "white",
                               # segment.color = "black",
                               # fontface = 'bold',
                               # box.padding = unit(0.2, "lines"),
                               # point.padding = unit(0.5, "lines")
                               ) +
    theme(
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),

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


# Heatmap dataframe cleaning

heatmap_data <- function(input_hm, df) {
  #Delete IL4 data
  if(input_hm == "DMSO"){
    df <- TOFA[TOFA$Condition != "M-DMSO-IL4", ]
  }else{
    df <- TOFA[TOFA$Condition != "M-TOFA-IL4", ]
  }

}


#Heatmap dataframe preparation

dataframe_function <- function(dataframe){
  names <- dataframe$Donor
  dataframe <- dataframe[,-2]
  dataframe <- dataframe[,-1]
  dataframe <- as.data.frame(dataframe)
  rownames(dataframe) <- names

  return(dataframe)
}


# Heatmap function

heatmap_plot <- function(input_hm, df) {

  if(input_hm == "DMSO"){
    #Color
    col_fun = colorRamp2(c(-2,0,1,3), c("green", "black", "darkred","red"))
  }else{
    #Color
    col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), c("green", "darkgreen","black","darkred","red"))
  }

  ht_opt$TITLE_PADDING = unit(c(10, 10), "points")

  #Gene vector
  gene_vector <- c("TNF", "IDO1", "SOCS3", "IRF1", "OAS1", "MX1", "CXCL1", "CXCL5", "CXCL8", "CXCL10", "CCL5", "IL1B", "IL6", "IL23A", "IL10", "CD209", "MMP9", "CLEC5A", "SPP1", "ACOD1", "INHBA")

  #Group rows
  row_groups <- factor(
    c(rep("IFNg", 6), rep("Inflammatory cytokines", 6), rep("JAK dependent", 3), rep("M2", 2), rep("M1", 4))
  )

  #Colnames
  col_names <- c("LPS1", "LPS2", "LPS3", "LPS4", "TNF1", "TNF2", "TNF3", "TNF4", "IFNG1", "IFNG2", "IFNG3", "IFNG4")

  #Create col groups
  col_groups <- factor(substring(col_names, 1, 3), levels = c("LPS", "TNF", "IFNG"),
                       labels = c("LPS Group", "TNF Group", "IFNG Group"))

  # note how we set the width of this empty annotation
  ha = rowAnnotation(foo = anno_empty(border = FALSE,
                                      width =unit(0.5, "mm")))


  # Identify numeric columns in the elisa dataframe
  numeric_cols <- sapply(df, is.numeric)

  df[numeric_cols] <- log10(df[numeric_cols])

  df <- df[1:16,]

  if(input_hm == "DMSO"){
    #LPS
    LPS <- subset(df, Condition == "M-DMSO-LPS")
    LPS <- dataframe_function(LPS)
    LPS <- data.frame(LPS, row.names = NULL)
    rownames(LPS) <- c("LPS1","LPS2","LPS3","LPS4")

    #IFNg
    IFNG <- subset(df, Condition == "M-DMSO-IFNg")
    IFNG <- dataframe_function(IFNG)
    IFNG <- data.frame(IFNG, row.names = NULL)
    rownames(IFNG) <- c("IFNG1","IFNG2","IFNG3","IFNG4")


    #TNFa
    TNFa <- subset(df, Condition == "M-DMSO-TNFa")
    TNFa <- dataframe_function(TNFa)
    TNFa <- data.frame(TNFa, row.names = NULL)
    rownames(TNFa) <- c("TNF1","TNF2","TNF3","TNF4")
  }else{
    #LPS
    LPS <- subset(df, Condition == "M-TOFA-LPS")
    LPS <- dataframe_function(LPS)
    LPS <- data.frame(LPS, row.names = NULL)
    rownames(LPS) <- c("LPS1","LPS2","LPS3","LPS4")

    #IFNg
    IFNG <- subset(df, Condition == "M-TOFA-IFNg")
    IFNG <- dataframe_function(IFNG)
    IFNG <- data.frame(IFNG, row.names = NULL)
    rownames(IFNG) <- c("IFNG1","IFNG2","IFNG3","IFNG4")


    #TNFa
    TNFa <- subset(df, Condition == "M-TOFA-TNFa")
    TNFa <- dataframe_function(TNFa)
    TNFa <- data.frame(TNFa, row.names = NULL)
    rownames(TNFa) <- c("TNF1","TNF2","TNF3","TNF4")
  }

  t_df <- rbind(LPS,TNFa,IFNG)
  t_df <- t_df[gene_vector]
  t_df <- t(t_df)

  desired_order <- c("LPS1","LPS2","LPS3","LPS4","TNF1","TNF2","TNF3","TNF4","IFNG1","IFNG2","IFNG3","IFNG4")
  ## Obtain database

  heatmap <- Heatmap(t_df,
                     na_col = "white",
                     name = "Legend",
                     col = col_fun,
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     show_heatmap_legend = FALSE,
                     column_split = col_groups,
                     row_split = row_groups,
                     column_title = c("","",""),
                     show_column_names = FALSE,
                     right_annotation = ha,
                     show_row_names = FALSE,
                     row_title=NULL

  )

  draw(heatmap)

  #Row groups colors
  heatmap_colors <- c("#B4DC7F", "#FF6392", "#F9A03F", "#8BAEC7", "#516C7B")

  for(i in 1:5) {
    decorate_annotation("foo", slice = i, {
      grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = heatmap_colors[i], col =NA), just = "left")
    })
  }

  dev.off()
}

#Legend function

heatmap_lgd <- function(input_hm) {

  if(input_hm == "DMSO"){
    #Color
    col_fun = colorRamp2(c(-2,0,1,3), c("green", "black", "darkred","red"))
    #At
    at_vec <- c(-2,0,1,3)
    #Label
    lab_vec <- c("","","","")
  }else{
    #Color
    col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), c("green", "darkgreen","black","darkred","red"))
    #At
    at_vec <- c(-1,-0.5,0,0.5,1)
    #Label
    lab_vec <- c("","","","","")
  }

  lgd = Legend(col_fun = col_fun, direction = "horizontal", legend_width = unit(7, "cm"),at = at_vec,labels = lab_vec)
  draw(lgd)

  dev.off()
}
