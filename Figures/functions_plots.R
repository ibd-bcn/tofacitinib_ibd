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

#Heatmaps of progeny

heatmap_progeny <- function(stimuli, anot = "reduced_anot"){

  #Open dfs
  todas_intermediate_w0_NR <- read_csv(paste0("~/tofacitinib_ibd/Analysis/PROGENy/data/todas_",anot,"_w0_NR.csv", sep = ""))
  cn_w0_NR <- todas_intermediate_w0_NR$row_names
  jak_w0_NR <- todas_intermediate_w0_NR[[stimuli]]

  todas_intermediate_w0_R <- read_csv(paste0("~/tofacitinib_ibd/Analysis/PROGENy/data/todas_",anot,"_w0_R.csv", sep = ""))
  cn_w0_R <- todas_intermediate_w0_R$row_names
  jak_w0_R <- todas_intermediate_w0_R[[stimuli]]

  todas_intermediate_w8_R <- read_csv(paste0("~/tofacitinib_ibd/Analysis/PROGENy/data/todas_",anot,"_w8_R.csv", sep = ""))
  cn_w8_R <- todas_intermediate_w8_R$row_names
  jak_w8_R <- todas_intermediate_w8_R[[stimuli]]

  todas_intermediate_w8_NR <- read_csv(paste0("~/tofacitinib_ibd/Analysis/PROGENy/data/todas_",anot,"_w8_NR.csv", sep = ""))
  cn_w8_NR <- todas_intermediate_w8_NR$row_names
  jak_w8_NR <- todas_intermediate_w8_NR[[stimuli]]

  #Modify df
  df_w0_NR <- data.frame(cn = cn_w0_NR, jak_w0_NR = jak_w0_NR)
  df_w0_R <- data.frame(cn = cn_w0_R, jak_w0_R = jak_w0_R)
  df_w8_R <- data.frame(cn = cn_w8_R, jak_w8_R = jak_w8_R)
  df_w8_NR <- data.frame(cn = cn_w8_NR, jak_w8_NR = jak_w8_NR)

  combined_df <- df_w0_NR %>%
    full_join(df_w0_R, by = "cn") %>%
    full_join(df_w8_R, by = "cn") %>%
    full_join(df_w8_NR, by = "cn")

  colnames(combined_df) <- c("cn", "w0_NR", "w0_R", "w8_R", "w8_NR")
  cells <- combined_df$cn
  rownames(combined_df) <- combined_df$cn
  combined_df <- combined_df[,-c(1)]
  data_t <- t(combined_df)
  colnames(data_t) <- cells
  data_t <- data_t[c("w0_NR","w0_R","w8_NR","w8_R"),]
  data <- data_t[,c("T cell","Plasma cell","B cell","Macrophages","Neutrophils","Mast cells","Inflammatory monocytes","DCs","Eosinophils","Fibroblasts","Endothelium","Glia","Epithelium")]
  paletteLength = 100
  myColor <- colorRamp2(range(na.omit(data)), hcl_palette = "Reds", reverse = TRUE)
  data <- data[c("w8_NR", "w8_R", "w0_NR", "w0_R"),]
  colnames(data) <- gsub(pattern = "Inflammatory monocytes",replacement = "Inf mono",x = colnames(data))

  #Heatmap
  progeny_hmap <- Heatmap(data,
                          name = "PROGENy (500)",
                          col = myColor,
                          show_row_names = TRUE,
                          show_column_names = TRUE,
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          rect_gp = gpar(col = NA),
                          row_title = NULL,
                          column_title = NULL,
                          row_names_gp = gpar(fontsize = 12),
                          column_names_gp = gpar(fontsize = 18)
  )

  # Draw the heatmap
  draw(progeny_hmap, heatmap_legend_side = "right", annotation_legend_side = "left")
}
