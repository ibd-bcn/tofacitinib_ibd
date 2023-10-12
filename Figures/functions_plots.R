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

feature_plot <- function(object, features, size = 3, cols = c('lightgray', 'red'), split.by = NULL, ncol = NULL){
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

