
seurat_to_pca <- function(object){
  print('NORMALIZATION')
  object <- NormalizeData(object)
  print('VARIABLE FEATURES')
  object <- FindVariableFeatures(object)
  print('SCALE DATA')
  object <- ScaleData(object)
  print('RUN PCA')
  object <- RunPCA(object, npcs = 100, ndims.print = 1, nfeatures.print = 5)
  return(object)
}

funfun <- function(x,y,z){ # x = pval, y = log2foldchange & z = fdr
  column <- rep(0, length(x))
  UPS <- intersect(which(x < 0.05), which(y > log2(1.2)))
  DWS <- intersect(which(x < 0.05), which(y < -log2(1.2)))
  UPP <- intersect(which(y > log2(1.2)), which(z< 0.05))
  DWW <- intersect(which(y < -log2(1.2)), which(z< 0.05))
  column[UPS] <- "UP"
  column[DWS] <- "DW"
  column[UPP] <- "UPP"
  column[DWW] <- "DWW"
  return(column)
}

resolutions <- function(object,
                        resolutions = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5),
                        workingdir,
                        title = object$orig.ident[1]){
  require(clustree)
  oldworkingdir <- getwd()
  setwd(workingdir)
  for(i in 1:length(resolutions)){
    object <- FindClusters(object, resolution = resolutions[i])

    markers <- FindAllMarkers(object = object,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              thresh.use = 0.25)

    write.table(x = markers, file = paste0(title,'_markers_resolution_',resolutions[i],'.csv'),
                row.names = F, sep = ';', dec= ',', col.names = T)

    png(filename = paste0(title,'_resolution_',resolutions[i],'.png'), res = 300, width = 16, height = 12, units = 'cm')
    p <- DimPlot(object, label = T) +
      theme(legend.position = 'right') + labs(title = title)
    print(p)
    dev.off()

    png(filename = paste0(title,'_Violin_mtper_resolution_',resolutions[i],'.png'),
        res = 300, width = 25, height = 20, units = 'cm')
    p <- VlnPlot(object, features = 'percent.mt') +
      theme(legend.position = 'none') + labs(title = title)
    print(p)
    dev.off()
  }

  png(filename = paste0(title,'_clustree_resolutions.png'), res = 300, width = 25, height = 20, units = 'cm')
  p <- clustree(object, prefix = 'RNA_snn_res.') + labs(subtitle = title )
  print(p)
  dev.off()

  setwd(oldworkingdir)

  rm(markers, oldworkingdir, p)
  return(object)
}


resolutions_integrated <- function(object,
                                   resolutions = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5),
                                   workingdir,
                                   title = object$orig.ident[1]){
  require(clustree)
  oldworkingdir <- getwd()
  setwd(workingdir)
  for(i in 1:length(resolutions)){
    object <- FindClusters(object, resolution = resolutions[i])

    markers <- FindAllMarkers(object = object,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              thresh.use = 0.25)

    write.table(x = markers, file = paste0(title,'_markers_resolution_',resolutions[i],'.csv'),
                row.names = F, sep = ';', dec= ',', col.names = T)

    png(filename = paste0(title,'_resolution_',resolutions[i],'.png'), res = 300, width = 16, height = 12, units = 'cm')
    p <- DimPlot(object, label = T) +
      theme(legend.position = 'right') + labs(title = title)
    print(p)
    dev.off()

    # png(filename = paste0(title,'_Violin_mtper_resolution_',resolutions[i],'.png'),res = 300, width = 25, height = 20, units = 'cm')
    # p <- VlnPlot(object, features = 'percent.mt') +
    #   theme(legend.position = 'none') + labs(title = title)
    # print(p)
    # dev.off()
  }

  png(filename = paste0(title,'_clustree_resolutions.png'), res = 300, width = 25, height = 20, units = 'cm')
  p <- clustree(object, prefix = 'integrated_snn_res.') + labs(subtitle = title )
  print(p)
  dev.off()

  setwd(oldworkingdir)

  rm(markers, oldworkingdir, p)
  return(object)
}

resolutions_2 <- function(object,
                          resolutions = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5),
                          workingdir,
                          title = object$orig.ident[1]){
  require(clustree)
  oldworkingdir <- getwd()
  setwd(workingdir)
  for(i in 1:length(resolutions)){
    object <- FindClusters(object, resolution = resolutions[i])

    markers <- FindAllMarkers(object = object,
                              only.pos = TRUE)
    # min.pct = 0.25,
    # thresh.use = 0.25)

    write.table(x = markers, file = paste0(title,'_markers_resolution_',resolutions[i],'.csv'),
                row.names = F, sep = ';', dec= ',', col.names = T)

    png(filename = paste0(title,'_resolution_',resolutions[i],'.png'), res = 300, width = 16, height = 12, units = 'cm')
    p <- DimPlot(object, label = T) +
      theme(legend.position = 'right') + labs(title = title)
    print(p)
    dev.off()

    png(filename = paste0(title,'_Violin_mtper_resolution_',resolutions[i],'.png'),res = 300, width = 25, height = 20, units = 'cm')
    p <- VlnPlot(object, features = 'percent.mt') +
      theme(legend.position = 'none') + labs(title = title)
    print(p)
    dev.off()
  }

  png(filename = paste0(title,'_clustree_resolutions.png'), res = 300, width = 25, height = 20, units = 'cm')
  p <- clustree(object, prefix = 'RNA_snn_res.') + labs(subtitle = title )
  print(p)
  dev.off()

  setwd(oldworkingdir)

  rm(markers, oldworkingdir, p)
  return(object)
}


select_pcs <- function(object, sdevi){
  require(tidyr)
  require(dplyr)
  if(!is.null(object@reductions$pca)){
    df <- data.frame(variance = object@reductions$pca@global,
                     sdev =object@reductions$pca@stdev)
    df <- df %>%
      mutate(diferencia = sdev - lag(sdev, 1))

    df$diferencia[nrow(df)] <- NA

    PCS_1 <- max(which(df$diferencia < -0.01 & df$sdev > sdevi))
    rm(df)
    return(PCS_1)
  }else{
    simpleError(message = 'No PCA reduction in this object')
  }
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


barplot_figure <- function(seurat, genes, group, cols = NULL){
  require(Seurat)
  require(gridExtra)
  require(ggpubr)
  library(formulaic)
  data <- FetchData(seurat, vars = c(genes, group))
  colnames(data)[ncol(data)] <- 'class'
  data <- data[order(data$class),]
  data$cells <- rownames(data)
  data$cells <- factor(data[,ncol(data)], levels =data[,ncol(data)])
  head(data)
  plot_list <- list()
  if(is.null(cols)){
    for(i in 1:length(genes)){
      p = ggplot(data, aes_string(x='cells',
                                  y=add.backtick(genes[i],
                                                 include.backtick = 'all'),
                                  fill = 'class', color = 'class')) +
        geom_bar(stat = "identity") +
        labs(x = '')+
        theme_classic() +
        scale_y_continuous(
          labels = scales::number_format(accuracy = 0.1))+
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'none')
      plot_list[[i]] = p
    }
    j <- ggplot(data, aes_string(x='cells',
                                 y=add.backtick(genes[i],
                                                include.backtick = 'all'),
                                 fill = 'class', color = 'class')) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(legend.position = 'bottom')
    p = get_legend(j)

    plot_list[[i+1]] = p
  }else{
    for(i in 1:length(genes)){
      p = ggplot(data, aes_string(x='cells',
                                  y=add.backtick(genes[i],
                                                 include.backtick = 'all'),
                                  fill = 'class', color = 'class')) +
        geom_bar(stat = "identity") +
        labs(x = '')+
        theme_classic() +
        scale_color_manual(values = cols)+
        scale_y_continuous(
          labels = scales::number_format(accuracy = 0.1))+
        scale_fill_manual(values = cols)+
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'none')
      plot_list[[i]] = p
    }
    j <- ggplot(data, aes_string(x='cells',
                                 y=add.backtick(genes[i],
                                                include.backtick = 'all'),
                                 fill = 'class', color = 'class')) +
      geom_bar(stat = "identity") +
      theme_classic() +
      scale_color_manual(values = cols)+
      scale_fill_manual(values = cols)+
      theme(legend.position = 'bottom')
    p = get_legend(j)

    plot_list[[i+1]] = p
  }

  return(do.call("ggarrange", c(plot_list, ncol=1)))
}

Col2Hex <- function(...) {
  colors <- as.character(x = c(...))
  alpha <- rep.int(x = 255, times = length(x = colors))
  if (sum(sapply(X = colors, FUN = grepl, pattern = '^#')) != 0) {
    hex <- colors[which(x = grepl(pattern = '^#', x = colors))]
    hex.length <- sapply(X = hex, FUN = nchar)
    if (9 %in% hex.length) {
      hex.alpha <- hex[which(x = hex.length == 9)]
      hex.vals <- sapply(X = hex.alpha, FUN = substr, start = 8, stop = 9)
      dec.vals <- sapply(X = hex.vals, FUN = strtoi, base = 16)
      alpha[match(x = hex[which(x = hex.length == 9)], table = colors)] <- dec.vals
    }
  }
  colors <- t(x = col2rgb(col = colors))
  colors <- mapply(
    FUN = function(i, alpha) {
      return(rgb(colors[i, , drop = FALSE], alpha = alpha, maxColorValue = 255))
    },
    i = 1:nrow(x = colors),
    alpha = alpha
  )
  return(colors)
}

violin_plot <- function(object, features,
                        split.by = NULL,
                        idents = NULL,
                        colors=NULL,
                        group.by=NULL) {
  require(ggplot2)
  require(Seurat)
  require(patchwork)
  require(formulaic)
  require(gplots)
  library(scales)
  library(ggridges)
  require(cowplot)

  if (sum(features %in% rownames(object)) != length(features)){
    print(paste('Feature(s)', features[!(features %in% rownames(object))], 'not present in the object'))
    features <- features[features %in% rownames(object)]
  }

  if(!is.null(idents)){
    if (!is.null(x = group.by)) {
      object <- SetIdent(object, value = group.by)
      object <- subset(object, idents = idents)
    }}

  #
  # select cells
  #
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }

  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- droplevels(factor(x = idents))
  }

  idents <- droplevels(factor(x = idents, levels = levels(idents)[order(levels(idents))]))
  if (is.null(x = colors)) {
    colors <- scales::hue_pal()(length(x = levels(x = idents)))
    colors <- alpha(colors, alpha = 0.5)
  } else {
    colors <- Col2Hex(colors)
  }

  y <- 'ident'
  xlab <- 'Expression Level'
  ylab <- 'Identity'

  data <- FetchData(object, vars = c(features, split.by, group.by),
                    cells = cells, slot = 'data')

  data[,group.by] <- factor( data[,group.by] , levels = levels(idents))
  if(is.null(group.by)){
    data <- cbind(data, object@active.ident)
    group.by <- 'object@active.ident'
  }
  if(is.null(split.by)){
    if(length(features) == 1){
      plot <- ggplot(data, aes_string(y=add.backtick(features),
                                      x=group.by,
                                      fill = group.by)) +
        geom_violin() +
        theme_classic() +
        theme(axis.text.x = element_text(angle=90, vjust = 0.5),
              axis.title.x = element_blank())
    }else{
      plot <- NULL
      list_plot <- vector(mode = "list", length = length(features))
      names(list_plot) <- features
      for(feature in features){
        if(feature != features[length(features)]){
          a <- ggplot(data, aes_string(y=add.backtick(feature),
                                       x=group.by,
                                       fill = group.by))  +
            geom_violin() +
            theme_classic() +
            theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank())
          list_plot[[feature]] <- a}
        else{
          a <- ggplot(data, aes_string(y=add.backtick(feature),
                                       x=group.by,
                                       fill = group.by))  +
            geom_violin() +
            theme_classic() +
            theme(axis.text.x = element_text(angle=90, vjust = 0.5),
                  axis.title.x = element_blank())
          list_plot[[feature]] <- a

        }
      }
      plot <- wrap_plots(list_plot, nrow = length(features),
                         guides = 'collect')
    }

  }
  #
  # split by!
  #
  if(!is.null(split.by)){
    if(length(split.by) == 1){
      #
      # ONE SPLIT.BY
      #
      colnames(data)[colnames(data) %in% split.by] <- make.unique(rep('split.by', length(split.by)))
      if(length(features) == 1){
        plot <- ggplot(data, aes_string(y=add.backtick(features),
                                        x=group.by, fill = group.by)) +
          geom_violin() +
          theme_classic() +
          theme(axis.text.x = element_text(angle=90, vjust = 0.5),
                axis.title.x = element_blank())+
          facet_grid(cols = vars(split.by), scales = 'free')
      }
      if(length(features) > 1){
        plot <- NULL
        list_plot <- vector(mode = "list", length = length(features))
        names(list_plot) <- features
        for(feature in features){
          if(feature == features[1]){
            a <-ggplot(data, aes_string(y=add.backtick(feature),
                                        x=group.by,
                                        fill = group.by)) +
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by), scales = 'free') +
              theme(axis.text.x = element_blank(),
                    axis.title.x = element_blank())

            list_plot[[feature]] <- a
          }
          if(feature != features[length(features)] &
             feature != features[1]){
            a <-ggplot(data, aes_string(y=add.backtick(feature),
                                        x=group.by,
                                        fill = group.by)) +
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by), scales = 'free')+
              theme(axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    strip.background = element_blank(),
                    strip.text = element_blank()
              )

            list_plot[[feature]] <- a

          }
          if(feature == features[length(features)]){
            a <-ggplot(data, aes_string(y=add.backtick(feature),
                                        x=group.by,
                                        fill = group.by)) +
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by), scales = 'free')+
              theme(axis.text.x = element_text(angle=90,
                                               vjust = 0.5),
                    axis.title.x = element_blank(),
                    strip.background = element_blank(),
                    strip.text = element_blank()
              )
            list_plot[[feature]] <- a
          }
        }
        plot <- wrap_plots(list_plot, nrow = length(features),
                           guides = 'collect')
      }
    }
    if(length(split.by) > 1){
      #
      # MULTIPLE SPLIT.BY
      #
      colnames(data)[colnames(data) %in% split.by] <- make.unique(rep('split.by', length(split.by)))
      if(length(features) == 1){
        plot <- ggplot(data, aes_string(y=add.backtick(features),
                                        x=group.by, fill = group.by)) +
          geom_violin() +
          theme_classic() +
          theme(axis.text.x = element_text(angle=90, vjust = 0.5),
                axis.title.x = element_blank())+
          facet_grid(cols = vars(split.by),
                     rows = vars(split.by.1),
                     scales = 'free')
      }
      if(length(features) > 1){
        plot <- NULL
        list_plot <- vector(mode = "list", length = length(features))
        names(list_plot) <- features
        for(feature in features){
          if(feature == features[1]){
            a <-ggplot(data, aes_string(y=add.backtick(feature),
                                        x=group.by,
                                        fill = group.by)) +
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by),
                         rows = vars(split.by.1),
                         scales = 'free')+
              theme(axis.text.x = element_blank(),
                    axis.title.x = element_blank())

            list_plot[[feature]] <- a
          }
          if(feature != features[length(features)] &
             feature != features[1]){
            a <-ggplot(data, aes_string(y=add.backtick(feature),
                                        x=group.by,
                                        fill = group.by)) +
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by),
                         rows = vars(split.by.1),
                         scales = 'free')+
              theme(axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    strip.background.x = element_blank(),
                    strip.text.x = element_blank()
              )

            list_plot[[feature]] <- a

          }
          if(feature == features[length(features)]){
            a <-ggplot(data, aes_string(y=add.backtick(feature),
                                        x=group.by,
                                        fill = group.by)) +
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by),
                         rows = vars(split.by.1),
                         scales = 'free')+
              theme(axis.text.x = element_text(angle=90,
                                               vjust = 0.5),
                    axis.title.x = element_blank(),
                    strip.background.x = element_blank(),
                    strip.text.x = element_blank()
              )
            list_plot[[feature]] <- a
          }
        }
        plot <- wrap_plots(list_plot, nrow = length(features),
                           guides = 'collect')
      }
    }

  }


  if(group.by == 'object@active.ident'){
    names(colors) <- NULL
    plot <- plot+
      scale_fill_manual(values = colors) +
      labs(y = 'Idents')+
      NoLegend()
  }
  return(plot)

}



ridge_plot <- function(object, features,
                       split.by = NULL,
                       idents = NULL,
                       colors=NULL,
                       group.by=NULL) {
  require(ggplot2)
  require(Seurat)
  require(patchwork)
  require(formulaic)
  require(gplots)
  library(scales)
  library(ggridges)
  require(cowplot)
  if (sum(features %in% rownames(object)) != length(features)){
    print(paste('Feature(s)', features[!(features %in% rownames(object))], 'not present in the object'))

  }
  if (length(features) > 1){
    print(paste('Only one gene please!'))
  }
  if(!is.null(idents)){
    if (!is.null(x = group.by)) {
      object <- SetIdent(object, value = group.by)
      object <- subset(object, idents = idents)
    }}

  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }


  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  idents <- factor(x = idents, levels = levels(idents)[order(levels(idents))])
  if (is.null(x = split.by)) {
    split <- NULL
  } else {
    split <- object[[split.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = split)) {
    split <- factor(x = split)
  }
  if (is.null(x = colors)) {
    colors <- hue_pal()(length(x = levels(x = idents)))
    colors <- alpha(colors, alpha = 0.5)
  } else if (length(x = colors) == 1 && colors == 'interaction') {
    split <- interaction(idents, split)
    colors <- hue_pal()(length(x = levels(x = idents)))
  } else {
    colors <- Col2Hex(colors)
  }
  if (length(x = colors) < length(x = levels(x = split))) {
    colors <- Interleave(colors, InvertHex(hexadecimal = colors))
  }
  if(length(split)>0){
    colors <- rep_len(x = colors, length.out = length(x = levels(x = split)))
    names(x = colors) <- levels(x = split)
  }
  if(length(split)==0){
    names(x = colors) <- levels(idents)
  }

  y <- 'ident'
  xlab <- 'Expression Level'
  ylab <- 'Identity'

  if(sum(features %in% rownames(object)) == length(features)){
    data <- FetchData(object, vars = c(features, split.by, group.by),
                      cells = cells, slot = 'data')
    data[,group.by] <- factor( data[,group.by] , levels = levels(idents))
    if(is.null(group.by)){
      data <- cbind(data, object@active.ident)
      group.by <- 'object@active.ident'
    }
    if(is.null(split.by)){
      plot <- ggplot(data, aes_string(x=features,
                                      y=group.by,
                                      fill = group.by)) +
        geom_density_ridges()+
        theme_cowplot() +
        labs(title = features)+
        scale_fill_manual(values = colors, labels = names(colors))
    }
    if(!is.null(split.by)){
      plot <- ggplot(data, aes_string(x=features,
                                      y=group.by, fill = split.by)) +
        geom_density_ridges()+
        theme_cowplot() +
        labs(title = features)  +
        scale_fill_manual(values = colors, labels = names(colors))
    }
  }
  if(group.by == 'object@active.ident'){
    names(colors) <- NULL
    plot <- plot+
      scale_fill_manual(values = colors) +
      labs(y = 'Idents')+
      NoLegend()
  }
  return(plot)

}


ridge_plot_2 <- function(object, features,
                         split.by = NULL,
                         idents = NULL,
                         colors=NULL,
                         group.by=NULL) {
  require(ggplot2)
  require(Seurat)
  require(patchwork)
  require(formulaic)
  require(gplots)
  library(scales)
  library(ggridges)
  require(cowplot)
  if (sum(features %in% rownames(object)) != length(features)){
    print(paste('Feature(s)', features[!(features %in% rownames(object))], 'not present in the object'))
    features <- features[features %in% rownames(object)]
  }

  if(!is.null(idents)){
    if (!is.null(x = group.by)) {
      object <- SetIdent(object, value = group.by)
      object <- subset(object, idents = idents)
    }}

  #
  # select cells
  #
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }

  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- droplevels(factor(x = idents))
  }

  idents <- droplevels(factor(x = idents, levels = levels(idents)[order(levels(idents))]))
  if (is.null(x = colors)) {
    colors <- scales::hue_pal()(length(x = levels(x = idents)))
    colors <- alpha(colors, alpha = 0.5)
  } else {
    colors <- Col2Hex(colors)
  }

  y <- 'ident'
  xlab <- 'Expression Level'
  ylab <- 'Identity'

  data <- FetchData(object, vars = c(features, split.by, group.by),
                    cells = cells, slot = 'data')

  data[,group.by] <- factor( data[,group.by] , levels = levels(idents))
  if(is.null(group.by)){
    data <- cbind(data, object@active.ident)
    group.by <- 'object@active.ident'
  }

  if(is.null(split.by)){
    if(length(features) == 1){
      plot <- ggplot(data, aes_string(x=features,
                                      y=group.by,
                                      fill = group.by)) +
        geom_density_ridges()+
        theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
        theme(axis.title.y = element_blank(),
              axis.title.x = element_blank())+
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        coord_cartesian(clip = "off") +
        labs(title = features)

    }else{
      plot <- NULL
      list_plot <- vector(mode = "list", length = length(features))
      names(list_plot) <- features
      for(feature in features){
        if(feature != features[length(features)]){
          a <- ggplot(data, aes_string(x=add.backtick(feature),
                                       y=group.by,
                                       fill = group.by))  +
            geom_density_ridges() +
            labs(title = feature) +
            theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.background = element_rect(fill = '#4287f500')) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            coord_cartesian(clip = "off")
          list_plot[[feature]] <- a}
        else{
          a <- ggplot(data, aes_string(x=add.backtick(feature),
                                       y=group.by,
                                       fill = group.by))  +
            geom_density_ridges() +
            labs(title = feature)+
            theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
            theme(
              axis.title.x = element_blank(),
              axis.title.y = element_blank())+
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            coord_cartesian(clip = "off")
          list_plot[[feature]] <- a

        }
      }
      plot <- wrap_plots(list_plot, ncol = length(features),
                         guides = 'collect')
    }

  }
  #
  # split by!
  #

  if(!is.null(split.by)){
    if(length(split.by) == 1){
      #
      # ONE SPLIT.BY
      #
      colnames(data)[colnames(data) %in% split.by] <- make.unique(rep('split.by', length(split.by)))
      if(length(features) == 1){
        plot <- ggplot(data, aes_string(x=add.backtick(features),
                                        y=group.by, fill = group.by)) +
          geom_density_ridges() +
          labs(title = features)+
          theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank())+
          facet_grid(cols = vars(split.by), scales = 'free')
      }
      if(length(features) > 1){
        plot <- NULL
        list_plot <- vector(mode = "list", length = length(features))
        names(list_plot) <- features
        for(feature in features){
          if(feature == features[1]){
            a <-ggplot(data, aes_string(x=add.backtick(feature),
                                        y=group.by,
                                        fill = group.by)) +
              geom_density_ridges() +
              labs(title = feature)+
              theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
              facet_grid(cols = vars(split.by), scales = 'free') +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill = '#4287f500'))

            list_plot[[feature]] <- a
          }
          if(feature != features[length(features)] &
             feature != features[1]){
            a <-ggplot(data, aes_string(x=add.backtick(feature),
                                        y=group.by,
                                        fill = group.by)) +
              geom_density_ridges() +
              labs(title = feature)+
              theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
              facet_grid(cols = vars(split.by), scales = 'free')+
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_blank(),
                    strip.text = element_blank()
              )

            list_plot[[feature]] <- a

          }
          if(feature == features[length(features)]){
            a <-ggplot(data, aes_string(x=add.backtick(feature),
                                        y=group.by,
                                        fill = group.by)) +
              geom_density_ridges() +
              labs(title = feature)+
              theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
              facet_grid(cols = vars(split.by), scales = 'free')+
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_blank(),
                    strip.text = element_blank()
              )
            list_plot[[feature]] <- a
          }
        }
        plot <- wrap_plots(list_plot, nrow = length(features),
                           guides = 'collect')
      }
    }
    if(length(split.by) > 1){
      #
      # MULTIPLE SPLIT.BY
      #
      colnames(data)[colnames(data) %in% split.by] <- make.unique(rep('split.by', length(split.by)))
      if(length(features) == 1){
        plot <- ggplot(data, aes_string(x=add.backtick(features),
                                        y=group.by, fill = group.by)) +
          geom_density_ridges() +
          labs(title = features)+
          theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank())+
          facet_grid(cols = vars(split.by),
                     rows = vars(split.by.1),
                     scales = 'free')
      }
      if(length(features) > 1){
        plot <- NULL
        list_plot <- vector(mode = "list", length = length(features))
        names(list_plot) <- features
        for(feature in features){
          if(feature == features[1]){
            a <-ggplot(data, aes_string(x=add.backtick(feature),
                                        y=group.by,
                                        fill = group.by)) +
              geom_density_ridges() +
              labs(title = feature)+
              theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
              facet_grid(cols = vars(split.by),
                         rows = vars(split.by.1),
                         scales = 'free')+
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill = '#4287f500'))

            list_plot[[feature]] <- a
          }
          if(feature != features[length(features)] &
             feature != features[1]){
            a <-ggplot(data, aes_string(x=add.backtick(feature),
                                        y=group.by,
                                        fill = group.by)) +
              geom_density_ridges() +
              labs(title = feature)+
              theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
              facet_grid(cols = vars(split.by),
                         rows = vars(split.by.1),
                         scales = 'free')+
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background.x = element_blank(),
                    strip.text.x = element_blank()
              )

            list_plot[[feature]] <- a

          }
          if(feature == features[length(features)]){
            a <-ggplot(data, aes_string(x=add.backtick(feature),
                                        y=group.by,
                                        fill = group.by)) +
              geom_density_ridges() +
              labs(title = feature)+
              theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
              facet_grid(cols = vars(split.by),
                         rows = vars(split.by.1),
                         scales = 'free')+
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background.x = element_blank(),
                    strip.text.x = element_blank()
              )
            list_plot[[feature]] <- a
          }
        }
        plot <- wrap_plots(list_plot, nrow = length(features),
                           guides = 'collect')
      }
    }

  }

  if(group.by == 'object@active.ident'){
    names(colors) <- NULL
    plot <- plot+
      scale_fill_manual(values = colors) +
      labs(y = 'Idents')+
      NoLegend()
  }

  plot <- plot &
    scale_x_continuous(expand = c(0, 0)) &
    scale_y_discrete(expand = c(0, 0)) &
    coord_cartesian(clip = "off") &
    theme(strip.background.y = element_rect(fill = '#4287f500'),
          strip.background.x = element_rect(fill = '#4287f500'))

  return(plot)

}

Interleave <- function(...) {
  return(as.vector(x = t(x = as.data.frame(x = list(...)))))
}

InvertHex <- function(hexadecimal) {
  return(vapply(
    X = toupper(x = hexadecimal),
    FUN = function(hex) {
      hex <- unlist(x = strsplit(
        x = gsub(pattern = '#', replacement = '', x = hex),
        split = ''
      ))
      key <- toupper(x = as.hexmode(x = 15:0))
      if (!all(hex %in% key)) {
        stop('All hexadecimal colors must be valid hexidecimal numbers from 0-9 and A-F')
      }
      if (length(x = hex) == 8) {
        alpha <- hex[7:8]
        hex <- hex[1:6]
      } else if (length(x = hex) == 6) {
        alpha <- NULL
      } else {
        stop("All hexidecimal colors must be either 6 or 8 characters in length, excluding the '#'")
      }
      value <- rev(x = key)
      inv.hex <- vapply(
        X = hex,
        FUN = function(x) {
          return(value[grep(pattern = x, x = key)])
        },
        FUN.VALUE = character(length = 1L)
      )
      inv.hex <- paste(inv.hex, collapse = '')
      return(paste0('#', inv.hex, paste(alpha, collapse = '')))
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = FALSE
  ))
}

