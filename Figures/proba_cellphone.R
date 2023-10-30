library(dplyr)

### DATA -----------------------------------------------------------------------
labels_values <- data.frame(
  stringsAsFactors = FALSE,
  Cluster.name = c("Cycling myeloids","IDA macrophages",
                   "DN EOMES","M2","M1",
                   "Inflammatory monocytes","IDA macrophages",
                   "DN EOMES","Inflammatory fibroblasts",
                   "Cycling myeloids","IDA macrophages","DN EOMES",
                   "M2"),
  X = c("RECEPTOR_R","RECEPTOR_R",
        "RECEPTOR_R","RECEPTOR_R","SENDER","SENDER",
        "SENDER","SENDER","SENDER",
        "RECEPTOR_NR","RECEPTOR_NR","RECEPTOR_NR",
        "RECEPTOR_NR"),
  Value = c(4,3,2,1,4,3.2,2.4,
            1.6,0.8,4,3,2,1)
)

mean_values <- data.frame(
  stringsAsFactors = FALSE,
  RECEPTOR_NR = c("Cycling myeloids",
                  "Cycling myeloids","IDA macrophages","DN EOMES","DN EOMES",
                  "DN EOMES","DN EOMES","DN EOMES","IDA macrophages",
                  "Cycling myeloids","M2","M2"),
  RECEPTOR_R = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  SENDER = c("M1","Inflammatory monocytes",
             "Inflammatory monocytes","Inflammatory monocytes","M1",
             "IDA macrophages","DN EOMES",
             "Inflammatory fibroblasts","IDA macrophages","IDA macrophages","M1",
             "DN EOMES"),
  MEAN = c(1.136,0.373,0.349,
           0.419,1.182,0.435,0.372,0.278,0.365,0.388,
           0.419,0.45),
  TX = c("PRE","PRE","PRE","PRE",
         "PRE","PRE","PRE","PRE","PRE","PRE","POST","POST")
)

# transform data ---------------------------------------------------------------

mean_values$RECEPTOR_NR <- plyr::mapvalues(
  x = mean_values$RECEPTOR_NR,
  from = labels_values$Cluster.name[labels_values$X == 'RECEPTOR_NR'],
  to =  labels_values$Value[labels_values$X == 'RECEPTOR_NR']
)
mean_values$RECEPTOR_R <- plyr::mapvalues(
  mean_values$RECEPTOR_R,
  from = labels_values$Cluster.name[labels_values$X == 'RECEPTOR_R'],
  to =  labels_values$Value[labels_values$X == 'RECEPTOR_R']
)
mean_values$SENDER <- plyr::mapvalues(
  mean_values$SENDER,
  from = labels_values$Cluster.name[labels_values$X == 'SENDER'],
  to =  labels_values$Value[labels_values$X == 'SENDER']
)

mean_values$interaction <- 1:nrow(mean_values)

mean_values2 <- reshape2::melt(mean_values[mean_values$TX == 'PRE',] %>%
                                 select(-TX, -MEAN) %>%
                                 type.convert(),
                               id = c('interaction'))

mean_values2$mean <- plyr::mapvalues(x = mean_values2$interaction,
                                     from = mean_values$interaction,
                                     to = mean_values$MEAN)

mean_values2$variable <- factor(mean_values2$variable,
                                levels = c('RECEPTOR_R', 'SENDER', 'RECEPTOR_NR')
                                )

# colors -----------------------------------------------------------------------
colores <- c('Cycling myeloids' = '#FEF1CA',
            'IDA macrophages' = '#F8C05F',
            'DN EOMES' = '#A28ED5',
            'M2' = '#269987',
            'M1' = '#6784A6',
            'Inflammatory monocytes' = '#EB7799',
            'Inflammatory fibroblasts' = '#ABC57B'
            )

## PLOT ------------------------------------------------------------------------


ggplot() +
  geom_line(data = mean_values2,
            inherit.aes = F,
            mapping = aes(x = variable,
                          y = value,
                          group = interaction,
                          color = mean), linewidth = 2)+
  geom_label(data = labels_values,
             family = 'Helvetica',
             inherit.aes = F,
             size = 6,
             mapping =  aes(x = X,
                            y = as.numeric(Value),
                            label = Cluster.name,
                            fill = Cluster.name),
             na.rm = F) +
  scale_fill_manual(values = colores) +
  labs(title = 'PRE-TX') +
  guides(fill =  "none") +
  theme(text = element_text(family = 'Helvetica'),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.title = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# proba2 -----------------------------------------------------------------------
cellinfo <- data.frame(
  stringsAsFactors = FALSE,
  RECEPTOR_NR = c("Cycling myeloids",
                  "Cycling myeloids","IDA macrophages","DN EOMES","DN EOMES",
                  "DN EOMES","DN EOMES","DN EOMES","IDA macrophages",
                  "Cycling myeloids","M2","M2",NA,NA,NA,NA,NA,NA,NA,NA,
                  NA,NA,NA,NA),
  RECEPTOR_R = c(NA,NA,NA,NA,NA,NA,NA,NA,
                 NA,NA,NA,NA,"Cycling myeloids","Cycling myeloids",
                 "IDA macrophages","DN EOMES","DN EOMES","DN EOMES",
                 "DN EOMES","DN EOMES","IDA macrophages",
                 "Cycling myeloids","M2","M2"),
  SENDER = c("M1","Inflammatory monocytes",
             "Inflammatory monocytes","Inflammatory monocytes","M1",
             "IDA macrophages","DN EOMES",
             "Inflammatory fibroblasts","IDA macrophages","IDA macrophages","M1",
             "DN EOMES","M1","Inflammatory monocytes",
             "Inflammatory monocytes","Inflammatory monocytes","M1","IDA macrophages",
             "DN EOMES","Inflammatory fibroblasts","IDA macrophages",
             "IDA macrophages","M1","DN EOMES"),
  MEAN = c(1.136,0.373,0.349,
           0.419,1.182,0.435,0.372,0.278,0.365,0.388,
           0.419,0.45,0,0,0,0,0,0,0,0,
           0,0,0,0),
  TX = c("PRE","PRE","PRE","PRE",
         "PRE","PRE","PRE","PRE","PRE","PRE","POST","POST",
         "PRE","PRE","PRE","PRE","PRE","PRE","PRE","PRE",
         "PRE","PRE","POST","POST"),
  INDIVIDUAL = c(1L,2L,3L,4L,5L,6L,7L,8L,
                 9L,10L,11L,12L,13L,14L,15L,16L,17L,18L,19L,
                 20L,21L,22L,23L,24L),
  PROP = c(rep(1,24))
)


cellinfo2 <- pivot_longer(cellinfo[cellinfo$TX == 'PRE',], cols = 1:3, values_to = 'poblacio', names_to = 'xaxis')
cellinfo2$xaxis <- factor(cellinfo2$xaxis, levels = c('RECEPTOR_R', 'SENDER', 'RECEPTOR_NR'))
cellinfo2 <- cellinfo2[!is.na(cellinfo2$poblacio),]
cellinfo2$PROP[cellinfo2$xaxis == 'SENDER' & cellinfo2$MEAN == 0] <- 0

ggplot(cellinfo2,
       aes(x = xaxis,
           stratum = poblacio,
           alluvium = INDIVIDUAL,
           y = PROP,
           fill = MEAN,
           label = poblacio)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5,
               fill = c("#F8C05F", # IDA macrophages
                        "#A28ED5", # DN
                        "#FEF1CA", # cy my
                        "#6784A6", # M1
                        "#EB7799", # Infl mono
                        "#ABC57B", # Infl fibro
                        "#F8C05F", # IDA macro
                        "#A28ED5", # DN
                        "#F8C05F", # IDA macrophages
                        "#A28ED5", # DN
                        "#FEF1CA" # cy my
               )
  ) +
  geom_text(stat = "stratum", size = 4) +
  scale_colour_gradient(
    low = "#EEFCFB00",
    high = "#21A59F",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  theme_void() +
  labs(title = 'PRE tx')

