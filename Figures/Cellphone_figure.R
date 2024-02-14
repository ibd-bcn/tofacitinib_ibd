
cellinfo <- read.csv('Figures/extra_data/cellinfo.csv')

library(tidyr)

cellinfo2 <- pivot_longer(cellinfo[cellinfo$TX == 'PRE',],
                          cols = 1:3, values_to = 'poblacio',
                          names_to = 'xaxis')
cellinfo2$xaxis <- factor(cellinfo2$xaxis, levels = c('RECEPTOR_R', 'SENDER', 'RECEPTOR_NR'))
cellinfo2 <- cellinfo2[!is.na(cellinfo2$poblacio),]
cellinfo2$PROP[cellinfo2$xaxis == 'SENDER' & cellinfo2$MEAN == 0 & cellinfo2$poblacio != 'M2'] <- 0



ggplot(cellinfo2,
       aes(x = xaxis,
           stratum = poblacio,
           alluvium = INDIVIDUAL,
           y = PROP,
           fill = MEAN#,
           # label = poblacio
           )) +
  # scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  # geom_text(stat = "stratum", size = 3) +
  scale_fill_gradient(
    low = "#132B4300",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  theme_void()

pre <- ggplot(cellinfo2,
       aes(x = xaxis,
           stratum = poblacio,
           alluvium = INDIVIDUAL,
           y = PROP,
           fill = MEAN#,
           # label = poblacio
           )) +
  # scale_x_discrete(expand = c(-.1, -.1)) +
  geom_flow() +
  geom_stratum(alpha = 1,
               width = 0.5,
               fill = c(
                 "#EBEBEB"
                 # 'gray', # M2
                 # "#F8C05F", # IDA macrophages
                 # "#A28ED5", # DN
                 # "#FEF1CA", # cy my
                 # "gray", #M2
                 # "#6784A6", # M1
                 # "#EB7799", # Infl mono
                 # "#ABC57B", # Infl fibro
                 # "#F8C05F", # IDA macro
                 # "#A28ED5", # DN
                 # 'gray', #M2
                 # "#F8C05F", # IDA macrophages
                 # "#A28ED5", # DN
                 # "#FEF1CA" # cy my
               )
  ) +
  # geom_text(stat = "stratum", size = 3) +
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

cellinfo_post <- read.csv('Figures/extra_data/cellinfo2.csv')
cellinfo_post$xaxis <- factor(cellinfo_post$xaxis, levels = c('RECEPTOR_R', 'SENDER', 'RECEPTOR_NR'))

post <- ggplot(cellinfo_post,
              aes(x = xaxis,
                  stratum = poblacio,
                  alluvium = INDIVIDUAL,
                  y = PROP,
                  fill = MEAN#,
                  # label = poblacio
                  )) +
  geom_flow() +
  geom_stratum(alpha = 1,
               width = 0.5,
               fill = c(
                 "#EBEBEB"
                 # 'gray',    # M2
                 # "#F8C05F", # IDA macrophages
                 # "#A28ED5", # DN
                 # "#FEF1CA", # cy my
                 # "gray",    # M2
                 # "#6784A6", # M1
                 # "#EB7799", # Infl mono
                 # "#ABC57B", # Infl fibro
                 # "#F8C05F", # IDA macro
                 # "#A28ED5", # DN
                 # 'gray',    # M2
                 # "#F8C05F", # IDA macrophages
                 # "#A28ED5", # DN
                 # "#FEF1CA"  # cy my
               )
  ) +
  # geom_text(stat = "stratum", size = 3) +
  scale_colour_gradient(
    low = "#EEFCFB00",
    high = "#21A59F",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  theme_void() +
  labs(title = 'POST tx') +
  theme(axis.text.x = element_text(family = 'helvetica', size = 12))

pre / post
