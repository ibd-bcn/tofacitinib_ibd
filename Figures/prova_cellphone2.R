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


cellinfo2 <- pivot_longer(cellinfo[cellinfo$TX == 'PRE',],
                          cols = 1:3, values_to = 'poblacio',
                          names_to = 'xaxis')
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
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  scale_fill_gradient(
    low = "#132B4300",
    high = "red",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  theme_void()

ggplot(cellinfo2,
       aes(x = xaxis,
           stratum = poblacio,
           alluvium = INDIVIDUAL,
           y = PROP,
           fill = MEAN,
           label = poblacio)) +
  # scale_x_discrete(expand = c(-.1, -.1)) +
  geom_flow() +
  geom_stratum(alpha = 1,
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
  geom_text(stat = "stratum", size = 3) +
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
