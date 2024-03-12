options(stringsAsFactors = FALSE)
library(openxlsx)
library(ggplot2)
library(readr)
library(nlcor)
library(ggpubr)
library(patchwork)

## Figure 1 C: Severity Index---------------------------------------------------

#Colors R & NR
colors_response <- c('R' = '#70ADE6',
                     'NR' = '#FF8E47')
#Dataframe severity
data <- data.frame(
  val = c(
    2,
    2.05,
    2.5,
    3,
    1,
    1.95,
    1.5,
    NA,
    NA,
    NA,

    0,
    0.05,
    0.1,
    1.05,
    -0.1,
    0.95,
    -0.05,
    NA,
    NA,
    NA,

    1.5,
    2.5,
    1.55,
    1.45,
    2,
    1,
    1.9,
    1.95,
    2.1,
    2.05,

    2,
    3.1,
    2.45,
    1.95,
    3,
    2.9,
    2.95,
    2.55,
    3.05,
    2.05
  ),
  condition = c(
    rep("Pre_tx_R", 10),
    rep("Post_tx_R", 10),
    rep("Pre_tx_NR", 10),
    rep("Post_tx_NR", 10)
  ),
  pat = c(rep(
    c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10"), 2
  ), rep(
    c("11", "12", "13", "14", "15", "16", "17", "18", "19", "20"), 2
  )),

  colors = c(rep("R", 20), rep("NR", 20))
)
#Factor condition
data$condition <-
  factor(data$condition,
         c("Pre_tx_R", "Post_tx_R", "Pre_tx_NR", "Post_tx_NR"))
#Plot
png(
  "Figures/output/severity.png",
  width = 15,
  height = 10,
  units = "in",
  res = 1200
)
ggplot(data, aes(
  x = condition,
  y = val,
  col = colors,
  group = pat
)) +
  geom_line(color = "black", size = 1) +
  geom_point(size = 10) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(size = 1.5)
  ) +
  scale_color_manual(values = colors_response)
dev.off()

## Figure 1 E: Correlation plot-------------------------------------------------
corr_tofa <- read_delim("Figures/extra_data/corr_tofa_sangre.csv",
                        delim = ";", escape_double = FALSE,
                        locale = locale(decimal_mark = ",",
                                        grouping_mark = "."),
                        na = c("ND","n.d"), trim_ws = TRUE)

SOCS1 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, SOCS1)) +
  geom_point(aes(color = Response)) +
  geom_smooth(formula = y ~ log(x), method = 'lm') +
  scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
  labs(x = 'Tofa conc (ng/ml)') +
  theme_classic()+
  theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))

SOCS3 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, SOCS3)) +
  geom_point(aes(color = Response)) +
  geom_smooth(formula = y ~ log(x), method = 'lm') +
  scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
  labs(x = 'Tofa conc (ng/ml)') +
  theme_classic()+
  theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))

IRF1 <- ggplot(corr_tofa, aes(`Conc tofa (ng/ml)`, IRF1)) +
  geom_point(aes(color = Response)) +
  geom_smooth(formula = y ~ log(x), method = 'lm') +
  scale_color_manual(values = c('R' = '#70ADE6', 'NR' = '#FF8E47')) +
  labs(x = 'Tofa conc (ng/ml)') +
  theme_classic()+
  theme(legend.position = 'top', text = element_text(family = 'Helvetica', size = 8))


patchwork::wrap_plots(SOCS1, SOCS3, IRF1, guides = 'collect') &
  theme(legend.position = 'top', axis.text = element_text(colour = 'black', size = 8))
