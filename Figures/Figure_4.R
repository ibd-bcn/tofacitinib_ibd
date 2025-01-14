options(stringsAsFactors = FALSE, bitmapType = "cairo")
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readr)
library(dplyr)
library(ggrepel)
library(readxl)
source('Figures/functions_plots.R')

#Figure_4C----------------------------------------------------------------------
#Heatmap Macrophages C

DMSO <-
  read_excel(
    "Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
    sheet = "Figure 4C",
    col_types = c(
      "text",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric"
    ),
    skip = 1
  )

DMSO <- DMSO %>%
  mutate(across(where(is.numeric), log2))

colnames(DMSO)[1] <- "Condition"
row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
gene_vector <-
  c(
    "TNF",
    "SOCS3",
    "MX1",
    "CXCL1",
    "CXCL10",
    "CXCL8",
    "IL1B",
    "IL23A",
    "IL6",
    "INHBA",
    "CLEC5A"
  )
numeric_cols <- sapply(DMSO, is.numeric)

#LPS
LPS <- subset(DMSO, Condition == "LPS")
LPS <- data.frame(LPS, row.names = NULL)
rownames(LPS) <- c("LPS1", "LPS2", "LPS3", "LPS4", "LPS5")
LPS <- LPS[, 2:length(colnames(LPS))]
LPS <- colMeans(LPS)
#TNFa
TNFa <- subset(DMSO, Condition == "TNF")
TNFa <- data.frame(TNFa, row.names = NULL)
rownames(TNFa) <- c("TNF1", "TNF2", "TNF3", "TNF4", "TNF5")
TNFa <- TNFa[, 2:length(colnames(TNFa))]
TNFa <- colMeans(TNFa)
#IFNg
IFNG <- subset(DMSO, Condition == "IFNg")
IFNG <- data.frame(IFNG, row.names = NULL)
rownames(IFNG) <- c("IFNG1", "IFNG2", "IFNG3", "IFNG4", "IFNG5")
IFNG <- IFNG[, 2:length(colnames(IFNG))]
IFNG <- colMeans(IFNG)

t_DMSO <- t(data.frame(LPS = LPS,
                       TNFa = TNFa,
                       IFNG = IFNG))

t_DMSO <- t_DMSO[, gene_vector]
col_names <- colnames(t_DMSO)

#Create groups
grupos <-
  c(
    "IFN",
    "IFN",
    "IFN",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "JAK",
    "IF",
    "IF"
  )
col_groups <-
  factor(grupos, levels = c("IFN", "Inflam_cyt", "JAK", "IF"))
colnames(t_DMSO) <- paste0(colnames(t_DMSO), " ")
rownames(t_DMSO) <- paste0(rownames(t_DMSO), " ")

#Statistics
macros_dmso_stats <-
  as.data.frame(read_excel("Figures/extra_data/macros_dmso_stats.xlsx"))
genes_adjusted <- paste(gene_vector, "_p_adjusted", sep = "")
macros_dmso_stats <-
  macros_dmso_stats[, c("Condition", genes_adjusted)]
macros_dmso_stats <- unique(macros_dmso_stats)
rownames(macros_dmso_stats) <- macros_dmso_stats$Condition
macros_dmso_stats <-
  macros_dmso_stats[c("M-DMSO-LPS", "M-DMSO-TNF?", "M-DMSO-IFN?"),]
macros_dmso_stats <- macros_dmso_stats[, 2:ncol(macros_dmso_stats)]
sig_mat <-
  ifelse(macros_dmso_stats < 0.05,
         ifelse(
           macros_dmso_stats < 0.005,
           ifelse(
             macros_dmso_stats < 0.001,
             ifelse(macros_dmso_stats < 0.0001, "4", "3"),
             "2"
           ),
           "1"
         ),
         "")

## Obtain database
png(
  "Figures/output/DMSO.png",
  width = 10,
  height = 6,
  units = "in",
  res = 600
)
heatmap <- Heatmap(
  t_DMSO,
  na_col = "white",
  name = "Legend",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  row_split = as.factor(c("IFNG", "LPS", "TNFa")),
  row_title = NULL,
  row_names_gp = gpar(fontsize = 25),
  column_names_gp =  gpar(fontsize = 25, fontface = "italic"),
  column_names_rot = 60,
  border_gp = gpar(col = "white", lwd = 2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] == "1") {
      grid.text("*", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "2") {
      grid.text("**", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "3") {
      grid.text("***", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "4") {
      grid.text("****", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
  }


)

draw(heatmap)
dev.off()


#Figure_4E----------------------------------------------------------------------
#Heatmap Macrophages C
col_fun = colorRamp2(c(-3, 0, 3), c("green", "black", "red"))


TOFA <-
  read_excel(
    "Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
    sheet = "Figure 4E",
    col_types = c(
      "text",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric"
    ),
    skip = 1
  )

TOFA <- TOFA %>%
  mutate(across(where(is.numeric), log2))
colnames(TOFA)[1] <- "Condition"
row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
gene_vector <-
  c(
    "TNF",
    "SOCS3",
    "MX1",
    "CXCL1",
    "CXCL10",
    "CXCL8",
    "IL1B",
    "IL23A",
    "IL6",
    "INHBA",
    "CLEC5A"
  )
numeric_cols <- sapply(TOFA, is.numeric)

#LPS
LPS <- subset(TOFA, Condition == "LPS")
LPS <- data.frame(LPS, row.names = NULL)
rownames(LPS) <- c("LPS1", "LPS2", "LPS3", "LPS4", "LPS5")
LPS <- LPS[, 2:length(colnames(LPS))]
LPS <- colMeans(LPS)
#TNFa
TNFa <- subset(TOFA, Condition == "TNF")
TNFa <- data.frame(TNFa, row.names = NULL)
rownames(TNFa) <- c("TNF1", "TNF2", "TNF3", "TNF4", "TNF5")
TNFa <- TNFa[, 2:length(colnames(TNFa))]
TNFa <- colMeans(TNFa)
#IFNg
IFNG <- subset(TOFA, Condition == "IFNg")
IFNG <- data.frame(IFNG, row.names = NULL)
rownames(IFNG) <- c("IFNG1", "IFNG2", "IFNG3", "IFNG4", "IFNG5")
IFNG <- IFNG[, 2:length(colnames(IFNG))]
IFNG <- colMeans(IFNG)

t_TOFA <- t(data.frame(LPS = LPS,
                       TNFa = TNFa,
                       IFNG = IFNG))


t_TOFA <- t_TOFA[, gene_vector]
col_names <- colnames(t_TOFA)
#Create groups
grupos <-
  c(
    "IFN",
    "IFN",
    "IFN",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "JAK",
    "IF",
    "IF"
  )
col_groups <-
  factor(grupos, levels = c("IFN", "Inflam_cyt", "JAK", "IF"))
colnames(t_TOFA) <- paste0(colnames(t_TOFA), " ")
rownames(t_TOFA) <- paste0(rownames(t_TOFA), " ")

#AStatistics
macros_tofa_stats <-
  as.data.frame(read_excel("Figures/extra_data/macros_tofa_stats.xlsx"))
genes_adjusted <- paste(gene_vector, "_p_adjusted", sep = "")
macros_tofa_stats <-
  macros_tofa_stats[, c("Condition", genes_adjusted)]
macros_tofa_stats <- unique(macros_tofa_stats)
rownames(macros_tofa_stats) <- macros_tofa_stats$Condition
macros_tofa_stats <-
  macros_tofa_stats[c("M-TOFA-LPS", "M-TOFA-TNF?", "M-TOFA-IFN?"),]
macros_tofa_stats <- macros_tofa_stats[, 2:ncol(macros_tofa_stats)]
sig_mat <-
  ifelse(macros_tofa_stats < 0.05,
         ifelse(
           macros_tofa_stats < 0.005,
           ifelse(
             macros_tofa_stats < 0.001,
             ifelse(macros_tofa_stats < 0.0001, "4", "3"),
             "2"
           ),
           "1"
         ),
         "")


## Obtain database
png(
  "Figures/output/TOFA.png",
  width = 10,
  height = 6,
  units = "in",
  res = 600
)
heatmap <- Heatmap(
  t_TOFA,
  na_col = "white",
  name = "Legend",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  row_split = as.factor(c("IFNG", "LPS", "TNFa")),
  row_title = NULL,
  row_names_gp = gpar(fontsize = 25),
  column_names_gp =  gpar(fontsize = 25, fontface = "italic"),
  column_names_rot = 60,
  border_gp = gpar(col = "white", lwd = 2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] == "1") {
      grid.text("*", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "2") {
      grid.text("**", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "3") {
      grid.text("***", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "4") {
      grid.text("****", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
  }


)

draw(heatmap)
dev.off()

#Figure_4D----------------------------------------------------------------------
#Heatmap fibros D

col_fun = colorRamp2(c(-5, 0, 5), c("green", "black", "red"))

DMSO <-
  read_excel(
    "Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
    sheet = "Figure 4D",
    col_types = c(
      "text",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric"
    ),
    skip = 1
  )

DMSO <- DMSO %>%
  mutate(across(where(is.numeric), log2))
colnames(DMSO)[1] <- "Sample"
gene_vector <-
  c("TNF",
    "SOCS3",
    "CXCL1",
    "CXCL10",
    "CXCL5",
    "IL1B",
    "IL6",
    "CHI3L1",
    "INHBA",
    "WNT5A")
row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
ha = columnAnnotation(foo = anno_empty(border = FALSE,
                                       width = unit(0.5, "mm")))
numeric_cols <- sapply(DMSO, is.numeric)

#LPS
LPS <- subset(DMSO, Sample == "LPS")
LPS <- data.frame(LPS, row.names = NULL)
LPS <- LPS[, 2:length(colnames(LPS))]
LPS <- colMeans(LPS)
#TNFa
TNFa <- subset(DMSO, Sample == "TNF")
TNFa <- data.frame(TNFa, row.names = NULL)
TNFa <- TNFa[, 2:length(colnames(TNFa))]
TNFa <- colMeans(TNFa)
#IFNg
IFNG <- subset(DMSO, Sample == "IFNg")
IFNG <- data.frame(IFNG, row.names = NULL)
IFNG <- IFNG[, 2:length(colnames(IFNG))]
IFNG <- colMeans(IFNG)

t_DMSO <- t(data.frame(LPS = LPS,
                       TNFa = TNFa,
                       IFNG = IFNG))

t_DMSO <- t_DMSO[, gene_vector]
col_names <- colnames(t_DMSO)

# Crear groups
grupos <-
  c(
    "IFN",
    "IFN",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "JAK",
    "IF",
    "IF",
    "IF"
  )
col_groups <-
  factor(grupos, levels = c("IFN", "Inflam_cyt", "JAK", "IF"))
colnames(t_DMSO) <- paste0(colnames(t_DMSO), " ")
rownames(t_DMSO) <- paste0(rownames(t_DMSO), " ")

#Statistics
fibros_dmso_stats <-
  as.data.frame(read_excel("Figures/extra_data/fibros_dmso_stats.xlsx"))
genes_adjusted <- paste(gene_vector, "_p_adjusted", sep = "")
fibros_dmso_stats <-
  fibros_dmso_stats[, c("Condition", genes_adjusted)]
fibros_dmso_stats <- unique(fibros_dmso_stats)
rownames(fibros_dmso_stats) <- fibros_dmso_stats$Condition
fibros_dmso_stats <- fibros_dmso_stats[c("LPS", "TNF", "IFN"),]
fibros_dmso_stats <- fibros_dmso_stats[, 2:ncol(fibros_dmso_stats)]
sig_mat <-
  ifelse(fibros_dmso_stats < 0.05,
         ifelse(
           fibros_dmso_stats < 0.005,
           ifelse(
             fibros_dmso_stats < 0.001,
             ifelse(fibros_dmso_stats < 0.0001, "4", "3"),
             "2"
           ),
           "1"
         ),
         "")


png(
  "Figures/output/DMSO.png",
  width = 10,
  height = 6,
  units = "in",
  res = 600
)
heatmap <- Heatmap(
  t_DMSO,
  na_col = "white",
  name = "Legend",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  row_split = as.factor(c("IFNG", "LPS", "TNFa")),
  row_title = NULL,
  row_names_gp = gpar(fontsize = 25),
  column_names_gp =  gpar(fontsize = 25, fontface = "italic"),
  column_names_rot = 60,
  border_gp = gpar(col = "white", lwd = 2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] == "1") {
      grid.text("*", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "2") {
      grid.text("**", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "3") {
      grid.text("***", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "4") {
      grid.text("****", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
  }


)

draw(heatmap)
dev.off()

#Figure_4F----------------------------------------------------------------------
#Heatmap fibros F
col_fun = colorRamp2(c(-2, 0, 2), c("green", "black", "red"))

TOFA <-
  read_excel(
    "Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
    sheet = "Figure 4F",
    col_types = c(
      "text",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric"
    ),
    skip = 1
  )

TOFA <- TOFA %>%
  mutate(across(where(is.numeric), log2))
colnames(TOFA)[1] <- "Sample"
gene_vector <-
  c("TNF",
    "SOCS3",
    "CXCL1",
    "CXCL10",
    "CXCL5",
    "IL1B",
    "IL6",
    "CHI3L1",
    "INHBA",
    "WNT5A")
row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
ha = columnAnnotation(foo = anno_empty(border = FALSE,
                                       width = unit(0.5, "mm")))
numeric_cols <- sapply(TOFA, is.numeric)

#LPS
LPS <- subset(TOFA, Sample == "LPS")
LPS <- data.frame(LPS, row.names = NULL)
LPS <- LPS[, 2:length(colnames(LPS))]
LPS <- colMeans(LPS)
#TNFa
TNFa <- subset(TOFA, Sample == "TNF")
TNFa <- data.frame(TNFa, row.names = NULL)
TNFa <- TNFa[, 2:length(colnames(TNFa))]
TNFa <- colMeans(TNFa)
#IFNg
IFNG <- subset(TOFA, Sample == "IFNg")
IFNG <- data.frame(IFNG, row.names = NULL)
IFNG <- IFNG[, 2:length(colnames(IFNG))]
IFNG <- colMeans(IFNG)

t_TOFA <- t(data.frame(LPS = LPS,
                       TNFa = TNFa,
                       IFNG = IFNG))

t_TOFA <- t_TOFA[, gene_vector]
col_names <- colnames(t_TOFA)

#Create groups
grupos <-
  c(
    "IFN",
    "IFN",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "Inflam_cyt",
    "JAK",
    "IF",
    "IF",
    "IF"
  )
col_groups <-
  factor(grupos, levels = c("IFN", "Inflam_cyt", "JAK", "IF"))

colnames(t_TOFA) <- paste0(colnames(t_TOFA), " ")
rownames(t_TOFA) <- paste0(rownames(t_TOFA), " ")

#Statistics
fibros_tofa_stats <-
  as.data.frame(read_excel("Figures/extra_data/fibros_tofa_stats.xlsx"))
genes_adjusted <- paste(gene_vector, "_p_adjusted", sep = "")
fibros_tofa_stats <-
  fibros_tofa_stats[, c("Condition", genes_adjusted)]
fibros_tofa_stats <- unique(fibros_tofa_stats)
rownames(fibros_tofa_stats) <- fibros_tofa_stats$Condition
fibros_tofa_stats <-
  fibros_tofa_stats[c("LPS+TOFA", "TNF+TOFA", "IFN+TOFA"),]
fibros_tofa_stats <- fibros_tofa_stats[, 2:ncol(fibros_tofa_stats)]
sig_mat <-
  ifelse(fibros_tofa_stats < 0.05,
         ifelse(
           fibros_tofa_stats < 0.005,
           ifelse(
             fibros_tofa_stats < 0.001,
             ifelse(fibros_tofa_stats < 0.0001, "4", "3"),
             "2"
           ),
           "1"
         ),
         "")


png(
  "Figures/output/TOFA.png",
  width = 10,
  height = 6,
  units = "in",
  res = 600
)
heatmap <- Heatmap(
  t_TOFA,
  na_col = "white",
  name = "Legend",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  row_split = as.factor(c("IFNG", "LPS", "TNFa")),
  row_title = NULL,
  row_names_gp = gpar(fontsize = 25),
  column_names_gp =  gpar(fontsize = 25, fontface = "italic"),
  column_names_rot = 60,
  border_gp = gpar(col = "white", lwd = 2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] == "1") {
      grid.text("*", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "2") {
      grid.text("**", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "3") {
      grid.text("***", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
  }

)

draw(heatmap)
dev.off()

#Legends of Figure 4C-4D-4E-4F -------------------------------------------------
#Legend MACROS TOFA
col_fun = colorRamp2(c(-3, 0, 3), c("green", "black", "red"))

png(
  "Figures/output/legend_macs_TOFA.png",
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
lgd = Legend(
  col_fun = col_fun,
  direction = "vertical",
  legend_width = unit(7, "cm"),
  at = c(-3, 0, 3),
  labels = c(-3, 0, 3)
)
draw(lgd)
dev.off()

#Legend MACROS DMSO
col_fun = colorRamp2(c(-7, 0, 7), c("green", "black", "red"))

png(
  "Figures/output/legend_macs_DMSO.png",
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
lgd = Legend(
  col_fun = col_fun,
  direction = "vertical",
  legend_width = unit(30, "cm"),
  at = c(-7, 0, 7),
  labels = c("-7", "0", "7")
)
draw(lgd)
dev.off()

#Legend FIBROS TOFA
col_fun = colorRamp2(c(-2, 0, 2), c("green", "black", "red"))

png(
  "Figures/output/legend_fibros_TOFA.png",
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
lgd = Legend(
  col_fun = col_fun,
  direction = "vertical",
  legend_width = unit(7, "cm"),
  at = c(-2, 0, 2),
  labels = c(-2, 0, 2)
)
draw(lgd)
dev.off()

#Legend FIBROS DMSO
col_fun = colorRamp2(c(-5, 0, 5), c("green", "black", "red"))

png(
  "Figures/output/legend_fibros_DMSO.png",
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
lgd = Legend(
  col_fun = col_fun,
  direction = "vertical",
  legend_width = unit(7, "cm"),
  at = c(-5, 0, 5),
  labels = c(-5, 0, 5)
)
draw(lgd)
dev.off()
