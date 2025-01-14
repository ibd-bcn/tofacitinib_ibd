library(readxl)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(dplyr)
options(bitmapType='cairo')

#Supplementary Figure 8A--------------------------------------------------------

#DMSO

#Color
col_fun = colorRamp2(c(-7,0,7), c("green", "black","red"))

DMSO <-
  read_excel(
    "Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
    sheet = "Sup figure 8A",
    skip = 1,
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
      "numeric"
    )
  )

colnames(DMSO)[1] <- "Condition"
row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
gene_vector <- c("IDO1", "IRF1", "OAS1","CCL5","IL10","SPP1","ACOD1","CD209","MMP9")

DMSO <- DMSO %>%
  mutate(across(where(is.numeric), log2))

# Identify numeric columns in the elisa dataframe
numeric_cols <- sapply(DMSO, is.numeric)

#LPS
LPS <- subset(DMSO, Condition == "LPS")
LPS <- data.frame(LPS, row.names = NULL)
rownames(LPS) <- c("LPS1","LPS2","LPS3","LPS4","LPS5")
LPS <- LPS[,2:length(colnames(LPS))]
LPS <- colMeans(LPS)
#TNFa
TNFa <- subset(DMSO, Condition == "TNF")
TNFa <- data.frame(TNFa, row.names = NULL)
rownames(TNFa) <- c("TNF1","TNF2","TNF3","TNF4","TNF5")
TNFa <- TNFa[,2:length(colnames(TNFa))]
TNFa <- colMeans(TNFa)
#IFNg
IFNG <- subset(DMSO, Condition == "IFNg")
IFNG <- data.frame(IFNG, row.names = NULL)
rownames(IFNG) <- c("IFNG1","IFNG2","IFNG3","IFNG4","IFNG5")
IFNG <- IFNG[,2:length(colnames(IFNG))]
IFNG <- colMeans(IFNG)

t_DMSO <- t(data.frame(
  LPS = LPS,
  TNFa = TNFa,
  IFNG = IFNG
))

t_DMSO <- t_DMSO[,gene_vector]
col_names <- colnames(t_DMSO)

#Create_Grups
colnames(t_DMSO) <- paste0(colnames(t_DMSO), " ")
rownames(t_DMSO) <- paste0(rownames(t_DMSO), " ")

#Statistics
macros_dmso_stats <- as.data.frame(read_excel("Figures/extra_data/macros_dmso_stats.xlsx"))
genes_adjusted <- paste(gene_vector,"_p_adjusted",sep = "")
macros_dmso_stats <- macros_dmso_stats[,c("Condition",genes_adjusted)]
macros_dmso_stats <- unique(macros_dmso_stats)
rownames(macros_dmso_stats) <- macros_dmso_stats$Condition
macros_dmso_stats <- macros_dmso_stats[c("M-DMSO-LPS","M-DMSO-TNF?","M-DMSO-IFN?"),]
macros_dmso_stats <- macros_dmso_stats[,2:ncol(macros_dmso_stats)]
sig_mat <- ifelse(macros_dmso_stats < 0.05, ifelse(macros_dmso_stats < 0.005,ifelse(macros_dmso_stats < 0.001, ifelse(macros_dmso_stats < 0.0001, "4", "3"), "2"),"1"), "")

## Obtain database
png("Figures/output/comp_DMSO.png",width=10,height=6,units="in",res=600)
heatmap <- Heatmap(t_DMSO,
                   na_col = "white",
                   name = "Legend",
                   col = col_fun,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_names_side = "left",
                   show_heatmap_legend = FALSE,
                   row_split = as.factor(c("IFNG","LPS","TNFa")),
                   row_title = NULL,
                   row_names_gp = gpar(fontsize = 25),
                   column_names_gp =  gpar(fontsize = 25, fontface ="italic"),
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

#Supplementary Figure 8B--------------------------------------------------------
# TOFA
col_fun = colorRamp2(c(-1.5,0,1.5), c("green","black","red"))
TOFA <-
  read_excel(
    "Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
    sheet = "Sup figure 8B",
    skip = 1,
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
      "numeric"
    )
  )
TOFA <- TOFA %>%
  mutate(across(where(is.numeric), log2))

colnames(TOFA)[1] <- "Condition"
row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
gene_vector <- c("IDO1", "IRF1", "OAS1","CCL5", "IL10","SPP1","ACOD1","CD209","MMP9")
numeric_cols <- sapply(TOFA, is.numeric)

#LPS
LPS <- subset(TOFA, Condition == "LPS")
LPS <- data.frame(LPS, row.names = NULL)
rownames(LPS) <- c("LPS1","LPS2","LPS3","LPS4","LPS5")
LPS <- LPS[,2:length(colnames(LPS))]
LPS <- colMeans(LPS)
#TNFa
TNFa <- subset(TOFA, Condition == "TNF")
TNFa <- data.frame(TNFa, row.names = NULL)
rownames(TNFa) <- c("TNF1","TNF2","TNF3","TNF4","TNF5")
TNFa <- TNFa[,2:length(colnames(TNFa))]
TNFa <- colMeans(TNFa)
#IFNg
IFNG <- subset(TOFA, Condition == "IFNg")
IFNG <- data.frame(IFNG, row.names = NULL)
rownames(IFNG) <- c("IFNG1","IFNG2","IFNG3","IFNG4","IFNG5")
IFNG <- IFNG[,2:length(colnames(IFNG))]
IFNG <- colMeans(IFNG)

t_TOFA <- t(data.frame(
  LPS = LPS,
  TNFa = TNFa,
  IFNG = IFNG
))

t_TOFA <- t_TOFA[,gene_vector]
col_names <- colnames(t_TOFA)

#Create Grups
colnames(t_TOFA) <- paste0(colnames(t_TOFA), " ")
rownames(t_TOFA) <- paste0(rownames(t_TOFA), " ")

#Statistics
macros_tofa_stats <- as.data.frame(read_excel("Figures/extra_data/macros_tofa_stats.xlsx"))
genes_adjusted <- paste(gene_vector,"_p_adjusted",sep = "")
macros_tofa_stats <- macros_tofa_stats[,c("Condition",genes_adjusted)]
macros_tofa_stats <- unique(macros_tofa_stats)
rownames(macros_tofa_stats) <- macros_tofa_stats$Condition
macros_tofa_stats <- macros_tofa_stats[c("M-TOFA-LPS","M-TOFA-TNF?","M-TOFA-IFN?"),]
macros_tofa_stats <- macros_tofa_stats[,2:ncol(macros_tofa_stats)]
sig_mat <- ifelse(macros_tofa_stats < 0.05, ifelse(macros_tofa_stats < 0.005,ifelse(macros_tofa_stats < 0.001, ifelse(macros_tofa_stats < 0.0001, "4", "3"), "2"),"1"), "")

png("Figures/output/comp_TOFA.png",width=10,height=6,units="in",res=600)
heatmap <- Heatmap(t_TOFA,
                   na_col = "white",
                   name = "Legend",
                   col = col_fun,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_names_side = "left",
                   show_heatmap_legend = FALSE,
                   row_split = as.factor(c("IFNG","LPS","TNFa")),
                   row_title = NULL,
                   row_names_gp = gpar(fontsize = 25),
                   column_names_gp =  gpar(fontsize = 25, fontface ="italic"),
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

#Legends------------------------------------------------------------------------
#Legend Figure 8B---------------------------------------------------------------
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("green", "black", "red"))

png(
  "Figures/output/supp_macs_TOFA.png",
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
lgd = Legend(
  col_fun = col_fun,
  direction = "vertical",
  legend_width = unit(7, "cm"),
  at = c(-1.5, 0, 1.5),
  labels = c(-1.5, 0, 1.5)
)
draw(lgd)
dev.off()

#Legend Figure 8A---------------------------------------------------------------
col_fun = colorRamp2(c(-7, 0, 7), c("green", "black", "red"))

png(
  "Figures/output/supp_macs_DMSO.png",
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
lgd = Legend(
  col_fun = col_fun,
  direction = "vertical",
  legend_width = unit(7, "cm"),
  at = c(-7, 0, 7),
  labels = c(-7, 0, 7)
)
draw(lgd)

dev.off()

