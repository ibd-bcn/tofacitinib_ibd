library(readxl)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(dplyr)
options(bitmapType='cairo')

#Supplementary Figure 8C--------------------------------------------------------

#DMSO

#Color
col_fun = colorRamp2(c(-7,0,7), c("green", "black","red"))

DMSO <-
  read_excel(
    "Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
    sheet = "Sup figure 9C",
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
df_dmso_stats <- as.data.frame(read_excel("Figures/extra_data/df_dmso_stats.xlsx"))
genes_adjusted <- paste(gene_vector,"_p_adjusted",sep = "")
df_dmso_stats <- df_dmso_stats[,c("Condition",genes_adjusted)]
df_dmso_stats <- unique(df_dmso_stats)
rownames(df_dmso_stats) <- df_dmso_stats$Condition
df_dmso_stats <- df_dmso_stats[c("M-DMSO-LPS","M-DMSO-TNF?","M-DMSO-IFN?"),]
df_dmso_stats <- df_dmso_stats[,2:ncol(df_dmso_stats)]
sig_mat <- ifelse(df_dmso_stats < 0.05, ifelse(df_dmso_stats < 0.005,ifelse(df_dmso_stats < 0.001, ifelse(df_dmso_stats < 0.0001, "4", "3"), "2"),"1"), "")

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

#Supplementary Figure 9E--------------------------------------------------------
# TOFA
col_fun = colorRamp2(c(-1.5,0,1.5), c("green","black","red"))
TOFA <-
  read_excel(
    "Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
    sheet = "Sup figure 9D",
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
df_tofa_stats <- as.data.frame(read_excel("Figures/extra_data/df_tofa_stats.xlsx"))
genes_adjusted <- paste(gene_vector,"_p_adjusted",sep = "")
df_tofa_stats <- df_tofa_stats[,c("Condition",genes_adjusted)]
df_tofa_stats <- unique(df_tofa_stats)
rownames(df_tofa_stats) <- df_tofa_stats$Condition
df_tofa_stats <- df_tofa_stats[c("M-TOFA-LPS","M-TOFA-TNF?","M-TOFA-IFN?"),]
df_tofa_stats <- df_tofa_stats[,2:ncol(df_tofa_stats)]
sig_mat <- ifelse(df_tofa_stats < 0.05, ifelse(df_tofa_stats < 0.005,ifelse(df_tofa_stats < 0.001, ifelse(df_tofa_stats < 0.0001, "4", "3"), "2"),"1"), "")

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
#Legend Figure 9E---------------------------------------------------------------
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

#Legend Figure 9C---------------------------------------------------------------
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

# Supplementary Figure 9 D


df <-read_excel("Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
                  sheet = "sup 9D")


df <- df %>%
  dplyr::mutate(across(where(is.numeric), log2))

colnames(df)[1] <- "Sample"
colnames(df)[2] <- "Condition"

df

# Identify numeric columns
numeric_cols <- sapply(df, is.numeric)

#LPS
LPS <- subset(df, Condition == "LPS")
LPS <- data.frame(LPS, row.names = NULL)

#TNFa
TNFa <- subset(df, Condition == "TNF")
TNFa <- data.frame(TNFa, row.names = NULL)

#IFNg
IFNG <- subset(df, Condition == "IFN")
IFNG <- data.frame(IFNG, row.names = NULL)


results <- IFNG %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  IFNG[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  IFNG[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- LPS %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  LPS[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  LPS[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- TNFa %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  TNFa[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  TNFa[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}


heatmap_tt <- rbind(LPS, IFNG, TNFa)

heatmap_tt <- heatmap_tt %>% relocate(Sample)

openxlsx::write.xlsx(heatmap_tt, "~/stats_df_control_2519.xlsx", rowNames=T)


df <-read_excel("Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
                  sheet = "sup 9D",
                  col_types = c(
                    "text",
                    "text",
                    "numeric",
                    "numeric",
                    "numeric",
                    "numeric",
                    "numeric",
                    "numeric",
                    "numeric",
                    "numeric"
                  ),
                  skip = 0
)

df <- df %>%
  dplyr::mutate(across(where(is.numeric), log2))
colnames(df)[1] <- "Sample"
colnames(df)[2] <- "Condition"
# colnames(df)[1] <- "Condition"
row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
gene_vector <-
  c(
    "SOCS3", "CXCL10", "CXCL1", "IL6", "IL1B", "TNF", "CHI3L1", "WNT5A"

  )
numeric_cols <- sapply(df, is.numeric)

#LPS
df$Condition <- as.character(df$Condition)
LPS <- subset(df, Condition == "LPS")
LPS <- data.frame(LPS, row.names = NULL)
rownames(LPS) <- c("LPS1", "LPS2", "LPS3", "LPS4","LPS5", "LPS6", "LPS7")
LPS <- LPS[, 3:length(colnames(LPS))]
LPS <- colMeans(LPS)
#TNFa
TNFa <- subset(df, Condition == "TNF")
TNFa <- data.frame(TNFa, row.names = NULL)
rownames(TNFa) <- c("TNF1", "TNF2", "TNF3", "TNF4", "TNF5", "TNF6", "TNF7")
TNFa <- TNFa[, 3:length(colnames(TNFa))]
TNFa <- colMeans(TNFa)
#IFNg
IFNG <- subset(df, Condition == "IFN")
IFNG <- data.frame(IFNG, row.names = NULL)
rownames(IFNG) <- c("IFNG1", "IFNG2", "IFNG3", "IFNG4", "IFNG5", "IFNG6", "IFNG7")
IFNG <- IFNG[, 3:length(colnames(IFNG))]
IFNG <- colMeans(IFNG)

t_df <- t(data.frame(LPS = LPS,
                       TNFa = TNFa,
                       IFNG = IFNG))

t_df <- t_df[, gene_vector]
col_names <- colnames(t_df)

valid_genes <- intersect(gene_vector, colnames(t_df))
t_df <- t_df[, valid_genes, drop = FALSE]



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
colnames(t_df) <- paste0(colnames(t_df), " ")
rownames(t_df) <- paste0(rownames(t_df), " ")

#Statistics
df_df_stats <-
  as.data.frame(read_excel("~/stats_df_control_2519.xlsx"))
genes_adjusted <- paste(valid_genes, "_p_adjusted", sep = "")
df_df_stats <-
  df_df_stats[, c("Condition", genes_adjusted)]
df_df_stats <- unique(df_df_stats)
rownames(df_df_stats) <- df_df_stats$Condition
df_df_stats <-
  df_df_stats[c("LPS", "TNF", "IFN"),]
df_df_stats <- df_df_stats[, 2:ncol(df_df_stats)]
sig_mat <-
  ifelse(df_df_stats < 0.05,
         ifelse(
           df_df_stats < 0.005,
           ifelse(
             df_df_stats < 0.001,
             ifelse(df_df_stats < 0.0001, "4", "3"),
             "2"
           ),
           "1"
         ),
         "")
col_fun = colorRamp2(c(-5, 0, 5), c("green", "black", "red"))
## Obtain database
png(
  "~/df_control_new.png",
  width = 10,
  height = 6,
  units = "in",
  res = 600
)
heatmap <- Heatmap(
  t_df,
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
  heatmap_legend_param = list(
    at = c(-7, 0, 7),  # Force the legend range
    color_bar = "continuous"
  ),
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

# Supplementary Figure 9F


## TOFA
TOFA <-read_excel("Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
                  sheet = "Sup 9F")

TOFA <- TOFA %>%
  dplyr::mutate(across(where(is.numeric), log2))

colnames(TOFA)[1] <- "Sample"
colnames(TOFA)[2] <- "Condition"
TOFA
# Identify numeric columns in the elisa dataframe
numeric_cols <- sapply(TOFA, is.numeric)

#LPS
LPS <- subset(TOFA, Condition == "LPS+TOFA")
LPS <- data.frame(LPS, row.names = NULL)

#TNFa
TNFa <- subset(TOFA, Condition == "TNF+TOFA")
TNFa <- data.frame(TNFa, row.names = NULL)

#IFNg
IFNG <- subset(TOFA, Condition == "IFN+TOFA")
IFNG <- data.frame(IFNG, row.names = NULL)



results <- IFNG %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  IFNG[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  IFNG[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- LPS %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  LPS[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  LPS[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- TNFa %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  TNFa[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  TNFa[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}


heatmap_tt <- rbind(LPS, IFNG, TNFa)

heatmap_tt <- heatmap_tt %>% relocate(Sample)

openxlsx::write.xlsx(heatmap_tt, "~/df_tofa_stats_2519.xlsx", rowNames=T)

col_fun = colorRamp2(c(-2, 0, 2), c("green", "black", "red"))


TOFA <-read_excel("Figures/extra_data/Supporting data values Melon-Ardanaz et al.xlsx",
                  sheet = "Sup 9F")
colnames(TOFA)[1] <- "Sample"
colnames(TOFA)[2] <- "Condition"
TOFA <- TOFA %>%
  dplyr::mutate(across(where(is.numeric), log2))
row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
gene_vector <-
  c(
    "SOCS3", "CXCL10", "CXCL1", "IL6", "IL1B", "TNF", "CHI3L1", "WNT5A"

  )
numeric_cols <- sapply(TOFA, is.numeric)

#LPS
LPS <- subset(TOFA, Condition == "LPS+TOFA")
LPS <- data.frame(LPS, row.names = NULL)
rownames(LPS) <- c("LPS1", "LPS2", "LPS3", "LPS4", "LPS5", "LPS6", "LPS7")
LPS <- LPS[, 3:length(colnames(LPS))]
LPS <- colMeans(LPS)
#TNFa
TNFa <- subset(TOFA, Condition == "TNF+TOFA")
TNFa <- data.frame(TNFa, row.names = NULL)
rownames(TNFa) <- c("TNF1", "TNF2", "TNF3", "TNF4", "TNF5", "TNF6", "TNF7")
TNFa <- TNFa[, 3:length(colnames(TNFa))]
TNFa <- colMeans(TNFa)
#IFNg
IFNG <- subset(TOFA, Condition == "IFN+TOFA")
IFNG <- data.frame(IFNG, row.names = NULL)
rownames(IFNG) <- c("IFNG1", "IFNG2", "IFNG3", "IFNG4", "IFNG5", "IFNG6", "IFNG7")
IFNG <- IFNG[, 3:length(colnames(IFNG))]
IFNG <- colMeans(IFNG)

t_TOFA <- t(data.frame(LPS = LPS,
                       TNFa = TNFa,
                       IFNG = IFNG))


t_TOFA <- t_TOFA[, gene_vector]
col_names <- colnames(t_TOFA)

valid_genes <- intersect(gene_vector, colnames(t_TOFA))
t_TOFA <- t_TOFA[, valid_genes, drop = FALSE]
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
df_tofa_stats <-
  as.data.frame(read_excel("~/df_tofa_stats_2519.xlsx"))
genes_adjusted <- paste(valid_genes, "_p_adjusted", sep = "")
df_tofa_stats <-
  df_tofa_stats[, c("Condition", genes_adjusted)]
df_tofa_stats <- unique(df_tofa_stats)
rownames(df_tofa_stats) <- df_tofa_stats$Condition
df_tofa_stats <-
  df_tofa_stats[c("LPS+TOFA", "TNF+TOFA", "IFN+TOFA"),]

df_tofa_stats <- df_tofa_stats[, 2:ncol(df_tofa_stats)]
sig_mat <-
  ifelse(df_tofa_stats < 0.05,
         ifelse(
           df_tofa_stats < 0.005,
           ifelse(
             df_tofa_stats < 0.001,
             ifelse(df_tofa_stats < 0.0001, "4", "3"),
             "2"
           ),
           "1"
         ),
         "")


## Obtain database
png(
  "~/TOFA_df_new.png",
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
  heatmap_legend_param = list(
    at = c(-2, 0, 2),  # Force the legend range
    color_bar = "continuous"
  ),
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




